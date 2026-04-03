#' Build a TMB data list for FLicc
#'
#' Converts a user-supplied input list into the numeric data structures required
#' by the FLicc TMB objective function. The function validates the main inputs,
#' constructs observed length-frequency matrices by gear, derives length-bin
#' midpoints, sets selectivity indexing, standardizes relative catch weights,
#' builds Gauss-Laguerre quadrature inputs, and prepares optional prior vectors.
#'
#' @param dl A named list containing the FLicc model inputs. Expected elements
#'   include:
#'   \describe{
#'     \item{\code{fq}}{A named list of length-frequency vectors, one per gear.}
#'     \item{\code{LLB}}{Numeric vector of lower length-bin boundaries.}
#'     \item{\code{sel_fun}}{Character vector of selectivity functions by gear.
#'       Currently supported: \code{"logistic"} and \code{"dsnormal"}.}
#'     \item{\code{Catch}}{Numeric vector of relative catch proportions by gear.}
#'     \item{\code{gear_names}}{Optional character vector of gear names.}
#'     \item{\code{model_name}}{Character string used to identify the mortality
#'       model. If it contains \code{"length-inverse"}, inverse-length natural
#'       mortality is assumed.}
#'     \item{\code{ref_length}}{Reference length for inverse-length natural
#'       mortality.}
#'     \item{\code{a}, \code{b}}{Length-weight parameters used in post hoc SPR
#'       calculations.}
#'     \item{\code{L50}}{Length at 50 percent maturity.}
#'     \item{\code{GL}}{Number of Gauss-Laguerre quadrature points.}
#'     \item{\code{Linf}}{Numeric vector or scalar used for initialization and
#'       optionally prior specification. If length 2, the second value is treated
#'       as the prior SD on the log scale.}
#'     \item{\code{Mk}}{Numeric vector or scalar used for initialization and
#'       optionally prior specification. If length 2, the second value is treated
#'       as the prior SD on the log scale.}
#'   }
#'
#' @details
#' The returned list is intended for direct use in \code{TMB::MakeADFun()}.
#' Selectivity functions are internally coded as integers, with
#' \code{1 = logistic} and \code{2 = dsnormal}. Selectivity parameters are packed
#' into a single vector in C++, with \code{sm_start} providing zero-based start
#' indices by gear.
#'
#' Gauss-Laguerre quadrature nodes and weights are generated using
#' \code{statmod::gauss.quad()} and are used by the fishblicc-style
#' population-at-length recursion.
#'
#' @return A named list containing validated and transformed data for TMB,
#'   including:
#'   \describe{
#'     \item{\code{obs}}{Observed counts matrix, length bins by gears.}
#'     \item{\code{Lmid}, \code{LLB}}{Length-bin midpoints and lower bounds.}
#'     \item{\code{sel_type}, \code{sm_start}, \code{sm_n}}{Selectivity
#'       bookkeeping objects for the TMB template.}
#'     \item{\code{node}, \code{quad_wt}}{Gauss-Laguerre quadrature nodes and
#'       weights.}
#'     \item{\code{prior_mu}, \code{prior_sd}, \code{prior_code}}{Optional prior
#'       specification vectors.}
#'     \item{\code{Linf_init}, \code{Mk_init}}{Initial values for fixed or
#'       estimated biological parameters.}
#'   }
#'
#' @examples
#' \dontrun{
#' tmb_data <- data_tmb_flicc(dl)
#' str(tmb_data)
#' }
#'
#' @export

data_tmb_flicc <- function(dl) {
  stopifnot(is.list(dl))

  # ---- observed frequencies ----
  if (is.null(dl$fq))
    stop("dl$fq not found")

  gear_names <- names(dl$fq)
  if (is.null(gear_names)) {
    if (!is.null(dl$gear_names)) {
      gear_names <- dl$gear_names
    } else {
      gear_names <- paste0("gear", seq_along(dl$fq))
    }
  }

  obs <- do.call(cbind, lapply(dl$fq, as.numeric))
  colnames(obs) <- gear_names

  nlen  <- nrow(obs)
  ngear <- ncol(obs)

  # ---- lengths ----
  LLB <- as.numeric(dl$LLB)
  if (length(LLB) != nlen) stop("length(LLB) must equal number of rows in obs")

  step <- c(diff(LLB), tail(diff(LLB), 1))
  Lmid <- LLB + step / 2

  # ---- biology for SPR ----
  a_wl    <- if (!is.null(dl$a)) as.numeric(dl$a) else 1
  b_wl    <- if (!is.null(dl$b)) as.numeric(dl$b) else 3
  L50_mat <- if (!is.null(dl$L50)) as.numeric(dl$L50) else stats::median(Lmid)

  # ---- selectivity ----
  sel_fun <- dl$sel_fun
  if (length(sel_fun) != ngear)
    stop("length(sel_fun) must equal number of gears")

  sel_type <- ifelse(sel_fun == "logistic", 1L,
                     ifelse(sel_fun == "dsnormal", 2L, NA_integer_))
  if (any(is.na(sel_type))) {
    stop("Only logistic and dsnormal supported in this TMB prototype")
  }

  sm_n <- ifelse(sel_fun == "logistic", 2L, 3L)
  sm_start <- cumsum(c(0L, head(sm_n, -1L)))  # already 0-based for C++

  stopifnot(all(sm_start >= 0))
  stopifnot(all(sm_start + sm_n <= sum(sm_n)))

  # ---- catch weights ----
  catch_wt <- as.numeric(dl$Catch)
  if (length(catch_wt) != ngear)
    stop("length(Catch) must equal number of gears")
  catch_wt <- catch_wt / sum(catch_wt)

  # ---- mortality model ----
  M_model <- if (grepl("length-inverse", tolower(dl$model_name))) 2L else 1L
  ref_length <- if (!is.null(dl$ref_length)) as.numeric(dl$ref_length) else mean(Lmid)

  # ---- quadrature inputs for pop_len_glq ----
  # Adjust these names if fishblicc stores them differently
  gl <- statmod::gauss.quad(dl$GL, kind = "laguerre")

  node <- gl$nodes
  quad_wt <- gl$weights

  if (is.null(node) || is.null(quad_wt)) {
    stop("dl must contain quadrature vectors 'node' and 'quad_wt' for pop_len_glq()")
  }

  node <- as.numeric(node)
  quad_wt <- as.numeric(quad_wt)

  if (length(node) != length(quad_wt)) {
    stop("node and quad_wt must have the same length")
  }

  # ---- fixed-value / init extraction ----
  # Linf can be c(mean, sd) or similar
  if (is.null(dl$Linf))
    stop("dl$Linf not found")

  Linf_init <- as.numeric(dl$Linf[1])

  # Mk can be scalar or c(mean, sd)
  if (is.null(dl$Mk))
    stop("dl$Mk not found")

  Mk_init <- as.numeric(dl$Mk[1])

  # ---- simple prior scaffold on transformed scale ----
  prior_mu <- numeric()
  prior_sd <- numeric()
  prior_code <- integer()

  # codes:
  # 1 log_Linf
  # 2 log_Galpha
  # 3 log_Mk
  # 4 log_phi
  # 100+g log_Fk[g]
  # 200+j Sm[j]

  # only add Linf prior if uncertainty element exists
  if (length(dl$Linf) >= 2 && !is.na(dl$Linf[2])) {
    prior_mu   <- c(prior_mu, log(Linf_init))
    prior_sd   <- c(prior_sd, as.numeric(dl$Linf[2]))
    prior_code <- c(prior_code, 1L)
  }

  # only add Mk prior if uncertainty element exists
  if (length(dl$Mk) >= 2 && !is.na(dl$Mk[2])) {
    prior_mu   <- c(prior_mu, log(Mk_init))
    prior_sd   <- c(prior_sd, as.numeric(dl$Mk[2]))
    prior_code <- c(prior_code, 3L)
  }

  list(
    nlen       = nlen,
    ngear      = ngear,
    Lmid       = as.numeric(Lmid),
    LLB        = as.numeric(LLB),
    obs        = unname(obs),
    catch_wt   = catch_wt,
    sel_type   = as.integer(sel_type),
    sm_start   = as.integer(sm_start),
    sm_n       = as.integer(sm_n),
    ref_length = as.numeric(ref_length),
    M_model    = as.integer(M_model),
    a_wl       = a_wl,
    b_wl       = b_wl,
    L50_mat    = L50_mat,
    node       = node,
    quad_wt    = quad_wt,
    Linf_init  = Linf_init,
    Mk_init    = Mk_init,
    prior_mu   = as.numeric(prior_mu),
    prior_sd   = as.numeric(prior_sd),
    prior_code = as.integer(prior_code),
    gear_names = gear_names
  )
}

#' Generate initial parameter values for FLicc TMB fits
#'
#' Builds a named parameter list suitable for \code{TMB::MakeADFun()} from a
#' TMB data list produced by \code{data_tmb_flicc()}. Initial values are derived
#' from the data where possible, using observed modal length by gear as the
#' starting point for selectivity location parameters.
#'
#' @param tmb_data A list produced by \code{data_tmb_flicc()}.
#'
#' @details
#' Initial values are currently set as follows:
#' \itemize{
#'   \item \code{log_Linf}: initialized from \code{tmb_data$Linf_init}
#'   \item \code{log_Mk}: initialized from \code{tmb_data$Mk_init}
#'   \item \code{log_Galpha}: initialized to \code{log(20)}
#'   \item \code{log_Fk}: initialized to \code{log(0.5)} for all gears
#'   \item \code{log_phi}: initialized to \code{log(20)}
#'   \item \code{Sm}: selectivity parameters initialized by gear using the
#'     observed modal length and default spread values
#' }
#'
#' For logistic selectivity, the packed parameter block is
#' \code{c(L50, log_steep)}.
#' For double-sided normal selectivity, the packed parameter block is
#' \code{c(mode, log_lsd, log_rsd)}.
#'
#' @return A named list of parameter vectors and scalars for use in
#'   \code{TMB::MakeADFun()}.
#'
#' @examples
#' \dontrun{
#' tmb_data <- data_tmb_flicc(dl)
#' par0 <- init_tmb_flicc(tmb_data)
#' str(par0)
#' }
#'
#' @export
init_tmb_flicc <- function(tmb_data) {
  ngear    <- tmb_data$ngear
  sel_type <- tmb_data$sel_type
  Lmid     <- tmb_data$Lmid
  obs      <- tmb_data$obs

  Sm <- numeric(sum(tmb_data$sm_n))
  pos <- 1

  for (g in seq_len(ngear)) {
    obs_g <- obs[, g]
    mode_g <- Lmid[which.max(obs_g)]

    if (sel_type[g] == 1L) {
      # logistic: L50, log_steep
      Sm[pos:(pos + 1)] <- c(mode_g, log(0.3))
      pos <- pos + 2

    } else if (sel_type[g] == 2L) {
      # dsnormal: mode, log_lsd, log_rsd
      Sm[pos:(pos + 2)] <- c(mode_g, log(5), log(5))
      pos <- pos + 3
    }
  }

  list(
    log_Linf   = log(tmb_data$Linf_init),
    log_Galpha = log(20),
    log_Mk     = log(tmb_data$Mk_init),
    log_Fk     = rep(log(0.5), ngear),
    Sm         = Sm,
    log_phi    = log(20)
  )
}


#' Fit the FLicc TMB model
#'
#' Compiles, loads, and fits the FLicc TMB objective function to a single
#' equilibrium multi-gear length dataset. The function prepares the TMB data and
#' parameter lists, optionally recompiles the C++ template, fixes selected
#' biological parameters via a TMB map, runs numerical optimization, and returns
#' fitted parameters together with reported model outputs.
#'
#' @param dl A named input list defining the FLicc model data and settings.
#'   Passed to \code{data_tmb_flicc()}.
#' @param compile Logical. If \code{TRUE}, the C++ template is recompiled before
#'   fitting. Default is \code{TRUE}.
#' @param silent Logical. Passed to \code{TMB::MakeADFun()}. Default is
#'   \code{TRUE}.
#' @param dll Character string giving the base name of the TMB dynamic library.
#'   Default is \code{"fishblicc_tmb"}.
#'
#' @details
#' The current implementation fixes \code{log_Linf} and \code{log_Mk} through a
#' TMB map. Optimization is performed with \code{stats::nlminb()}, and standard
#' errors are obtained from \code{TMB::sdreport()}.
#'
#' Returned reported quantities depend on the C++ template, but typically include
#' fitted fishing mortalities by gear, selectivity, mortality-at-length,
#' population-at-length, expected counts, and SPR-related outputs.
#'
#' @return An object of class \code{"blicc_tmb_fit"}, a list containing:
#'   \describe{
#'     \item{\code{opt}}{The \code{nlminb} optimization result.}
#'     \item{\code{obj}}{The TMB ADFun object.}
#'     \item{\code{rep}}{The \code{sdreport} object.}
#'     \item{\code{par}}{A tibble of fixed-effect estimates and standard errors.}
#'     \item{\code{report}}{A named list of reported quantities from the C++
#'       template.}
#'     \item{\code{tmb_data}}{The TMB data list used for the fit.}
#'   }
#'
#' @examples
#' \dontrun{
#' fit <- fiticc(dl, compile = FALSE)
#' fit$par
#' fit$report$Fk
#' fit$report$spr
#' }
#'
#' @seealso \code{\link{data_tmb_flicc}}, \code{\link{init_tmb_flicc}}
#'
#' @export
fiticc <- function(dl,
                   compile = FALSE,
                   silent = TRUE,
                   dll = "FLicc") {

  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Package 'TMB' is required")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required")
  }

  tmb_data <- data_tmb_flicc(dl)
  parameters <- init_tmb_flicc(tmb_data)

  dll_path <- normalizePath(file.path("src", paste0(dll, ".dll")),
                            winslash = "/", mustWork = FALSE)

  if (compile) {
    if (file.exists(dll_path)) {
      try(dyn.unload(dll_path), silent = TRUE)
    }
    if (file.exists(file.path("src", paste0(dll, ".o")))) {
      unlink(file.path("src", paste0(dll, ".o")))
    }
    if (file.exists(file.path("src", paste0(dll, ".dll")))) {
      unlink(file.path("src", paste0(dll, ".dll")))
    }

    TMB::compile(file.path("src", paste0(dll, ".cpp")))
  }

  dll_path <- normalizePath(file.path("src", paste0(dll, ".dll")),
                            winslash = "/")

  if (!is.loaded(paste0("R_init_", dll))) {
    dyn.load(dll_path)
  }

  map <- list(
    log_Linf = factor(NA),
    log_Mk   = factor(NA)
  )

  obj <- TMB::MakeADFun(
    data = tmb_data[setdiff(names(tmb_data), c("gear_names", "Linf_init", "Mk_init"))],
    parameters = parameters,
    map = map,
    DLL = dll,
    silent = silent
  )

  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(iter.max = 1000, eval.max = 1000)
  )

  obj$env$last.par.best <- opt$par

  rep <- TMB::sdreport(obj)
  adrep <- summary(rep, "fixed")

  par_tab <- data.frame(
    par = rownames(adrep),
    mpd = adrep[, 1],
    se  = adrep[, 2],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  rep_list <- obj$report()

  structure(
    list(
      opt = opt,
      obj = obj,
      rep = rep,
      par = par_tab,
      report = rep_list,
      tmb_data = tmb_data
    ),
    class = "flicc_tmb_fit"
  )
}
