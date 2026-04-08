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
#'     \item{\code{catch_sd}}{Optional scalar controlling the strength of the
#'       log-ratio penalty linking modeled gear proportions to relative catch
#'       proportions. Smaller values imply stronger anchoring. Default is
#'       \code{0.01}.}
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
#'     \item{\code{prior_mu}, \code{prior_sd}, \code{prior_code},\code{catch_sd}}{Optional prior
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

  # ---- catch-weight penalty SD ----
  catch_sd <- if (!is.null(dl$catch.sd)) as.numeric(dl$catch.sd) else 0.05

  if (length(catch_sd) != 1 || !is.finite(catch_sd) || catch_sd <= 0) {
    stop("dl$catch_sd must be a single positive numeric value")
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



  # ---- selectivity ----
  sel_fun <- dl$sel_fun
  if (length(sel_fun) != ngear)
    stop("length(sel_fun) must equal number of gears")

  sel_type <- ifelse(sel_fun == "logistic", 1L,
                     ifelse(sel_fun == "dsnormal", 2L,
                            ifelse(sel_fun %in% c("normal"), 3L,
                                   NA_integer_)))

  if (any(is.na(sel_type))) {
    stop("Supported selectivity options are: logistic, dsnormal, normal")
  }

  sm_n <- ifelse(sel_fun == "logistic", 2L,
                 ifelse(sel_fun == "dsnormal", 3L,
                        ifelse(sel_fun %in% c("normal"), 2L,
                               NA_integer_)))
  sm_start <- cumsum(c(0L, head(sm_n, -1L)))  # already 0-based for C++

  stopifnot(all(sm_start >= 0))
  stopifnot(all(sm_start + sm_n <= sum(sm_n)))

  # ---- catch weights ----
  catch_wt <- as.numeric(dl$catch_by_gear)
  if (length(catch_wt) != ngear)
    stop("length(Catch) must equal number of gears")
  catch_wt <- catch_wt / sum(catch_wt)


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

  # ---- biological inputs ----

  if (is.null(dl$Linf)) stop("dl$Linf not found")
  Linf_init <- as.numeric(dl$Linf[1])

  if (is.null(dl$Mk)) stop("dl$Mk not found")
  Mk_init <- as.numeric(dl$Mk[1])

  if (is.null(dl$CVL)) {
    CVL_init <- 0.1
    CVL_input <- 0.1
  } else {
    CVL_input <- dl$CVL
    CVL_init <- as.numeric(CVL_input[1])
  }

  if (Linf_init <= 0) stop("Linf must be > 0")
  if (Mk_init <= 0) stop("Mk must be > 0")
  if (CVL_init <= 0) stop("CVL must be > 0")

  Galpha_init <- 1 / CVL_init^2


  prior_mu <- numeric()
  prior_sd <- numeric()
  prior_code <- integer()

  # 1 log_Linf
  if (!is.null(dl$linf.sd)) {
    prior_mu   <- c(prior_mu, log(Linf_init))
    prior_sd   <- c(prior_sd, as.numeric(dl$linf.sd))
    prior_code <- c(prior_code, 1L)
  }

  # 3 log_Mk
  if (!is.null(dl$Mk.sd)) {
    prior_mu   <- c(prior_mu, log(Mk_init))
    prior_sd   <- c(prior_sd, as.numeric(dl$Mk.sd))
    prior_code <- c(prior_code, 3L)
  }

  # 2 log_Galpha, derived from CVL prior
  if (!is.null(dl$CVL.sd)){
    prior_mu   <- c(prior_mu, -2 * log(CVL_init))
    prior_sd   <- c(prior_sd,  2 * as.numeric(dl$CVL.sd))
    prior_code <- c(prior_code, 2L)
  }

  Mscaler <- dl$mL/mean(dl$mL)

  list(
    nlen       = nlen,
    ngear      = ngear,
    Lmid       = as.numeric(Lmid),
    LLB        = as.numeric(LLB),
    obs        = unname(obs),
    catch_wt   = catch_wt,
    catch_sd   = catch_sd,
    sel_type   = as.integer(sel_type),
    sm_start   = as.integer(sm_start),
    sm_n       = as.integer(sm_n),
    Mscaler   =  as.numeric(Mscaler),
    mat         = as.numeric(dl$mat),
    wt         = as.numeric(dl$wt),
    node       = node,
    quad_wt    = quad_wt,
    Linf_init   = Linf_init,
    Mk_init     = Mk_init,
    CVL_init    = CVL_init,
    Galpha_init = Galpha_init,
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
#' Selectivity parameters are stored in a packed vector \code{Sm}, with the
#' meaning depending on selectivity type.
#'
#' For logistic selectivity, the packed parameter block is
#' \code{c(log_SL50_rel, log_dSL_rel)}, where \code{SL50 = exp(log_SL50_rel) * Linf}
#' and \code{SL95 = SL50 + exp(log_dSL_rel) * Linf}.
#'
#' For normal selectivity, the packed parameter block is
#' \code{c(log_mode_rel, log_sd_rel)}, where \code{mode = exp(log_mode_rel) * Linf}
#' and \code{sd = exp(log_sd_rel) * Linf}.
#'
#' For double-sided normal selectivity, the packed parameter block is
#' \code{c(log_mode_rel, log_lsd_rel, log_rsd_rel)}, where
#' \code{mode = exp(log_mode_rel) * Linf}, \code{lsd = exp(log_lsd_rel) * Linf},
#' and \code{rsd = exp(log_rsd_rel) * Linf}.
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
  Linf     <- tmb_data$Linf_init

  Sm <- numeric(sum(tmb_data$sm_n))
  pos <- 1

  for (g in seq_len(ngear)) {
    obs_g  <- obs[, g]
    mode_g <- Lmid[which.max(obs_g)]
    mode_rel <- max(mode_g / Linf, 1e-4)

    if (sel_type[g] == 1L) {
      # logistic: c(log_SL50_rel, log_dSL_rel)
      SL50_rel <- mode_rel
      dSL_rel  <- 0.10
      Sm[pos:(pos + 1)] <- c(log(SL50_rel), log(dSL_rel))
      pos <- pos + 2

    } else if (sel_type[g] == 2L) {
      # dsnormal: c(log_mode_rel, log_lsd_rel, log_rsd_rel)
      lsd_rel <- 0.10
      rsd_rel <- 0.10
      Sm[pos:(pos + 2)] <- c(log(mode_rel), log(lsd_rel), log(rsd_rel))
      pos <- pos + 3

    } else if (sel_type[g] == 3L) {
      # normal: c(log_mode_rel, log_sd_rel)
      sd_rel <- 0.10
      Sm[pos:(pos + 1)] <- c(log(mode_rel), log(sd_rel))
      pos <- pos + 2
    }
  }

  list(
    log_Linf   = log(tmb_data$Linf_init),
    log_Galpha = log(tmb_data$Galpha_init),
    log_Mk     = log(tmb_data$Mk_init),
    log_Fk     = rep(log(0.5), ngear),
    Sm         = Sm,
    log_phi    = log(10)
  )
}

#' Fit the FLicc TMB model
#'
#' Compiles, loads, and fits the FLicc TMB objective function to a single
#' equilibrium multi-gear length dataset.
#'
#' @param dl A named input list defining the FLicc model data and settings.
#' @param compile Logical. If \code{TRUE}, the C++ template is recompiled before
#'   fitting. Default is \code{FALSE}.
#' @param silent Logical. Passed to \code{TMB::MakeADFun()}. Default is
#'   \code{TRUE}.
#' @param dll Character string giving the base name of the TMB dynamic library.
#'   Default is \code{"FLicc"}.
#' @param sel_fixed Optional named list of fixed selectivity values passed to
#'   \code{apply_sel_fixed_flicc()}.
#'
#' @return An object of class \code{"flicc_tmb_fit"}.
#'
#' @export
fiticc <- function(lfd, stklen,
                   sel_fun = c("logistic","normal","dsnormal"),
                   catch_by_gear,
                   settings = list(CVL=0.1,GL=50,catch.sd=0.05),
                   compile = FALSE,
                   silent = TRUE,
                   dll = "FLicc",
                   sel_fixed = NULL){

  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Package 'TMB' is required")
  }



  lens<- an(dimnames(lfd[[1]])$len)
  yrs = an(dimnames(lfd[[1]])$year)
  nyrs = length(yrs)
  lhpar <- stklen@lhpar
  syrs = an(dimnames(stklen)$year)

  i = nyrs

  stkl <- window(stklen,start=syrs[i],end=syrs[i])

  fq <- lapply(lfd,function(x){
    (as.data.frame(yearMeans(x[,ac(yrs[i])]))$data)
  })
  pars <- dimnames(lhpar)$params

  dl <- list(
    LLB = lens,
    fq = fq,
    sel_fun = sel_fun,
    catch_by_gear= catch_by_gear,
    gear_names = names(fq),
    Linf = an(lhpar["linf"]),
    Mk = an(lhpar["Mk"]),
    wt = as.data.frame(catch.wt(stkl))$data,
    mat= as.data.frame(mat(stkl))$data,
    mL = as.data.frame(m(stkl))$data,
    CVL=settings$CVL,
    GL = settings$GL,
    catch.sd = settings$catch.sd,
    linf.sd =  settings$linf.sd,
    Mk.sd =  settings$Mk.sd
  )

  tmb_data <- data_tmb_flicc(dl)
  parameters <- init_tmb_flicc(tmb_data)

  # Apply optional fixed selectivity settings
  fix_res <- apply_sel_fixed_flicc(
    parameters = parameters,
    tmb_data   = tmb_data,
    sel_fixed  = sel_fixed
  )
  parameters <- fix_res$parameters

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

  map <- list()

  # Fixed by default; estimated only if a log.sd is supplied
  if (length(dl$Linf) < 2 || is.na(dl$Linf[2])) {
    map$log_Linf <- factor(NA)
  }

  if (length(dl$Mk) < 2 || is.na(dl$Mk[2])) {
    map$log_Mk <- factor(NA)
  }

  if (is.null(dl$CVL) || length(dl$CVL) < 2 || is.na(dl$CVL[2])) {
    map$log_Galpha <- factor(NA)
  }

  # Optional fixed selectivity parameters
  if (!is.null(fix_res$map_Sm)) {
    map$Sm <- fix_res$map_Sm
  }

  if (length(map) == 0) map <- NULL

  obj <- TMB::MakeADFun(
    data = tmb_data[setdiff(names(tmb_data),
                            c("gear_names", "Linf_init", "Mk_init", "CVL_init", "Galpha_init"))],
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


#' Apply fixed selectivity settings to FLicc parameters
#'
#' Updates the packed selectivity parameter vector \code{Sm} and builds a TMB
#' mapping object for fixing selected selectivity parameters by gear. This
#' provides a user-facing interface for constraining selectivity parameters on
#' the natural length scale while retaining the internal relative-to-\code{Linf}
#' parameterization used by the TMB model.
#'
#' @param parameters A named parameter list, typically produced by
#'   \code{init_tmb_flicc()}, containing at least the packed selectivity vector
#'   \code{Sm}.
#' @param tmb_data A TMB data list produced by \code{data_tmb_flicc()}.
#' @param sel_fixed An optional named list specifying fixed selectivity values by
#'   gear. Names may be either gear names or character representations of gear
#'   indices. Each element should itself be a list of fixed parameter values on
#'   the natural length scale.
#'
#' @details
#' The function allows users to fix selectivity parameters by gear without
#' needing to work directly with the internal packed \code{Sm} vector.
#'
#' Supported user-facing selectivity inputs are:
#'
#' \strong{Logistic selectivity} (\code{sel_type == 1})
#' \itemize{
#'   \item \code{SL50}: length at 50 percent selectivity
#'   \item \code{SL95}: length at 95 percent selectivity
#'   \item \code{dSL}: difference \code{SL95 - SL50}
#' }
#'
#' \strong{Normal selectivity} (\code{sel_type == 3})
#' \itemize{
#'   \item \code{mode}: modal length of selectivity
#'   \item \code{sd}: standard deviation of the symmetric dome
#' }
#'
#' \strong{Double-sided normal selectivity} (\code{sel_type == 2})
#' \itemize{
#'   \item \code{mode}: modal length of selectivity
#'   \item \code{lsd}: left-side standard deviation
#'   \item \code{rsd}: right-side standard deviation
#' }
#'
#' Internally, all selectivity parameters are stored relative to \code{Linf} on
#' the log scale. The function converts natural-scale user inputs to this
#' internal parameterization and updates the TMB map so that selected values are
#' held fixed during optimization.
#'
#' This makes it possible, for example, to:
#' \itemize{
#'   \item fix a logistic selectivity curve using \code{SL50} and \code{SL95},
#'   \item fix a symmetric dome using \code{mode} and \code{sd},
#'   \item or approximate a one-sided shoulder by fixing a very large
#'     \code{rsd} for a \code{dsnormal} gear.
#' }
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{parameters}}{The updated parameter list with modified
#'       selectivity values in \code{Sm}.}
#'     \item{\code{map_Sm}}{A factor vector suitable for inclusion in the TMB
#'       \code{map} argument. Elements set to \code{NA} are fixed during
#'       optimization. Returns \code{NULL} if no selectivity parameters are
#'       fixed.}
#'   }
#'
#' @examples
#' \dontrun{
#' tmb_data <- data_tmb_flicc(dl)
#' parameters <- init_tmb_flicc(tmb_data)
#'
#' # Fix logistic selectivity using SL50 and SL95
#' sel_fixed <- list(
#'   "Gill net" = list(SL50 = 22, SL95 = 27)
#' )
#'
#' fix_res <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed)
#'
#' # Fix right tail of a double-sided normal to be very broad
#' sel_fixed2 <- list(
#'   "Gill net" = list(rsd = 1000)
#' )
#'
#' fix_res2 <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed2)
#'
#' # Fix a symmetric normal dome
#' sel_fixed3 <- list(
#'   "Marine set bagnet" = list(mode = 24, sd = 4)
#' )
#'
#' fix_res3 <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed3)
#' }
#'
#' @export

apply_sel_fixed_flicc <- function(parameters, tmb_data, sel_fixed = NULL) {
  if (is.null(sel_fixed) || length(sel_fixed) == 0) {
    return(list(parameters = parameters, map_Sm = NULL))
  }

  Sm_map <- factor(seq_along(parameters$Sm))

  gear_names <- tmb_data$gear_names
  sel_type   <- tmb_data$sel_type
  sm_start   <- tmb_data$sm_start
  sm_n       <- tmb_data$sm_n
  Linf       <- tmb_data$Linf_init

  for (nm in names(sel_fixed)) {
    spec <- sel_fixed[[nm]]
    if (is.null(spec)) next

    if (nm %in% gear_names) {
      g <- match(nm, gear_names)
    } else if (grepl("^[0-9]+$", nm)) {
      g <- as.integer(nm)
      if (g < 1 || g > length(gear_names)) {
        stop("Gear index out of range in sel_fixed: ", nm)
      }
    } else {
      stop("Unknown gear in sel_fixed: ", nm)
    }

    start <- sm_start[g] + 1L

    if (sel_type[g] == 1L) {
      # logistic: c(log_SL50_rel, log_dSL_rel)

      current_SL50 <- exp(parameters$Sm[start]) * Linf
      current_dSL  <- exp(parameters$Sm[start + 1L]) * Linf

      if (!is.null(spec$SL50)) {
        current_SL50 <- as.numeric(spec$SL50)
        if (current_SL50 <= 0) stop("SL50 must be > 0")
        parameters$Sm[start] <- log(current_SL50 / Linf)
        Sm_map[start] <- NA
      }

      if (!is.null(spec$SL95)) {
        dSL <- as.numeric(spec$SL95) - current_SL50
        if (dSL <= 0) stop("For logistic selectivity, SL95 must be > SL50")
        parameters$Sm[start + 1L] <- log(dSL / Linf)
        Sm_map[start + 1L] <- NA
      }

      if (!is.null(spec$dSL)) {
        dSL <- as.numeric(spec$dSL)
        if (dSL <= 0) stop("For logistic selectivity, dSL must be > 0")
        parameters$Sm[start + 1L] <- log(dSL / Linf)
        Sm_map[start + 1L] <- NA
      }

    } else if (sel_type[g] == 2L) {
      # dsnormal: c(log_mode_rel, log_lsd_rel, log_rsd_rel)

      if (!is.null(spec$mode)) {
        mode <- as.numeric(spec$mode)
        if (mode <= 0) stop("mode must be > 0")
        parameters$Sm[start] <- log(mode / Linf)
        Sm_map[start] <- NA
      }

      if (!is.null(spec$lsd)) {
        lsd <- as.numeric(spec$lsd)
        if (lsd <= 0) stop("lsd must be > 0")
        parameters$Sm[start + 1L] <- log(lsd / Linf)
        Sm_map[start + 1L] <- NA
      }

      if (!is.null(spec$rsd)) {
        rsd <- as.numeric(spec$rsd)
        if (rsd <= 0) stop("rsd must be > 0")
        parameters$Sm[start + 2L] <- log(rsd / Linf)
        Sm_map[start + 2L] <- NA
      }

    } else if (sel_type[g] == 3L) {
      # normal: c(log_mode_rel, log_sd_rel)

      if (!is.null(spec$mode)) {
        mode <- as.numeric(spec$mode)
        if (mode <= 0) stop("mode must be > 0")
        parameters$Sm[start] <- log(mode / Linf)
        Sm_map[start] <- NA
      }

      if (!is.null(spec$sd)) {
        sd <- as.numeric(spec$sd)
        if (sd <= 0) stop("sd must be > 0")
        parameters$Sm[start + 1L] <- log(sd / Linf)
        Sm_map[start + 1L] <- NA
      }

    } else {
      stop("Unsupported selectivity type for gear ", g)
    }

    if ((start + sm_n[g] - 1L) > length(parameters$Sm)) {
      stop("Selectivity indexing exceeded Sm length for gear ", g)
    }
  }

  list(parameters = parameters, map_Sm = Sm_map)
}
