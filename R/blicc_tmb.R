#' Build multi-year TMB data list for FLicc
#'
#' @param lfd named list of FLQuant objects, one per gear
#' @param stklen FLStockLen-like object with matching years
#' @param sel_fun selectivity function per gear
#' @param catch_by_gear numeric vector of relative catch by gear
#' @param settings list with CVL, GL, catch.sd, optional linf.sd, Mk.sd, CVL.sd
#' @param years optional vector of years to include
#'
#' @return list for TMB::MakeADFun
data_tmb_flicc <- function(lfd, stklen, sel_fun, catch_by_gear,
                             settings = list(CVL = 0.1, GL = 50, catch.sd = 0.05),
                             years = NULL) {

  stopifnot(is.list(lfd))
  if (is.null(names(lfd))) {
    names(lfd) <- paste0("gear", seq_along(lfd))
  }
  gear_names <- names(lfd)
  ngear <- length(lfd)

  lens <- as.numeric(dimnames(lfd[[1]])$len)
  yrs_all <- as.numeric(dimnames(lfd[[1]])$year)

  if (is.null(years)) {
    years <- yrs_all
  }
  years <- as.numeric(years)
  nyear <- length(years)
  nlen  <- length(lens)

  # observed length frequencies: obs[len, year, gear]
  obs <- array(NA_real_,
               dim = c(nlen, nyear, ngear),
               dimnames = list(len = lens, year = years, gear = gear_names))

  for (g in seq_len(ngear)) {
    x <- lfd[[g]]
    xyrs <- as.numeric(dimnames(x)$year)
    if (!all(years %in% xyrs)) {
      stop("Not all requested years found in lfd gear: ", gear_names[g])
    }

    for (iy in seq_along(years)) {
      y <- years[iy]
      obs[, iy, g] <- as.numeric(as.data.frame(yearMeans(x[, ac(y)]))$data)
    }
  }

  # stock years
  stk_years <- as.numeric(dimnames(stklen)$year)
  if (!all(years %in% stk_years)) {
    stop("Not all requested years found in stklen")
  }

  stkw <- window(stklen, start = min(years), end = max(years))

  # extract year-varying biology as len x year matrices
  wt_df  <- as.data.frame(catch.wt(stkw))
  mat_df <- as.data.frame(mat(stkw))
  m_df   <- as.data.frame(m(stkw))

  # assume ordering len-fast, then year
  wt_mat  <- matrix(wt_df$data,  nrow = nlen, ncol = nyear)
  mat_mat <- matrix(mat_df$data, nrow = nlen, ncol = nyear)
  m_mat   <- matrix(m_df$data,   nrow = nlen, ncol = nyear)

  # scale M around annual mean so Mk remains the average level
  Mscaler <- sweep(m_mat, 2, colMeans(m_mat, na.rm = TRUE), "/")

  # lengths
  LLB <- lens
  step <- c(diff(LLB), tail(diff(LLB), 1))
  Lmid <- LLB + step / 2

  # selectivity
  if (length(sel_fun) != ngear) {
    stop("length(sel_fun) must equal number of gears")
  }

  sel_type <- ifelse(sel_fun == "logistic", 1L,
                     ifelse(sel_fun == "dsnormal", 2L,
                            ifelse(sel_fun == "normal", 3L, NA_integer_)))
  if (any(is.na(sel_type))) {
    stop("Supported selectivity options are logistic, dsnormal, normal")
  }

  sm_n <- ifelse(sel_fun == "logistic", 2L,
                 ifelse(sel_fun == "dsnormal", 3L,
                        ifelse(sel_fun == "normal", 2L, NA_integer_)))
  sm_start <- cumsum(c(0L, head(sm_n, -1L)))  # 0-based for C++





  # catch weights by gear
  catch_wt <- do.call(
    cbind,
    lapply(catch_by_gear, function(x) as.data.frame(x)$data)
  )
  rownames(catch_wt) <- dimnames(catch_by_gear[[1]])$year

  catch_wt <- catch_wt / rowSums(catch_wt)


  if (ncol(catch_wt) != ngear) {
    stop("length(catch_by_gear) must equal number of gears")
  }
  if (nrow(catch_wt) != nyear) {
    stop("catch_by_gear must have one value per year")
  }

  catch_sd <- if (!is.null(settings$catch.sd)) as.numeric(settings$catch.sd) else 0.05

  # quadrature
  gl <- statmod::gauss.quad(settings$GL, kind = "laguerre")
  node <- as.numeric(gl$nodes)
  quad_wt <- as.numeric(gl$weights)

  # life history initial values
  lhpar <- stklen@lhpar
  Linf_init <- as.numeric(lhpar["linf"])
  Mk_init   <- as.numeric(lhpar["Mk"])

  CVL_init <- if (!is.null(settings$CVL)) as.numeric(settings$CVL) else 0.1
  Galpha_init <- 1 / CVL_init^2

  prior_mu <- numeric()
  prior_sd <- numeric()
  prior_code <- integer()

  if (!is.null(settings$linf.sd)) {
    prior_mu   <- c(prior_mu, log(Linf_init))
    prior_sd   <- c(prior_sd, as.numeric(settings$linf.sd))
    prior_code <- c(prior_code, 1L)
  }

  if (!is.null(settings$Mk.sd)) {
    prior_mu   <- c(prior_mu, log(Mk_init))
    prior_sd   <- c(prior_sd, as.numeric(settings$Mk.sd))
    prior_code <- c(prior_code, 3L)
  }

  if (!is.null(settings$CVL.sd)) {
    prior_mu   <- c(prior_mu, -2 * log(CVL_init))
    prior_sd   <- c(prior_sd, 2 * as.numeric(settings$CVL.sd))
    prior_code <- c(prior_code, 2L)
  }

  list(
    nlen       = nlen,
    nyear      = nyear,
    ngear      = ngear,
    Lmid       = as.numeric(Lmid),
    LLB        = as.numeric(LLB),
    obs        = unname(obs),
    Mscaler    = unname(Mscaler),
    wt         = unname(wt_mat),
    mat        = unname(mat_mat),
    catch_wt   = unname(catch_wt),
    catch_sd   = catch_sd,
    sel_type   = as.integer(sel_type),
    sm_start   = as.integer(sm_start),
    sm_n       = as.integer(sm_n),
    node       = node,
    quad_wt    = quad_wt,
    prior_mu   = as.numeric(prior_mu),
    prior_sd   = as.numeric(prior_sd),
    prior_code = as.integer(prior_code),
    Linf_init   = Linf_init,
    Mk_init     = Mk_init,
    CVL_init    = CVL_init,
    Galpha_init = Galpha_init,
    gear_names  = gear_names,
    year_names  = years
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
#' data(alfonsino)
#' tmb_data <- data_tmb_flicc(lfd,stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear = c(0.7,0.3))
#' par0 <- init_tmb_flicc(tmb_data)
#' str(par0)
#' }
#'
#' @export
#' Initial values for multi-year FLicc TMB
init_tmb_flicc <- function(tmb_data) {

  nlen     <- tmb_data$nlen
  nyear    <- tmb_data$nyear
  ngear    <- tmb_data$ngear
  sel_type <- tmb_data$sel_type
  Lmid     <- tmb_data$Lmid
  obs      <- tmb_data$obs
  Linf     <- tmb_data$Linf_init

  Sm <- numeric(sum(tmb_data$sm_n))
  pos <- 1L

  for (g in seq_len(ngear)) {
    obs_g <- apply(obs[, , g, drop = FALSE], 1, sum, na.rm = TRUE)
    mode_g <- Lmid[which.max(obs_g)]
    mode_rel <- max(mode_g / Linf, 1e-4)

    if (sel_type[g] == 1L) {
      SL50_rel <- mode_rel
      dSL_rel  <- 0.10
      Sm[pos:(pos + 1L)] <- c(log(SL50_rel), log(dSL_rel))
      pos <- pos + 2L

    } else if (sel_type[g] == 2L) {
      lsd_rel <- 0.10
      rsd_rel <- 0.10
      Sm[pos:(pos + 2L)] <- c(log(mode_rel), log(lsd_rel), log(rsd_rel))
      pos <- pos + 3L

    } else if (sel_type[g] == 3L) {
      sd_rel <- 0.10
      Sm[pos:(pos + 1L)] <- c(log(mode_rel), log(sd_rel))
      pos <- pos + 2L
    }
  }

  list(
    log_Linf   = log(tmb_data$Linf_init),
    log_Galpha = log(tmb_data$Galpha_init),
    log_Mk     = log(tmb_data$Mk_init),
    log_Fk     = matrix(log(0.5), nrow = nyear, ncol = ngear,
                        dimnames = list(tmb_data$year_names, tmb_data$gear_names)),
    Sm         = Sm,
    log_phi    = log(10),
    log_sigmaF = tmb_data$prior_sigmaF_mean
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
#' Fit multi-year FLicc TMB model
fiticc <- function(lfd, stklen,
                     sel_fun = c("logistic", "normal", "dsnormal"),
                     catch_by_gear,
                     settings = list(
                       CVL = 0.1,
                       GL = 50,
                       catch.sd = 0.05,
                       prior_sigmaF = c(log(0.7), 0.1, 1)
                     ),
                     years = NULL,
                     iter = 1,
                     compile = FALSE,
                     silent = TRUE,
                     dll = "FLicc",
                     sel_fixed = NULL,FLRreport=TRUE ) {

  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Package 'TMB' is required")
  }


  default_settings <- list(
    CVL = 0.1,
    GL = 50,
    catch.sd = 0.05,
    prior_sigmaF = c(log(0.3), 0.3, 1)   # optional
  )

  # ---- merge user settings with defaults ----
  if (is.null(settings)) settings <- list()

  if(any(settings$prior_sigmaF == FALSE)) settings$prior_sigmaF <- c(NA_real_, NA_real_, 0)

  # fill missing elements
  for (nm in names(default_settings)) {
    if (is.null(settings[[nm]])) {
      settings[[nm]] <- default_settings[[nm]]
    }
  }

  unknown <- setdiff(names(settings), names(default_settings))
  if (length(unknown) > 0) {
    stop("Unknown settings: ", paste(unknown, collapse = ", "))
  }

  # subset iter in R; do not loop over iter in TMB
  lfd_i <- lapply(lfd, function(x) x[, , iter])
  stk_i <- stklen[, , , iter]

  if (is.null(years)) {
    years <- as.numeric(dimnames(lfd_i[[1]])$year)
  }
  years <- as.numeric(years)

  if(!inherits(catch_by_gear, "FLQuants")){

 catch_by_gear  <-  FLCore::FLQuants(Map(function(x,y){
    flq = FLCore::quantSums(x)
    flq[] <- y
    flq
  },x=lfd_i,y=catch_by_gear))

 }

  tmb_data <- data_tmb_flicc(
    lfd = lfd_i,
    stklen = stk_i,
    sel_fun = sel_fun,
    catch_by_gear = catch_by_gear,
    settings = settings,
    years = years
  )

  # Set priors for random walk on F (if nyear > 0)
  tmb_data$prior_sigmaF_mean <- as.numeric(settings$prior_sigmaF[1])
  tmb_data$prior_sigmaF_sd   <- as.numeric(settings$prior_sigmaF[2])
  tmb_data$prior_sigmaF_use  <- as.integer(settings$prior_sigmaF[3])


  parameters <- init_tmb_flicc(tmb_data)

  fix_res <- apply_sel_fixed_flicc(
    parameters = parameters,
    tmb_data   = tmb_data,
    sel_fixed  = sel_fixed
  )
  parameters <- fix_res$parameters

  map <- list()

  if (is.null(settings$linf.sd)) {
    map$log_Linf <- factor(NA)
  }
  if (is.null(settings$Mk.sd)) {
    map$log_Mk <- factor(NA)
  }
  if (is.null(settings$CVL.sd)) {
    map$log_Galpha <- factor(NA)
  }
  if (!is.null(fix_res$map_Sm)) {
    map$Sm <- fix_res$map_Sm
  }
  if (length(map) == 0) {
    map <- NULL
  }

  if (compile) {
    TMB::compile("src/FLicc.cpp")
    dyn.load(dynlib("src/FLicc"))
  }

  obj <- TMB::MakeADFun(
    data = tmb_data[setdiff(names(tmb_data),
                            c("gear_names", "year_names",
                              "Linf_init", "Mk_init", "CVL_init", "Galpha_init"))],
    parameters = parameters,
    map = map,
    DLL = dll,
    silent = silent
  )

  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(iter.max = 1000, eval.max = 1000)
  )

  obj$env$last.par.best <- opt$par

  rep <- TMB::sdreport(obj)
  adrep <- summary(rep, "fixed")
  report <- obj$report()
  par_tab <- data.frame(
    par = rownames(adrep),
    mpd = adrep[, 1],
    se  = adrep[, 2],
    row.names = NULL,
    stringsAsFactors = FALSE
  )


  fit <- structure(
    list(
      opt = opt,
      obj = obj,
      rep = rep,
      par = par_tab,
      report = report,
      stklen= stklen,
      tmb_data = tmb_data
    ),
    class = "flicc_tmb_fit"
  )
  if(FLRreport){fit$report <- as_FLQuants(fit,stklen)} else {fit$report$FLReport=FALSE}
  return(fit)

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
