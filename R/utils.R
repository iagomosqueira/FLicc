#' Convert FLicc TMB report output to FLQuant/FLQuants objects
#'
#' Converts selected elements from a fitted FLicc TMB report and its
#' associated `tmb_data` input list into a structured list of FLR objects.
#' The returned list contains `FLQuant` and `FLQuants` objects for
#' abundance, observed and predicted length compositions, selectivity,
#' fishing mortality by gear, and spawning potential ratio.
#'
#' This helper is intended to provide a more FLR-friendly representation of
#' model output than the raw arrays returned by `obj$report()`.
#'
#' @param report A list of reported model outputs, typically returned by
#'   `obj$report()` from a fitted TMB model. Expected elements include
#'   `N_y`, `Sel`, `plen`, `plen_all_y`, `Fk`, and `spr_y`.
#' @param tmb_data A TMB data list used to fit the model. Expected elements
#'   include `year_names`, `gear_names`, `LLB`, and `obs`.
#'
#' @details
#' The returned object is a named list with the following components:
#' \describe{
#'   \item{`N`}{An `FLQuant` containing estimated numbers-at-length by year.}
#'   \item{`sel`}{An `FLQuants` object with one `FLQuant` per gear,
#'   containing selectivity-at-length.}
#'   \item{`obslen`}{An `FLQuants` object with one `FLQuant` per gear,
#'   containing observed length frequencies by year.}
#'   \item{`predlen_gear`}{An `FLQuants` object with one `FLQuant` per gear,
#'   containing predicted length compositions by year and gear.}
#'   \item{`predlen`}{An `FLQuants` object containing predicted length
#'   compositions by year for each gear-specific component returned in
#'   `report$plen_all_y`.}
#'   \item{`Fk`}{An `FLQuants` object with one `FLQuant` per gear,
#'   containing fishing mortality by year.}
#'   \item{`spr`}{An `FLQuant` containing spawning potential ratio by year.}
#' }
#'
#' Length classes are mapped to the `quant`/`len` dimension, years to the
#' `year` dimension, and gears are represented as separate elements in
#' `FLQuants` objects.
#'
#' @return A named list containing `FLQuant` and `FLQuants` objects.
#'
#' @seealso [FLCore::FLQuant()], [FLCore::FLQuants()]
#'
#' @examples
#' \dontrun{
#' rep <- fit$obj$report()
#' flqs <- as_FLQuants(rep, fit$tmb_data)
#'
#' flqs$N
#' flqs$obslen
#' flqs$predlen_gear
#' flqs$predlen
#' flqs$spr
#' names(flqs$obslen)
#' names(flqs$predlen_gear)
#' }
#'
#' @export

as_FLQuants <- function(fit,stklen) {

  pop_model <- if (!is.null(fit$pop_model)) fit$pop_model else
    if (!is.null(fit$settings$pop_model)) fit$settings$pop_model else "gamma"

  obs_model <- if (!is.null(fit$obs_model)) fit$obs_model else
    if (!is.null(fit$settings$obs_model)) fit$settings$obs_model else "nb"


  tmb_data <- fit$tmb_data
  report <- fit$report
  yrs   <- an(tmb_data$year_names)
  gears <- tmb_data$gear_names
  lens  <- an(tmb_data$LLB)
  k <- an(stklen@lhpar["k"])


  out <- list()

  out$N <- FLCore::FLQuant(
    report$N_y,
    dimnames = list(
      lens = lens,
      year = yrs,
      unit = "unique",
      season = "all",
      area = "unique",
      iter = "1"
    )
  )
  units(out$N) <- "1000"


  out$w_at_len <- stklen@catch.wt


  # FLQuants
  out$sel_gear<- FLCore::FLQuants(lapply(1:length(gears),function(x){
    dat <- report$Sel[,x]
    FLCore::FLQuant(
      dat,
      dimnames = list(
        len = lens,
        year = "all",
        unit = "unique",
        season = "all",
        area = "unique",
        iter = "1"
      ))
  }))
  names(out$sel_gear) <- gears

  # FLQuants
  out$sel <-FLCore::FLQuant(
    report$sel_joint,
      dimnames = list(
        len = lens,
        year = yrs,
        unit = "unique",
        season = "all",
        area = "unique",
        iter = "1"
    ))
  out$sel <-  out$sel%/%apply( out$sel,2,max)



  # FLQuants
  out$catch_by_gear <- FLCore::FLQuants(lapply(1:length(gears),function(x){
    dat <- tmb_data$catch_wt[,x]
    FLCore::FLQuant(
      dat,
      dimnames = list(
        len = "all",
        year = yrs,
        unit = "unique",
        season = "all",
        area = "unique",
        iter = "1"
      ))
  }))
  names(out$catch_by_gear) <-  gears


  # FLQuants
  out$predlen_gear <- FLCore::FLQuants(lapply(1:length(gears),function(x){
    dat <- report$plen[,,x]
    FLCore::FLQuant(
      dat,
      dimnames = list(
        len = lens,
        year = yrs,
        unit = "unique",
        season = "all",
        area = "unique",
        iter = "1"
      ),unit="cm")
  }))
  names(out$predlen_gear) <- gears



    dat <- report$plen_all_y
    out$predlen <- FLCore::FLQuant(
      dat,
      dimnames = list(
        len = lens,
        year = yrs,
        unit = "unique",
        season = "all",
        area = "unique",
        iter = "1"
      ),unit="cm")
    #units(flq) <- "proportion"

  names(out$predlen) <- gears

  out$obslen <- FLCore::FLQuants(lapply(1:length(gears),function(x){
    dat <- tmb_data$obs[,,x]
    FLCore::FLQuant(
      dat,
      dimnames = list(
        len = lens,
        year = yrs,
        unit = "unique",
        season = "all",
        area = "unique",
        iter = "1"
      ),unit="cm")
  }))
  names(out$obslen) <- gears

  out$Mk <- FLCore::FLQuant(
    report$M,
    dimnames = list(
      len = lens,
      year = yrs,
      unit = "unique",
      season = "all",
      area = "unique",
      iter = "1"
    ))
  units(out$Mk) <- "Mk"

  out$Fk <- FLQuants(lapply(1:length(gears),function(x){
    dat <- report$Fk[,x]
    flq = FLCore::FLQuant(
      dat,
      dimnames = list(
        len = "all",
        year = yrs,
        unit = "unique",
        season = "all",
        area = "unique",
        iter = "1"
      ))
    units(flq) <- "Fk"
    flq
  }))
  names(out$Fk) <- gears

  out$F <-FLCore::FLQuant(
    report$Fk_l*k,
    dimnames = list(
      len = lens,
      year = yrs,
      unit = "unique",
      season = "all",
      area = "unique",
      iter = "1"
    ))

  out$Fap <- apply(out$F,2:6,max,na.rm=TRUE)

  out$M <- FLCore::FLQuant(
    report$Mk*k,
    dimnames = list(
      len = lens,
      year = yrs,
      unit = "unique",
      season = "all",
      area = "unique",
      iter = "1"
    ))

  out$mat <- FLCore::mat(stklen)

  out$spr <- FLCore::FLQuant(
    matrix(report$spr_y, nrow = 1),
    dimnames = list(
      quant = "all",
      year = yrs,
      unit = "unique",
      season = "all",
      area = "unique",
      iter = "1"
    )
  )
  units(out$spr) <- "spr"

  out$lhpar <- stklen@lhpar
  out$lhpar["linf"] <- report$Linf
  out$lhpar["Mk"] <- report$Mk
  out$lhpar["M"] <- out$lhpar["k"]*report$Mk

  out$selpars <- selpars_flicc(fit)
  out$pars <- FLPar(Linf=report$Linf, Galpha= report$Galpha,
                    Gbeta = report$Gbeta,phi=report$phi)

  par_list <- list(Linf = report$Linf)

  if (pop_model == "gamma") {
    par_list$Galpha <- report$Galpha
    par_list$Gbeta  <- report$Gbeta
  }

  if (pop_model == "gtg") {
    if (!is.null(tmb_data$ngtg))  par_list$ngtg  <- tmb_data$ngtg
    if (!is.null(fit$settings$maxsd)) par_list$maxsd <- fit$settings$maxsd
    if (!is.null(fit$settings$Mpow))  par_list$Mpow  <- fit$settings$Mpow
  }

  if (obs_model == "nb") {
    par_list$phi <- report$phi
  }

  out$pars <- do.call(FLCore::FLPar, par_list)


  out$Lmid = tmb_data$Lmid
  if (pop_model == "gamma") {
    out$node <- tmb_data$node
    out$quad_wt <- tmb_data$quad_wt
  }
  out$logLik <- -fit$opt$objective
  out$pop_model <- pop_model
  out$obs_model <- obs_model

  out$FLReport <- TRUE
  return(out)

}




#' Lorenzen natural mortality at length
#'
#' Computes natural mortality at length using a Lorenzen-style allometric
#' relationship referenced to a chosen length.
#'
#' @param len Numeric vector of lengths, usually bin midpoints.
#' @param Mref Natural mortality at the reference length.
#' @param reflen Reference length corresponding to \code{Mref}.
#' @param a Length-weight coefficient.
#' @param b Length-weight exponent.
#' @param exponent Lorenzen exponent on weight. Default is \code{-0.288}.
#'
#' @return Numeric vector of natural mortality values at length.
#' @export
m_lorenzen <- function(len, Mref, reflen, a, b, exponent = -0.288) {
  stopifnot(all(len > 0), Mref > 0, reflen > 0, a > 0, b > 0)

  w <- a * len^b
  wref <- a * reflen^b

  Mref * (w / wref)^exponent
}

#' Gislason natural mortality at length
#'
#' Computes natural mortality at length using the Gislason et al. length-based
#' empirical formula.
#'
#' @param len Numeric vector of lengths, usually bin midpoints.
#' @param linf Asymptotic length.
#' @param k Von Bertalanffy growth coefficient.
#'
#' @return Numeric vector of natural mortality values at length.
#' @export
m_gislason <- function(len, linf, k) {
  stopifnot(all(len > 0), linf > 0, k > 0)

  exp(0.55 - 1.61 * log(len) + 1.44 * log(linf) + log(k))
}

#' Selectivity parameter names by selectivity function
#'
#' Returns the expected user-facing selectivity parameter names for each
#' supported FLicc selectivity function. These names are intended for
#' reporting, helper constructors, and user input on the natural length scale.
#'
#' For logistic selectivity, the user-facing parameterization is
#' \code{SL50} and \code{SL95}, where \code{SL50} is the length at 50 percent
#' selectivity and \code{SL95} is the length at 95 percent selectivity.
#' Internally, the TMB model uses \code{SL50} and
#' \code{dSL = SL95 - SL50}, but this internal representation is hidden from
#' the user.
#'
#' For double-sided normal selectivity, parameters are
#' \code{mode}, \code{lsd}, and \code{rsd}, corresponding to the modal length,
#' left-side standard deviation, and right-side standard deviation.
#'
#' For symmetric normal selectivity, parameters are
#' \code{mode} and \code{sd}, corresponding to the modal length and standard
#' deviation.
#'
#' @return A named list with one character vector per supported selectivity
#'   function:
#'   \describe{
#'     \item{\code{logistic}}{\code{c("SL50", "SL95")}}
#'     \item{\code{dsnormal}}{\code{c("mode", "lsd", "rsd")}}
#'     \item{\code{normal}}{\code{c("mode", "sd")}}
#'   }
#'
#' @details
#' These parameter names are intended as the standard user-facing convention
#' throughout FLicc. They are especially useful when:
#' \itemize{
#'   \item constructing \code{FLPar} or \code{FLPars} objects for selectivity,
#'   \item extracting fitted selectivity parameters for reporting,
#'   \item validating user input in helper functions,
#'   \item mapping between natural-scale selectivity parameters and the internal
#'     packed TMB parameter vector \code{Sm}.
#' }
#'
#' Internally, the TMB objective function parameterizes selectivity relative to
#' \code{Linf} on the log scale. For example, logistic selectivity is fitted
#' internally as \code{log(SL50 / Linf)} and \code{log((SL95 - SL50) / Linf)}.
#' However, reporting and helper functions should generally use the more
#' interpretable natural-scale names returned here.
#'
#' @examples
#' sel_parnames <- list(
#'   logistic = c("SL50", "SL95"),
#'   dsnormal = c("mode", "lsd", "rsd"),
#'   normal   = c("mode", "sd")
#' )
#'
#' sel_parnames$logistic
#' sel_parnames$dsnormal
#' sel_parnames$normal
#'
#' # Example FLPar objects on the natural length scale
#' trawl_par <- FLCore::FLPar(setNames(c(24, 30), c("SL50", "SL95")))
#' gillnet_par <- FLCore::FLPar(setNames(c(28, 4, 8), c("mode", "lsd", "rsd")))
#'
#' @seealso \code{\link[FLCore]{FLPar}}, \code{\link[FLCore]{FLPars}}
#'
#'@export

selpars_flicc <- function(fit) {

  parlist <- fit$obj$env$parList()
  Sm <- parlist$Sm
  Linf <- exp(parlist$log_Linf)

  sel_type   <- fit$tmb_data$sel_type
  sm_start   <- fit$tmb_data$sm_start
  sm_n       <- fit$tmb_data$sm_n
  gear_names <- fit$tmb_data$gear_names

  out <- vector("list", length(gear_names))
  names(out) <- gear_names

  for (g in seq_along(gear_names)) {

    start <- sm_start[g] + 1L
    idx <- start:(start + sm_n[g] - 1L)
    pars <- Sm[idx]

    if (sel_type[g] == 1L) {
      # internal: SL50, dSL
      SL50 <- exp(pars[1]) * Linf
      dSL  <- exp(pars[2]) * Linf
      SL95 <- SL50 + dSL

      out[[g]] <- FLCore::FLPar(
        setNames(c(SL50 , SL95), c("SL50", "SL95"))
      )

    } else if (sel_type[g] == 2L) {
      mode <- exp(pars[1]) * Linf
      lsd  <- exp(pars[2]) * Linf
      rsd  <- exp(pars[3]) * Linf

      out[[g]] <- FLCore::FLPar(
        setNames(c(mode, lsd, rsd), c("mode", "lsd", "rsd"))
      )

    } else if (sel_type[g] == 3L) {
      mode <- exp(pars[1]) * Linf
      sd   <- exp(pars[2]) * Linf

      out[[g]] <- FLCore::FLPar(
        setNames(c(mode, sd), c("mode", "sd"))
      )

    } else {
      stop("Unknown sel_type for gear ", gear_names[g])
    }
    FLCore::units(out[[g]]) <- "cm"
  }

  FLPars(out)


}

#' Extract log-likelihood from a fitted FLicc model
#'
#' Returns the fitted log-likelihood from a \code{"flicc_tmb_fit"} object as a
#' \code{logLik}-class object. The value is computed as the negative of the
#' TMB objective function at the optimum.
#'
#' The returned object includes \code{df}, taken as the number of estimated
#' parameters in \code{object$opt$par}, and \code{nobs}, taken as the number of
#' rows in \code{as.data.frame(object$report$obslen)}.
#'
#' @param object A fitted \code{"flicc_tmb_fit"} object.
#' @param ... Additional arguments, ignored.
#'
#' @return An object of class \code{logLik}.
#'
#' @details
#' This helper assumes that \code{object$opt$objective} is the negative
#' log-likelihood returned by TMB. If penalties or priors were included in the
#' objective function, the returned value should be interpreted as the
#' corresponding penalized log-likelihood or objective-based criterion rather
#' than a pure data likelihood.
#'
#' @examples
#' \dontrun{
#' LL <- LLflicc(fit)
#' LL
#' AIC(LL)
#' BIC(LL)
#' }
#'
#' @export

LLflicc <- function(object, ...) {
  structure(
    -object$opt$objective,
    df = length(object$opt$par),
    nobs = nrow(as.data.frame(object$report$obslen)),
    class = "logLik"
  )
}



#' Rescale length frequencies to effective sample size by gear
#'
#' Rescales an \code{FLQuant} of length frequencies so that, within each
#' year-gear-iter combination, the total over length equals a user-supplied
#' effective sample size (ESS). This is useful for composition-based likelihoods
#' such as multinomial or Dirichlet-multinomial, where raw sample sizes are not
#' intended to determine the weight of the fit.
#'
#' @param lfd An \code{FLQuant} of length frequencies with at least dimensions
#'   \code{len}, \code{year}, \code{unit}, and \code{iter}.
#' @param ess.g Numeric scalar or vector of effective sample sizes by gear. If a
#'   scalar is supplied, it is recycled across gears. If a vector is supplied,
#'   its length must match the number of gears in the \code{unit} dimension.
#'
#' @return An \code{FLQuant} with the same dimensions as \code{lfd}, but scaled
#'   so that totals over \code{len} equal the requested ESS for each gear.
#'
#' @examples
#' \dontrun{
#' lfd_ess <- lfdess(lfd, ess.g = c(Trawl = 150, Gillnet = 100))
#' }
#'
#' @export
lfdess <- function(lfd, ess.g=100) {
  if (!inherits(lfd, "FLQuants")) {
    stop("lfd must be an FLQuant")
  }

  dn <- dimnames(lfd[[1]])
  gears <- dn$unit
  years <- dn$year
  seasons <- dn$season
  areas <- dn$area
  iters <- dn$iter

  ng <- length(gears)

  if (length(ess.g) == 1 && length(lfd)>1) {
    ess.g <- rep(as.numeric(ess.g), ng)
    names(ess.g) <- gears
  } else {
    ess.g <- as.numeric(ess.g)

  }

  out <- FLQuants(Map(function(x,y){
    res <- x%/%apply(x,2:6,sum,na.rm=T)
    res*y
  },x=lfd,y=ess.g))


  out
}


#' Convert wide length-frequency tables to FLQuants
#'
#' Converts a wide \code{data.frame} of length-frequency observations into
#' an \code{FLQuants} object suitable as FLicc length-frequency input.
#'
#' The input table is expected to contain:
#' \itemize{
#'   \item column 1: length class values,
#'   \item column 2: gear names,
#'   \item remaining columns: annual length-frequency data by year.
#' }
#'
#' One \code{FLQuant} is created per gear, and the resulting collection is
#' returned as an \code{FLQuants} object.
#'
#' @param lfd A \code{data.frame} with length in the first column, gear in the
#'   second column, and yearly observations in the remaining columns.
#' @param midL Logical. If \code{TRUE}, the first column is interpreted as
#'   length midpoints and converted to lower bin bounds by subtracting half
#'   the bin width. If \code{FALSE} (default), the first column is assumed
#'   to already contain lower bin bounds.
#' @param unit Character string giving the length unit of the first column.
#'   If \code{"mm"}, lengths are converted to cm by dividing by 10.
#'   Default is \code{"mm"}.
#'
#' @details
#' The returned \code{FLQuants} object contains one \code{FLQuant} per gear.
#' Length classes are stored in the \code{len} dimension and years in the
#' \code{year} dimension. Output units are set to \code{"cm"}.
#'
#' If \code{midL = TRUE}, the bin width is calculated from the first two
#' length values in the input and half a bin width is subtracted to obtain
#' lower bin bounds.
#'
#' @return An \code{FLQuants} object with one element per gear.
#'
#' @examples
#' \dontrun{
#' lfd_df <- data.frame(
#'   Length = c(25, 35, 45, 25, 35, 45),
#'   gear   = c("Trawl", "Trawl", "Trawl", "Gillnet", "Gillnet", "Gillnet"),
#'   `2020` = c(10, 20, 15, 5, 8, 6),
#'   `2021` = c(12, 18, 17, 4, 9, 7)
#' )
#'
#' lfd_flq <- FLQuantLen(lfd_df, midL = TRUE, unit = "mm")
#' lfd_flq
#' }
#'
#' @export

FLQuantLen <- function(lfd,midL=FALSE,unit="mm") {

  gear <- unique(lfd[,2])
  half <- 0
  div <- 1
  if(unit=="mm") lfd[,1] <- lfd[,1]/10
  if(midL) half <- (lfd[2,1]-lfd[1,1])/2
  names(lfd)[2] <- "gear"

  flqs <- FLQuants(lapply(gear,function(x){
  FLQuant(as.matrix(lfd[lfd$gear%in%x,-c(1:2)]),dimnames=list(len=ac((lfd[lfd$gear%in%x,1]-half)/div), year=colnames(lfd)[-c(1:2)]),unit="cm")

  }))
  names(flqs) <- gear
  return(flqs)
}



#' Convert long-format LFD to FLQuantLen input format
#'
#' Converts a long-format length-frequency dataset into the wide format
#' required by \code{FLQuantLen()}, with one row per length class and one
#' column per year for each gear.
#'
#' The input data must contain the columns:
#' \itemize{
#'   \item \code{len}: length class values
#'   \item \code{year}: year
#'   \item \code{data}: observed counts or frequencies
#'   \item \code{qname}: gear identifier
#' }
#'
#' The output is a wide \code{data.frame} with columns:
#' \itemize{
#'   \item \code{len}: length class
#'   \item \code{gear}: gear name
#'   \item one column per year containing observations
#' }
#'
#' @param df A \code{data.frame} in long format with columns
#'   \code{len}, \code{year}, \code{data}, and \code{qname}.
#'
#' @return A wide \code{data.frame} suitable as input to
#'   \code{FLQuantLen()}.
#'
#' @details
#' The function splits the data by gear, reshapes each subset from long to
#' wide format using \code{stats::reshape()}, and then combines the results.
#' Missing combinations of length and year are filled with zeros.
#'
#' The resulting table can be passed directly to \code{FLQuantLen()} to
#' construct an \code{FLQuants} object for FLicc.
#'
#' @examples
#' \dontrun{
#' df_long <- as.data.frame(lfd_alfonsino)
#'
#' df_wide <- lfd_long_to_wide(df_long)
#'
#' lfd_flq <- FLQuantLen(df_wide, midL = FALSE,unit="cm")
#' }
#'
#' @seealso \code{\link{FLQuantLen}}
#'
#' @export

lfd_long_to_wide <- function(df) {
  df <- df[, c("len", "year", "data", "qname")]

  gears <- split(df, df$qname)

  out <- do.call(
    rbind,
    lapply(names(gears), function(g) {
      d <- gears[[g]]

      # keep only len, year, data for reshape
      d <- d[, c("len", "year", "data")]

      w <- reshape(
        d,
        idvar = "len",
        timevar = "year",
        direction = "wide"
      )

      names(w) <- sub("^data\\.", "", names(w))
      w <- w[order(w$len), , drop = FALSE]
      w$gear <- g

      w <- w[, c("len", "gear", setdiff(names(w), c("len", "gear"))), drop = FALSE]
      w
    })
  )

  rownames(out) <- NULL
  out
}
