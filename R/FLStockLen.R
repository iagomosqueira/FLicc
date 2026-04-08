

#' Construct an FLStockLen object from length bins and life-history parameters
#'
#' Builds an \code{FLStockLen} object from the length bins in a length-frequency
#' \code{FLQuants} together with scalar life-history inputs in an \code{FLPar}
#' object. The constructor derives midpoint lengths, weight-at-length,
#' maturity-at-length, and natural mortality-at-length under a chosen mortality
#' model.
#'
#' @param lfds An \code{FLQuants} containing length-frequency data by gear. The
#'   \code{len} dimension is used to define the length bins.
#' @param lhpar An \code{FLPar} object containing life-history parameters.
#'   Expected entries include \code{linf}, \code{L50}, \code{a}, and \code{b},
#'   and optionally \code{k}, \code{M}, \code{MK}, \code{L95}, and \code{s}.
#'   If \code{s} is not supplied, a default value of \code{0.7} is used.
#' @param m_model Character string specifying the natural mortality-at-length
#'   model. Supported options are \code{"constant"}, \code{"inverse"},
#'   \code{"Lorenzen"}, and \code{"Gislason"}.
#' @param reflen Optional reference length used by mortality models that require
#'   it, such as \code{"inverse"} and \code{"Lorenzen"}.
#'
#' @details
#' The \code{len} dimension of \code{lfd} is interpreted as lower length-bin
#' boundaries. Biological schedules are evaluated at bin midpoints derived from
#' these lower bounds.
#'
#' Weight at length is calculated as:
#' \deqn{W_L = aL^b}
#'
#' Maturity at length is represented by a logistic ogive parameterized by
#' \code{L50} and \code{L95}. If \code{L95} is not supplied, a default value is
#' constructed internally. If \code{s} is not supplied in \code{lhpar}, it
#' defaults to \code{0.7}.
#'
#' Natural mortality at length depends on \code{m_model}:
#' \itemize{
#'   \item \code{"constant"}: constant natural mortality at all lengths
#'   \item \code{"inverse"}: inverse-length mortality referenced to
#'     \code{reflen}
#'   \item \code{"Lorenzen"}: Lorenzen-style allometric mortality
#'   \item \code{"Gislason"}: Gislason length-based empirical mortality
#' }
#'
#' If \code{M} is not supplied directly but both \code{MK} and \code{k} are
#' available, then \code{M} is derived as \code{M = MK * k}. Likewise,
#' \code{MK} may be derived from \code{M / k} when needed.
#'
#' @return An \code{FLStockLen} object containing at least:
#'   \describe{
#'     \item{\code{stock.wt}}{Weight-at-length schedule}
#'     \item{\code{mat}}{Maturity-at-length schedule}
#'     \item{\code{m}}{Natural mortality-at-length schedule}
#'   }
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#' lhpar <- FLPar(
#'   linf = 55.7,
#'   k    = 0.08,
#'   M    = 0.162,
#'   L50  = 31.1859,
#'   a    = 0.004721956,
#'   b    = 3.146168
#' )
#'
#' stklen <- stocklen(
#'   lfd_alfonsino,
#'   lhpar = lhpar,
#'   m_model = "constant"
#' )
#'  stklen@lhpar
#' }
#'
#' @export

stocklen <- function(lfds, lhpar, m_model = c("constant", "inverse", "Lorenzen", "Gislason"),reflen=NULL) {
  m_model <- match.arg(m_model)

  pars <- dimnames(lhpar)$params


  if (!inherits(lfds, "FLQuants")) {
    stop("lfds must be class FLQuants")
  }
  if (!inherits(lhpar, "FLPar")) {
    stop("lhpar must be an FLPar")
  }


  getpar <- function(x, default = NA_real_) {
    if (x %in% pars) as.numeric(lhpar[x]) else default
  }

  linf   <- getpar("linf")
  k      <- getpar("k")
  M      <- getpar("M")
  Mk     <- getpar("Mk")
  L50    <- getpar("L50")
  L95    <- getpar("L95")
  a      <- getpar("a")
  b      <- getpar("b")
  s      <- getpar("s")
  if(is.null(reflen)) reflen <- getpar("reflen", L50)

  if (is.na(M) && !is.na(Mk) && !is.na(k)) M <- Mk * k
  if (is.na(Mk) && !is.na(M) && !is.na(k)) Mk <- M / k

  if(is.na(s)) s <- 0.7
  if (is.na(linf)) stop("lhpar must contain 'linf'")
  if (is.na(L50)) stop("lhpar must contain 'L50'")
  if (is.na(a) || is.na(b)) stop("lhpar must contain 'a' and 'b'")

  lens <- an(dimnames(lfds[[1]])$len)
  year <- an(dimnames(lfds[[1]])$year)
  its <- an(dimnames(lfds[[1]])$iter)
  if (anyNA(lens)) stop("lfds must have numeric 'len' dimnames")

  step <- c(diff(lens), tail(diff(lens), 1))
  midL <- lens + step / 2



  if (is.na(L95)) {
    L95 <- L50 + 0.2 * linf
  }
  if (L95 <= L50) stop("L95 must be > L50")

  # weight at length
  wt <- a * midL^b

  # maturity at length using L50/L95 parameterization
  dL <- L95 - L50
  mat <- 1 / (1 + exp(-log(19) * (midL - L50) / dL))

  # natural mortality at length
  if (is.na(M) && m_model %in% c("constant", "inverse")) {
    stop("Need M (or Mk and k) for selected m_model")
  }



  if (m_model %in% c("inverse", "Lorenzen") && is.null(reflen)) {
    stop("reflen must be supplied for m_model = '", m_model, "'")
  }
  if (m_model == "Gislason" && is.na(k)) {
    stop("k must be supplied in lhpar for m_model = 'Gislason'")
  }

  m <- switch(
    m_model,
    constant = rep(M, length(midL)),
    inverse  = M * reflen / midL,
    Lorenzen = m_lorenzen(midL, Mref = M, reflen = reflen, a = a, b = b),
    Gislason = m_gislason(midL, linf = linf, k = k)
  )

  # make FLQuants
  dn <- list(
    len = ac(lens),
    year = ac(year),
    unit = "unique",
    season = "all",
    area = "unique",
    iter = ac(its)
  )

  stock.n <- FLCore::FLQuant(rep(NA_real_, length(lens)), dimnames = dn)
  stock.wt <- FLCore::FLQuant(wt, dimnames = dn)
  mat.q <- FLCore::FLQuant(mat, dimnames = dn)
  m.q <- FLCore::FLQuant(m, dimnames = dn)
  catch.n <- FLCore::FLQuant(rep(NA_real_, length(lens)), dimnames = dn)
  catch.wt <- stock.wt
  flq0 <-   FLCore::FLQuant(rep(0, length(lens)), dimnames = dn)

  out <- FLCore::FLStockLen(
    stock.n = stock.n,
    catch.n = catch.n,
    catch.wt = catch.wt,
    stock.wt = stock.wt,
    mat = mat.q,
    m = m.q,
    harvest.spwn = flq0,
    m.spwn = flq0,
    discards.n = flq0,
    discards.wt = catch.wt,
    name = "FLicc stocklen"
  )

  lhpar.out <- FLPar(linf=linf,k=k,M=M,Mk=Mk,L50=L50,L95=L95,a=a,b=b,s=0.7)
  units(out) <- standardUnits(out)
  attr(out,"lhpar") <- lhpar.out
  return(out)
}

#' Spawning stock biomass for FLStockLen
#'
#' @param stk An \code{FLStockLen} object.
#'
#' @return An \code{FLQuant}.
#' @export
ssbl <- function(stk) {
  FLCore::quantSums(stock.n(stk) * stock.wt(stk) * mat(stk))
}

#' Compute catch numbers from stock numbers and mortalities
#'
#' Computes catch-at-length or catch-at-age in numbers using the Baranov
#' catch equation.
#'
#' @param stock.n An \code{FLQuant} of stock numbers.
#' @param harvest An \code{FLQuant} of fishing mortality.
#' @param m An \code{FLQuant} of natural mortality.
#'
#' @return An \code{FLQuant} of catch numbers.
#' @export
compute_catch.n <- function(stock.n, harvest, m) {

  z <- harvest + m
  cn <- stock.n * harvest / z * (1 - exp(-z))
  cn[z == 0] <- 0

  return(cn)
}


#' Mean fishing mortality over a length range
#'
#' @param stk An \code{FLStockLen} object.
#' @param lrange Optional numeric vector of length 2 giving min and max length.
#'
#' @return An \code{FLQuant}.
#' @export
fbarl <- function(stk, lrange = NULL) {

  L <- as.numeric(dimnames(harvest(stk))$len)

  idx <- if (is.null(lrange)) {
    seq_along(L)
  } else {
    which(L >= lrange[1] & L <= lrange[2])
  }

  FLCore::quantMeans(harvest(stk)[idx,,,,,])
}


#' Update an FLStockLen object with FLicc fitted outputs
#'
#' Fills fitted values into an existing \code{FLStockLen} template, updating
#' \code{stock.n}, \code{catch.n}, and \code{harvest} for a selected year.
#'
#' @param fit A fitted \code{flicc_tmb_fit} object returned by \code{fiticc()}.
#' @param stklen An \code{FLStockLen} object used as template.
#' @param year Optional year to update. Defaults to the last year in \code{stklen},
#'   consistent with the current \code{fiticc()} implementation.
#' @param R0 Recruitment scaling factor used to scale per-recruit outputs
#'   (default = 1000).
#' @param report_names Optional named list giving the names of the elements in
#'   \code{fit$report} corresponding to fitted vectors. Defaults are
#'   \code{N_F}, \code{catch_n}, and \code{FM}.
#'
#' @return An updated \code{FLStockLen} object.
#' @export
upd_stklen <- function(fit, stklen, year = NULL, R0 = 1000) {

  stopifnot(inherits(fit, "flicc_tmb_fit"))
  stopifnot(methods::is(stklen, "FLStockLen"))

  if (is.null(year)) {
    year <- tail(dimnames(stklen)$year, 1)
  } else {
    year <- as.character(year)
  }

  rep <- fit$report

  stky <- FLCore::window(stklen, start = year, end = year)

  dn <- dimnames(stock.n(stky))
  d  <- dim(stock.n(stky))

  make_FLQuant <- function(x) {
    FLCore::FLQuant(array(as.numeric(x), dim = d, dimnames = dn))
  }

  stock.n(stky) <- make_FLQuant(rep$N * R0)
  harvest(stky) <- make_FLQuant(rep$F_len)

  catch.n(stky)    <- compute_catch.n(stock.n(stky), harvest(stky), m(stky))
  landings.n(stky) <- catch.n(stky)
  discards.n(stky) <- catch.n(stky) * 0

  stock(stky)    <- FLCore::computeStock(stky)
  catch(stky)    <- FLCore::computeCatch(stky)
  landings(stky) <- catch(stky)
  discards(stky) <- catch(stky) * 0

  stklen[, year] <- stky

  return(stklen)
}


