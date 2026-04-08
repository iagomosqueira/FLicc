
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




#' Construct an FLStockLen object from length bins and life-history parameters
#'
#' Builds an \code{FLStockLen} object from the length bins in a length-frequency
#' \code{FLQuant} together with scalar life-history inputs in an \code{FLPar}
#' object. The constructor derives midpoint lengths, weight-at-length,
#' maturity-at-length, and natural mortality-at-length under a chosen mortality
#' model.
#'
#' @param lfd An \code{FLQuant} containing length-frequency data. The
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
#'   lfd = lfds,
#'   lhpar = lhpar,
#'   m_model = "constant"
#' )
#'  stklen@lhpar
#' }
#'
#' @export

stocklen <- function(lfd, lhpar, m_model = c("constant", "inverse", "Lorenzen", "Gislason"),reflen=NULL) {
  m_model <- match.arg(m_model)

  pars <- dimnames(lhpar)$params


  if (!inherits(lfd, "FLQuants")) {
    stop("lfd must be class FLQuants")
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
  if (anyNA(lens)) stop("lfd must have numeric 'len' dimnames")

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

  attr(out,"lhpar") <- lhpar.out
  return(out)
}
