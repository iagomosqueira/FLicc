
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




