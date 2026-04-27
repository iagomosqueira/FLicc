
#' plot_hcrflicc
#'
#' Plot an SPR-based harvest control rule (HCR) for FLicc.
#'
#' The x-axis shows current SPR and the y-axis shows the relative change
#' in TAC or effort. The horizontal zero line indicates no change.
#'
#' Stock status zones are defined independently from the HCR triggers:
#'   - below blim: critical
#'   - blim to bthr: cautious / rebuilding
#'   - above bthr: healthyW
#'
#' The HCR is defined by four trigger points:
#'   - below b1: fixed low change (dlow)
#'   - b1 to b2: linear increase from dlow to dopt
#'   - b2 to b3: constant at dopt
#'   - b3 to b4: linear increase from dopt to dup
#'   - above b4: constant at dup
#'
#' @param b1 Lower HCR trigger for reduced effort/TAC.
#' @param b2 HCR trigger where no-change plateau starts.
#' @param b3 HCR trigger where upper slope starts.
#' @param b4 HCR trigger where maximum increase is reached.
#' @param blim Lower biological limit for background zone.
#' @param bthr Biological threshold separating cautious and healthy zones.
#' @param dlow Relative change below b1, e.g. -0.20 for -20%.
#' @param dopt Relative change in the target range, usually 0.
#' @param dup Relative change above b4, e.g. 0.10 for +10%.
#' @param metric Character string for y-axis label: "Effort" or "TAC".
#' @param spr_target Optional target SPR shown as a vertical dashed line.
#' @param show_spr Logical, if TRUE annotate the SPR target.
#' @param obs Optional numeric vector of observed SPR values.
#' @param years Optional year labels for obs.
#' @param xmax Maximum x-axis value.
#' @param ymin Minimum y-axis value.
#' @param ymax Maximum y-axis value.
#' @param alpha Transparency for background shading.
#' @param triggers Logical, if TRUE annotate trigger points and levels.
#' @param refpts Logical, if TRUE annotate trigger points and levels.
#' @param status.text Logical, if TRUE add stock-zone labels.
#' @param line_colour Colour of HCR curve.
#' @param line_size Line width of HCR curve.
#'
#' @return A ggplot object.
#' @examples
#'
#' plot_hcrflicc(
#'   spr_target = 0.4,
#'   show_target = TRUE,
#'   spr_target_label = "SPR target"
#' )
#' # example code
#'
#' @export
plot_hcrspr <- function(
    b1 = 0.1,
    b2 = 0.3,
    b3 = 0.5,
    b4 = 0.8,
    blim = 0.10,
    bthr = 0.3,
    dlow = -0.20,
    dopt =  0.00,
    dup  =  0.15,
    metric = c("TAC","Effort"),
    spr_target = 0.4,
    show_target = TRUE,
    spr_target_label = "SPR target",
    obs = NULL,
    years = NULL,
    xmax = 0.9,
    ymin = -0.30,
    ymax =  0.30,
    alpha = 0.18,
    triggers = TRUE,
    refpts = TRUE,
    status.text = FALSE,
    line_colour = "blue",
    line_size = 1.2
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  metric <- match.arg(metric)

  stopifnot(
    is.numeric(c(b1, b2, b3, b4, blim, bthr, dlow, dopt, dup)),
    b1 < b2,
    b2 <= b3,
    b3 < b4,
    blim < bthr
  )

  # background zones based on stock status, not HCR triggers
  zone_dat <- data.frame(
    xmin = c(0, blim, bthr),
    xmax = c(blim, bthr, xmax),
    ymin = ymin,
    ymax = ymax,
    zone = factor(
      c("Critical", "Cautious", "Healthy"),
      levels = c("Critical", "Cautious", "Healthy")
    )
  )

  # piecewise HCR curve
  hcr_fun <- function(x) {
    ifelse(
      x <= b1, dlow,
      ifelse(
        x < b2, dlow + (dopt - dlow) * (x - b1) / (b2 - b1),
        ifelse(
          x <= b3, dopt,
          ifelse(
            x < b4, dopt + (dup - dopt) * (x - b3) / (b4 - b3),
            dup
          )
        )
      )
    )
  }

  xseq <- seq(0, xmax, length.out = 500)
  hcr_dat <- data.frame(
    spr = xseq,
    change = hcr_fun(xseq)
  )

  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_rect(
      data = zone_dat,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = zone),
      alpha = alpha,
      colour = NA
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Critical" = "indianred2",
        "Cautious" = "wheat3",
        "Healthy"  = "lightblue3"
      ),
      guide = "none"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6) +
    ggplot2::geom_line(
      data = hcr_dat,
      ggplot2::aes(x = spr, y = change),
      colour = line_colour,
      linewidth = line_size
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, xmax),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(ymin, ymax),
      expand = c(0, 0),
      labels = function(z) paste0(round(z * 100), "%")
    ) +
    ggplot2::xlab(expression(SPR[curr])) +
    ggplot2::ylab(paste("Relative change in", metric)) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank()
    )

  p <- p + ggplot2::annotate("text", x = xmax * 0.98, y = dlow, label = paste0(round(100 * dlow), "%"), hjust = 1, vjust = -0.2) +
    ggplot2::annotate("text", x = xmax * 0.98, y = dopt, label = paste0(round(100 * dopt), "%"), hjust = 1, vjust = -0.2) +
    ggplot2::annotate("text", x = xmax * 0.98, y = dup,  label = paste0("+", round(100 * dup), "%"), hjust = 1, vjust = -0.2)

  if (!is.null(spr_target)) {
    p <- p +
      ggplot2::geom_segment(
        aes(
          x = spr_target, xend = spr_target,
          y = dopt, yend = dup*0.75
        ),
        linetype = 5,
        linewidth = 0.7,
        colour = "blue"
      )


      p <- p +
        ggplot2::annotate(
          "text",
          x = spr_target,
          y = dup*0.75 + 0.01 * (ymax - ymin),
          label = spr_target_label,
          colour = "blue",
          hjust = 0.5,
          vjust = 0,
          size = 3.5
        )

  }

  if (triggers) {
    p <- p +
      ggplot2::annotate("text", x = b1, y = ymin + 0.03 * (ymax - ymin), label = "b1", hjust = 1.1) +
      ggplot2::annotate("text", x = b2, y = ymin + 0.03 * (ymax - ymin), label = "b2", hjust = 1.1) +
      ggplot2::annotate("text", x = b3, y = ymin + 0.03 * (ymax - ymin), label = "b3", hjust = 1.1) +
      ggplot2::annotate("text", x = b4, y = ymin + 0.03 * (ymax - ymin), label = "b4", hjust = 1.1) +
      ggplot2::geom_vline(xintercept = c(b1, b2, b3, b4), linewidth = 0.4, colour = "blue",linetype=2)
  }
  if(refpts) {
    p <- p + ggplot2::geom_vline(xintercept = c(blim, bthr), linewidth = 0.6, colour = "black")
  }

  if (status.text) {
    p <- p +
      ggplot2::annotate("text", x = blim / 2, y = ymax * 0.80, label = "Critical", angle = 90) +
      ggplot2::annotate("text", x = (blim + bthr) / 2, y = ymax * 0.80, label = "Cautious", angle = 90) +
      ggplot2::annotate("text", x = (bthr + xmax) / 2, y = ymax * 0.80, label = "Healthy", angle = 90)
  }

  if (!is.null(obs)) {
    obs <- as.numeric(obs)
    obs_dat <- data.frame(
      spr = obs,
      change = hcr_fun(obs),
      idx = seq_along(obs)
    )

    p <- p +
      ggplot2::geom_path(
        data = obs_dat,
        ggplot2::aes(x = spr, y = change, group = 1),
        colour = "grey40",
        alpha = 0.7
      ) +
      ggplot2::geom_point(
        data = obs_dat,
        ggplot2::aes(x = spr, y = change),
        shape = 21,
        fill = c(rep("grey80", max(0, nrow(obs_dat) - 1)), "blue"),
        colour = "black",
        size = c(rep(2, max(0, nrow(obs_dat) - 1)), 3)
      )

    if (!is.null(years) && length(years) == length(obs)) {
      obs_dat$year <- years
      p <- p +
        ggplot2::geom_text(
          data = obs_dat,
          ggplot2::aes(x = spr, y = change, label = year),
          nudge_y = 0.02,
          size = 3
        )
    }
  }

  return(p)
}


#' hcr_spr (flicc)
#'
#' SPR-based harvest control rule for FLicc.
#'
#' Returns a multiplicative TAC/effort adjustment factor based on current SPR.
#' The rule follows the same piecewise structure as plot_hcrflicc():
#'   - spr <= b1: dlow
#'   - b1 < spr < b2: linear slope from dlow to dopt
#'   - b2 <= spr <= b3: dopt
#'   - b3 < spr < b4: linear slope from dopt to dhi
#'   - spr >= b4: dhi
#'
#' Relative changes are converted to multipliers as:
#'   multiplier = 1 + change
#'
#' @param spr Numeric SPR value(s).
#' @param b1 Lower trigger point.
#' @param b2 Trigger where no-change plateau starts.
#' @param b3 Trigger where upper slope starts.
#' @param b4 Trigger where maximum increase is reached.
#' @param dlow Relative change below b1, e.g. -0.20.
#' @param dopt Relative change between b2 and b3, usually 0.
#' @param dhi Relative change above b4, e.g. 0.10.
#' @param multiplier Logical, if TRUE return multiplier (default),
#'   otherwise return relative change.
#'
#' @return Numeric vector of multipliers or relative changes.
#'
#' @examples
#' # example code
#'
#' spr <- seq(0,1,0.02)
#' change <- spr_rule(spr=seq(0,1,0.02),
#' b1 = 0.20,
#' b2 = 0.35,
#' b3 = 0.6,
#' b4 = 0.8,
#' dlow = -0.20,
#' dopt =  0.00,
#' dhi  =  0.15)
#' plot(spr,change-1,type="l",col=4,lwd=2,ylab=("Change TAC"),xlab="SPR")
#' abline(0,0,lty=2)
#' plot(spr,change,type="l",col=4,lwd=2,ylab=("Multiplier"),xlab="SPR")
#' abline(0,0,lty=2)
#'
#' @export
spr_rule <- function(
    spr,
    b1 = 0.10,
    b2 = 0.3,
    b3 = 0.50,
    b4 = 0.75,
    dlow = -0.20,
    dopt =  0.00,
    dhi  =  0.15,
    multiplier = TRUE
) {
  stopifnot(
    is.numeric(spr),
    is.numeric(c(b1, b2, b3, b4, dlow, dopt, dhi)),
    b1 < b2,
    b2 <= b3,
    b3 < b4
  )

  change <- ifelse(
    spr <= b1, dlow,
    ifelse(
      spr < b2,
      dlow + (dopt - dlow) * (spr - b1) / (b2 - b1),
      ifelse(
        spr <= b3, dopt,
        ifelse(
          spr < b4,
          dopt + (dhi - dopt) * (spr - b3) / (b4 - b3),
          dhi
        )
      )
    )
  )

  if (multiplier) {
    return(1 + change)
  } else {
    return(change)
  }
}
#' Biomass ratio trend rule (2 over 3)
#'
#' Calculate a simple trend ratio from a biomass-related index time series.
#'
#' This rule compares the mean of the most recent \code{n1} years of an index
#' with the mean of the preceding \code{n2} years, returning the ratio:
#'
#' \deqn{
#' \bar{I}_{recent} / \bar{I}_{previous}
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\bar{I}_{recent}} is the mean over the last \code{n1} years
#'   \item \eqn{\bar{I}_{previous}} is the mean over the \code{n2} years
#'         immediately before that
#' }
#'
#' Values greater than 1 indicate that the recent biomass index is higher than
#' the earlier period, while values below 1 indicate a declining trend.
#'
#' @param idx An FLQuant-like index object, for example one element of
#'   \code{LBIspr(fit)}.
#' @param n1 Number of most recent years used to calculate the current mean.
#'   Default is 2.
#' @param n2 Number of preceding years used as the comparison period.
#'   Default is 3.
#'
#' @return A numeric scalar returned via \code{an()}, giving the ratio of the
#'   recent mean index to the previous mean index.
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item extracting the index values with \code{index(idx)}
#'   \item taking the last \code{n1} years and calculating their mean
#'   \item taking the \code{n2} years immediately before those and calculating
#'         their mean
#'   \item returning the ratio of recent mean to previous mean
#' }
#'
#' This can be used as a simple decision rule component for biomass or
#' recruitment indicators.
#'
#' @examples
#' # Biomass ratio trend rule
#' #
#' # idx <- LBIspr(fit)[["Gillnet"]]
#' # r_rule(idx)
#' # r_rule(idx, n1 = 3, n2 = 5)
#'
#' @export
r_rule <- function(idx, n1 = 2, n2 = 3) {
  an(
    yearMeans(tail(index(idx), n1)) /
      yearMeans(head(tail(index(idx), n1 + n2), n2))
  )
}

#' Combined SPR and index trend harvest control rule
#'
#' Apply a combined harvest control rule using current SPR and a recent
#' biomass-ratio trend indicator.
#'
#' The rule first calculates a gear-specific biomass trend ratio using
#' \code{r_rule()} applied to \code{LBIspr(fit, thresh = thresh)[[gear]]}.
#' It then calculates an SPR-based change multiplier using \code{spr_rule()}.
#'
#' The SPR multiplier is applied only when the biomass trend and SPR rule point
#' in the same direction:
#' \itemize{
#'   \item if \code{r < 1} and \code{s < 1}, apply \code{s}
#'   \item if \code{r > 1} and \code{s > 1}, apply \code{s}
#'   \item otherwise, return \code{1}, meaning no change
#' }
#'
#' @param fit A fitted FLicc object.
#' @param gear Numeric or character index selecting the gear-specific indicator
#'   from \code{LBIspr(fit, thresh = thresh)}.
#' @param thresh Threshold passed to \code{LBIspr()}.
#' @param b1 Lower SPR trigger point.
#' @param b2 SPR trigger where the no-change plateau starts.
#' @param b3 SPR trigger where the upper slope starts.
#' @param b4 SPR trigger where the maximum increase is reached.
#' @param dlow Relative change below \code{b1}, e.g. \code{-0.20} for a
#'   20\% reduction.
#' @param dopt Relative change between \code{b2} and \code{b3}, usually
#'   \code{0}.
#' @param dhi Relative change above \code{b4}, e.g. \code{0.15} for a
#'   15\% increase.
#' @param nyrs Number of most recent years of \code{fit$report$spr} used to
#'   calculate mean current SPR. Default is 1.
#' @param n1 Number of recent years used by \code{r_rule()} for the numerator
#'   mean. Default is 2.
#' @param n2 Number of preceding years used by \code{r_rule()} for the
#'   denominator mean. Default is 3.
#'
#' @return A numeric multiplier. Values below 1 imply a reduction, 1 implies
#'   no change, and values above 1 imply an increase.
#'
#' @examples
#' # Default combined SPR-ratio rule
#' # hcr_sprr(fit, gear = 1)
#'
#' # Custom SPR control points
#' # hcr_sprlbi(
#' #   fit, gear = 1,
#' #   b1 = 0.10, b2 = 0.35, b3 = 0.50, b4 = 0.75,
#' #   dlow = -0.20, dopt = 0, dhi = 0.15
#' # )
#'
#' # Use a three-year mean SPR and alternative trend windows
#' # hcr_sprr(fit, gear = 1, nyrs = 3, n1 = 3, n2 = 5)
#'
#' @export
hcr_sprlbi <- function(fit, gear,
                     thresh = 0.75,
                     b1 = 0.10,
                     b2 = 0.35,
                     b3 = 0.50,
                     b4 = 0.75,
                     dlow = -0.20,
                     dopt =  0.00,
                     dhi  =  0.15,
                     nyrs = 1,
                     n1 = 2,
                     n2 = 3) {

  idx <- LBIspr(fit, thresh = thresh)[[gear]]

  r <- r_rule(idx, n1 = n1, n2 = n2)

  spr <- mean(tail(fit$report$spr, nyrs))

  s <- spr_rule(
    spr,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    b4 = b4,
    dlow = dlow,
    dopt = dopt,
    dhi = dhi
  )

  mult <- ifelse((r < 1 && s < 1) | (r > 1 && s > 1), s, 1)

  return(mult)
}

