
#' Plot observed and fitted length compositions by gear using ggplot2
#'
#' Produces a faceted ggplot comparing observed length-frequency counts with
#' fitted expected counts for each gear in a fitted FLicc/TMB model.
#'
#' @param fit A fitted model object returned by \code{fiticc()}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_fiticc <- function(fit) {
  obs <- fit$tmb_data$obs
  mu  <- fit$report$mu
  L   <- fit$tmb_data$Lmid
  gn  <- fit$tmb_data$gear_names

  nlen  <- nrow(obs)
  ngear <- ncol(obs)

  df_obs <- data.frame(
    length = rep(L, times = ngear),
    gear   = rep(gn, each = nlen),
    count  = as.vector(obs),
    stringsAsFactors = FALSE
  )

  df_fit <- data.frame(
    length = rep(L, times = ngear),
    gear   = rep(gn, each = nlen),
    count  = as.vector(mu),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = df_obs,
      ggplot2::aes(x = length, xend = length, y = 0, yend = count),
      linewidth = 0.8
    ) +
    ggplot2::geom_line(
      data = df_fit,
      ggplot2::aes(x = length, y = count),
      linewidth = 1
    ) +
    ggplot2::facet_wrap(~ gear, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = "Length", y = "Count") +
    ggplot2::theme_bw()
}

#' Plot FLicc size compositions in an LBSPR-style layout
#'
#' Creates a ggplot showing observed and fitted size compositions as
#' grey bars plus a fitted line. By default the plot uses the catch-weighted
#' overall proportions reported by TMB. Optionally, gear-specific panels can
#' be drawn using the within-gear observed and predicted proportions.
#'
#' @param fit A fitted FLicc object returned by \code{fiticc()}.
#' @param by_gear Logical. If \code{FALSE} (default), plot the overall
#'   catch-weighted observed and predicted proportions at length. If \code{TRUE},
#'   plot separate panels by gear.
#' @param title Optional plot title.
#' @param observed_fill Fill colour for observed bars.
#' @param observed_colour Outline colour for observed bars.
#' @param fitted_colour Line colour for fitted proportions.
#' @param fitted_linewidth Line width for fitted proportions.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_plen <-  function(fit,
                       by_gear = FALSE,
                       title = NULL,
                       observed_fill = "grey75",
                       observed_colour = "black",
                       fitted_colour = "blue",
                       fitted_linewidth = 1,
                       len_by = 2) {

  L <- fit$tmb_data$Lmid
  binwidth <- if (length(L) > 1) min(diff(L)) * 0.9 else 0.9

  x_breaks <- seq(
    from = floor(min(L) / len_by) * len_by,
    to   = ceiling(max(L) / len_by) * len_by,
    by   = len_by
  )

  if (!by_gear) {
    df_obs <- data.frame(
      length = L,
      proportion = fit$report$obs_p_l,
      stringsAsFactors = FALSE
    )

    df_fit <- data.frame(
      length = L,
      proportion = fit$report$pred_p_l,
      stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot() +
      ggplot2::geom_col(
        data = df_obs,
        ggplot2::aes(x = length, y = proportion),
        width = binwidth,
        fill = observed_fill,
        colour = observed_colour
      ) +
      ggplot2::geom_line(
        data = df_fit,
        ggplot2::aes(x = length, y = proportion),
        colour = fitted_colour,
        linewidth = fitted_linewidth
      ) +
      ggplot2::scale_x_continuous(breaks = x_breaks) +
      ggplot2::labs(
        x = "Length",
        y = "Proportion at length",
        title = title
      ) +
      ggplot2::theme_bw()

    return(p)
  }

  gn <- fit$tmb_data$gear_names
  nlen <- length(L)
  ngear <- length(gn)

  df_obs <- data.frame(
    length = rep(L, times = ngear),
    gear = rep(gn, each = nlen),
    proportion = as.vector(fit$report$obs_p_lg),
    stringsAsFactors = FALSE
  )

  df_fit <- data.frame(
    length = rep(L, times = ngear),
    gear = rep(gn, each = nlen),
    proportion = as.vector(fit$report$pred_p_lg),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = df_obs,
      ggplot2::aes(x = length, y = proportion),
      width = binwidth,
      fill = observed_fill,
      colour = observed_colour
    ) +
    ggplot2::geom_line(
      data = df_fit,
      ggplot2::aes(x = length, y = proportion),
      colour = fitted_colour,
      linewidth = fitted_linewidth
    ) +
    ggplot2::facet_wrap(~ gear, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::labs(
      x = "Length",
      y = "Proportion at length",
      title = title
    ) +
    ggplot2::theme_bw()
}

#' Plot observed length-frequency data
#'
#' Plots raw or scaled length-frequency distributions (LFDs) from an
#' \code{FLQuants} object, faceted by year and coloured by fishery or gear.
#'
#' The function supports plotting raw frequencies or values scaled relative
#' to the maximum within each \code{FLQuant} (e.g. per gear), allowing
#' comparison of distribution shapes independent of absolute magnitude.
#'
#' @param lfd An \code{FLQuants} object containing observed length-frequency data.
#' @param iter Numeric. Iteration to plot (default = 1).
#' @param len_by Numeric. Step size used to generate x-axis breaks (default = 2).
#' @param yticks Optional numeric vector of y-axis breaks.
#' @param scales Character. Passed to \code{facet_wrap()} (e.g. "fixed", "free_y").
#' @param type Character. Either:
#' \describe{
#'   \item{"raw"}{Plot observed frequencies (default).}
#'   \item{"relmax"}{Scale each \code{FLQuant} by its maximum value.}
#' }
#' @param colours Optional named vector of colours corresponding to \code{qname}.
#' If \code{NULL}, colours are assigned automatically.
#' @param legend_title Character. Title for the legend.
#' @param ncol Optional number of columns in facet layout.
#' @param lwd Numeric. Line width (default = 0.5).
#'
#' @details
#' Length bins are extracted from the \code{FLQuant} dimension names, ensuring
#' consistency with the underlying FLR structure. The function converts the data
#' to a \code{data.frame} for plotting with \pkg{ggplot2}.
#'
#' When \code{type = "relmax"}, each \code{FLQuant} is normalised by its maximum
#' value across the length dimension, preserving relative distribution shapes.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{FLQuants}}, \code{\link{FLQuant}}
#'
#' @author Henning Winker
#'
#' @examples
#' \dontrun{
#' library(FLicc)
#'
#' data(alfonsino)
#'
#' # Raw observed LFDs
#' plot_lfd(lfd_alfonsino)
#'
#' # Relative distributions (scaled to max per gear)
#' plot_lfd(lfd_alfonsino, type = "relmax")
#'
#' # Custom axis ticks
#' plot_lfd(
#'   lfd_alfonsino,
#'   len_by = 20,
#'   yticks = seq(0, 80, by = 20)
#' )
#'
#' # Free y-scale for better visual comparison
#' plot_lfd(lfd_alfonsino, scales = "free_y")
#'
#' # Custom colours
#' plot_lfd(
#'   lfd_alfonsino,
#'   colours = c(Trawl = "steelblue", Gillnet = "tomato")
#' )
#' }
#'
#' @export
plot_lfd <- function(lfd,
                     type=c("raw","relmax"),
                     iter = 1,
                     len_by = 10,
                     yticks = NULL,
                     scales = "fixed",
                     colours = NULL,
                     legend_title = "Fishery",
                     ncol = NULL,
                     lwd = 0.5) {

  type <- match.arg(type)

  L <- as.numeric(dimnames(lfd[[1]])$len)

  lfd <- FLCore::iter(lfd, iter)
  if(type=="relmax")
  lfd <- lapply(lfd, function(x){
    x%/%apply(x,2:5,max)
  })
  lfd <- as.data.frame(lfd)

  x_breaks <- seq(
    from = floor(min(L, na.rm = TRUE) / len_by) * len_by,
    to   = ceiling(max(L, na.rm = TRUE) / len_by) * len_by,
    by   = len_by
  )

  qlevels <- unique(as.character(lfd$qname))

  if (is.null(colours)) {
    base_cols <- c(
      "#D55E00", "#0072B2", "#009E73", "#CC79A7",
      "#E69F00", "#56B4E9", "#F0E442", "#999999"
    )
    colours <- stats::setNames(
      rep(base_cols, length.out = length(qlevels)),
      qlevels
    )
  } else {
    if (is.null(names(colours))) {
      stop("'colours' must be a named vector with names matching qname levels")
    }
    colours <- colours[qlevels]
  }

  p <- ggplot2::ggplot(lfd) +
    ggplot2::geom_line(
      ggplot2::aes(x = len, y = data, colour = qname),
      linewidth = lwd
    ) +
    ggplot2::facet_wrap(~year, scales = scales, ncol = ncol) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Length", y = "Frequency", colour = legend_title) +
    ggplot2::theme(
      legend.position = "right",
      strip.background = ggplot2::element_rect(fill = "grey90"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(yticks)) {
    p <- p + ggplot2::scale_y_continuous(breaks = yticks)
  }

  p
}

