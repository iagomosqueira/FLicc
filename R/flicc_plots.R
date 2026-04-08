
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
