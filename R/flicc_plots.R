
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

