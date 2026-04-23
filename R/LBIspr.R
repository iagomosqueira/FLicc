#' Gear-specific length-based indicator relative to Fspr reference
#'
#' Computes a gear-specific length-based indicator by comparing the observed
#' proportion above a reference threshold length to the expected proportion
#' above that threshold under equilibrium numbers-at-length at Fsprx,
#' filtered through gear-specific selectivity.
#'
#' The threshold is defined separately for each gear from:
#'   Nf(Fsprx) * sel_gear
#'
#' @param fit fitted flicc_tmb_fit object
#' @param gear optional character vector of gears to include; defaults to all
#' @param spr target SPR percentage, default 40
#' @param thresh cumulative threshold used to define Lref, default 0.75
#' @param nyears number of terminal years to average
#' @param scale_sel logical; passed to nf_flicc
#'
#' @return *FLIndices* with *FLIndexBiomass*
#' @export
LBIspr<- function(fit, gear = NULL, spr = 40, thresh = 0.75,
                         nyears = 1, scale_sel = TRUE) {


  gear_names <- names(fit$report$sel_gear)
  if (is.null(gear)) {
    gear <- gear_names
  }
  gear <- as.character(gear)
  yrs <- dimnames(fit$report$obslen[[1]])$year
  lens <- dimnames(fit$report$obslen[[1]])$len

  sel.pattern <- FLQuant(
    dimnames = list(
      len = lens,
      year = yrs,
      unit = "unique",
      season = "all",
      area = "unique",
      iter = "1"
    )
  )

  flqs = FLQuants(lapply(fit$report$Fk,function(x){
      x <-  FLQuant(x,quant='age')
      x[] <- NA
      units(x) = "lbi"
      x
  }))

  if (!all(gear %in% gear_names)) {
    stop("Unknown gear name(s): ", paste(setdiff(gear, gear_names), collapse = ", "))
  }

  # reference F on F scale
  Ftgt <- fspr_flicc(fit, spr = spr, nyears = nyears, input = "F")

  # equilibrium numbers-at-length at Fsprx
  Nref <- nf_flicc(fit, nyears = nyears, F = Ftgt, scale_sel = scale_sel)
  Len  <- as.numeric(dimnames(Nref)$len)



  vals <- numeric(length(gear))
  Lref_out <- numeric(length(gear))
  pref_out <- numeric(length(gear))
  pobs_out <- numeric(length(gear))

  LFDobs <- fit$report$obslen
  LFDref <- fit$report$obslen
  Lref <- NULL

  for (i in seq_along(gear)) {
    g <- gear[i]
    # gear-specific selectivity
    sg <- fit$report$sel_gear[[g]]


    # reference vulnerable numbers for this gear
    vref <- Nref * sg

    if (sum(vref, na.rm = TRUE) <= 0) {
      vals[i] <- NA_real_
      Lref_out[i] <- NA_real_
      pref_out[i] <- NA_real_
      pobs_out[i] <- NA_real_
      next
    }

    # define threshold length from cumulative vulnerable numbers
    vref2 <- vref[-1]
    Len2  <- Len[-1]

    cums = apply(vref2, 2:6, cumsum)
    n_thresh <- sum(vref2, na.rm = TRUE) * thresh
    Lthresh<- Len2[which.min((n_thresh - cums)^2)]
    Li <- ac(Len2[Len2>= Lthresh])

    LFDref[[g]][] <- vref
    Lref <- c(Lref,Lthresh)


    # expected proportion above threshold
    pref <- sum(vref[Li,], na.rm = TRUE) / sum(vref, na.rm = TRUE)

    for(y in seq(yrs)){
    # observed gear-specific length composition
    obs_g <- fit$report$obslen[[g]][,y]

    if (sum(obs_g, na.rm = TRUE) <= 0) {
      next
    }

    pobs <- sum(obs_g[Li,], na.rm = TRUE) / sum(obs_g, na.rm = TRUE)

    flqs[[g]][,y] <- pobs / pref
  }
  }

  out <- FLIndices(Map(function(x,y){
    idx <-FLIndexBiomass(index=x)
    idx@range[c("startf","endf")] =c(0.4,0.6)
    idx
    },x=flqs,y= fit$report$sel_gear))

 names(Lref) <- gear
 attr(out,"LFDref") <- LFDref
 attr(out,"Lref") <- Lref

 return(out)

}

