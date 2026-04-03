
library(FLicc)

dl <- list(
  model_name = "Base: Mk ~ Length-inverse",
  LLB = 3:45,
  fq = list(
    `Estuarine set bagnet` = c(2,19,19,25,95,106,105,220,246,268,266,185,274,
                               213,240,206,165,137,119,122,88,65,43,29,27,15,5,
                               15,5,6,6,14,25,1,1,5,0,9,0,2,0,0,3),
    `Gill net` = c(0,0,0,0,0,0,1,14,11,42,38,49,125,133,160,222,226,341,195,
                   251,283,318,322,183,203,282,142,154,42,62,42,33,28,15,12,
                   11,6,9,1,4,0,0,0),
    `Marine set bagnet` = c(1,0,0,1,1,3,14,54,63,118,130,228,303,344,532,649,
                            692,954,895,1032,828,793,755,621,581,472,295,362,
                            208,203,94,95,35,38,13,33,9,18,7,11,0,2,2)
  ),
  Linf = c(40, 0.01),
  sel_fun = c("dsnormal", "dsnormal", "dsnormal"),
  Catch = c(0.1802070, 0.2101353, 0.6096577),
  gear_names = c("Estuarine set bagnet", "Gill net", "Marine set bagnet"),
  Mk = c(2,0.01),
  ref_length = 21.93152,
  a = 0.004721956,
  b = 3.146168,
  L50 = 23.2335,
  GL = 20
)


fit <- fiticc(dl)

fit$par
fit$opt$convergence
fit$report$Fk
fit$report$Sel
fit$report$mu
fit$report$spr
fit$report$Mk
fit$report$Galpha

plot_fiticc(fit)
