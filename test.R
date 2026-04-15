
library(FLicc)

data("alfonsino")

lhpar <- FLPar(
  linf = 55.7,
  k    = 0.08,
  M    = 0.162,
  L50  = 31.1859,
  a    = 0.004721956/1000,
  b    = 3.146168
)

lfd<- lfd_alfonsino
stklen <- stocklen(lfd,lhpar)

plot_lfd(lfd,type="relmax")

fit <- fiticc(lfd_alfonsino,stklen_alfonsino,sel_fun=c("dsnormal","logistic"),catch_by_gear = c(0.7,0.3))

# Option to create FLStockLen
stkl <- flicc_stklen(fit,stklen_alfonsino)

# Plotting
plot_len(fit)
plot_len(fit,by_gear = T,year=2020:2024)
# Status
plot_spr(fit)
plot_lbfao(fit)
# LBIspr
plot_LBIspr(fit,thresh = 0.7)

# Equilibrium dynamics
eqstk <- eqstklen(fit,s=0.75)
plot_eqcurves(eqstk)

