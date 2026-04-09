
library(FLicc)
data("alfonsino")

lhpar <- FLPar(
  linf = 55.7,
  k    = 0.08,
  M    = 0.162,
  L50  = 31.1859,
  a    = 0.004721956,
  b    = 3.146168,
  s    = 0.7
)

lfds<- lfd_alfonsino
stklen <- stocklen(lfds,lhpar)

plot_lfd(lfds,type="relmax")


fit <- fiticc(lfds,stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear = c(0.78,0.22),settings=list(CVL=0.1,GL=50,catch.sd=0.05))

plot_plen(fit)
plot_plen(fit,by_gear=T)
matplot(fit$tmb_data$Lmid,fit$report$Sel,type="l",lwd=2)
fit$report$spr





