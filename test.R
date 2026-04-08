
library(FLicc)




lhpar <- FLPar(
  linf = 55.7,
  k    = 0.08,
  M    = 0.162,
  L50  = 31.1859,
  a    = 0.004721956,
  b    = 3.146168,
  s    = 0.7
)

stklen <- stocklen(lfds,lhpar)

lfd_mu = FLQuants(lapply(lfds,function(x){
  yearSums(x[,ac(2022)])
}))


stklen <- stocklen(lfds,lhpar,m_model = "constant")
lfd=lfds
fit <- fiticc(lfd,stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear = c(0.78,0.22),settings=list(CVL=0.1,GL=50,catch.sd=0.05))

fit <- fiticc(lfd,stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear = c(0.75,0.25),settings=list(CVL=0.1,GL=50,catch.sd=0.05))

plot_fiticc(fit)
plot_plen(fit)
plot_plen(fit,by_gear=T)
matplot(fit$tmb_data$Lmid,fit$report$Sel,type="l",lwd=2)
fit$report$spr



plotAdvice(iter(om,1))






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
  Linf = c(40),
  CVL  = 0.10,
  Mk = c(2),
  sel_fun = c("dsnormal", "dsnormal", "logistic"),
  Catch = c(0.1802070, 0.2101353, 0.6096577),
  catch_sd = 0.05,
  gear_names = c("Estuarine set bagnet", "Gill net", "Marine set bagnet"),
  ref_length = 21.93152,
  a = 0.004721956,
  b = 3.146168,
  L50 = 23.2335,
  GL = 50
)


fit <- fiticc(dl,compile = F)


fit$report$spr

plot_fiticc(fit)
plot_plen(fit)
plot_plen(fit,by_gear=T)
matplot(dl$LLB,fit$report$Sel,type="l",lwd=2)



fit$par
fit$opt$convergence
fit$report$Fk
fit$report$Sel
fit$report$mu
fit$report$Mk
fit$report$Galpha # Change this to CVL

fit$report$spr


# SIM data


LLB <- an(dimnames(lfds[[1]])$len)
yrs = an(dimnames(lfds[[1]])$year)
fq <- lapply(lfds,function(x){
  (as.data.frame(yearMeans(x[,ac(2024)]))$data)
})


lhpar <- FLCore::FLPar(
  linf = 55.7,
  k    = 0.08,
  M    = 0.162,
  L50  = 31.1859,
  a    = 0.004721956,
  b    = 3.146168,
  s    = 0.7
)


dl <- list(
  model_name = "Base: M",
  LLB = LLB,
  fq = fq,
  Linf = c(55.7),
  sel_fun = c("dsnormal","logistic"),
  Catch = c(0.7,0.3),
  catch_sd = 0.05,
  gear_names = names(fq),
  Mk = c(0.1620/0.0800),
  CVL=0.05,
  ref_length = 50,
  a = 0.004721956,
  b = 3.146168,
  L50 =  31.1859,
  GL = 100
)


fit <- fiticc(dl,compile=F)
fit$report$spr

plot_fiticc(fit)
plot_plen(fit)
plot_plen(fit,by_gear=T)
matplot(dl$LLB,fit$report$Sel,type="l",lwd=2)

fit$par
fit$opt$convergence
fit$report$Fk
fit$report$Sel
fit$report$mu
fit$report$Mk
fit$report$Galpha # Change this to CVL

fit$report$spr


tmb_data <- data_tmb_flicc(dl)
parameters <- init_tmb_flicc(tmb_data)
#' # Fix logistic selectivity using SL50 and SL95
 sel_fixed <- list(
   "Trawl"  = list(rsd = 0.2)
)
 fix_res <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed)


 #'
#' fix_res <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed)
#'
#' # Fix right tail of a double-sided normal to be very broad
#' sel_fixed2 <- list(
#'   "Gill net" = list(rsd = 1000)
#' )
#'
#' fix_res2 <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed2)
