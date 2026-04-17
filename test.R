
library(FLicc)

data("alfonsino")

# Specify Life History
lhpar <- FLPar(
  linf = 55.7,
  k    = 0.08,
  M    = 0.162,
  L50  = 31.1859,
  a    = 0.004721956/1000,
  b    = 3.146168
)

# LFD observations
lfd<- lfd_alfonsino
plot_lfd(lfd,type="relmax")

# Build FLStockLen input
stklen <- stocklen(lfd,lhpar,m_model="constant")

# check life history input
plot_lw(stklen)
plot_mat(stklen)
plot_m(stklen)

# Check m model shapes
m_models = c("constant", "inverse", "Lorenzen", "Gislason")
# loop through

m_stks <- FLStocks(lapply(m_models,function(x){
  stocklen(lfd,lhpar,m_model=x)
}))
names(m_stks) <- m_models
plot_m(m_stks)

# Fit model
fit <- fiticc(lfd, stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear = c(0.7,0.3))


# Plotting observed vs predicted
plot_len(fit)
plot_len(fit,by_gear = T,year=2020:2024)

# Plot Selectivity by gear
plot_sel(fit)

# Option to create FLStockLen
stkl <- flicc_stklen(fit)
# Plot, e.g., fishery selectivity weighted by the ratio of catches
plot_sel(stkl)


# Equilibrium dynamics
eqstk <- eqstklen(fit,s=0.75)
eqstk@refpts
plot_eqcurves(eqstk)



# LBIspr indicator
plot_LBIspr(fit,thresh = 0.75)

# Status
plot_spr(fit)
plot_lbfao(fit)

# Convert to simplified FLStockR
stk <- flicc2FLStockR(fit)
stk@refpts
plot_LBAdvice(stk)
# relative to MSY proxy
stkr <- flicc2FLStockR(fit,rel=T)
plot_LBAdvice(stkr)


# LBIspr indicator
plot_LBIspr(fit,thresh = 0.75)

