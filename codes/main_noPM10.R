#######################################################
## Influenza model of Sweden during season 2017-2018 ##
#######################################################
folder = ### set working directory
setwd(folder)

library(dplyr)
library(ggplot2)
library(gridExtra)
#library(rgdal) # will be retired during 2023
library(sp)
library(gstat) # Inverse Distance Weighting interpolation
library(rgeos)
library(raster)
library(INLA)
library(TMB)
# ------------------------------------------------------------------------------
south_flag = FALSE
plot_flag = FALSE
unifPopDensity_flag = TRUE

## Load helper function
#source('TMB/createU_EB_sweden_noPM10.R')
source('TMB/createU_EB_sweden_noPM10_priors.R')

## Setup Sweden
swe.shp = shapefile('SwedenData/Shapefile/SWE_adm1.shp')
flu = read.csv('SwedenData/df2_1718.csv')
cities = read.csv('SwedenData/coord_cities.csv')[-5,]
pop.dat = raster('SwedenData/popCount_swe2018_1km_image.tif')
pm10 = read.csv('SwedenData/data_1718.csv')[,-12]
border1 = read.csv('SwedenData/border1.csv')
border2 = read.csv('SwedenData/border3.csv')
if(south_flag){
  border1 = read.csv('SwedenData/bordersouth1.csv')
  border2 = read.csv('SwedenData/bordersouth2.csv')
}
colnames(border1) = colnames(border2) = c('X', 'Longitude', 'Latitude')

flu.weeks = unique(flu$Week)
flu.weeks[c(1:7, 9:21)]
flu.weeks[1:7] #param1.RData
flu.weeks[9:15] #param2.RData
flu.weeks[16:21] #param3.RData

## Getting the estimates -------------------------------------------------------
param = list() #9:21 done
for(i_week in 1:7){ #c(1:7, 9:15, 16:21)){
  tmp_par = list(week = NA, time = NA, mesh.df = NA, par = NA, hessian = NA)
  week = flu.weeks[i_week]
  tmp_par$week = week
  cat('Studying week',week,'\n')
  
  ## Get flu and population for the current week
  if(south_flag){ 
    cat('Focusing on the South of Sweden\n')
    source('get_flu_sweden_south.R')
  }else{
    source('get_flu_sweden.R')
  }
  tmp_par$mesh.df = mesh.df
  
  ## Prep for TMB & HMC
  spde = inla.spde2.matern(mesh.true)
  
  ## Influenza model with PM10_mean 
  Data.tmb = list('y_i'=flu.df.grouped$Counts,
                  'N_i'=flu.df.grouped$PopCounts,
                  'M0'=spde$param.inla$M0,
                  'M1'=spde$param.inla$M1,
                  'M2'=spde$param.inla$M2,
                  'A'=inla.as.sparse(D))
  
  Params.tmb = list('beta_0'=0,
                    'theta'=spde$param.inla$theta.initial,
                    'S_j'=rep(0, spde$n.spde))
  
  
  ## Empirical Bayes Approach
  cat('Computing objective function in week',week,'\n')
  Obj.eb = MakeADFun(data = Data.tmb,
                     parameters = Params.tmb,
                     DLL = 'U_EB_noPM10',
                     random = c('S_j'),
                     silent = T)
  # Obj.eb$fn(c(0,0,0))
  
  start = Sys.time()
  cat('Finding optimal parameters in week',week,'--- Started at',as.character(start),'\n')
  Opt.eb = optim(par = Obj.eb$par,
                 fn = Obj.eb$fn,
                 gr = Obj.eb$gr,
                 method = 'BFGS',
                 hessian = T)
  
  tmp_par$par = Opt.eb$par
  tmp_par$hessian = Opt.eb$hessian
  tmp_par$time = Sys.time() - start
  
  cat('*** Week',week,'is done! --- *** Time:',Sys.time() - start,'\n\n')
  param = c(param, list(tmp_par))
}
param[[1]]
signif(param[[1]]$par, 3)
signif(exp(param[[1]]$par), 3)
# saveRDS(param, 'SwedenData/par_noPM10_new_unifpop_prior.RData')

param.short = readRDS('SwedenData/par_noPM10_new_unifpop_prior.RData')
param.short.df = as.data.frame(
  t(sapply(1:length(param.short), function(x){
    c(param.short[[x]]$week, param.short[[x]]$par)
  }))
)
param.short.df = as.data.frame(apply(param.short.df, 2, as.numeric)) %>%
  dplyr::add_row(V1 = c(51,52), .before = 8)
param.short.df$ID.order = 1:nrow(param.short.df)
param.short.df$beta_0.exp = signif(exp(param.short.df$beta_0), 3)

param.long = readRDS('SwedenData/par_noPM10_new_unifpop.RData')
param.df = as.data.frame(
  t(sapply(1:length(param.long), function(x){
    c(param.long[[x]]$week, param.long[[x]]$par)
  }))
)
param.df = as.data.frame(apply(param.df, 2, as.numeric)) %>%
  dplyr::add_row(V1 = c(51,52), .before = 8)
param.df$ID.order = 1:nrow(param.df)
param.df$beta_0.exp = signif(exp(param.df$beta_0), 3)

plot(param.df$ID.order, param.df$beta_0.exp, type = 'b', lwd = 2, pch = 19,
     xaxt = 'n', ylim = c(0,1.3), col = 1,
     xlab = 'Weeks', ylab = 'exp(beta_0)',
     cex.axis = 1.5, cex.lab=1.5,
     main = 'Estimation of exp(beta_0) with uniform pop density')
lines(param.short.df$ID.order, param.short.df$beta_0.exp, type = 'b', lwd = 2, pch = 19,
      col = 6)
axis(1, at = param.df$ID.order, labels = param.df$V1)
legend('topright', legend = c('No prior', 'Prior'), lty = 1, 
       col = c(1,6), lwd = 2)

# ------------------------------------------------------------------------------


## Plotting --------------------------------------------------------------------
param.short = readRDS('SwedenData/par_noPM10_simple_unifpop.RData')
param.short.df = as.data.frame(
  t(sapply(1:length(param.short), function(x){
    c(param.short[[x]]$week, param.short[[x]]$par)
  }))
)
param.short.df = as.data.frame(apply(param.short.df, 2, as.numeric)) %>% 
  dplyr::add_row(V1 = c(51,52), .before = 8)
param.short.df$ID.order = 1:nrow(param.short.df)
param.short.df$beta_0.exp = signif(exp(param.short.df$beta_0), 3)

param.long = readRDS('SwedenData/par_noPM10_new_unifpop.RData')
param.df = as.data.frame(
  t(sapply(1:length(param.long), function(x){
    c(param.long[[x]]$week, param.long[[x]]$par)
  }))
)
param.df = as.data.frame(apply(param.df, 2, as.numeric))%>%
  dplyr::add_row(V1 = c(51,52), .before = 8)
param.df$ID.order = 1:nrow(param.df)
param.df$beta_0.exp = signif(exp(param.df$beta_0), 3)

plot(param.df$ID.order, param.df$beta_0.exp, type = 'b', lwd = 2, pch = 19, 
     xaxt = 'n', xlab = 'Weeks', ylab = 'exp(beta_0)', ylim = c(0,0.8),
     main = 'Estimation of exp(beta_0) with simplified uniform pop density')
lines(param.short.df$ID.order, param.short.df$beta_0.exp, type = 'b', lwd = 2, pch = 19,
      col = 4)
axis(1, at = param.df$ID.order, labels = param.df$V1)
legend('topright', legend = c('Unif Pop', 'Simple Unif Pop'), lty = 1, cex = 0.8, 
       col = c(1,4), lwd = 2)


tmp = readRDS('../Influenza/par_noPM10_unifpop.RData')
param.df.old = as.data.frame(
  t(sapply(c(1,3:length(tmp)), function(x){
    c(tmp[[x]]$week, tmp[[x]]$par)
  }))
)
param.df.old = as.data.frame(apply(param.df.old, 2, as.numeric)) %>%
  dplyr::add_row(V1 = 45, .after = 1) %>%
  dplyr::add_row(V1 = c(51,52), .before = 8) 
param.df.old$ID.order = 1:nrow(param.df.old)
param.df.old$beta_0.exp = signif(exp(param.df.old$beta_0), 3)

lines(param.df.old$ID.order, param.df.old$beta_0.exp, type = 'b', lwd = 2, pch = 19,
      col = 2)
legend('topright', legend = c('Unif Pop', 'Simple Unif Pop','Old Unif Pop'), lty = 1, 
       cex = 0.8, col = c(1,4,2), lwd = 2)
# ------------------------------------------------------------------------------


## Summary stats with delta method ---------------------------------------------
cov.eb = solve(Opt.eb$hessian)
# Estimate for exp(beta_0)
signif(exp(Opt.eb$par[1]), 3)
signif(exp(Opt.eb$par[1] + qnorm(c(0.025, 0.975))* sqrt(cov.eb[1, 1])), 3)


# On log scale for phi and lambda2
cov.eb.log.sp.parms <- matrix(c(-2, 0, -2, -1), nrow=2) %*% 
  cov.eb[2:3, 2:3] %*% t(matrix(c(-2, 0, -2, -1), nrow=2))

# Estimate for phi
sqrt(8)/exp(Opt.eb$par[3])
exp(log(sqrt(8)) - Opt.eb$par[3] + qnorm(c(0.025, 0.975)) * 
      sqrt(cov.eb.log.sp.parms[2,2]))

# Estimate for lambda^2
1/(4*pi*exp(2*Opt.eb$par[3] + 2*Opt.eb$par[2]))
exp(-log(4*pi) - 2 * Opt.eb$par[3] - 2 * Opt.eb$par[2] + qnorm(c(0.025, 0.975)) * 
      sqrt(cov.eb.log.sp.parms[1,1]))
# ------------------------------------------------------------------------------


