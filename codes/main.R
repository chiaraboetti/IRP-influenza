#######################################################
## Influenza model of Sweden during season 2017-2018 ##
#######################################################
folder = ### set working directory
setwd(folder)

library(dplyr)
library(ggplot2)
library(gridExtra)
library(sp)
library(gstat) # Inverse Distance Weighting interpolation
library(rgeos)
library(raster)
library(INLA)
library(TMB)
# ------------------------------------------------------------------------------
south_flag = FALSE
plot_flag = FALSE
unifPopDensity_flag = FALSE

## Load helper function
source('TMB/createU_EB_sweden.R')
#source('TMB/createU_EB_sweden_priors.R')

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

## Influenza in specific weeks
flu.weeks = unique(flu$Week)
flu.weeks[1:7] #par_mean1.RData
flu.weeks[9:15] #par_mean2.RData
flu.weeks[16:21] #par_mean3.RData

## Getting mesh with flu and population during the whole season
if(south_flag){
  cat('Focusing on the South of Sweden\n')
  source('get_mesh_sweden_south.R')
}else{
  source('get_mesh_sweden.R')
}
flu.df.grouped.original = flu.df.grouped

#  -----------------------------------------------------------------------------


## Computing logPM10 and estimating influenza using mean(logPM10) --------------
# Pollution during the whole season
if(south_flag){ ### TO BE CHECKED
  cat('Focusing on the South of Sweden\n')
  source('est_pollution_sweden_south.R')
}else{
  source('est_pollution_sweden.R')
}
mesh.df.original = mesh.df

# Actual for loop
par_mean = list()
i_week = 15
for(i_week in 1:7){ #c(1:7, 9:15, 16:21)){
  tmp_par = list(week = NA, time = NA, mesh.df = NA, par = NA, hessian = NA, beta_table = NA)
  week = flu.weeks[i_week] #i_week is index, whereas week is the actual week
  tmp_par$week = week
  cat('Studying week',week,'\n')
  
  ## Resetting variables at every iteration
  flu.df.grouped = flu.df.grouped.original[flu.df.grouped.original$Week == week, ]
  mesh.df = mesh.df.original
  
  start = Sys.time()
  ## Predicting weekly pollution
  source('pred_pollution_sweden.R')
  tmp_par$beta_table = beta_table
  cat('\nPredicting weekly pollution took:',Sys.time() - start,'\n')
  
  ## We get the following message for the whole Sweden:
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=172
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=612
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=2680
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=3384
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=4088
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=4220
  
  ## and for the South ofSweden:
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=136
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=486
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=2131
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=2691
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=3251
  # GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=3356
  
  ## Taking the mean
  mesh.df$PM10_mean = apply(mesh.df[,8:14], 1, mean)
  mesh.df$exp_PM10_mean = exp(mesh.df$PM10_mean)
  tmp_par$mesh.df = mesh.df
  
  ## Prep for TMB & HMC
  spde = inla.spde2.matern(mesh.true)
  
  ## Influenza model with PM10_mean 
  Data.tmb = list('y_i'=flu.df.grouped$Counts,
                  'N_i'=flu.df.grouped$PopCounts,
                  'Z_j'=mesh.df$PM10_mean,
                  #'Z_j'=mesh.df$exp_PM10_mean,
                  'M0'=spde$param.inla$M0,
                  'M1'=spde$param.inla$M1,
                  'M2'=spde$param.inla$M2,
                  'A'=inla.as.sparse(D))
  
  Params.tmb = list('beta_0'=0,
                    'beta_1'=1,
                    'theta'=spde$param.inla$theta.initial,
                    'S_j'=rep(0, spde$n.spde))
  
  
  ## Empirical Bayes Approach
  cat('Computing objective function in week',week,'\n')
  Obj.eb = MakeADFun(data = Data.tmb,
                     parameters = Params.tmb,
                     DLL = 'U_EB',
                     random = c('S_j'),
                     silent = T)
  # Obj.eb$fn(c(0,0,0,0))
  
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
  par_mean = c(par_mean, list(tmp_par))
}
tmp_par

signif(exp(par_mean[[1]]$par[1]), 3)
signif(exp(par_mean[[1]]$par[2]), 3) 

head(par_mean[[1]]$mesh.df, 10)

saveRDS(par_mean, 'SwedenData/all_par_mean_south1.RData')
# readRDS('SwedenData/all_par_mean2.RData')

##RMK: when using exp, weeks 2,14-21 (45,6-13) give an error:
## Error in optim(par = Obj.eb$par, fn = Obj.eb$fn, gr = Obj.eb$gr, method = "BFGS",  : 
##                  initial value in 'vmmin' is not finite


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



#  -----------------------------------------------------------------------------


## Estimating influenza using max(logPM10) -------------------------------------
par_max = list()
i_week = 1
for(i_week in 16:21){ #c(1:7, 9:15, 16:21)){
  tmp_par = list(week = NA, time = NA, mesh.df = NA, par = NA, hessian = NA)
  week = flu.weeks[i_week] #i_week is index, whereas week is the actual week
  tmp_par$week = week
  cat('Studying week',week,'\n')
  
  ## Resetting variables at every iteration
  flu.df.grouped = flu.df.grouped.original[flu.df.grouped.original$Week == week, ]

  ## Getting weekly pollution
  if(i_week < 8){
    if(south_flag){ ### TBD
      cat('Focusing on the South of Sweden\n')
      mesh.df = readRDS('SwedenData/all_par_mean_south.RData')[[i_week]]$mesh.df
    }else{
      mesh.df = readRDS('SwedenData/all_par_mean.RData')[[i_week]]$mesh.df
    }
  }else{
    if(south_flag){ ### TBD
      cat('Focusing on the South of Sweden\n')
      mesh.df = readRDS('SwedenData/all_par_mean_south.RData')[[i_week-1]]$mesh.df
    }else{
      mesh.df = readRDS('SwedenData/all_par_mean.RData')[[i_week-1]]$mesh.df
    }
  }
  
  ## Taking the max
  mesh.df$PM10_max = apply(mesh.df[,8:14], 1, max)
  mesh.df$exp_PM10_max = exp(mesh.df$PM10_max)
  tmp_par$mesh.df = mesh.df
  
  ## Prep for TMB & HMC
  spde = inla.spde2.matern(mesh.true)
  
  ## Influenza model with PM10_mean 
  Data.tmb = list('y_i'=flu.df.grouped$Counts,
                  'N_i'=flu.df.grouped$PopCounts,
                  'Z_j'=mesh.df$PM10_max,
                  #'Z_j'=mesh.df$exp_PM10_max,
                  'M0'=spde$param.inla$M0,
                  'M1'=spde$param.inla$M1,
                  'M2'=spde$param.inla$M2,
                  'A'=inla.as.sparse(D))
  
  Params.tmb = list('beta_0'=0,
                    'beta_1'=1,
                    'theta'=spde$param.inla$theta.initial,
                    'S_j'=rep(0, spde$n.spde))
  
  
  ## Empirical Bayes Approach
  cat('Computing objective function in week',week,'\n')
  Obj.eb = MakeADFun(data = Data.tmb,
                     parameters = Params.tmb,
                     DLL = 'U_EB',
                     random = c('S_j'),
                     silent = T)
  # Obj.eb$fn(c(0,0,0,0))
  
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
  par_max = c(par_max, list(tmp_par))
}
par_max

signif(exp(par_max[[1]]$par[1]), 3)
signif(exp(par_max[[1]]$par[2]), 3) 

saveRDS(par_max, 'all_par_max.RData')
# readRDS('all_par_max.RData')

# ------------------------------------------------------------------------------


## Summary stats with delta method ---------------------------------------------
cov.eb = solve(Opt.eb$hessian)
# Estimate for exp(beta_0)
signif(exp(Opt.eb$par[1]), 3) #0.000109
signif(exp(Opt.eb$par[1] + qnorm(c(0.025, 0.975))* sqrt(cov.eb[1, 1])), 3) #3.24e-05 3.66e-04

# Estimate for exp(beta_1)
signif(exp(Opt.eb$par[2]), 3) #1.17
signif(exp(Opt.eb$par[2] + qnorm(c(0.025, 0.975))* sqrt(cov.eb[2, 2])), 3) #0.76 1.81

# On log scale for phi and lambda2  #NOT SURE
cov.eb.log.sp.parms <- matrix(c(-2, 0, -2, -1), nrow=2) %*%
  cov.eb[3:4, 3:4] %*% t(matrix(c(-2, 0, -2, -1), nrow=2))

# Estimate for phi  #NOT SURE
sqrt(8)/exp(Opt.eb$par[3]) #576.2849
exp(log(sqrt(8)) - Opt.eb$par[4] + qnorm(c(0.025, 0.975)) *
      sqrt(cov.eb.log.sp.parms[2,2])) #0.004015636 1.868412635

# Estimate for lambda^2 #NOT SURE
1/(4*pi*exp(2*Opt.eb$par[4] + 2*Opt.eb$par[3])) #3.098217
exp(-log(4*pi) - 2 * Opt.eb$par[4] - 2 * Opt.eb$par[3] + qnorm(c(0.025, 0.975)) *
      sqrt(cov.eb.log.sp.parms[1,1])) #0.0218717 438.8753255
# ------------------------------------------------------------------------------

