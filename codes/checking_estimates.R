#######################################################
## Checking influenza estimates                      ##
#######################################################
folder = '~/Documents/SAMBa/IRP/AAAFinal'
setwd(folder)

library(dplyr)
library(ggplot2)
##

get_params = function(tmp){
  params = as.data.frame(
    t(sapply(1:length(tmp), function(x){
      c(tmp[[x]]$week, tmp[[x]]$par)#, tmp[[x]]$time)
    }))
  )
  colnames(params) = c('week', 'beta_0', 'beta_1', 'theta_1', 'theta_2')#, 'time')
  params = as.data.frame(apply(params, 2, as.numeric)) %>%
    dplyr::add_row(week = c(51,52), .before = 8)
  
  params$ID.order = 1:nrow(params)
  params$beta_0.exp = signif(exp(params$beta_0), 3)
  params$beta_1.exp = signif(exp(params$beta_1), 3)
  return(params)
}

get_params_noPM10 = function(tmp, prior_flag = FALSE){
  params = as.data.frame(
    t(sapply(1:length(tmp), function(x){
      c(tmp[[x]]$week, tmp[[x]]$par, tmp[[x]]$time)
    }))
  )
  if(prior_flag){
    colnames(params) = c('week', 'theta_1', 'theta_2', 'beta_0', 'time')
  }else{
    colnames(params) = c('week', 'beta_0', 'theta_1', 'theta_2', 'time')
  }
  params = as.data.frame(apply(params, 2, as.numeric)) %>%
    dplyr::add_row(week = c(51,52), .before = 8)
  
  params$ID.order = 1:nrow(params)
  params$beta_0.exp = signif(exp(params$beta_0), 3)
  return(params)
}


## PM10 during the whole season and mean and max of logPM10 --------------------
tmp = readRDS('SwedenData/all_par_mean.RData')
params = get_params(tmp)
head(params, 5)
params_max = get_params(readRDS('SwedenData/all_par_max.RData'))
head(params_max, 5)

par(mfrow=c(2,1))
plot(params$ID.order, params$beta_0.exp, type = 'b', lwd = 2, pch = 19, xaxt = 'n',
     ylim = c(0,0.8), col = 1, cex.axis = 1.5, cex.lab = 1.5,
     xlab = 'Weeks', ylab = 'exp(beta_0)', main = 'Estimation of exp(beta_0)')
lines(params_max$ID.order, params_max$beta_0.exp, type = 'b', lwd = 2, pch = 19,
      col = 4)
axis(1, at = params$ID.order, labels = params$week)
legend('topright', legend = c('Mean', 'Max'), lty = 1, lwd = 2, col = c(1,4))
legend('bottomright', legend = c('Mean', 'Max'), lty = 1, lwd = 2, cex = 0.5,
       col = c(1,4), inset=c(0,1),horiz =TRUE, xpd = TRUE)

plot(params$ID.order, params$beta_1.exp, type = 'b', lwd = 2, pch = 19, xaxt = 'n',
     ylim = c(0,2.5), col = 1, cex.axis = 1.5, cex.lab = 1.5,
     xlab = 'Weeks', ylab = 'exp(beta_1)', main = 'Estimation of exp(beta_1)')
abline(h = 1, lty = 2, col = 'red', lwd = 2)
lines(params$ID.order, params$beta_1.exp, type = 'b', lwd = 2, pch = 19,
      col = 1)
lines(params_max$ID.order, params_max$beta_1.exp, type = 'b', lwd = 2, pch = 19,
      col = 4)
axis(1, at = params$ID.order, labels = params$week, cex.axis = 1.5, cex.lab = 1.5)
legend('topright', legend = c('Mean', 'Max'), lty = 1, lwd = 2, col = c(1,4))
# legend('bottomright', legend = c('Mean', 'Max'), lty = 1, lwd = 2, cex = 0.5,
#        col = c(1,4), inset=c(0,1), horiz =TRUE, xpd = TRUE)

par(mfrow=c(1,1))
# plot(params$ID.order, params$time, type = 'l', lwd = 2, xaxt = 'n',
#      xlab = 'Weeks', ylab = 'Minutes', main = 'Computatuional time')
# axis(1, at = params$ID.order, labels = params$week)

# ------------------------------------------------------------------------------


## Weekly PM10 estimations and mean and max of logPM10 -------------------------
params = get_params(readRDS('../Influenza/par_mean.RData'))
head(params, 5)
params_max = get_params(readRDS('../Influenza/par_max.RData'))
head(params_max, 5)

plot(params$ID.order, params$beta_0.exp, type = 'b', lwd = 2, pch = 19, xaxt = 'n',
     ylim = c(0,0.8), col = 1, cex.axis = 1.5, cex.lab = 1.5,
     xlab = 'Weeks', ylab = 'exp(beta_0)', main = 'Estimation of exp(beta_0)')
axis(1, at = params$ID.order, labels = params$week)
legend('topright', legend = c('Mean', 'Max'), lty = 1, cex = 0.8,
       col = c(1,4))

plot(params$ID.order, params$beta_1.exp, type = 'b', lwd = 2, pch = 19, xaxt = 'n',
     ylim = c(0.3, 3.3), cex.axis = 1.5, cex.lab = 1.5,
     xlab = 'Weeks', ylab = 'exp(beta_1)', main = 'Estimation of exp(beta_1)')
abline(h = 1, lty = 2, col = 'red', lwd = 2)
lines(params$ID.order, params$beta_1.exp, type = 'b', lwd = 2, pch = 19,
      col = 1)
lines(params_max$ID.order, params_max$beta_1.exp, type = 'b', lwd = 2, pch = 19, 
      col = 4)
axis(1, at = params$ID.order, labels = params$week, cex.axis = 1.5, cex.lab = 1.5)
legend('topright', legend = c('Mean', 'Max'), lty = 1, lwd = 2,col = c(1,4))

# ------------------------------------------------------------------------------


## Weekly PM10 estimations and mean of PM10 ------------------------------------
#REM: not every week because no output
tmp = readRDS('../Influenza/par_mean_exp.RData')
params = as.data.frame(
  t(sapply(1:length(tmp), function(x){
    c(tmp[[x]]$week, tmp[[x]]$par, tmp[[x]]$time)
  }))
)
colnames(params) = c('week', 'beta_0', 'beta_1', 'theta_1', 'theta_2', 'time')
params = as.data.frame(apply(params, 2, as.numeric)) %>%
  dplyr::add_row(week = 45, .after = 1) %>%
  dplyr::add_row(week = c(51,52), .before = 8) %>%
  dplyr::add_row(week = 6:13, .after = 14)
params$ID.order = 1:nrow(params)
params$beta_0.exp = signif(exp(params$beta_0), 3)
params$beta_1.exp = signif(exp(params$beta_1), 3)
head(params, 10)

plot(params$ID.order, params$beta_0.exp, type = 'b', lwd = 2, pch = 19, xaxt = 'n',
     xlab = 'Weeks', ylab = 'exp(beta_0)', main = 'Estimation of exp(beta_0)')
axis(1, at = params$ID.order, labels = params$week)

plot(params$ID.order, params$beta_1.exp, type = 'b', lwd = 2, pch = 19, xaxt = 'n',
     xlab = 'Weeks', ylab = 'exp(beta_1)', main = 'Estimation of exp(beta_1)')
axis(1, at = params$ID.order, labels = params$week)

# ------------------------------------------------------------------------------


## NEW: No PM10 estimations and unifPop and with/without priors ----------------
param.df = get_params_noPM10(readRDS('SwedenData/par_noPM10_new_unifpop.RData'))
#RMK: when there's a prior, beta_0 is in the last position
param_pr.df = get_params_noPM10(readRDS('SwedenData/par_noPM10_new_unifpop_prior.RData'),
                                 prior_flag = TRUE)

plot(param.df$ID.order, param.df$beta_0.exp, type = 'b', lwd = 2, pch = 19,
     xaxt = 'n', ylim = c(0,1.3), col = 1,
     xlab = 'Weeks', ylab = 'exp(beta_0)',
     cex.axis = 1.5, cex.lab=1.5,
     main = 'Estimation of exp(beta_0) with uniform pop density')
lines(param_pr.df$ID.order, param_pr.df$beta_0.exp, type = 'b', lwd = 2, pch = 19,
      col = 6)
axis(1, at = param.df$ID.order, labels = param.df$week, cex.axis = 1.5)
legend('topright', legend = c('No prior', 'Prior'), lty = 1, 
       col = c(1,6), lwd = 2)
# ------------------------------------------------------------------------------

## OLD: No PM10 estimations and with/without priors or unifPop -----------------
params = get_params_noPM10(readRDS('../Influenza/par_noPM10.RData'))
params_pr = get_params_noPM10(readRDS('../Influenza/par_noPM10_priors.RData'),
                              prior_flag = TRUE)
params_prup = get_params_noPM10(readRDS('../Influenza/par_noPM10_priors_unifpop.RData'),
                                prior_flag = TRUE)

#REM: when using unif pop and no prior, week 45 (i_week 2) required too much time,
#     so it was forced to stop
tmp = readRDS('../Influenza/par_noPM10_unifpop.RData')
tmp[[1]]$par
for(x in 1:20){
  print(tmp[[x]]$week)
  print(tmp[[x]]$par)
  print(tmp[[x]]$time)
}

params_up = as.data.frame(
  t(sapply(c(1, 3:length(tmp)), function(x){
    c(tmp[[x]]$week, tmp[[x]]$par, tmp[[x]]$time)
  }))
)
colnames(params_up) = c('week', 'beta_0', 'theta_1', 'theta_2', 'time')
params_up = as.data.frame(apply(params_up, 2, as.numeric)) %>%
  dplyr::add_row(week = 45, .after = 1) %>%
  dplyr::add_row(week = c(51,52), .before = 8)
params_up$ID.order = 1:nrow(params_up)
params_up$beta_0.exp = signif(exp(params_up$beta_0), 3)
params_up$time[2] = 3
head(params_up, 10) #REM: small times are actually hours


# plot(params_up$ID.order, params_up$beta_0.exp, type = 'b', lwd = 2, pch = 19, col = 1,
#      xaxt = 'n', ylim = c(0,1.5),
#      xlab = 'Weeks', ylab = 'exp(beta_0)', main = 'Estimation of exp(beta_0)')
plot(param.df$ID.order, param.df$beta_0.exp, type = 'b', lwd = 2, pch = 19, col = 1,
     xaxt = 'n', ylim = c(0,1.5), cex.axis = 1.5, cex.lab = 1.5,
     xlab = 'Weeks', ylab = 'exp(beta_0)', main = 'Estimation of exp(beta_0)')
# lines(params_prup$ID.order, params_prup$beta_0.exp, type = 'b', lwd = 2, pch = 19, 
#       col = 7)
lines(param_pr.df$ID.order, param_pr.df$beta_0.exp, type = 'b', lwd = 2, pch = 19, 
      col = 7)
lines(params_pr$ID.order, params_pr$beta_0.exp, type = 'b', lwd = 2, pch = 19, 
      col = 4)
lines(params$ID.order, params$beta_0.exp, type = 'b', lwd = 2, pch = 19,
      col = 6)
axis(1, at = params$ID.order, labels = params$week, cex.axis = 1.5)
legend('topright', legend = c('Unif Pop', 'Unif Pop, Prior', 'Prior','Nothing'),
       lty = 1, cex = 0.8, lwd = 2, col = c(1,7,4,6))

plot(params_prup$ID.order, params_prup$beta_0.exp, type = 'b', lwd = 2, pch = 19, xaxt = 'n',
     col = 7, ylim = c(0, 1.5),
     xlab = 'Weeks', ylab = 'exp(beta_0)', main = 'Estimation of exp(beta_0)')
lines(params_up$ID.order, params_up$beta_0.exp, type = 'b', lwd = 2, pch = 19,
      col = 6)
axis(1, at = params_up$ID.order, labels = params_up$week)
legend('topright', legend = c('Unif Pop', 'Prior, Unif Pop'), lty = 1, cex = 0.8,
       col = c(6,7))

plot(params_prup$ID.order, params_prup$beta_0.exp, type = 'b', lwd = 2, pch = 19, xaxt = 'n',
     col = 7, ylim = c(0, 1.5),
     xlab = 'Weeks', ylab = 'exp(beta_0)', main = 'Estimation of exp(beta_0)')
lines(params_pr$ID.order, params_pr$beta_0.exp, type = 'b', lwd = 2, pch = 19,
      col = 4)
axis(1, at = params_up$ID.order, labels = params_up$week)
legend('topright', legend = c('Prior', 'Prior, Unif Pop'), lty = 1, cex = 0.8,
       col = c(4,7))

# ------------------------------------------------------------------------------


## Comparing exp(beta_0) when there is also PM10 [TO DO] -----------------------

# ------------------------------------------------------------------------------


## Getting hessian -------------------------------------------------------------
par_mean = readRDS('../Influenza/par_mean.RData')
par_max = readRDS('../Influenza/par_max.RData')

hes_mean = lapply(1:length(par_mean), function(x){
    list(week = par_mean[[x]]$week, hessian = signif(par_mean[[x]]$hessian,3))
  })
hes_max = lapply(1:length(par_max), function(x){
  list(week = par_max[[x]]$week, hessian = signif(par_max[[x]]$hessian,3))
})

# for(i in 1:20){ #2nd, 3rd obs and weeks 1, 6, 8, 12 
#   cat('Week',hes_mean[[i]]$week,'\nMean:\n')
#   print(signif(hes_mean[[i]]$hessian), 3)
#   cat('\nMax:\n')
#   print(signif(hes_max[[i]]$hessian), 3)
#   cat('\n----------------------------------------------------\n')
# }

# Weeks where Hessian matrices have NaN:
# - Mean: 44, 48, 49, 50, 1, 6, 11 (all)
# - Max: 46, 47, 49, 50, 1, 6, 9, 12

#RMK: sometimes matrices of the same week have values of different order (big/small)
#RMK: Weeks 2, 3, 4, and 10 have big values in both matrices, 
#     Weeks 48, 5, and 8 have big values only in Max
#     Week 6 have big values only in Mean
# Overall, Mean matrices seem to be more stable

for(i in c(8,13,15,19)){ #weeks 1, 6, 8, and 12, when difference in beta_1
  cat('Week',hes_mean[[i]]$week,'\nMean:\n')
  print(par_mean[i,7:8])
  cat('\n')
  print(hes_mean[[i]]$hessian)
  cat('\nMax:\n')
  print(par_max[i,7:8])
  cat('\n')
  print(hes_max[[i]]$hessian)
  cat('\n----------------------------------------------------\n')
}


# ------------------------------------------------------------------------------



##### Plotting PM10
## Comparing PM10_mean with true values ----------------------------------------
myweeks = c('44', '45', '46', '47', '48', '49', '50', 
            '1', '2', '3', '4', '5', '6', '7', 
            '8', '9', '10', '11', '12', '13')
pm10 = read.csv('SwedenData/data_1718.csv')[,-12] %>%
  dplyr::filter(Week %in% myweeks)
pm10$logPM10 = log(pm10$PM10)#column 14
pm10$WeeklylogPM10 = log(pm10$WeeklyPM10) #column 15
head(pm10)

# par_mean = readRDS('../Influenza/par_mean.RData')
par_mean = readRDS('SwedenData/all_par_max.RData')
head(par_mean[[1]]$mesh.df) #mean(logPM10) at col 15, max(logPM10) at col 17

logPM10.range = t(sapply(1:length(myweeks), function(i){
  c(myweeks[i],
    range(pm10[pm10$Week == myweeks[i],15], na.rm = TRUE),
    range(par_mean[[i]]$mesh.df$PM10_mean),
    as.numeric(sum(is.na(pm10[pm10$Week == myweeks[i],15])))/7)
}))
colnames(logPM10.range) = c('Week', 'TrueMin', 'TrueMax', 'EstMin', 'EstMax', 'HowManyNa')
logPM10.range = as.data.frame(apply(logPM10.range, 2, as.numeric))
logPM10.range = signif(logPM10.range, 3)
logPM10.range$ID.order = 1:nrow(logPM10.range)
head(logPM10.range)

range(logPM10.range$EstMin)
range(logPM10.range$EstMax)

## Plotting
plot(logPM10.range$ID.order, logPM10.range$TrueMin, type = 'l', lwd = 2, col = 1,
     ylim = c(-3.5,6), xaxt = 'n', #ylim = c(-1,5),
     xlab = 'Weeks', ylab = 'logPM10', main = 'Ranges of PM10_mean (log)')
lines(logPM10.range$ID.order, logPM10.range$TrueMax, type = 'l', lwd = 2, col = 1)
lines(logPM10.range$ID.order, logPM10.range$EstMin, type = 'l', lwd = 2, col = 4)
lines(logPM10.range$ID.order, logPM10.range$EstMax, type = 'l', lwd = 2, col = 4)
axis(1, at = logPM10.range$ID.order, labels = logPM10.range$Week)
legend('bottomright', legend = c('True', 'Estimated'),
       col = c(1,4), lty = 1, lwd = 2, cex = 0.8)

# ------------------------------------------------------------------------------


## Plotting PM10 ---------------------------------------------------------------
pm10_stat = read.csv('SwedenData/data_1718.csv')[1:44, c(1,4,5)]
# par_mean = readRDS('../Influenza/par_mean.RData')
par_mean = readRDS('SwedenData/all_par_max.RData')
border1 = read.csv('SwedenData/border1.csv')
border2 = read.csv('SwedenData/border3.csv')
colnames(border1) = colnames(border2) = c('X', 'Longitude', 'Latitude')

# PM10 levels
lower_lim = 0
upper_lim = 40
for(i in 1:length(par_mean)){
  p = dplyr::filter(par_mean[[i]]$mesh.df, !is.na(ID.Region)) %>%
    ggplot(aes(x = Longitude, y = Latitude)) +
    geom_tile(aes(fill = exp(PM10_mean))) + 
    scale_fill_continuous(limits = c(lower_lim, upper_lim)) + 
    scale_fill_gradient(low = 'green', high = 'red', limits = c(lower_lim, upper_lim)) +
    geom_polygon(data = border1, aes(x = Longitude, y = Latitude),
                 colour = 'black', fill = 'black', alpha = 0.2) +
    geom_polygon(data = border2, aes(x = Longitude, y = Latitude),
                 colour = 'black', fill = 'black', alpha = 0.2) +
    geom_point(data = pm10_stat, aes(x = Longitude, y = Latitude),
               colour = 'black', fill ='orange', pch = 21, size = 3) +
    labs(title = paste('Mean of PM10 (exp) in week', par_mean[[i]]$week))
  print(p)
  # breaks = seq(lower_lim, upper_lim, by = 10)
}

# logPM10 levels
tmp = dplyr::filter(par_mean[[i]]$mesh.df, !is.na(ID.Region),
                    Longitude > 12.5, Longitude < 16,
                    Latitude > 60, Latitude < 64) 
head(tmp)
range(tmp$PM10_mean)
range(tmp$PM10_max)

lower_lim = 0
upper_lim = 3.5
for(i in 1:length(par_mean)){
  p = dplyr::filter(par_mean[[i]]$mesh.df, !is.na(ID.Region)) %>%
    ggplot(aes(x = Longitude, y = Latitude)) +
    geom_tile(aes(fill = PM10_mean)) + 
    scale_fill_continuous(limits = c(lower_lim, upper_lim)) + 
    scale_fill_gradient(low = 'green', high = 'red', limits = c(lower_lim, upper_lim)) +
    geom_polygon(data = border1, aes(x = Longitude, y = Latitude),
                 colour = 'black', fill = 'black', alpha = 0.2) +
    geom_polygon(data = border2, aes(x = Longitude, y = Latitude),
                 colour = 'black', fill = 'black', alpha = 0.2) +
    geom_point(data = pm10_stat, aes(x = Longitude, y = Latitude),
               colour = 'black', fill ='orange', pch = 21, size = 3) +
    labs(title = paste('Mean of logPM10 in week', par_mean[[i]]$week))
  print(p)
}
# Weeks 45, 46, 48 and 1 (i.e. 2, 3, 5 and 8) have negative values
# Week 13 (i.e. 20) have much bigger values than usual

range(dplyr::filter(par_mean[[19]]$mesh.df, !is.na(ID.Region))$PM10_mean)
range(dplyr::filter(par_mean[[8]]$mesh.df, !is.na(ID.Region))$PM10_mean)


# library(animation)
# saveGIF({
#   for (i in 1:length(par_mean)){
#     p = ...
#     print(p)
#     }
#   }, interval = .2, movie.name='test.gif')


# ------------------------------------------------------------------------------
