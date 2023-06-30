## Get pollution in Swedish during season 2017-2018
folder = ### set working directory
setwd(folder)

###################################
## Estimating pollution          ##
###################################
# pm10.df = pm10[pm10$Week == week,]
pm10.df = pm10
stations = unique(pm10.df$SamplingPoint) # 44 and 36 stations under max_lat
dates = unique(pm10.df$Date) # 7 days if weekly based, else 224 days
coord_all = pm10.df[,c(1,4,5)]
coord = coord_all[1:length(stations),]

pm10.df$Time = rep(1:length(dates), each = length(stations)) 
pm10.df$logPM10 = log(pm10.df$PM10)
head(pm10.df)
pm10.df_not_scaled = pm10.df

## Standardizing covariates
mean_cov = apply(pm10.df[,4:11], 2, mean)
sd_cov = apply(pm10.df[,4:11], 2, sd)
pm10.df[,4:11] = scale(pm10.df[,4:11], center = mean_cov, scale = sd_cov)
head(pm10.df)

mesh.pm10 = inla.mesh.2d(
  loc = coord[,2:3],
  loc.domain = cbind(c(border1$Longitude, border2$Longitude), c(border1$Latitude, border2$Latitude)),
  max.edge = c(4, 20),
  offset = c(6, 10),
  cutoff = 0.03)
mesh.pm10$n # 228 vertices in the triangulation

if(plot_flag){
  plot(mesh.pm10, main ='Triangulation of South of Sweden')
  lines(border1$Longitude, border1$Latitude, lwd = 3, col = 'black')
  lines(border2$Longitude, border2$Latitude, lwd = 3, col = 'black')
  points(coord[,2:3], pch = 21, bg = 'orange', col = 'black')
}

## Projector matrix A 
A.est.pm10 = inla.spde.make.A(mesh = mesh.pm10,
                              loc = as.matrix(coord_all[,2:3]),
                              group = pm10.df$Time, #time grouping = n_station per day
                              n.group = length(dates)) #number of groups = n_days

## Checking
# table(apply(A.est.pm10, 1 , nnzero)) # 3 non-zero values per each of the 9865 rows
# dim(A.est.pm10)
# nrow(A.est.pm10) == nrow(pm10.df)
# ncol(A.est.pm10) == mesh.pm10$n*length(dates)

## Matern SPDE 
spde.pm10 = inla.spde2.matern(mesh = mesh.pm10, alpha = 2)
s_index = inla.spde.make.index(name = 'spatial.field',
                               n.spde = spde.pm10$n.spde,
                               n.group = length(dates))

## Creating the stack object for estimating the spatio-temporal model
stack.est.pm10 = inla.stack(data = list(logPM10 = pm10.df$logPM10),
                            A = list(A.est.pm10, 1),
                            effects = list(c(s_index, list(Intercept = 1)),
                                           list(pm10.df[,4:11])),
                            tag = 'est') # my tag for calling it later

cat('\nEstimating logPM10: done!\n')


############################################################
## Function for Predicting pollution at the mesh.df level ##
############################################################
grid.pred.pm10 = mesh.df[,1:2]

get_pollution = function(i_day){
  A.pred.pm10 = inla.spde.make.A(mesh = mesh.pm10,
                                 loc = as.matrix(grid.pred.pm10),
                                 group = i_day,
                                 n.group = length(dates))
  # table(apply(A.pred.pm10, 1 , nnzero)) # 3 non-zero values per each of the 2264 rows
  
  pm10.df.day = pm10.df_not_scaled[pm10.df_not_scaled$Time == i_day, ]
  cov.interp = data.frame(Longitude = grid.pred.pm10$Longitude, Latitude = grid.pred.pm10$Latitude)
  
  # Inverse distance weighted interpolation:
  cov.interp$Altitude = idw(pm10.df.day$Altitude ~ 1,
                            SpatialPoints(pm10.df.day[,4:5]),
                            SpatialPoints(grid.pred.pm10),
                            idp = 2)$var1.pred
  cov.interp$Temperature = idw(pm10.df.day$Temperature ~ 1,
                               SpatialPoints(pm10.df.day[,4:5]),
                               SpatialPoints(grid.pred.pm10),
                               idp = 2)$var1.pred
  cov.interp$AirPressure = idw(pm10.df.day$AirPressure ~ 1,
                               SpatialPoints(pm10.df.day[,4:5]),
                               SpatialPoints(grid.pred.pm10),
                               idp = 2)$var1.pred
  cov.interp$Humidity = idw(pm10.df.day$Humidity ~ 1,
                            SpatialPoints(pm10.df.day[,4:5]),
                            SpatialPoints(grid.pred.pm10),
                            idp = 2)$var1.pred
  cov.interp$Precipitation = idw(pm10.df.day$Precipitation ~ 1,
                                 SpatialPoints(pm10.df.day[,4:5]),
                                 SpatialPoints(grid.pred.pm10),
                                 idp = 2)$var1.pred
  cov.interp$WindSpeed = idw(pm10.df.day$WindSpeed ~ 1,
                             SpatialPoints(pm10.df.day[,4:5]),
                             SpatialPoints(grid.pred.pm10),
                             idp = 2)$var1.pred
  
  cov.interp_not_scaled = cov.interp
  
  ## Plotting idw (NOTE: we cannot plot points outside Sweden)
  if(plot_flag){
    ggplot(cov.interp_not_scaled, aes(x = Longitude, y = Latitude, fill = Altitude)) +
      geom_tile() +
      scale_fill_gradient(low = 'white', high = 'blue') +
      geom_point(data = pm10.df_not_scaled[1:length(stations),], aes(x = Longitude, y = Latitude),
                 colour = 'black', fill ='orange', pch = 21, size = 2) +
      geom_polygon(data = border1, aes(x = Longitude, y = Latitude),
                   colour = 'black', fill = 'black', alpha = 0.2) +
      geom_polygon(data = border2, aes(x = Longitude, y = Latitude),
                   colour = 'black', fill = 'black', alpha = 0.2) +
      labs(title = 'IDW interpolation of Altitude')
  }
  
  ## Standardizing covariates, i.e. Longitude, Latitude, Altitude, Temperature,
  # AirPressure, Humidity, Precipitation, WindSpeed
  mean_cov = apply(cov.interp, 2, mean)
  sd_cov = apply(cov.interp, 2, sd)
  cov.interp = scale(cov.interp, center = mean_cov, scale = sd_cov) %>%
    as.data.frame()
  # head(cov.interp, 10)
  # head(pm10.df[4:11])
  # range(cov.interp$Latitude)
  # range(pm10.df$Latitude)
  
  # Actually predicting
  stack.pred.pm10 = inla.stack(data = list(logPM10 = NA),
                               A = list(A.pred.pm10, 1),
                               effects = list(c(s_index, list(Intercept = 1)),
                                              list(cov.interp)),
                               tag = 'pred')
  
  stack.pm10 = inla.stack(stack.est.pm10, stack.pred.pm10)
  
  formula = logPM10 ~ -1 + Intercept + #intercept
    Longitude + Latitude + Altitude + Temperature + AirPressure + Humidity + Precipitation + WindSpeed + #predictors
    f(spatial.field, model = spde.pm10,
      group = spatial.field.group, control.group = list(model = 'ar1'))
  
  output = inla(formula,
                data = inla.stack.data(stack.pm10, spde = spde.pm10),
                family = 'gaussian',
                control.predictor = list(A = inla.stack.A(stack.pm10), compute = TRUE),
                control.compute = list(return.marginals.predictor = TRUE))
  
  cat('PM10 in day',i_day,'of week',i_week,': done!\n\n')
  return(list(cov.interp_not_scaled = cov.interp_not_scaled, 
              stack.pm10 = stack.pm10,
              output = output))
}
