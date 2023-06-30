## Get pollution in Swedish during season 2017-2018
folder = ### set working directory
setwd(folder)

#######################################
## Computing daily PM10              ##
#######################################

#first day
result = get_pollution(1*i_week)
output = result$output
stack.pm10 = result$stack.pm10

beta_j = round(output$summary.fixed, 3)
beta_j$i_day = 1*i_week
beta_table = beta_j

index.pred = inla.stack.index(stack.pm10, 'pred')$data
pm10.days.pred = data.frame(Longitude = grid.pred.pm10$Longitude,
                            Latitude = grid.pred.pm10$Latitude,
                            PM10_1 = output$summary.linear.predictor[index.pred,'mean'])

if(plot_flag){
  tmp = pm10.days.pred
  colnames(tmp)[colnames(tmp) == 'PM10_1'] = 'logPM10'
  tmp$`St. Dev.` = output$summary.linear.predictor[index.pred,'sd']
  tmp = left_join(x = mesh.df, y = tmp) %>%
    dplyr::filter(!is.na(ID.Region))
  # saveRDS(tmp, 'SwedenData/pm10_iweek15-day1.RData')
  # tmp = readRDS('SwedenData/pm10_iweek15-day1.RData')
  
  p_mean = ggplot(tmp, aes(Longitude, Latitude, fill = logPM10)) +
    geom_tile() +
    scale_fill_gradient(low = 'green2', high = 'red2') +
    geom_polygon(data = border1, aes(x = Longitude, y = Latitude),
                 colour = 'black', fill = 'black', alpha = 0.2) +
    geom_polygon(data = border2, aes(x = Longitude, y = Latitude),
                 colour = 'black', fill = 'black', alpha = 0.2) +
    geom_point(data = pm10[1:44, c(1,4,5)], aes(x = Longitude, y = Latitude),
               colour = 'black', fill ='orange', pch = 21, size = 3) +
    coord_fixed(ratio = 1) +
    theme(legend.title = element_text(size = 12)) +
    geom_label(aes(x = 14.1, y = 64.5, label = "Bredkälen"), 
               hjust = 0, 
               vjust = 0.5, 
               colour = "#000000", 
               fill = "white", 
               label.size = NA, 
               family="Helvetica", 
               size = 6) +
    scale_x_continuous(name = 'Longitude', breaks = seq(10, 25, 5),
                       labels = c('10°E','15°E','20°E','25°E')) + 
    scale_y_continuous(name = 'Latitude', breaks = seq(55, 70, 5),
                       labels = c('55°N','60°N','65°N','70°N'))
  #+ labs(title = 'Posterior mean of the spatial latent field')
  
  p_sd = ggplot(tmp, aes(Longitude, Latitude, fill = `St. Dev.`)) +
    geom_tile() +
    scale_fill_gradient(low = 'white', high = 'blue') +
    geom_polygon(data = border1, aes(x = Longitude, y = Latitude),
                 colour = 'black', fill = 'black', alpha = 0.2) +
    geom_polygon(data = border2, aes(x = Longitude, y = Latitude),
                 colour = 'black', fill = 'black', alpha = 0.2) +
    geom_point(data = pm10[1:44, c(1,4,5)], aes(x = Longitude, y = Latitude),
               colour = 'black', fill ='orange', pch = 21, size = 3) +
    coord_fixed(ratio = 1) +
    theme(legend.title = element_text(size = 12)) +
    geom_label(aes(x = 14.1, y = 64.5, label = "Bredkälen"), 
               hjust = 0, 
               vjust = 0.5, 
               colour = "#000000", 
               fill = "white", 
               label.size = NA, 
               family="Helvetica", 
               size = 6) +
    scale_x_continuous(name = 'Longitude', breaks = seq(10, 25, 5),
                       labels = c('10°E','15°E','20°E','25°E')) + 
    scale_y_continuous(name = 'Latitude', breaks = seq(55, 70, 5),
                       labels = c('55°N','60°N','65°N','70°N'))
  # + labs(title = 'Standard deviation of the spatial latent field')
  
  grid.arrange(p_mean, p_sd, ncol = 2)
}

#second day
result = get_pollution(2*i_week)
output = result$output
stack.pm10 = result$stack.pm10
beta_j = round(output$summary.fixed, 3)
beta_j$i_day = 2*i_week
beta_table = rbind(beta_table, beta_j)

index.pred = inla.stack.index(stack.pm10, 'pred')$data
pm10.days.pred$PM10_2 = output$summary.linear.predictor[index.pred,'mean']

#third day
result = get_pollution(3*i_week)
output = result$output
stack.pm10 = result$stack.pm10
beta_j = round(output$summary.fixed, 3)
beta_j$i_day = 3*i_week
beta_table = rbind(beta_table, beta_j)

index.pred = inla.stack.index(stack.pm10, 'pred')$data
pm10.days.pred$PM10_3 = output$summary.linear.predictor[index.pred,'mean']

#fourth day
result = get_pollution(4*i_week)
output = result$output
stack.pm10 = result$stack.pm10
beta_j = round(output$summary.fixed, 3)
beta_j$i_day = 4*i_week
beta_table = rbind(beta_table, beta_j)

index.pred = inla.stack.index(stack.pm10, 'pred')$data
pm10.days.pred$PM10_4 = output$summary.linear.predictor[index.pred,'mean']

#fifth day
result = get_pollution(5*i_week)
output = result$output
stack.pm10 = result$stack.pm10
beta_j = round(output$summary.fixed, 3)
beta_j$i_day = 5*i_week
beta_table = rbind(beta_table, beta_j)

index.pred = inla.stack.index(stack.pm10, 'pred')$data
pm10.days.pred$PM10_5 = output$summary.linear.predictor[index.pred,'mean']

#sixth day
result = get_pollution(6*i_week)
output = result$output
stack.pm10 = result$stack.pm10
beta_j = round(output$summary.fixed, 3)
beta_j$i_day = 6*i_week
beta_table = rbind(beta_table, beta_j)

index.pred = inla.stack.index(stack.pm10, 'pred')$data
pm10.days.pred$PM10_6 = output$summary.linear.predictor[index.pred,'mean']

#seventh day
result = get_pollution(7*i_week)
output = result$output
stack.pm10 = result$stack.pm10
beta_j = round(output$summary.fixed, 3)
beta_j$i_day = 7*i_week
beta_table = rbind(beta_table, beta_j)

index.pred = inla.stack.index(stack.pm10, 'pred')$data
pm10.days.pred$PM10_7 = output$summary.linear.predictor[index.pred,'mean']

# write.csv(pm10.days.pred, 'SwedenInfluenzaData/PM10_week7.csv')


#######################################
## Joining daily PM10 to mesh.df     ##
#######################################
mesh.df.noPM10 = mesh.df
mesh.df = left_join(x = mesh.df, y = pm10.days.pred)
head(mesh.df, 10)

sum(is.na(mesh.df$PM10_mean))  == 0
nrow(mesh.df[!is.na(mesh.df$PM10_mean),]) == nrow(mesh.df)
cat('\nGetting pollution of week', week, ': done!\n')

# write.csv(mesh.df, 'SwedenInfluenzaData/mesh.df_week7.csv')
