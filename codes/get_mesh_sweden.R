## Setup of Swedish during season 2017-2018
folder = ### set working directory
setwd(folder)

###################################
## Setting                       ##
###################################
flu$Laboratory[flu$Laboratory == 'Karlskrona1'] = 'Karlskrona'
flu$Laboratory[flu$Laboratory == 'Växjö1'] = 'Växjö'
flu$Laboratory[flu$Laboratory == 'Skåne'] = 'Malmö'
flu$Laboratory[flu$Laboratory == 'Eskilstuna/Unilabs'] = 'Eskilstuna'
flu$Laboratory[flu$Laboratory == 'Skövde/Unilabs'] = 'Skövde'
flu$Laboratory[flu$Laboratory == 'Stockholm/Unilabs'] = 'Stockholm'
flu$Laboratory[flu$Laboratory == 'Göteborg'] = 'Gothenburg'

flu$County[flu$County == 'Jönköpings län'] = 'Jönköping'
flu$County[flu$County == 'Kalmar län'] = 'Kalmar'
flu$County[flu$County == 'Örebro län'] = 'Orebro'
flu$County[flu$County == 'Region Skåne'] = 'Skåne'
flu$County[flu$County == 'Sörmland'] = 'Södermanland'
flu$County[flu$County == 'Uppsala län'] = 'Uppsala'


# Adding population density per county 
# WikiPopCounts = c(2377081, 1725881, 1377827, 465495, 383713, 363599, 333848, 
#                   304805, 29754, 287966, 287382, 282414, 275845, 271736, 250093,
#                   245446, 245347, 201469, 159606, 13081,59686)
# WikiArea = c(6519.3, 23948.8, 11034.5, 10602.0, 8207.2, 10495.1, 5460.7, 8545.6, 
#              6102.3, 28188.8, 18198.9, 17591.0, 5145.8, 55186.2, 98244.8, 11217.8,
#              21683.8, 8466.0, 2946.4, 49341.2, 3151.4)
# WikiDens = c(360, 68, 120, 43, 45, 34, 59, 35, 52, 10, 16, 16, 37, 4.9, 2.6, 22,
#              11, 23, 52, 3.4, 19)

pop.county = data.frame(
  Region = c('Stockholm', 'Västra Götaland', 'Skåne', 'Östergötland', 'Uppsala',
             'Jönköping', 'Halland', 'Orebro', 'Södermanland', 'Dalarna', 'Gävleborg',
             'Värmland', 'Västmanland', 'Västerbotten', 'Norrbotten', 'Kalmar', 
             'Västernorrland', 'Kronoberg', 'Blekinge', 'Jämtland', 'Gotland'),
  PopCounts = c(2415139, 1744859, 1402425, 469704, 395026, 367064, 340243, 306792, 
                301801, 288387, 287767, 283196, 278967, 274563, 249693, 247175, 
                244193, 203340, 158937, 132054, 61001))

###################################
## Creating Swedish flu data     ##
###################################
flu.df = flu %>%
  dplyr::filter(County != 'WeeklyTotal') %>%
  mutate(Counts = Current.week.cases.A + Current.week.cases.B) %>%
  dplyr::select(Region = County, Laboratory, Counts, Samples, Week) %>%
  left_join(y = cities, by = c('Laboratory'='city')) %>%
  dplyr::select(Region, Laboratory, Counts, Samples, Latitude=lat, Longitude=lng, Week)
flu.df[flu.df$Region == 'Stockholm',5] = 59.3294
flu.df[flu.df$Region == 'Stockholm',6] = 18.0686
flu.df$ID.Order = 1:nrow(flu.df)
head(flu.df,10)

flu.df.grouped = flu.df %>%
  group_by(Region, Week) %>%
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude),
            Counts = sum(Counts), Samples = sum(Samples), ID.Order = mean(ID.Order)) %>%
  ungroup() %>%
  as.data.frame() %>%
  arrange(ID.Order)
flu.df.grouped$ID.Order = 1:nrow(flu.df.grouped)
head(flu.df.grouped, 21)

flu.df.grouped = dplyr::left_join(x = flu.df.grouped, y = swe.shp@data, by = c('Region'='NAME_1')) %>%
  dplyr::select(Region, ID.Region = ID_1, Longitude, Latitude, Week, Counts, Samples, ID.Order)
flu.df.grouped$ID.Region = as.numeric(flu.df.grouped$ID.Region)
flu.df.grouped = dplyr::left_join(x = flu.df.grouped, y = pop.county, by = c('Region'='Region'))
head(flu.df.grouped, 21)


flu.sp.grouped = swe.shp # flu.sp.grouped = readOGR('../Sweden_shp/SWE_adm1.shp')
flu.sp.grouped@data = flu.df.grouped

if(plot_flag){
  plot(flu.sp.grouped@data$Longitude, flu.sp.grouped@data$Latitude, pch = 21, bg = 'orange',
       main = 'Influenza labs of Sweden', xlab = 'Longitude', ylab = 'Latitude',
       xlim = c(10,25), ylim = c(55,70),)
  lines(border1[,2:3], lwd = 3, col = 'black')
  lines(border2[,2:3], lwd = 3, col = 'black')
}
cat('\nGetting regional flu counts: done!\n')

###################################
## Set-up variables for analysis ##
###################################
# - create mesh: "mesh.true"
# - create D matrix: "D"
#     - note: this is a 25 x m matrix
# - TMB objects

pop.dat = setMinMax(pop.dat)
swe.extent = extent(10.95708, 24.16542, 55.33292, 69.06625)
pop.dat.swe = crop(pop.dat, swe.extent)
pop.dat.swe = projectRaster(pop.dat.swe, crs = CRS("+proj=longlat +datum=WGS84 +no_defs")) 

if(plot_flag){
  plot(pop.dat.swe, main = 'Population density in Sweden')
}


## Creating Mesh/grid points from Sweden shp
set.seed(2311)
grid.sweden = spsample(swe.shp, 2000, 'regular')
dim(grid.sweden@coords) #2001 grid points

if(plot_flag){
  plot(grid.sweden, main = 'Grid of Sweden with influenza labs')
  lines(border1[,2:3], lwd = 3, col = 'black')
  lines(border2[,2:3], lwd = 3, col = 'black')
  points(flu.df.grouped$Longitude, flu.df.grouped$Latitude, axes = TRUE, pch = 21, bg = 'orange')
}

pop.at.grid = data.frame(grid.sweden@coords)#/1000 #coordinates of grid points
pop.at.gridX = raster::extract(pop.dat.swe, pop.at.grid) #pop at each grid point
pop.at.gridX[is.na(pop.at.gridX)] = 0
pop.at.grid$pop = pop.at.gridX #pop at each grid point
colnames(pop.at.grid) = c('Longitude', 'Latitude', 'PopDensity')
dim(pop.at.grid)

## Creating mesh of triangulation
# mesh.true = inla.mesh.2d(pop.at.grid[, 1:2],
#                          offset = c(1,2), 
#                           max.edge = 50*c(1,5))
mesh.true = inla.mesh.2d(
  loc = grid.sweden@coords, # loc = pop.at.grid[,1:2],
  loc.domain = cbind(c(border1$Longitude, border2$Longitude), c(border1$Latitude, border2$Latitude)),
  max.edge = c(25, 100),
  offset = c(1.5, 5),
  cutoff = 0.15)
mesh.true$n #2264 mesh points

# write.csv(data.frame(Longitude = mesh.true$loc[,1], Latitude = mesh.true$loc[,2]),
#           'SwedenInfluenzaData/mesh.true_coords.csv')

#NOTE: mesh.true should have more points than in the grid [to be coherent with them]
#      In this way, we can join with no NA in pop.at.grid.new

if(plot_flag){
  plot(mesh.true, main = 'Triangulation of Sweden with influenza labs')
  lines(border1[,2:3], lwd = 3, col = 'black')
  lines(border2[,2:3], lwd = 3, col = 'black')
  points(flu.df.grouped$Longitude, flu.df.grouped$Latitude, axes = TRUE, pch = 21, bg = 'orange')
}

# Mesh points of the triangulation  
#NOTE: we match them according to the Region, because we have political boundaries
#for regions only.
#--> we will incorporate the 4 labs in Västra Götaland with the distances trick 
meshgrid = data.frame(Longitude = mesh.true$loc[,1],
                      Latitude = mesh.true$loc[,2])
coordinates(meshgrid) = ~ Longitude + Latitude
proj4string(meshgrid) = proj4string(flu.sp.grouped)
meshlocs = over(meshgrid, flu.sp.grouped)[,2] #ID.Region
#it matches each triangulation point with the respective Region
#-->gives me NA in the triangulation points not belonging to regions
stopifnot(min(meshlocs, na.rm = TRUE) == 1)
stopifnot(max(meshlocs, na.rm = TRUE) == 21)
mesh.df = data.frame(Longitude = mesh.true$loc[,1], Latitude = mesh.true$loc[,2],
                     ID.Region = meshlocs)
mesh.df$id = 1:nrow(mesh.df)

## Creating matrix D with population density in each mesh point
pop.at.grid.new = dplyr::left_join(x = pop.at.grid, y = mesh.df)
dim(pop.at.grid.new)[1] == dim(pop.at.grid)[1]
head(pop.at.grid.new, 10) #it has also ID.Region and id
#Being the grid of Sweden, each point should be associated to an id [hence, no NA]
sum(is.na(pop.at.grid.new$ID.Region)) == 0

if(unifPopDensity_flag){
  # Uniform population density
  cat('\nUsing uniform population density\n')
  pop.at.grid.new$PopDensityUnif = NA
  for(i_region in 1:21){
    tmp = pop.at.grid.new[pop.at.grid.new$ID.Region == i_region,]
    # pop.at.grid.new[pop.at.grid.new$ID.Region == i_region,6] = sum(tmp$PopDensity)/nrow(tmp)
    pop.at.grid.new[pop.at.grid.new$ID.Region == i_region,6] = 1/nrow(tmp)
  }
  D = inla.spde.make.A(mesh = mesh.true,
                       loc = as.matrix(pop.at.grid.new[,1:2]),
                       block = pop.at.grid.new$ID.Region,
                       weights = pop.at.grid.new$PopDensityUnif)
}else{
  D = inla.spde.make.A(mesh = mesh.true,
                       loc = as.matrix(pop.at.grid.new[,1:2]),
                       block = pop.at.grid.new$ID.Region,
                       weights = pop.at.grid.new$PopDensity)
}
#NOTE: here we will need to change the block
#Now blocks are wrt to the region

#Cheking
# dim(D) #21 regions and 2264 mesh points in total
# table(apply(D, 1, nnzero)) #1 row with 11 nnzero, 1 row with 17 nnzeros,...
# table(apply(D, 1, sum))
# table(apply(D, 2, sum)) #NOTE: some columns are negative and very very small --> zero?

D.tmp = list()
D.tmp$D = D/rowSums(D)
D.tmp$mesh.weights = colSums(D)
mesh.df$weight.unscaled = D.tmp$mesh.weights
D = D.tmp$D
cat('\nGetting population density matrix D: done!\n')


## Adding mesh weights relative to population density
mesh.df$weight.scaled = colSums(D)
nmesh.area = 1 / tapply(rep(1, nrow(mesh.df)), mesh.df$ID.Region, sum)
mesh.df$weight.comp = nmesh.area[mesh.df$ID.Region + 1]  # comparison weight
#REM: NAs in ID.Region when mesh point outside Sweden
#     weight.comp is the uniform areal weight

# write.csv(mesh.df, 'SwedenInfluenzaData/triangulation_meshdf.csv')
