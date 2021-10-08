##### Case study - Copernicus Land Cover Map (Po' Valley): spatial analysis and kriging #####
# written by Lucia Clarotto #

rm(list=ls())

## Set "path" to your local path and uncomment
# path = ''
# setwd(path)

library(ggplot2)
library(compositions) # only for ilr
library(Ternary)
library(lattice)
library(latex2exp)
library(raster)
library(rgdal)
library(ggmap)
library(RGeostats)
library(reshape2)
library(wesanderson)
library(Compositional) # only for helm

source("functions_simulation.R")

set.seed(123)

comp <- read.csv("data/comp_data_reduced_4_parts.csv")

## For Po' Valley
comp <- comp[which(comp$lng > 7.4 & comp$lng < 11.9 & comp$lat > 44.7 & comp$lat < 45.5), ]

## Grid
grid_lat_lng <- comp[samples,c("lng", "lat")]
grid_change <- grid_lat_lng
coordinates(grid_change)=~lng + lat
proj4string(grid_change)=CRS("+init=epsg:4326") # set it to lat-long
grid_change = spTransform(grid_change,CRS("+init=epsg:27700"))
grid <- grid_change@coords/1000

compositions <- comp[samples,c("crops", "grass", "shrub", "tree")]

plot(grid_lat_lng)
plot(compositions)

## Arithmetic mean of components
apply(compositions, 2, mean)
cols <- wes_palette(n=5, name="FantasticFox1")[c(1,2,3,5)]

## Maps
myLocation <- c(min(grid_lat_lng[,1])-0.3, min(grid_lat_lng[,2])-0.2, max(grid_lat_lng[,1])+0.3, max(grid_lat_lng[,2])+0.2)
map_physics <- get_stamenmap(myLocation, zoom = 8, maptype = "terrain-background")
rivers <- readOGR(dsn=path.expand("data/wise_large_rivers"), layer="Large_rivers")
lat_lng_rivers <- spTransform(rivers, CRS("+init=epsg:4326"))
po <- fortify(lat_lng_rivers@lines[[18]])

map <- ggmap(map_physics) +
  geom_line(data = po, aes(long, lat), color="black") 
map

## Complete dataset
dataset <- data.frame(cbind(grid_lat_lng, compositions))
colnames(dataset) <- c("x", "y", "crops", "grass", "shrub", "tree")

## Samples with some zeros
index0z <- rowSums(dataset[,3:6] == 0) == 0; sum(index0z)
index1z <- rowSums(dataset[,3:6] == 0) == 1; sum(index1z)
index2z <- rowSums(dataset[,3:6] == 0) == 2; sum(index2z)
index3z <- rowSums(dataset[,3:6] == 0) == 3; sum(index3z)

## Dataset without zeros
# no_zero <- which(compositions$crops !=0 & compositions$grass !=0 & compositions$shrub !=0 & compositions$tree !=0 )
# dataset <- dataset[no_zero,]

## Plots
plot_crops <- ggmap(map_physics) +
  geom_point(data = dataset, aes(x, y, colour=crops)) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  geom_line(data = po, aes(long, lat), color="black")

plot_grass <- ggmap(map_physics) +
  geom_point(data = dataset, aes(x, y, colour=grass)) +
  scale_color_gradient(low = "yellow", high = "green", na.value = NA) +
  geom_line(data = po, aes(long, lat), color="black")

plot_shrub <- ggmap(map_physics) +
  geom_point(data = dataset, aes(x, y, colour=shrub)) +
  scale_color_gradient(low = "yellow", high = "blue", na.value = NA) +
  geom_line(data = po, aes(long, lat), color="black")

plot_tree <- ggmap(map_physics) +
  geom_point(data = dataset, aes(x, y, colour=tree)) +
  scale_color_gradient(low = "yellow", high = "black", na.value = NA) +
  geom_line(data = po, aes(long, lat), color="black")

plot_crops
plot_grass
plot_shrub
plot_tree

## Plot of components
df <- dataset[,3:6]
df <- df[order(df$tree),]
df$event <- seq.int(nrow(df))   #create a column to indicate which values happened on the same column for each variable
df <- melt(df, id='event')      #reshape dataframe to make it readable to ggplot2

percentage_plot <- ggplot(df, aes(x = event, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  ggplot2::theme_bw() +
  theme(axis.text = element_text(size=22),
        axis.title = element_text(size=22),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15),
        plot.title = element_text(hjust = 0.5, size = 22)) +
  labs(title = TeX('original data'),x=TeX('$x_i$'), y=TeX('$\\%$')) +
  scale_fill_manual(values = cols)
percentage_plot

## Parameters
N_try <- 5
N_train <- 500 
alpha_krig=c(0.01, 0.12,0.3, 0.5, 0.75,1) ## for dataset with zeros
# alpha_krig=c(0, 0.12, 0.3, 0.5, 0.75,1)   ## for dataset without zeros
N <- length(alpha_krig)

yy <- dataset[,3:6]
coord <-dataset[,1:2]

points <- list()
krig_data_total <- list()
krig_data_euc_total <- list()
plot_krig_data_total <- list()
plot_krig_data_euc_total <- list()

H=helm(4)
trans <- list()
for (i in 1:N) {
  if (alpha_krig[i] == 0) {
    trans[[i]] <- ilr(as.matrix(yy),  V=t(H))
  } else {
    trans[[i]] <- data.frame(alpha.IT(as.matrix(yy), alpha_krig[i]))
  }
}
plot(trans[[1]])
plot(trans[[6]])

## Variography
number = 500
ncell = 100
nlag  <- 15
ndir <- 1
ang0  <- 0

model  = model.create(2,range=150,sill=0.5)
xloc1  = as.vector(coord[,1])
xloc2  = as.vector(coord[,2])

test.neigh <- neigh.create(ndim = 2, type = 0)
sample_vars <- list()
var_models <- list()
for (i in 1:N) {
  db.data = db.create(x1=xloc1,x2=xloc2)
  db.data = db.add(db.data, z1=trans[[i]][,1], z2=trans[[i]][,2], z3=trans[[i]][,3])

  # Computes experimental variograms 
  test.vario <- vario.calc(db.data,nlag=nlag)
  plot(test.vario)
  sample_vars[[i]] <- test.vario
  # Automatic Model Fitting
  test.model <- model.auto(vario=test.vario, struct=c("Nugget Effect","Exponential"), draw=TRUE)
  print(test.model)
  var_models[[i]] <- test.model
}

## Prediction (KRIGING)
for (j in 1:N_try) {
  ## Select some data from the simulated data
  pts <- sample(nrow(yy), N_train)
  
  points[[j]] <- pts
  
  ## Dataset with selected points and coordinates
  y <- yy[pts,]
  
  krig_data <- list()
  krig_data_euc <- list()
  data_final <- list()
  error_alfa_sim <- NULL
  
  for (i in 1:N) {
    
    ## Create spatial database
    db.data_sample <- db.create(x1=coord[pts,1],x2=coord[pts,2],z1=trans[[i]][pts,1], z2=trans[[i]][pts,2], z3=trans[[i]][pts,3])
    db.grid <- db.create(x1=coord[,1],x2=coord[,2])
    
    ## Apply kriging
    db.krig <- kriging(db.data_sample,db.grid, var_models[[i]],test.neigh)
    print(db.krig)
    
    new_pred <- data.frame("Z1"=db.krig$Kriging.z1.estim, "Z2"=db.krig$Kriging.z2.estim, "Z3"=db.krig$Kriging.z3.estim)
    krig_data_euc[[i]] <- new_pred
    
    ## Back-transform with inverse ilr for alpha=0 and with inverse alpha-IT in all the other cases
    if (alpha_krig[i] == 0) {
      y_new <- data.frame(ilrInv(new_pred, V=t(H)))
    }else {
      y_new <- data.frame(alpha.IT.inv(new_pred,alpha_krig[i]))
    }
    
    y_new[y_new < 10^(-5)] <- 0 # rounding small values to 0
    krig_data[[i]] <- y_new
    
  }
  
  krig_data_euc_total[[j]] <- krig_data_euc
  krig_data_total[[j]] <- krig_data
  
  print(j)
}

## Save results
results <- list(comp_data = yy,
                trans = trans,
                grid = grid,
                var_models = var_models,
                lat_long = grid_lat_lng,
                points = points, 
                krig_data_total = krig_data_total, 
                krig_data_euc_total = krig_data_euc_total,
                alpha_krig = alpha_krig,
                N_try = N_try)

## for dataset with zeros
name = paste("data/copernicus_zeros_cokriging_N_try",toString(N_try),"N_train",toString(N_train),"krig",sep="_")
## for dataset without zeros
# name = paste("data/copernicus_cokriging_N_try",toString(N_try),"N_train",toString(N_train),"krig",sep="_") 
for (i in 1:N) {
  name = paste(name, toString(alpha_krig[i]), sep="_")
}
name = paste(name, ".rdata", sep="")
name
save(results, file=name)
