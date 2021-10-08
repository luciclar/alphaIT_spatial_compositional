##### Simulation study: spatial analysis and kriging #####
# written by Lucia Clarotto #

rm(list=ls())

## Set "path" to your local path and uncomment
# path = ''
# setwd(path)

## Libraries
library(RandomFields)
library(Ternary)
library(lattice)
library(wesanderson)
library(latex2exp)
library(RGeostats)
library(ggplot2)
library(compositions) # only for ilr
library(viridis)
library(Compositional) # only for helm

source("functions_simulation.R")

## Parameters
alpha_krig = c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1) # values of alpha for kriging
alpha_0 = 0 # for the inverse alpha-IT in the simulation
N <- length(alpha_krig)

## Dimension of the simplex
D = 3
H = helm(D)

## N_samples is the number of observations
N_samples = 500
## N_try is the number of simulations MC
N_try = 10 # change to the number of MC simualations that you want to do (we use with N_try=100 in the paper)


####### START 1 ######
## Run this section if no simulation has been run yet ##

## Choice of the model for the bivariate random field
model <- RMbiwm(nudiag=c(1, 1), c=c(1,0.8,1)) # Whittle-Matérn model

x <- y <- seq(0, 15, length.out=200)
try <- RFsimulate(model, x, y)
plot(try@data)

## Run N_try realizations of the Matern 
N_locations <- 2000
locations <- sample(nrow(try), N_locations)

grid <- coordinates(try)[locations,]
plot(grid)
 
z <- list()
tot <- NULL
for (i in 1:N_try) {
  total <- RFsimulate(model, x, y)
  z[[i]] <- total@data[locations,]
  tot <- rbind(tot, z[[i]])
}

name_real = "data/realizations"
save(z, file=name_real)
name_tot = "data/tot_realizations"
save(tot, file=name_tot)
name_grid = "data/grid"
save(grid, file=name_grid)

## If you have already the simulations ##
load("data/realizations")
load("data/tot_realizations")
load("data/grid")

z_ce <- z
z_bo <- z
z_co <- z

## Shifting and scaling of data

## centre
data_ce <- data.frame(cbind(tot[,1], tot[,2]))

s <- seq(0,1,by=.0001)
z1 <- cbind(s,1-s,0)
z2 <- cbind(0,s,1-s)
z3 <- cbind(s,0,1-s)
H = helm(3)

ps_1 = t(H%*%(t(z1)^alpha_0-1))/alpha_0
ps_2 = t(H%*%(t(z2)^alpha_0-1))/alpha_0
ps_3 = t(H%*%(t(z3)^alpha_0-1))/alpha_0

### alpha=0.2
s_ce = 0.5
scaled_ce <- data.frame(scale(data_ce, center = FALSE, scale=c(1/s_ce, 1/s_ce)))
s_bo = 0.5
scaled_bo <- data.frame(scale(data_ce, center = c(-2.3, +1), scale=c(1/s_bo, 1/s_bo)))
s_co = 0.38
scaled_co <- data.frame(scale(data_ce, center = c(4, -3), scale=c(1/s_co, 1/s_co)))

### alpha=0.6
# s_ce=0.16
# scaled_ce <- data.frame(scale(data_ce, center = FALSE, scale=c(1/s_ce, 1/s_ce)))
# s_bo = 0.16
# scaled_bo <- data.frame(scale(data_ce, center =c(-2.3, +1), scale=c(1/s_bo, 1/s_bo)))
# s_co = 0.11
# scaled_co <- data.frame(scale(data_ce, center = c(4, -3), scale=c(1/s_co, 1/s_co)))

### alpha=1
# s_ce=0.065
# scaled_ce <- data.frame(scale(data_ce, center = FALSE, scale=c(1/s_ce, 1/s_ce)))
# s_bo = 0.065
# scaled_bo <- data.frame(scale(data_ce, center = c(-2.3, +1), scale=c(1/s_bo, 1/s_bo)))
# s_co = 0.045
# scaled_co <- data.frame(scale(data_ce, center = c(4, -3), scale=c(1/s_co, 1/s_co)))

### alpha=0
# scaled_ce <- data.frame(data_ce)
# scaled_bo <- data.frame(scale(data_ce, center = c(-2.3, +1)))
# scaled_co <- data.frame(scale(data_ce, center = c(4, -3)))


plot_euc_complete <- ggplot() +
  theme_bw() +
  geom_point(data = scaled_co, aes(X1, X2, colour = X1)) +
  geom_line(data = data.frame(ps_1), aes(X1, X2)) +
  geom_line(data = data.frame(ps_2), aes(X1, X2)) +
  geom_line(data = data.frame(ps_3), aes(X1, X2)) +
  scale_color_viridis(discrete = FALSE, option = "D") +
  scale_fill_viridis(discrete = FALSE) +
  labs(x=TeX('$z_1$'),y=TeX('$z_2$')) +
  xlim(-3,3) +
  ylim(-3,3) +
  theme(axis.text = element_text(size=22),
        axis.title = element_text(size=22),
        legend.position = "none")
plot_euc_complete

## Update the values of the realization with the right scaling
for (i in 1:N_try) {
  z_ce[[i]] <- scaled_ce[((i-1)*nrow(z[[1]])+1):(i*nrow(z[[1]])),]
  z_bo[[i]] <- scaled_bo[((i-1)*nrow(z[[1]])+1):(i*nrow(z[[1]])),]
  z_co[[i]] <- scaled_co[((i-1)*nrow(z[[1]])+1):(i*nrow(z[[1]])),]
}

## Back transform into the simplex with inverse alpha-IT
## (we use the inverse alpha-transformation by Tsagris for alpha=0, since it is equal to the ILR)

## CENTER data in simplex
yy_ce <- list()
for (j in 1:N_try) {
  if(alpha_0 == 0) {
    yy_ce[[j]] <- ilrInv(as.matrix(z_ce[[j]]),alpha_0)
  }
  else{
    yy_ce[[j]] <- alpha.IT.inv(z_ce[[j]],alpha_0)
  }
}

## BORDER data in simplex
yy_bo <- list()
for (j in 1:N_try) {
  if(alpha_0 == 0) {
    yy_bo[[j]] <- ilrInv(as.matrix(z_bo[[j]]),alpha_0)
  }
  else{
    yy_bo[[j]] <- alpha.IT.inv(z_bo[[j]],alpha_0)
  }
}

## CORNER data in simplex
yy_co <- list()
for (j in 1:N_try) {
  if(alpha_0 == 0) {
    yy_co[[j]] <- ilrInv(as.matrix(z_co[[j]]),alpha_0)
  }
  else{
    yy_co[[j]] <- alpha.IT.inv(z_co[[j]],alpha_0)
  }
}

## Plots of the compositional data
ternary(as.matrix(yy_ce[[1]]))
ternary(as.matrix(yy_bo[[1]]))
ternary(as.matrix(yy_co[[1]]))

## Define the training set
data_02 <- list(yy_ce = yy_ce,
                yy_bo = yy_bo,
                yy_co = yy_co)

# data_06 <- list(yy_ce = yy_ce,
#                 yy_bo = yy_bo,
#                 yy_co = yy_co)

# data_1 <- list(yy_ce = yy_ce,
#                 yy_bo = yy_bo,
#                 yy_co = yy_co)

# data_0 <- list(yy_ce = yy_ce,
#                 yy_bo = yy_bo,
#                 yy_co = yy_co)

name_data_02 = "data/official_02"
save(data_02, file=name_data_02)
# name_data_06 = "data/official_06"
# save(data_06, file=name_data_06)
# name_data_1 = "data/official_1"
# save(data_1, file=name_data_1)
# name_data_0 = "data/official_0"
# save(data_0, file=name_data_0)


####### START 2 ######
## Run from now on if you have already the compositional simulated realizations ##

load("data/official_02")
# load("data/official_06")
# load("data/official_1")
# load("data/official_0")

## Compositional data and observation points for each of the N_try simulations
load("data/grid")

## Change the name each time (ce, bo, co)
yy_ce <- data_02$yy_ce
yy_bo <- data_02$yy_bo
yy_co <- data_02$yy_co
alpha_0=0.2

# yy_ce <- data_06$yy_ce
# yy_bo <- data_06$yy_bo
# yy_co <- data_06$yy_co
# alpha_0=0.6

# yy_ce <- data_1$yy_ce
# yy_bo <- data_1$yy_bo
# yy_co <- data_1$yy_co
# alpha_0=1

# yy_ce <- data_0$yy_ce
# yy_bo <- data_0$yy_bo
# yy_co <- data_0$yy_co
# alpha_0=0


## Choose of one of the dataset (center, border or corner): we choose border)
yy = yy_bo
# yy = yy_ce
# yy = yy_co

## Number of observation points 
points <- sample(nrow(yy_ce[[1]]), N_samples)

plot(grid[points,])
name_file = paste("data/points_", toString(N_samples), sep="") 
save(points, file=name_file)

load("data/points_500")

## Log-likelihood
alpha.hat = NULL
for (i in 1:N_try){
  LL = NULL
  for (alpha in seq(0.01,1.5,by=0.01)){
    LL = c(LL,LogLik.alphaIT(as.matrix(yy[[i]][points,]),alpha))
  }
  plot(seq(0.01,1.5,by=0.01),LL)
  abline(v=alpha_0)
  alpha.hat[i] = seq(0.01,1.5,by=0.01)[which(LL == min(LL))]
}
hist(alpha.hat)
mean(alpha.hat)
sd(alpha.hat)

## Transform compositional data into Euclidean space with range of alphas
trans <- list()
for (j in 1:N_try) {
  transformed <- list()
  for (i in 1:N) {
    transformed[[i]] <- alpha.IT(as.matrix(yy[[j]]), alpha_krig[i])
  }
  trans[[j]] <- transformed
}

## Example plots
plot(trans[[1]][[1]])
plot(trans[[1]][[8]])

## Variography
nlag <- 30
test.neigh <- neigh.create(ndim = 2, type = 0)

var_models <- list()
for (j in 1:N_try) {
  var_models_i <- list()
  for (i in 1:N) {
    db.data = db.create(x1=grid[,1],x2=grid[,2])
    db.data = db.add(db.data, z1=trans[[j]][[i]][,1], z2=trans[[j]][[i]][,2])
    test.vario <- vario.calc(db.data,nlag=nlag)
    test.model <- model.auto(vario=test.vario, struct=8, draw=TRUE)
    print(test.model)

    var_models_i[[i]] <- test.model
  }
  var_models[[j]] <- var_models_i
}

## Prediction (KRIGING)
krig_data_total <- list()     # kriged data in Simplex space
krig_data_euc_total <- list() # kriged data in Euclidean space

for (j in 1:N_try) {

  ## Dataset with selected points and coordinates
  y <- yy[[j]][points,]
  coord <- grid[points,]
  
  krig_data <- list()
  krig_data_euc <- list()
  
  for (i in 1:N) {
    
    data_i = as.data.frame(cbind(coord, trans[[j]][[i]][points,]), c=c("x.coords", "y.coords", "V3", "V4"))
    db.grid <- db.create(x1=grid[,1],x2=grid[,2])
    
    ## Cokriging
    db.data_sample <- db.create(x1=data_i[,1],x2=data_i[,2],z1=data_i[,3], z2=data_i[,4])
    db.krig <- kriging(db.data_sample,db.grid, var_models[[j]][[i]], test.neigh)
    print(db.krig)
    new_pred <- data.frame(cbind("Z1"=db.krig$Kriging.z1.estim, "Z2"=db.krig$Kriging.z2.estim))
    
    ## Inverse alpha-IT to obtain compositional data (inverse ilr for alpha=0)
    if (alpha_krig[i] == 0) {
      y_new <- data.frame(ilrInv(as.matrix(new_pred), V=t(helm(3))))
    }
    else {
      y_new <- data.frame(alpha.IT.inv(new_pred,alpha_krig[i]))
    }
    
    krig_data_euc[[i]] <- new_pred
    krig_data[[i]] <- y_new
  }
  
  krig_data_total[[j]] <- krig_data
  krig_data_euc_total[[j]] <- krig_data_euc
  print(j)
}

## Save results
results <- list(comp_data = yy,
                points = points, 
                krig_data_total = krig_data_total, 
                krig_data_euc_total = krig_data_euc_total,
                alpha_krig = alpha_krig,
                alpha_0 = alpha_0,
                var_models = var_models, 
                N_try = N_try,
                trans = trans)

# Change the name centre, border, corner for the corresponding realization
name = paste("data/cokriging_border_Ntry_", toString(N_try), "_Nsamples_500_dom_", toString(alpha_0),sep="")
for (i in 1:N) {
  name = paste(name, toString(alpha_krig[i]), sep="_")
}
name = paste(name, ".rdata", sep="")
name
save(results, file=name)
