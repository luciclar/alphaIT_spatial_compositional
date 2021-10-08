##### Analysis of the errors of the case study: Copernicus Land Cover Map (Po' Valley) #####
# written by Lucia Clarotto #

rm(list=ls())

## Libraries
library(Ternary)
library(ggplot2)
library(lattice)
library(reshape2)
library(wesanderson)
library(latex2exp)
library(viridis)

## The data come from the outputs of the R_code "casestudy_kriging.R"
## We present here the results of the analysis of the case study with N_try=5, N_train=500 and 6 different values of alpha

## Dataset without zeros
load("data/copernicus_cokriging_N_try_5_N_train_500_krig_0_0.12_0.3_0.5_0.75_1.rdata")
## Dataset with zeros
# load("data/copernicus_zeros_cokriging_N_try_5_N_train_500_krig_0.01_0.12_0.3_0.5_0.75_1.rdata")

data <- results

alpha_krig=results$alpha_krig
alpha_names=as.character(alpha_krig)

N_try = results$N_try
N = length(alpha_krig)

### ERROR ANALYSIS

## Total Variation
support <- matrix(c(1,2,3,4), ncol=1)
error_tv <- NULL
for (k in 1:N){
  error_tv_k <- NULL
  for (j in 1:N_try) {
    e1 = data$krig_data_total[[j]][[k]][-data$points[[j]],]
    e2 = data$comp_data[-data$points[[j]],]
    error_tv_k = c(error_tv_k, mean(0.5*rowSums(abs(e1-e2))))
  }
  error_tv<- cbind(error_tv, error_tv_k)
}
error_tv <- data.frame(error_tv)
colnames(error_tv) <- alpha_names

mins_tv <- apply(error_tv,2,min)
maxs_tv <- apply(error_tv,2,max)

mean_tv <- data.frame(cbind(alpha_krig, colMeans(error_tv)))
mean_tv

## Hellinger
error_hel <- NULL
for (k in 1:N){
  error_hel_k <- NULL
  for (j in 1:N_try) {
    e1 = data$krig_data_total[[j]][[k]][-data$points[[j]],]
    e2 = data$comp_data[-data$points[[j]],]
    error_hel_k = c(error_hel_k, mean(0.5*rowSums((sqrt(e1)-sqrt(e2))^2)))
  }
  error_hel<- cbind(error_hel, error_hel_k)
}
error_hel <- data.frame(error_hel)
colnames(error_hel) <- alpha_names

mins_hel <- apply(error_hel,2,min)
maxs_hel <- apply(error_hel,2,max)

mean_hel <- data.frame(cbind(alpha_krig, colMeans(error_hel)))
mean_hel

mean <- round(cbind(mean_tv, mins_tv, maxs_tv, mean_hel[,2], mins_hel, maxs_hel),7)
colnames(mean) <- c("alpha","totvar","totvar_min", "totvar_max", "hell", "hell_min", "hell_max")
coeff <- mean[1,5]/mean[1,2]

colours <- wes_palette(n=5, name="Zissou1")
coltv <- colours[1]
colhel <- colours[5]

mean_plot <- ggplot(mean, aes(x=alpha)) + 
  geom_line( aes(y=totvar), size=1.5, color=coltv) + 
  geom_point( aes(y=totvar), size=3, color=coltv) +
  geom_line( aes(y=hell/coeff), size=1.5, color=colhel) +
  geom_point( aes(y=hell/coeff), size=3, color=colhel) +
  scale_y_continuous(
    # Features of the first axis
    name = TeX('$\\delta_{TV}$'),
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name=TeX('$\\delta_{H}$'))) + 
  theme_bw() +
  labs(title = TeX(''), x=TeX('$\\alpha$')) +
  theme(axis.title.y = element_text(color = coltv,size=22),
        axis.title.y.right = element_text(color = colhel,size=22),
        plot.title = element_text(hjust = 0.5, size = 22),
        axis.text = element_text(size=15),
        axis.title.x = element_text(size=22))
mean_plot

## RMSE before inverse alpha-IT
rmse <- NULL
for (k in 1:N){
  rmse_k <- NULL
  for (j in 1:N_try) {
    e1 = data.frame(data$krig_data_euc_total[[j]][[k]][-data$points[[j]],])
    e2 = data.frame(data$trans[[k]][-data$points[[j]],])
    mse = colSums((e1 - e2)^2)/(nrow(e2)*diag(var(as.matrix(data$trans[[k]]))))
    rmse_k = c(rmse_k, sqrt(sum(mse)))
  }
  rmse <- cbind(rmse, rmse_k)
}
rmse <- data.frame(rmse)
colnames(rmse) <- alpha_names
mean_rmse <- data.frame(cbind(alpha_krig, colMeans(rmse)))
mean_rmse


