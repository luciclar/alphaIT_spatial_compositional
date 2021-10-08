##### Analysis of the errors of the simulation study #####
# written by Lucia Clarotto #

rm(list=ls())

## Set "path" to your local path and uncomment
# path = ''
# setwd(path)

## Libraries
library(Ternary)
library(ggplot2)
library(reshape2)
library(wesanderson)
library(latex2exp)
library(viridis)
library(tidyr)

# The data come from the outputs of the R code "1_simulationcase_kriging.R"
# We present here the results of the analysis of the simulation with alpha = 0.2

## Analysis for alpha=0.2
load("data/cokriging_border_Ntry_10_Nsamples_500_dom_0.2_0_0.1_0.2_0.4_0.6_0.8_0.9_1.rdata")
data_bo <- results

# load("data/cokriging_centre_Ntry_10_Nsamples_500_dom_0.2_0_0.1_0.2_0.4_0.6_0.8_0.9_1.rdata")
# data_ce <- results
# load("data/cokriging_corner_Ntry_10_Nsamples_500_dom_0.2_0_0.1_0.2_0.4_0.6_0.8_0.9_1.rdata")
# data_co <- results

ternary(data_bo$comp_data[[1]])
# ternary(data_ce$comp_data[[1]])
# ternary(data_co$comp_data[[1]])

alpha_names=as.character(data_bo$alpha_krig)
N_try = data_bo$N_try
N = length(data_bo$alpha_krig)


### ERROR ANALYSIS

## CENTER data

## Total Variation
# support <- matrix(c(1,2,3), ncol=1)
# error_tv <- NULL
# for (k in 1:length(data_ce$alpha_krig)){
#   error_tv_k <- NULL
#   for (j in 1:N_try) {
#     e1 = data_ce$krig_data_total[[j]][[k]][-data_ce$points,1:3]
#     e2 = data_ce$comp_data[[j]][-data_ce$points,]
#     error_tv_k = c(error_tv_k, mean(0.5*rowSums(abs(e1-e2))))
#   }
#   error_tv<- cbind(error_tv, error_tv_k)
# }
# error_tv <- data.frame(error_tv)
# colnames(error_tv) <- alpha_names
# mean_tv_ce <- data.frame(cbind(data_ce$alpha_krig, colMeans(error_tv)))
# mean_tv_ce
# std_tv_ce <- data.frame(cbind(data_ce$alpha_krig, 2*sqrt(apply(error_tv,2,var))))
# std_tv_ce
# 
# mins_tv <- mean_tv_ce[,2] - std_tv_ce[,2]
# maxs_tv <- mean_tv_ce[,2] + std_tv_ce[,2]
# 
# ## Hellinger
# error_hel <- NULL
# for (k in 1:length(data_ce$alpha_krig)){
#   error_hel_k <- NULL
#   for (j in 1:N_try) {
#     e1 = data_ce$krig_data_total[[j]][[k]][-data_ce$points,1:3]
#     e2 = data_ce$comp_data[[j]][-data_ce$points,]
#     error_hel_k = c(error_hel_k, mean(0.5*rowSums((sqrt(e1)-sqrt(e2))^2)))
#   }
#   error_hel<- cbind(error_hel, error_hel_k)
# }
# error_hel <- data.frame(error_hel)
# colnames(error_hel) <- alpha_names
# mean_hel_ce <- data.frame(cbind(data_ce$alpha_krig, colMeans(error_hel)))
# 
# std_hel_ce <- data.frame(cbind(data_ce$alpha_krig, 2*sqrt(apply(error_hel,2,var))))
# std_hel_ce
# 
# mins_hel <- mean_hel_ce[,2] - std_hel_ce[,2]
# maxs_hel <- mean_hel_ce[,2] + std_hel_ce[,2]
# 
# coeff_ce <- mean_hel_ce[1,2]/mean_tv_ce[1,2]
# mean_ce <- round(cbind(mean_tv_ce, mean_hel_ce[,2]/coeff_ce),7)
# colnames(mean_ce) <- c("alpha","Total Variation", "Hellinger")
# df_ce <- gather(mean_ce, key = metric, value = Rate, 
#              c("Total Variation", "Hellinger"))
# 

## BORDER data

## Total Variation
support <- matrix(c(1,2,3), ncol=1)
error_tv <- NULL
for (k in 1:length(data_bo$alpha_krig)){
  error_tv_k <- NULL
  for (j in 1:N_try) {
    e1 = data_bo$krig_data_total[[j]][[k]][-data_bo$points,1:3]
    e2 = data_bo$comp_data[[j]][-data_bo$points,]
    error_tv_k = c(error_tv_k, mean(0.5*rowSums(abs(e1-e2))))
  }
  error_tv<- cbind(error_tv, error_tv_k)
}
error_tv <- data.frame(error_tv)
colnames(error_tv) <- alpha_names
mean_tv_bo <- data.frame(cbind(data_bo$alpha_krig, colMeans(error_tv)))

## Hellinger
error_hel <- NULL
for (k in 1:length(data_bo$alpha_krig)){
  error_hel_k <- NULL
  for (j in 1:N_try) {
    e1 = data_bo$krig_data_total[[j]][[k]][-data_bo$points,1:3]
    e2 = data_bo$comp_data[[j]][-data_bo$points,]
    error_hel_k = c(error_hel_k, mean(0.5*rowSums((sqrt(e1)-sqrt(e2))^2)))
  }
  error_hel<- cbind(error_hel, error_hel_k)
}
error_hel <- data.frame(error_hel)
colnames(error_hel) <- alpha_names
mean_hel_bo <- data.frame(cbind(data_bo$alpha_krig, colMeans(error_hel)))

coeff_bo <- mean_hel_bo[1,2]/mean_tv_bo[1,2]
mean_bo <- round(cbind(mean_tv_bo, mean_hel_bo[,2]/coeff_bo),6)
colnames(mean_bo) <- c("alpha","Total Variation", "Hellinger")
df_bo <- gather(mean_bo, key = metric, value = Rate, 
                c("Total Variation", "Hellinger"))

## CORNER data

## Total Variation
# support <- matrix(c(1,2,3), ncol=1)
# error_tv <- NULL
# for (k in 1:length(data_co$alpha_krig)){
#   error_tv_k <- NULL
#   for (j in 1:N_try) {
#     e1 = data_co$krig_data_total[[j]][[k]][-data_co$points,1:3]
#     e2 = data_co$comp_data[[j]][-data_co$points,]
#     error_tv_k = c(error_tv_k, mean(0.5*rowSums(abs(e1-e2))))
#   }
#   error_tv<- cbind(error_tv, error_tv_k)
# }
# error_tv <- data.frame(error_tv)
# colnames(error_tv) <- alpha_names
# mean_tv_co <- data.frame(cbind(data_co$alpha_krig, colMeans(error_tv)))
# 
# ## Hellinger
# error_hel <- NULL
# for (k in 1:length(data_co$alpha_krig)){
#   error_hel_k <- NULL
#   for (j in 1:N_try) {
#     e1 = data_co$krig_data_total[[j]][[k]][-data_co$points,1:3]
#     e2 = data_co$comp_data[[j]][-data_co$points,]
#     error_hel_k = c(error_hel_k, mean(0.5*rowSums((sqrt(e1)-sqrt(e2))^2)))
#   }
#   error_hel<- cbind(error_hel, error_hel_k)
# }
# error_hel <- data.frame(error_hel)
# colnames(error_hel) <- alpha_names
# mean_hel_co <- data.frame(cbind(data_co$alpha_krig, colMeans(error_hel)))
# 
# coeff_co <- mean_hel_co[1,2]/mean_tv_co[1,2]
# mean_co <- round(cbind(mean_tv_co, mean_hel_co[,2]/coeff_co),6)
# colnames(mean_co) <- c("alpha","Total Variation", "Hellinger")
# df_co <- gather(mean_co, key = metric, value = Rate, 
#                 c("Total Variation", "Hellinger"))

## Plots
colours <- wes_palette(n=5, name="Zissou1")
coltv <- colours[1]
colhel <- colours[5]

# mean_ce_plot <- ggplot(df_ce) + 
#   geom_line(aes(x=alpha, y = Rate, group = metric, colour = metric, linetype = metric), size=2) + 
#   geom_point(aes(x=alpha, y = Rate, group = metric, colour = metric, shape=metric), size=5) +
#   scale_y_continuous(
#     # Features of the first axis
#     name = TeX('$\\delta_{TV}$'),
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~.*coeff_ce, name=TeX('$\\delta_{H}$'))) + 
#   theme_bw() +
#   labs(title = TeX('center data, $\\alpha_0=1$'), x=TeX('$\\alpha$')) +
#   theme(axis.title.y = element_text(color = coltv,size=30),
#         axis.title.y.right = element_text(color = colhel,size=30),
#         axis.text = element_text(size=25),
#         axis.title.x = element_text(size=30),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 25),
#         legend.position = "bottom",
#         legend.key.width = unit(2, 'cm'),
#         plot.title = element_text(hjust = 0.5, size = 30))
# mean_ce_plot

mean_bo_plot <- ggplot(df_bo) + 
  geom_line(aes(x=alpha, y = Rate, group = metric, colour = metric, linetype = metric), size=2) + 
  geom_point(aes(x=alpha, y = Rate, group = metric, colour = metric, shape=metric), size=5) +
  scale_y_continuous(
    # Features of the first axis
    name = TeX('$\\delta_{TV}$'),
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff_bo, name=TeX('$\\delta_{H}$'))) + 
  theme_bw() +
  labs(title = TeX('border data, $\\alpha_0=0.2$'), x=TeX('$\\alpha$')) +
  theme(axis.title.y = element_text(color = coltv,size=30),
        axis.title.y.right = element_text(color = colhel,size=30),
        axis.text = element_text(size=25),
        axis.title.x = element_text(size=30),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position = "bottom",
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 30))
mean_bo_plot

# mean_co_plot <- ggplot(df_co) + 
#   geom_line(aes(x=alpha, y = Rate, group = metric, colour = metric, linetype = metric), size=2) + 
#   geom_point(aes(x=alpha, y = Rate, group = metric, colour = metric, shape=metric), size=5) +
#   scale_y_continuous(
#     # Features of the first axis
#     name = TeX('$\\delta_{TV}$'),
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~.*coeff_co, name=TeX('$\\delta_{H}$'))) + 
#   theme_bw() +
#   labs(title = TeX('corner data, $\\alpha_0=0.2$'), x=TeX('$\\alpha$')) +
#   theme(axis.title.y = element_text(color = coltv,size=30),
#         axis.title.y.right = element_text(color = colhel,size=30),
#         axis.text = element_text(size=25),
#         axis.title.x = element_text(size=30),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 25),
#         legend.position = "bottom",
#         legend.key.width = unit(2, 'cm'),
#         plot.title = element_text(hjust = 0.5, size = 30))
# mean_co_plot


### RMSE before inverse alpha-it

## CENTER data
# rmse <- NULL
# for (k in 1:length(data_ce$alpha_krig)){
#   rmse_k <- NULL
#   for (j in 1:N_try) {
#     e1 = data.frame(data_ce$krig_data_euc_total[[j]][[k]][-data_ce$points,])
#     e2 = data.frame(data_ce$trans[[j]][[k]][-data_ce$points,])
#     mse = colSums((e1 - e2)^2)/(nrow(e2)*diag(var(as.matrix(data_ce$trans[[j]][[k]]))))
#     rmse_k = c(rmse_k, sqrt(sum(mse)))
#   }
#   rmse <- cbind(rmse, rmse_k)
# }
# rmse <- data.frame(rmse)
# colnames(rmse) <- alpha_names
# mean_rmse_ce <- data.frame(cbind(data_ce$alpha_krig, colMeans(rmse)))

## BORDER data
rmse <- NULL
for (k in 1:length(data_bo$alpha_krig)){
  rmse_k <- NULL
  for (j in 1:N_try) {
    e1 = data.frame(data_bo$krig_data_euc_total[[j]][[k]][-data_bo$points,])
    e2 = data.frame(data_bo$trans[[j]][[k]][-data_bo$points,])
    mse = colSums((e1 - e2)^2)/(nrow(e2)*diag(var(as.matrix(data_bo$trans[[j]][[k]]))))
    rmse_k = c(rmse_k, sqrt(sum(mse)))
  }
  rmse <- cbind(rmse, rmse_k)
}
rmse <- data.frame(rmse)
colnames(rmse) <- alpha_names
mean_rmse_bo <- data.frame(cbind(data_bo$alpha_krig, colMeans(rmse)))

## CORNER data
# rmse <- NULL
# for (k in 1:length(data_co$alpha_krig)){
#   rmse_k <- NULL
#   for (j in 1:N_try) {
#     e1 = data.frame(data_co$krig_data_euc_total[[j]][[k]][-data_co$points,])
#     e2 = data.frame(data_co$trans[[j]][[k]][-data_co$points,])
#     mse = colSums((e1 - e2)^2)/(nrow(e2)*diag(var(as.matrix(data_co$trans[[j]][[k]]))))
#     rmse_k = c(rmse_k, sqrt(sum(mse)))
#   }
#   rmse <- cbind(rmse, rmse_k)
# }
# rmse <- data.frame(rmse)
# colnames(rmse) <- alpha_names
# mean_rmse_co <- data.frame(cbind(data_co$alpha_krig, colMeans(rmse)))

## Plots (if you have centre, border and corner data)

# means_rmse <- data.frame(cbind(mean_rmse_ce[,2], mean_rmse_bo[,2], mean_rmse_co[,2]))
# colnames(means_rmse) <- c("center","border","corner")
# stack_rmse <- cbind(stack(means_rmse),rep(data_ce$alpha_krig,3))
# colnames(stack_rmse) <- c("values","ind","alpha")

# lines <- c(1,2,3)
# cols <- c('mediumturquoise','darkmagenta','gold1')
# means_rmse_plot <- ggplot(stack_rmse, aes(x = alpha, y = values, col = ind, linetype = ind, shape=ind)) + 
#   geom_line( size=2) + geom_point( size=5) +
#   ggplot2::theme_bw() +
#   theme(axis.title = element_text(size=30),
#         axis.text = element_text(size=25),
#         axis.title.x = element_text(size=30),
#         plot.title = element_text(hjust = 0.5, size = 30),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 25),
#         legend.position = "bottom",
#         legend.key.width = unit(2, 'cm')) +
#   scale_fill_manual(values=cols) +
#   scale_linetype_manual(values=lines) +
#   labs(title = TeX('RMSE $\\alpha_0=0.2$'), x=TeX('$\\alpha$'), y=TeX('$\\bar{\\delta}_{\\alpha}/\\sigma$'))
# 
# means_rmse_plot
