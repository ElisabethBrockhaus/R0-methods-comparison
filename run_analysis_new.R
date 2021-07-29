library(dplyr)
library(bbmle)
library(EpiEstim)
library(R0)
library(ggplot2)
library(deSolve)
library(truncnorm)
library(ggridges)

# specify repository path 
# repo_path <- 'D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code/R0_methods_comparison_Elisabeth'
repo_path <- '/home/johannes/Documents/Projects/R0-methods-comparison-1'
setwd(repo_path)

# source functions for estimating R0 & data simulation
source('methods.R')
source('fit_R0_seq.R')
source('data_simulation.R')
source('assess_performance.R')


# helper function to plot:
scatterplot_results <- function(metrics, ...){
  plot(metrics$metrics_indiv$R0_true, metrics$metrics_indiv$R0,
       col = metrics$metrics_indiv$Nweeks/3, ylim = c(0, 8), xlim = c(0, 8),
       xlab = "true R0", ylab = "estimated R0", ...)
  legend("topleft", legend = paste(unique(metrics$metrics_indiv$Nweeks), "obs"),
         col = unique(metrics$metrics_indiv$Nweeks/3), pch = 1)
  abline(0:1, lty = "dotted")
}

# fix population size parameters:


#===== 1. Original version =====#

# simulate 250 datasets, assuming average infectious and incubation periods of 6 and 14 days
set.seed(125) # JB: setting seed
simulations <- simulate_data(N=250, gamma=1/6, sigma=1/14, lower_Npv = lower_Npv,
                             upper_Npv = upper_Npv, n_aggr = 7)
plot(simulations[[1]]$sim2[1:150])


# define the generation time distribution
GTd <- R0::generation.time("gamma", c(20/7, 7.4/7))

# list to store results from simulations
Sim_results <- list()

# fit each method at sequential time points of epidemic growth phase
# note: running this locally in a loop will take a long time.
# where possible, we recommend using parallel computing for this.

# for each of the 3 noise levels
for(i in 1:2){ # run only without noise for now
  
  # list for storage
  results_i <- list()
  
  # for each of the 250 simulations
  for(j in 1:ncol(simulations$sims)){ 
    
    # data for fitting - simulation j with noise level i
    sim_ij <- data.frame(week=1:nrow(simulations[[i]]), cases=simulations[[i]][ ,j], 
                         country=paste("sim", j, sep='_'))
    
    # fit & store results
    # run only for EpiEstim for now
    results_i[[j]] <- fit_R0_seq(data=sim_ij, mean_GT=20/7, sd_GT=7.4/7, GTd=GTd, GT_week=3, methods = "EpiEstim")
    print(j)
  }
  
  # store results from each noise level
  Sim_results[[i]] <- results_i
}


# calculate performance metrics
metrics <- list()
for(i in seq_along(Sim_results)){
  metrics[[i]] <- performance_metrics(trueR0=simulations$pars, estimates=Sim_results[[i]], max_weeks=15, min_peak=15)
}


# choose metrics to plot
get_metrics()
want <- c('Bias','Coverage','Uncertainty','RMSE')

# metrics summary plot & bias plot
bias_plot(metrics[[1]]$metrics_indiv)
bias_plot(metrics[[2]]$metrics_indiv)

# plot all results:
scatterplot_results(metrics[[1]], main = "No noise, original setting")
scatterplot_results(metrics[[2]], main = "Poisson noise, original setting")


#===== 2. With sd of generation time set to 15.25 =====#

# define the generation time distribution
GTd_15 <- R0::generation.time("gamma", c(20/7, 15.25/7))

# list to store results from simulations
Sim_results_15 <- list()

# fit each method at sequential time points of epidemic growth phase
# note: running this locally in a loop will take a long time.
# where possible, we recommend using parallel computing for this.

# for each of the 3 noise levels
for(i in 1:2){ # run only without noise for now
  
  # list for storage
  results_i <- list()
  
  # for each of the 250 simulations
  for(j in 1:ncol(simulations$sims)){ 
    
    # data for fitting - simulation j with noise level i
    sim_ij <- data.frame(week=1:nrow(simulations[[i]]), cases=simulations[[i]][ ,j], 
                         country=paste("sim", j, sep='_'))
    
    # fit & store results
    # run only for EpiEstim for now
    results_i[[j]] <- fit_R0_seq(data=sim_ij, mean_GT=20/7, sd_GT=15.25/7, GTd=GTd, GT_week=3, methods = "EpiEstim")
    print(j)
  }
  
  # store results from each noise level
  Sim_results_15[[i]] <- results_i
}


# calculate performance metrics
metrics_15 <- list()
for(i in seq_along(Sim_results_15)){
  metrics_15[[i]] <- performance_metrics(trueR0=simulations$pars, estimates=Sim_results_15[[i]], max_weeks=15, min_peak=15)
}

# choose metrics to plot
get_metrics()
want <- c('Bias','Coverage','Uncertainty','RMSE')

# metrics summary plot & bias plot
bias_plot(metrics_15[[1]]$metrics_indiv)
bias_plot(metrics_15[[2]]$metrics_indiv)

# plot all results:
scatterplot_results(metrics_15[[1]], main = "No noise, sd of gen time = 15.25")
scatterplot_results(metrics_15[[2]], main = "Poisson noise, sd of gen time = 15.25")


#===== 3. With sd of generation time set to 15.25 and daily data =====#

set.seed(125) # JB: setting seed
simulations_daily <- simulate_data(N=250, gamma=1/6, sigma=1/14, lower_Npv = lower_Npv, 
                                   upper_Npv = upper_Npv, n_aggr = 1)

# define the generation time distribution
GTd_15_daily <- R0::generation.time("gamma", c(20, 15.25))

# list to store results from simulations
Sim_results_15_daily <- list()

# fit each method at sequential time points of epidemic growth phase
# note: running this locally in a loop will take a long time.
# where possible, we recommend using parallel computing for this.

# for each of the 3 noise levels
for(i in 1:2){ # run only without noise for now
  
  # list for storage
  results_i <- list()
  
  # for each of the 250 simulations
  for(j in 1:ncol(simulations$sims)){ 
    
    # data for fitting - simulation j with noise level i
    sim_ij <- data.frame(week=1:nrow(simulations_daily[[i]]), cases=simulations_daily[[i]][ ,j], 
                         country=paste("sim", j, sep='_'))
    
    # fit & store results
    # run only for EpiEstim for now
    # undebug(fit_R0_seq)
    results_i[[j]] <- fit_R0_seq(data=sim_ij, mean_GT=20, sd_GT=15.25, GTd=GTd, GT_week=3*7, methods = "EpiEstim", min_peak_time = 7*15)
    print(j)
  }
  
  # store results from each noise level
  Sim_results_15_daily[[i]] <- results_i
}


# calculate performance metrics
metrics_15_daily <- list()
for(i in seq_along(Sim_results_15_daily)){
  metrics_15_daily[[i]] <- performance_metrics(trueR0=simulations_daily$pars, estimates=Sim_results_15_daily[[i]], max_weeks=7*15, min_peak=7*15)
}

# choose metrics to plot
get_metrics()
want <- c('Bias','Coverage','Uncertainty','RMSE')

# metrics summary plot & bias plot
bias_plot(metrics_15_daily[[1]]$metrics_indiv)
bias_plot(metrics_15_daily[[2]]$metrics_indiv)

# plot all results:
scatterplot_results(metrics_15_daily[[1]], main = "No noise, sd of gen time = 15.25, daily data")
scatterplot_results(metrics_15_daily[[2]], main = "Poisson noise, sd of gen time = 15.25, daily data")

png("scatterplots.png", width = 1200, height = 900)
par(mfcol = c(2, 3), cex = 1.3)

scatterplot_results(metrics[[1]], main = "No noise, original setting")
scatterplot_results(metrics[[2]], main = "Poisson noise, original setting")

scatterplot_results(metrics_15[[1]], main = "No noise, sd of gen time = 15.25")
scatterplot_results(metrics_15[[2]], main = "Poisson noise, sd of gen time = 15.25")

scatterplot_results(metrics_15_daily[[1]], main = "No noise, sd of gen time = 15.25,\n daily data")
scatterplot_results(metrics_15_daily[[2]], main = "Poisson noise, sd of gen time = 15.25,\n daily data")
dev.off()

# Could potentially be added:
#===== 4. With sd of generation time set to 15.25, daily data and a discrete-valued simulation process =====#
