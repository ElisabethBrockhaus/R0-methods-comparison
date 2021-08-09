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
source('assess_performance.R')


# fix population size parameters: fix population size at 10^6
lower_Npv <- 6
upper_Npv <- 6

# number of simulation runs:
N_sim <- 500

# limits for plots:
xlim_left <- -4
xlim_right <- 4

# number of pixels for plots:
wdt_bias_plot <- 600
hgt_bias_plot <- 250

#===== 1. Original version =====#

# simulate 250 datasets, assuming average infectious and incubation periods of 6 and 14 days
set.seed(123) # JB: setting seed
simulations <- simulate_data(N=N_sim, gamma=1/6, sigma=1/14, lower_Npv = lower_Npv,
                             upper_Npv = upper_Npv, n_aggr = 7)
plot(simulations[[1]]$sim2[1:150])


# define the generation time distribution
GTd <- R0::generation.time("gamma", c(20/7, 7.4/7))

# list to store results from simulations
Sim_results <- list()

# fit each method at sequential time points of epidemic growth phase
for(i in 1:1){ # run only without noise as shown in manuscript
  
  # list for storage
  results_i <- list()
  
  # for each of the 250 simulations
  for(j in 1:ncol(simulations$sims)){ 
    
    # data for fitting - simulation j with noise level i
    sim_ij <- data.frame(week=1:nrow(simulations[[i]]), cases=simulations[[i]][ ,j], 
                         country=paste("sim", j, sep='_'))
    
    # fit & store results
    # run only for EpiEstim
    results_i[[j]] <- fit_R0_seq(data=sim_ij, mean_GT=20/7, sd_GT=7.4/7, GTd=GTd, GT_week=3,
                                 methods = c("EpiEstim"))
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
png("text/figures/bias_7.png", width = wdt_bias_plot, height = hgt_bias_plot)
bias_plot(metrics[[1]]$metrics_indiv) + xlim(xlim_left, xlim_right)
dev.off()

# plot all results:
scatterplot_results(metrics[[1]], main = "No noise, original setting")


#===== 2. With sd of generation time set to 15.25 =====#

# define the generation time distribution
GTd_15 <- R0::generation.time("gamma", c(20/7, 15.25/7))

# list to store results from simulations
Sim_results_15 <- list()

# fit each method at sequential time points of epidemic growth phase
for(i in 1:1){ # run only without noise for now
  
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

# metrics summary plot & bias plot
png("text/figures/bias_15.png", width = wdt_bias_plot, height = hgt_bias_plot)
bias_plot(metrics_15[[1]]$metrics_indiv) + xlim(xlim_left, xlim_right)
dev.off()

# plot all results:
scatterplot_results(metrics_15[[1]], main = "No noise, sd of gen time = 15.25")


#===== 3. With sd of generation time set to 15.25 and daily data =====#

set.seed(123) # JB: setting seed
simulations_daily <- simulate_data(N=N_sim, gamma=1/6, sigma=1/14, lower_Npv = lower_Npv, 
                                   upper_Npv = upper_Npv, n_aggr = 1)

# define the generation time distribution
GTd_15_daily <- R0::generation.time("gamma", c(20, 15.25))

# list to store results from simulations
Sim_results_15_daily <- list()

# fit each method at sequential time points of epidemic growth phase
for(i in 1:1){ # run only without noise for now
  
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

# metrics summary plot & bias plot
png("text/figures/bias_15_daily.png", width = wdt_bias_plot, height = hgt_bias_plot)
bias_plot(metrics_15_daily[[1]]$metrics_indiv) + xlim(xlim_left, xlim_right)
dev.off()

# plot all results:
scatterplot_results(metrics_15_daily[[1]], main = "No noise, sd of gen time = 15.25, daily data")
scatterplot_results(metrics_15_daily[[2]], main = "Poisson noise, sd of gen time = 15.25, daily data")

#===== 4. With sd of generation time set to 15.25, daily data and a discrete-valued simulation process =====#

set.seed(123)
simulations_daily_stoch <- simulate_data_stoch(N=N_sim, gamma=1/6, sigma=1/14, lower_Npv = lower_Npv, 
                                               upper_Npv = upper_Npv, n_aggr = 1, min_Iv = 5, max_Iv = 5)

plot(simulations_daily_stoch$sims$sim12)

# list to store results from simulations
Sim_results_15_daily_stoch <- list()

# for each of the 3 noise levels
for(i in 1:1){ # only one slot exists
  
  # list for storage
  results_i <- list()
  
  # for each of the 250 simulations
  for(j in 1:ncol(simulations$sims)){ 
    
    # data for fitting - simulation j with noise level i
    sim_ij <- data.frame(week=1:nrow(simulations_daily_stoch[[i]]), cases=simulations_daily_stoch[[i]][ ,j], 
                         country=paste("sim", j, sep='_'))
    
    # fit & store results
    # run only for EpiEstim for now
    # undebug(fit_R0_seq)
    results_i[[j]] <- fit_R0_seq(data=sim_ij, mean_GT=20, sd_GT=15.25, GTd=GTd, GT_week=3*7, methods = "EpiEstim", min_peak_time = 7*15)
    print(j)
  }
  
  # store results from each noise level
  Sim_results_15_daily_stoch[[i]] <- results_i
}


# calculate performance metrics
metrics_15_daily_stoch <- list()
for(i in seq_along(Sim_results_15_daily_stoch)){
  metrics_15_daily_stoch[[i]] <- performance_metrics(trueR0=simulations_daily_stoch$pars,
                                                     estimates=Sim_results_15_daily_stoch[[i]], max_weeks=7*15, min_peak=7*15)
}

# metrics summary plot & bias plot
png("text/figures/bias_15_daily_stoch.png", width = wdt_bias_plot, height = hgt_bias_plot)
bias_plot(metrics_15_daily_stoch[[1]]$metrics_indiv) + xlim(xlim_left, xlim_right)
dev.off()

# plot all results:
scatterplot_results(metrics_15_daily_stoch[[1]], main = "Sd of gen time = 15.25, daily data, stochastic simulation")



png("text/figures/scatterplots.png", width = 1200, height = 800)
par(mfrow = c(2, 2), cex = 1.3)

scatterplot_results(metrics[[1]], main = "No noise, sd of gen time = 7.4 days")
# scatterplot_results(metrics[[2]], main = "Poisson noise, original setting")

scatterplot_results(metrics_15[[1]], main = "No noise, sd of gen time = 15.25 days")
# scatterplot_results(metrics_15[[2]], main = "Poisson noise, sd of gen time = 15.25")

scatterplot_results(metrics_15_daily[[1]], 
                    main = "No noise, sd of gen time = 15.25 days,\n daily data",
                    tag_legend = "days")
# scatterplot_results(metrics_15_daily[[2]], main = "Poisson noise, sd of gen time = 15.25,\n daily data")

# plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = 0:1, ylim = 0:1)
scatterplot_results(metrics_15_daily_stoch[[1]],
                    main = "Sd of gen time = 15.25 days, daily\n data, stochastic simulation",
                    tag_legend = "days")
dev.off()
