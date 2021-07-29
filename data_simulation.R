#=== Function to simulate N epidemic curves with 3 levels of noise:
# - no noise
# - Poisson random errors
# - Negative Binomial random errors

# new arguments:
# lower_Npv: lower end for sampling of population size
# upper_Npv: upper end
# n_aggr: how many days to aggregate into one observation? defaults to 7 = one week
simulate_data <- function(N, gamma, sigma, lower_Npv = 3, upper_Npv = 6, n_aggr = 7,
                          min_Iv = 3, max_Iv = 8){
  
  #=== Define the SEIR model
  
  # gamma = 1 / average infectious period
  # sigma = 1 / average incubation period
  
  model <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      
      dS <- -(R0*gamma) * S * I/(S + E + I + R)
      dE <- (R0*gamma) * S * I/(S + E + I + R) - sigma * E 
      dI <- sigma * E - gamma * I 
      dR <- gamma * I
      
      return(list(c(dS, dE, dI, dR)))
      
    })
  }
  
  
  
  #=== Functions to generate Poisson and Negative Binomial random errors
  
  # Poisson 
  add_noiseP <- function(inc){ 
    res <- round(rpois(length(inc), lambda=inc))
    return(res)
  }
  
  # Negative binomial (over-dispersed)
  add_noiseNB <- function(inc){ 
    res <- round(rnbinom(length(inc), size=2, mu=inc))
    return(res)
  }
  
  
  
  #=== Generate N simulated epidemic curves with randomly selected parameter values
  
  # randomly select values of R0, population size (Np), and inital number infected (I)
  # JB changed to min R = 2, mean = 3
  R0v <- round(rtruncnorm(N, a=2, b=8, mean=3, sd=2),2) 
  # JB: set up population size by a factor of 10^3
  Npv <- if(lower_Npv == upper_Npv) rep(10^lower_Npv, N) else 10^(round(runif(N, min=lower_Npv, max=upper_Npv),0))
  #Npv <- 10^(round(runif(N, min=3, max=6),0))  
  Iv <- round(runif(N, min=min_Iv, max=max_Iv),0)
  
  # daily time step for 2 years
  timestep <- seq(0, 2*365) 
  
  # assuming that 20% of infections result in reported cases 
  Reporting_fraction <- 0.2
  
  # dataframes to store simulations
  sim_data <- data.frame(matrix(ncol=N, nrow=105))
  colnames(sim_data) <- paste("sim", seq(1:N), sep="")
  sim_data_P <- sim_data
  sim_data_NB <- sim_data
  
  # store parameter values
  pars_data <- data.frame(sim=1, R0=R0v[1], Np=Npv[1], I0=Iv[1])
  
  
  
  #=== Simulate SEIR model & aggregate to weekly incidence data
  for(i in 1:N){
    
    # define parameters & initial states
    pars <- c(R0=R0v[i], gamma=gamma, sigma=sigma)
    state_init <- c(S=Npv[i], E=0, I=Iv[i], R=0) 
    
    # run the SEIR model
    out <- ode(y=state_init, times=timestep, func=model, parms=pars, method="rk4") 
    
    # extract the incidence of infections
    incidence <- round(pars[which(names(pars)=="sigma")] * out[,which(colnames(out)=="E")]) 
    
    # aggregate to weekly data
    Inc_weekly <- round(unname(tapply(incidence, (seq_along(incidence)-1) %/% n_aggr, sum))* Reporting_fraction) 
    
    # cut the time series to begin at time of 1st case
    indx <- which(Inc_weekly>0)
    Inc_weekly <- Inc_weekly[indx[1]:indx[length(indx)]]
    
    # add reporting noise to the data (after week 1)
    inc_P <- c(Inc_weekly[1], add_noiseP(Inc_weekly[2:length(Inc_weekly)]))
    inc_NB <- c(Inc_weekly[1], add_noiseNB(Inc_weekly[2:length(Inc_weekly)]))
    
    # store weekly incidence data and true parameter values
    sim_data[1:length(Inc_weekly) ,i] <- Inc_weekly
    sim_data_P[1:length(Inc_weekly),i] <- inc_P 
    sim_data_NB[1:length(Inc_weekly),i] <- inc_NB
    pars_data[i, ] <- c(i, R0v[i], Npv[i], Iv[i])
    
  }
  
  # return a list of simulated data, with varying noise levels, & parameter values used
  return(list(sims=sim_data, simsP=sim_data_P, simsNB=sim_data_NB, pars=pars_data))
}
  

# same for a stochastic model:
simulate_data_stoch <- function(N, gamma, sigma, lower_Npv = 3, upper_Npv = 6, n_aggr = 7,
                                min_Iv = 3, max_Iv = 8){
  
  #=== Define the SEIR model
  
  stochastic_sim <- function(gamma, sigma, R0, Npv, Iv, lgt = 2*365){
    S <- E <- I <- R <- inc <- rep(NA, lgt)
    S[1] <- Npv - Iv
    E[1] <- 0
    I[1] <- Iv
    R[1] <- 0
    inc[1] <- 0
    
    for(t in 2:lgt){
      S_to_E <- rbinom(1, S[t - 1], R0*gamma*I[t - 1]/Npv)
      E_to_I <- rbinom(1, E[t - 1], sigma)
      I_to_R <- rbinom(1, I[t - 1], gamma)
      
      S[t] <- S[t - 1] - S_to_E
      E[t] <- E[t - 1] + S_to_E - E_to_I
      I[t] <- I[t - 1] + E_to_I - I_to_R
      R[t] <- R[t - 1] + I_to_R
      
      inc[t] <- E_to_I
    }
    
    # reporting:
    reported_inc <- rbinom(lgt, inc, 0.2)
    
    return(reported_inc)
  }
  
  
  
  #=== Generate N simulated epidemic curves with randomly selected parameter values
  
  # randomly select values of R0, population size (Np), and inital number infected (I)
  R0v <- round(rtruncnorm(N, a=2, b=8, mean=3, sd=2),2) 
  # JB: set up population size by a factor of 10^3
  Npv <- if(lower_Npv == upper_Npv) rep(10^lower_Npv, N) else 10^(round(runif(N, min=lower_Npv, max=upper_Npv),0))
  #Npv <- 10^(round(runif(N, min=3, max=6),0))  
  Iv <- round(runif(N, min=min_Iv, max=max_Iv),0)
  
  
  # dataframes to store simulations
  sim_data <- data.frame(matrix(ncol=N, nrow=105))
  colnames(sim_data) <- paste("sim", seq(1:N), sep="")
  
  # store parameter values
  pars_data <- data.frame(sim=1, R0=R0v[1], Np=Npv[1], I0=Iv[1])
  
  #=== Simulate SEIR model & aggregate to weekly incidence data
  for(i in 1:N){
    
    # run the SEIR model
    incidence <- stochastic_sim(gamma = gamma, sigma = sigma, R0 = R0v[i], Npv = Npv[i], Iv = Iv[i], lgt = 150)
    
    # aggregate to weekly data
    Inc_weekly <- round(unname(tapply(incidence, (seq_along(incidence)-1) %/% n_aggr, sum))) 
    
    # cut the time series to begin at time of 1st case
    first_index_cases <- min(which(Inc_weekly>0))
    Inc_weekly <- Inc_weekly[first_index_cases:length(Inc_weekly)]
    Inc_weekly <- Inc_weekly[1:105]
    
    # store weekly incidence data and true parameter values
    sim_data[ ,i] <- Inc_weekly
    pars_data[i, ] <- c(i, R0v[i], Npv[i], Iv[i])
  }
  
  # return a list of simulated data, with varying noise levels, & parameter values used
  return(list(sims=sim_data, pars=pars_data))
}

