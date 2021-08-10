#=== Function to estimate time-constant R0 with each method at sequential ===#
#=== time points in the epidemic growth phase and output results as csv's ===#

# JB: made methods an argument
fit_R0_seq <- function(data, mean_GT, sd_GT, GTd, GT_week, min_peak_time = 15,
                       methods = c("EG_Lin", "EG_P", "EG_MLE", "EpiEstim", "WP", "WT")){
  
  # list to store results
  store <- list()
  
  # method names
  # EB: excluded to BR
  # methods <- c("EG_Lin", "EG_P", "EG_MLE", "EpiEstim", "WP", "WT")
  
  # define peak - max time point of highest reported cases
  peak <- which.max(data$cases)[length(which.max(data$cases))] 
  
  if(peak>=min_peak_time){
    
    # define sections: from 2 generation times (approximated in terms of weeks) on, up to peak
    sections <- seq(from=GT_week*2, to=peak, by=GT_week)
    
    # Loop to fit each method to each section with increasing number of time points
    # EB: only run over sections 6, 9, 12, 15
    for(t in 1:4){
      
      # extract section of the epidemic curve for fitting
      # JB: add leading zeros to avoid cutting of serial interval distribution
      s_t <- data.frame(week = 1:(sections[t] + 50),
                        cases = c(rep(0, 50), data$cases[1:sections[t]]),
                        country = data$country[1])
      
      # JB: fit selected methods
      # EB: excluded to BR
      if("EG_Lin" %in% methods) eglin <- EG_Lin(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      if("EG_P" %in% methods) egp <- EG_P(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      if("EG_MLE" %in% methods) egmle <- EG_MLE(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      if("EpiEstim" %in% methods) epiest <- epiestim(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      if("WP" %in% methods) wp <- WP(data=s_t, GTd=GTd)
      if("WT" %in% methods) wt <- WT(data=s_t, GTd=GTd)
      if("BR" %in% methods) br <- BR(data=s_t, GTd=GTd)
      
      # store results
      # JB: changed to data.frame rather than matrix (not sure why this problem occurs)
      #store[[t]] <- cbind(paste(data$country[1]), sections[t], methods, rbind(eglin, egp, egmle, epiest, wp, wt, br), peak)
      store[[t]] <- data.frame(paste(data$country[1]), sections[t], methods, 
                               rbind(if("EG_Lin" %in% methods) eglin,
                                     if("EG_P" %in% methods) egp,
                                     if("EG_MLE" %in% methods) egmle,
                                     if("EpiEstim" %in% methods) epiest,
                                     if("WP" %in% methods) wp,
                                     if("WT" %in% methods) wt,
                                     if("BR" %in% methods) br), peak)
      colnames(store[[t]]) <- c("Country", "Nweeks", "method", "R0", "CI_L", "CI_U", "peak")
    }
    
    # bind results by country and output as csv
    results_all <- do.call("rbind", store)
    return(results_all)
  }
}
