estimate_interval_size <- function(sampling_effort, Ppresence){
  # this function estimates the interval size that should be used 
  # in the Royle-Nichols occupancy model, given 
  # 1. the sampling effort (= total camera days)
  # 2. the proportion of cameras with detections of the species
  
  x1 <- sampling_effort
  x2 <- log(Ppresence/(1 - Ppresence))
  
  b <- read.table(
    './results/output/simulations/to_estimate_time_interval_sizes_II.txt', 
    header=TRUE)$Estimate
  
  y <- b[1] + 
    x1*b[2] + 
    x2*b[3] +  
    x1^2*b[4] + 
    x2^2*b[5] + 
    x1^3*b[6] + 
    x2^3*b[7] + 
    x1^4*b[8] + 
    x2^4*b[9] + 
    x1^5*b[10] + 
    x2^5*b[11] + 
    x1^6*b[12] + 
    x2^6*b[13] + 
    x1*x2*b[14] +
    x1^2*x2^2*b[15] + 
    x1^3*x2^3*b[16] +
    x1^4*x2^4*b[17] + 
    x1^5*x2^5*b[18] +
    x1^6*x2^6*b[19]
  dT <- ceiling(exp(y))
  
  return(dT)
}

