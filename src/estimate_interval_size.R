estimate_interval_size <- function(sampling_effort, Ppresence, max_interval_size){
  # this function estimates the interval size that should be used 
  # in the Royle-Nichols occupancy model, given 
  # 1. the sampling effort (= total camera days)
  # 2. the proportion of cameras with detections of the species
  # 3. the maximum interval size that can be used
  
  x1 <- 1:max_interval_size 
  x2 <- sampling_effort
  x3 <- log(Ppresence/(1 - Ppresence))
  
  b <- read.table(
    './results/output/simulations/to_estimate_z_values.txt', 
    header=TRUE)$Estimate
  
  z <- b[1] + 
    x1*b[2] + 
    x2*b[3] + 
    x3*b[4] + 
    
    x1^2*b[5] + 
    x2^2*b[6] + 
    x3^2*b[7] + 
    
    x1^3*b[8] + 
    x2^3*b[9] + 
    x3^3*b[10] + 
    
    x1^4*b[11] + 
    x2^4*b[12] + 
    x3^4*b[13] + 
    
    x1*x2*b[14] +
    x1*x3*b[15] +
    x2*x3*b[16] +

    x1^2*x3^2*b[17] + 
    x2^2*x3^2*b[18] + 
    
    x1^3*x2^3*b[19] +
    x1^3*x3^3*b[20] +
    x2^3*x3^3*b[21] +
    
    x1^4*x2^4*b[22] +
    x1^4*x3^4*b[23] +
    x2^4*x3^4*b[24] +
    
    x1*x2*x3*b[25] + 
    x1^2*x2^2*x3^2*b[26] + 
    x1^3*x2^3*x3^3*b[27] + 
    x1^4*x2^4*x3^4*b[28]
    
  plot(x1, z)
  dT <- x1[z == max(z)]
  
  return(dT)
}

