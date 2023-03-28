estimate_interval_size <- function(sampling_effort, Ppresence, max_interval_size, n_cams){
  # this function estimates the interval size that should be used 
  # in the Royle-Nichols occupancy model, given 
  # 1. the sampling effort (= total camera days)
  # 2. the proportion of cameras with detections of the species
  # 3. the maximum interval size that can be used
  
  x1 <- 1:max_interval_size 
  x2 <- 1/sampling_effort
  x3 <- log(Ppresence/(1 - Ppresence))
  x4 <- 1/n_cams
  
  b <- read.table(
    './results/output/simulations/to_estimate_z_values.txt', 
    header=TRUE)$Estimate
  
  z <- b[1] + 
    x1*b[2] + 
    x2*b[3] + 
    x3*b[4] + 
    x4*b[5] + 
    
    x1^2*b[6] + 
    x3^2*b[7] + 
    
    x1^3*b[8] + 
    x3^3*b[9] + 
    
    x1*x2*b[10] +
    x1*x3*b[11] +
    x2*x3*b[12] +
    x1*x4*b[13] +
    x2*x4*b[14] +
    x3*x4*b[15] +
    
    x2*x1^2*b[16] + 
    x1^2*x3^2*b[17] + 
    x2*x3^2*b[18] + 
    x4*x1^2*b[19] + 
    x4*x3^2*b[20] +
    
    x2*x1^3*b[21] + 
    x1^3*x3^3*b[22] +
    x2*x3^3*b[23] +
    x4*x1^3*b[24] +
    x4*x3^3*b[25] +
    
    x1*x2*x3*b[26] + 
    x1*x2*x4*b[27] + 
    x1*x3*x4*b[28] +
    x2*x3*x4*b[29] + 
    
    x2*x1^2*x3^2*b[30] + 
    x2*x4*x1^2*b[31] + 
    x4*x1^2*x3^2*b[32] + 
    x2*x4*x3^2*b[33] + 
    
    x2*x1^3*x3^3*b[34] + 
    x2*x4*x1^3*b[35] + 
    x4*x1^3*x3^3*b[36] + 
    x2*x4*x3^3*b[37] + 
    
    x1*x2*x3*x4*b[38] + 
    x2*x4*x1^2*x3^2*b[39] + 
    x2*x4*x1^3*x3^3*b[40]
    
  plot(x1, z)
  dT <- x1[z == max(z)]
  maxz <- max(z)
  
  return(c(dT, maxz))
}


