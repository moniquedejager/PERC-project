royle_nichols_stats <- function (sim_nr) {
  #' This function estimates relative occupancies between two locations
  #' and writes the z- and p-values to a file
  
  #' With the presence/absence matrices per density, 
  #' we can use the Royle-Nichols occupany model to 
  #' examine whether it provides a significant difference between
  #' pairs of simulations with different densities.
  
  library(unmarked)
  source('./src/functions/get_presence_absence_data.R')
  source('./src/functions/create_outputfile.R')
  source('./src/functions/adjust_matrices.R')
  create_outputfile(sim_nr)
  
  values_all <- read.table('./data/temp/parameter_combinations_royle_nichols_stats.txt', header=TRUE)
  
  for (i in 1:length(values_all$i1)){
    
    # Get the parameter values from the table.
    values <- values_all[i,]
    study_duration <- values$study_duration
    time_interval <- values$time_interval
    n_cams <- values$n_cams
    i1 <- values$i1
    i2 <- values$i2
    
    # get the presence/absence matrices for the given simulation runs. 
    data <- get_pres_abs_data(i1, i2, sim_nr)
    cov <- data[,1:2]
    dets <- data[,3:367]
      
    # We need to adjust the detection and covariates matrices to the
    # parameters used in the i'th row of the values table
    data <- adjust_matrices(study_duration, time_interval, dets, cov)
    cov2 <- data[,1:2]
    dets3 <- data[,3:(length(data$location))]
    
    # What is the proportion of cameras with sightings?
    p_presence <- mean(rowSums(dets3) > 0)
    
    # This code does the actual analysis on the adjusted matrices:
    umf <- unmarkedFrameOccu(y = dets3, siteCovs = cov2)
    m1 <- occuRN(~1 ~ location, umf)
    s1 <- summary(m1)
       
    # The data needs to be written to an output file,
    # so we can use this in the next analyses. 
    data <- data.frame(
      dens1 = dens1,
      dens2 = dens2,
      time_interval = 1,
      study_duration = study_duration,
      n_cams = n_cams,
      z = s1$state$z[2],
      P = s1$state$`P(>|z|)`[2],
      p_presence = p_presence)
    
    write.table(
      data, outputfile,
      append = TRUE,
      col.names = FALSE,
      row.names = FALSE)
  }
}