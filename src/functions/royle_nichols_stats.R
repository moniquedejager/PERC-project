# this function estimates relative occupancies between two locations
# and writes the z- and p-values to a file:
royle_nichols_stats <- function(sim_nr) {
  library(unmarked)
  
  # 1. create an outputfile:
  source('./src/functions/create_outputfile.R')
  
  # get the parameter values from the table:
  values_all <- read.table('./data/temp/parameter_combinations_royle_nichols_stats.txt', header=TRUE)
  
  # for each parameter value combination in values_all, 
  # do the analysis:
  for (i in 1:length(values_all$d1)){
    # 2. Get the parameter values and the presence/absence data:
    source('./src/functions/get_param_values.R')
    source('./src/functions/get_presence_absence_data.R')
      
    # 3. adjust the detection and covariates matrices:
    source('./src/functions/adjust_matrices.R')
    
    # 4. What is the proportion of cameras with sightings?
    p_presence <- mean(rowSums(dets3) > 0)
    
    # 5. Do the actual analysis:
    umf <- unmarkedFrameOccu(y = dets3, siteCovs = cov2)
    m1 <- occuRN(~1 ~ location, umf)
    s1 <- summary(m1)
       
    # 6. Write the results to the outputfile:          
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