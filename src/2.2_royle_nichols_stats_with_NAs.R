royle_nichols_stats <- function(sim_nr) {
  # We need to find out if camera traps with fewer survey days than the 
  # study duration need to be disregarded, or simply filled up with NAs. 
  # To do this, we are going to focus on the difference in abundance 
  # between dens 8 and dens 11. The 25 cameras per density will each get
  # a different sampling period. z- and p-values of the RN-models will
  # be recorded. 
  
  # with the presence/absence matrices per density, 
  # we can use the RN occupany model to examine whether 
  # it provides a significant difference between pairs of 
  # simulations with different densities. 
  
  # We can do this with different parameter settings:
  # - density-combinations (also same densties against each 
  # other to examine false significant effects)
  # - time intervals
  # - study durations
  # - number of camera traps
  # do this for 10 simulation runs, using 10 parallel sessions...
  library(unmarked)
  
  x <- vector(length=0)
  data <- data.frame(dens1 = x, 
                     dens2 = x, 
                     time_interval = x, 
                     original_study_duration = x,
                     study_duration = x, 
                     n_cams = x, 
                     z = x, 
                     P = x, 
                     type = x)
  outputfile <- paste('./results/output/simulations/royle_nichols_stats_with_NAs', 
                      sim_nr, '.txt')
  if (!file.exists(outputfile)) { 
    write.table(data, outputfile, col.names = TRUE, row.names = FALSE)
  }
  
  i1 <- 9
  i2 <- 13

  df1 <- as.matrix(read.table(paste(
    './data/raw/IBM results/PresAbs', i1, sim_nr, '.txt')))
  dens1 <- df1[1,1] 
  df1 <- df1[,2:366]
  
  df2 <- as.matrix(read.table(paste(
    './data/raw/IBM results/PresAbs', i2, sim_nr, '.txt')))
  dens2 <- df2[1,1] 
  df2 <- df2[,2:366]
  
  cov <- data.frame(location=factor(c(rep('A', 25), rep('B', 25))),
                    density=c(rep(dens1, 25), rep(dens2, 25)))
  dets <- rbind(df1, df2) 
  
  n_cams <- 25
  
  # generate different sampling periods per camera:
  # what if the sampling periods are shorter in location A than in B? 
  n_sampling_days <- round(c(runif(25, 60, 365), runif(25, 30, 150)))
  for (ni in 1:50){
    dets[ni, n_sampling_days[ni]:365] = NA
  }
  
  time_interval <- 1
  for (study_duration in 10:200) {
    stu_dur <- study_duration
    dets2 <- dets[,1:stu_dur]
    
    # without NAs:
    dets3 <- dets2[!is.na(rowMeans(dets2)),]
    cov3 <- cov[!is.na(rowMeans(dets2)),]
    
    # first the RN model with NAs:
    umf <- unmarkedFrameOccu(y=dets2, siteCovs=cov)
    m1 <- occuRN(~1 ~location, umf)
    s1 <- summary(m1)
    
    data <- data.frame(dens1 = dens1, 
                       dens2 = dens2, 
                       time_interval = time_interval, 
                       original_study_duration = study_duration, 
                       study_duration = stu_dur, 
                       n_cams = n_cams, 
                       z = s1$state$z[2], 
                       P = s1$state$`P(>|z|)`[2], 
                       type = 'With NAs')
    write.table(data, outputfile, 
                append=TRUE, col.names = FALSE, row.names = FALSE)

    # second the RN model without NAs:
    if (length(unique(cov3$location)) == 2){
      umf <- unmarkedFrameOccu(y=dets3, siteCovs=cov3)
      m1 <- occuRN(~1 ~location, umf)
      s1 <- summary(m1)
      
      data <- data.frame(dens1 = dens1, 
                         dens2 = dens2, 
                         time_interval = time_interval, 
                         original_study_duration = study_duration, 
                         study_duration = stu_dur, 
                         n_cams = n_cams, 
                         z = s1$state$z[2], 
                         P = s1$state$`P(>|z|)`[2], 
                         type = 'Without NAs')
      write.table(data, outputfile, 
                  append=TRUE, col.names = FALSE, row.names = FALSE)
    }
  }
}

# run in parallel:
library(parallel)
cl   = makeCluster(10)
results = parSapply(cl, 1:10, royle_nichols_stats)
stopCluster(cl)

