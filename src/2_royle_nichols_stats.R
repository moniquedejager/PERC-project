royle_nichols_stats <- function(sim_nr) {
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
                     p_presence = x)
  outputfile <- paste('./results/PresenceAbsence/royle_nichols_stats', 
                      sim_nr, '.txt')
  if (!file.exists(outputfile)) { 
    write.table(data, outputfile, col.names = TRUE, row.names = FALSE)
  }
  
  for (i1 in 1:18) 
  {
    for (i2 in (i1+1):19)
    {
      print(paste(i1, i2))
      
      df1 <- as.matrix(read.table(paste(
        './data/raw/IBM results/presence_absence', i1, sim_nr, '.txt')))
      dens1 <- df1[1,1] 
      df1 <- df1[,2:366]
      
      df2 <- as.matrix(read.table(paste(
        './results/PresenceAbsence/PresAbs', i2, sim_nr, '.txt')))
      dens2 <- df2[1,1] 
      df2 <- df2[,2:366]
      
      cov <- data.frame(location=factor(c(rep('A', 25), rep('B', 25))),
                        density=c(rep(dens1, 25), rep(dens2, 25)))
      dets <- rbind(df1, df2) 
      
      for (n_cams in 25:5) {
        for (study_duration in round(exp(log(365)*(10:30/30)))) {
          for (time_interval in unique(round(1.1^(4:54)))) {
            if (time_interval <= study_duration/2){
              stu_dur <- study_duration
              n_intervals <- floor(stu_dur/time_interval)
              stu_dur <- n_intervals*time_interval
              
              dets2 <- dets[,1:stu_dur]
              dets3 <- matrix(0, 50, n_intervals)
              # voor iedere rij in dets2, 
              # presence/absence samenvoegen per tijdsinterval...
              x <- sort(rep(1:n_intervals, time_interval))
              for (i3 in 1:50) {
                dets3[i3,] = tapply(dets2[i3,], x, max)
              }
              
              # randomly select cameras to use...
              cams_used <- sort(rank(runif(25))[1:n_cams]) 
              dets4 <- dets3[c(cams_used, cams_used+25),]
              cov2 <- cov[c(cams_used, cams_used+25),]
              
              p_presence <- mean(rowSums(dets4) > 0)
              
              umf <- unmarkedFrameOccu(y=dets4, siteCovs=cov2)
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
                                 p_presence=p_presence)
              write.table(data, outputfile, 
                          append=TRUE, col.names = FALSE, row.names = FALSE)
            }
          }
        }
      }
    }
  }
}

# run in parallel:
library(parallel)
cl   = makeCluster(10)
results = parSapply(cl, 1:10, royle_nichols_stats)
stopCluster(cl)

