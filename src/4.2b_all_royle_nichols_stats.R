# Calculate the royle-nichols stats per species, study duration, 
# and time interval (dT = 5 of optimal time interval):

royle_nichols_stats = function(i){
  library(unmarked)
  library(ggplot2)
  library(ggpubr)
  
  # i is the species number:
  spec <- read.table('./data/processed/FSC and nonFSC data/species.txt')$V1[i]
  filename <- paste('./data/processed/FSC and nonFSC data/pres_abs', 
                    spec, '.txt')
  dets <- as.matrix(read.table(filename))
  cov <- read.table('./data/raw/FSC and nonFSC data/EffortPerCam.csv', sep=';', header=T)
  cov$Cluster <- factor(cov$Cluster)
  cov$Visibility <- 3 - cov$Visibility
  
  # What is the survey effort per camera?
  dets2 <- dets
  dets2[dets2 == 0] <- 1
  dets2[is.na(dets2)] <- 0
  survey_effort_per_cam <- rowSums(dets2)
  
  outputfile <- paste('./results/output/FSC and nonFSC data/all_royle_nichols_stats', 
                      spec, '.txt')
  x <- vector(length=0)
  data <- data.frame(species = x, 
                     time_interval = x, 
                     min_number_of_intervals = x,
                     survey_effort = x, 
                     n_cams = x, 
                     z = x, 
                     P = x, 
                     p_presence=x)
  write.table(data, outputfile, 
              append=FALSE, col.names = TRUE, row.names = FALSE)
  
  # try out a range of time interval sizes and minimum numbers of intervals:
  for (nT in unique(round(exp(3:12*0.29)))){
    max_interval_size <- 
      quantile(survey_effort_per_cam, 0.95)/nT
    all_interval_sizes <- unique(round(exp((1:24*0.196))))
    
    for (time_interval in all_interval_sizes[all_interval_sizes <= max_interval_size]){
      # create new dets-matrix with the right time intervals:
      n_intervals <- floor(max(survey_effort_per_cam)/time_interval)
      stu_dur <- n_intervals*time_interval
      dets2 <- dets[,1:stu_dur]
      dets3 <- matrix(0, length(survey_effort_per_cam), n_intervals)
      x <- sort(rep(1:n_intervals, time_interval))
      for (i3 in 1:length(survey_effort_per_cam)) {
        dets3[i3,] <- tapply(dets2[i3,], x, max)
      }
  
      min_study_duration <- time_interval*nT  
  
      # we need to remove all camera traps with effort < study_duration:
      dets4 <- dets3[cov$Effort >= min_study_duration,]
      cov2 <- cov[cov$Effort >= min_study_duration,]
  
      survey_effort <- length(dets4[!is.na(dets4)])*time_interval
      dets5 <- dets4
      dets5[is.na(dets5)] <- 0
      Ppresence <- sum(rowSums(dets5) > 0) / length(dets5[,1])
      n_cams <- length(dets5[,1])
  
      umf <- unmarkedFrameOccu(y=dets4, siteCovs=cov2)
  
      if (length(unique(cov$Cluster)) == 1) {
        mod <- occuRN(~Visibility - 1 ~SiteType, umf)
      } else {
        mod <- occuRN(~Visibility - 1 ~Cluster + SiteType, umf)
      }
  
      s1 <- summary(mod)
  
      # sometimes the model has too many covariates and cannot compute all 
      # the parameters. We should then fall back on a simpler model: 
      if (is.na(s1$state[length(s1$state$Estimate),3])){
        mod <- occuRN(~Visibility - 1 ~SiteType, umf)
        s1 <- summary(mod)
      }
      
      data <- data.frame(species = spec, 
                         time_interval = time_interval, 
                         min_number_of_intervals = nT,
                         survey_effort = survey_effort, 
                         n_cams = n_cams, 
                         z = s1$state[length(s1$state$Estimate),3], 
                         P = s1$state[length(s1$state$Estimate),4], 
                         p_presence=Ppresence,
                         type='optimal dT')
      write.table(data, outputfile, 
                  append=TRUE, col.names = FALSE, row.names = FALSE)
    }
  }
}


i = 1:35 # there are 35 species that we want to analyse

# run in parallel:
library(parallel)
cl   = makeCluster(10)
results = parSapply(cl, i, royle_nichols_stats)
stopCluster(cl)
