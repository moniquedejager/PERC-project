# Calculate the royle-nichols stats per species, study duration, 
# and time interval (dT = 5 of optimal time interval):

royle_nichols_stats = function(i){
  library(unmarked)
  library(ggplot2)
  library(ggpubr)
  source("./src/estimate_interval_size.R")
  
  # i is the species number:
  spec <- read.table('./data/processed/FSC and nonFSC data/species.txt')$V1[i]
  
  filename <- paste('./data/processed/FSC and nonFSC data/pres_abs', 
                    spec, '.txt')
  dets <- as.matrix(read.table(filename))
  cov <- read.table('./data/raw/FSC and nonFSC data/EffortPerCam.csv', sep=';', header=T)
  cov$Cluster <- factor(cov$Cluster)
  cov$Visibility = 3 - cov$Visibility
  
  # What is the survey effort per camera?
  dets2 <- dets
  dets2[dets2 == 0] <- 1
  dets2[is.na(dets2)] <- 0
  survey_effort_per_cam <- rowSums(dets2)
  
  # depending on the chosen maximum interval size,
  # the number of cameras used in the analysis changes,
  # and thus the survey effort and Ppresence change as well. 
  # Thus, for a range of maximum interval sizes, 
  # we can estimate the optimal interval size, and find out which 
  # one we should use:
  # start at a maximum time interval size that is half of the median
  # survey effort per camera.
  
  df <- data.frame(max_interval_size = 
                     round(quantile(survey_effort_per_cam, 0.05)/2):
                     round(quantile(survey_effort_per_cam, 0.75)/2),
                   survey_effort = 0,
                   Ppresence = 0,
                   optimal_interval_size = 0,
                   z = 0)
  
  for (j in df$max_interval_size){
    dets3 <- dets[survey_effort_per_cam >= (j*2),]
    survey_effort <- length(dets3[!is.na(dets3)])
    df$survey_effort[df$max_interval_size == j] <- survey_effort
    dets3[is.na(dets3)] <- 0
    Ppresence <- sum(rowSums(dets3) > 0) / length(dets3[,1])
    df$Ppresence[df$max_interval_size == j] <- Ppresence
    n_cams <- length(dets3[,1])
      
    if ((Ppresence == 0)|(Ppresence == 1)){
      df$optimal_interval_size[df$max_interval_size == j] <- NA
      df$z[df$max_interval_size == j] <- NA
    }else{
      df$optimal_interval_size[df$max_interval_size == j] <- 
        estimate_interval_size(survey_effort, Ppresence, j, n_cams)[1]
      df$z[df$max_interval_size == j] <- 
        estimate_interval_size(survey_effort, Ppresence, j, n_cams)[2]
    }
  }
  
  df <- df[!is.na(df$optimal_interval_size),]
  optimal_interval_size <- df$optimal_interval_size[(df$survey_effort*df$z) == 
                                                      max(df$survey_effort*df$z)]
  coeff <- 1/4
  p1 <- ggplot(df, aes(x=max_interval_size, y=survey_effort/10000)) + 
    geom_point(color=' darkgreen') + 
    geom_line(color=' darkgreen') + 
    xlab('Maximum interval size (days)') + 
    ylab('Survey effort (10,000 days) ') +
    ggtitle(spec) + 
    geom_point(aes(x=max_interval_size, y=Ppresence/coeff), color=' darkblue') + 
    geom_line(aes(x=max_interval_size, y=Ppresence/coeff), color=' darkblue') +
    scale_y_continuous(sec.axis = 
                         sec_axis(~.*coeff, name="Proportion of cameras with detections")) +
    theme(plot.title = element_text(face='bold.italic'),
          axis.title.y.left = element_text(color='darkgreen'),
          axis.title.y.right = element_text(color='darkblue')) +
    geom_vline(xintercept=df$max_interval_size[df$z == max(df$z)], linetype='dashed')
  
  coeff2 <- 5
  p2 <- ggplot(df, aes(x=max_interval_size, y=z)) + 
    geom_point(color=' darkgreen') + 
    geom_line(color=' darkgreen') +
    xlab('Maximum interval size (days)') + 
    ylab('z-value') +
    geom_point(aes(x=max_interval_size, y=optimal_interval_size/coeff2), color=' darkblue') + 
    geom_line(aes(x=max_interval_size, y=optimal_interval_size/coeff2), color=' darkblue') +
    scale_y_continuous(sec.axis = 
                         sec_axis(~.*coeff2, name="Optimal interval size (days)")) +
    theme(plot.title = element_text(face='bold.italic'),
          axis.title.y.left = element_text(color='darkgreen'),
          axis.title.y.right = element_text(color='darkblue')) + 
    geom_vline(xintercept=df$max_interval_size[df$z == max(df$z)], linetype='dashed')
  
  tiff(filename=paste('./results/figures/FSC and nonFSC/optimal interval size', spec, '.tiff'), 
       height=5, width=5, units='in', res=300)
  plot(ggarrange(p1 + rremove("xlab"), p2, labels=c('A', 'B'), nrow = 2))
  dev.off()
  
  time_interval <- optimal_interval_size
  # create new dets-matrix with the right time intervals:
  n_intervals <- floor(200/time_interval)
  stu_dur <- n_intervals*time_interval
  dets2 <- dets[,1:stu_dur]
  dets3 <- matrix(0, 455, n_intervals)
  x <- sort(rep(1:n_intervals, time_interval))
  for (i3 in 1:455) {
    dets3[i3,] <- tapply(dets2[i3,], x, max)
  }
  
  study_duration <- optimal_interval_size*2  
    
  # we need to remove all camera traps with effort < study_duration:
  n_intervals <- (study_duration/time_interval)
  dets4 <- dets3[cov$Effort >= study_duration,]
  cov2 <- cov[cov$Effort >= study_duration,]
  
  Ppresence <- df$Ppresence[(df$survey_effort*df$z) == 
                             max(df$survey_effort*df$z)]
  survey_effort <- df$survey_effort[(df$survey_effort*df$z) == 
                                  max(df$survey_effort*df$z)]
  
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
      
  outputfile <- paste('./results/output/FSC and nonFSC data/royle_nichols_stats', 
                      spec, '.txt')
  
  data <- data.frame(species = spec, 
                     time_interval = time_interval, 
                     survey_effort = survey_effort, 
                     n_cams = length(cov2$Cam), 
                     z = s1$state[length(s1$state$Estimate),3], 
                     P = s1$state[length(s1$state$Estimate),4], 
                     p_presence=Ppresence,
                     type='optimal dT')
  write.table(data, outputfile, 
              append=FALSE, col.names = TRUE, row.names = FALSE)
  
  # with dT = 5:
  time_interval <- 5
  n_intervals <- floor(200/time_interval)
  stu_dur <- n_intervals*time_interval
  dets2 <- dets[,1:stu_dur]
  dets3 <- matrix(0, 455, n_intervals)
  x <- sort(rep(1:n_intervals, time_interval))
  for (i3 in 1:455) {
    dets3[i3,] <- tapply(dets2[i3,], x, max)
  }
  
  study_duration <- 10  
  
  # we need to remove all camera traps with effort < study_duration:
  n_intervals <- (study_duration/time_interval)
  dets4 <- dets3[cov$Effort >= study_duration,]
  cov2 <- cov[cov$Effort >= study_duration,]
  
  Ppresence <- df$Ppresence[df$max_interval_size == 10]
  survey_effort <- df$survey_effort[df$max_interval_size == 10]
  
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

  data2 <- data.frame(species = spec, 
                     time_interval = time_interval, 
                     survey_effort = survey_effort, 
                     n_cams = length(cov2$Cam), 
                     z = s1$state[length(s1$state$Estimate),3], 
                     P = s1$state[length(s1$state$Estimate),4], 
                     p_presence=Ppresence,
                     type='5 day interval')
  write.table(data2, outputfile, 
              append=TRUE, col.names = FALSE, row.names = FALSE)
}

royle_nichols_stats(1)

i = 1:35 # there are 35 species that we want to analyse

# run in parallel:
library(parallel)
cl   = makeCluster(10)
results = parSapply(cl, i, royle_nichols_stats)
stopCluster(cl)
