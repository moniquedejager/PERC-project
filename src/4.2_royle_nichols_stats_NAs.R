# this code is used to get rid of as many NAs as possible. 
# it should not be used later on (only once) as the previous
# script has now been adjusted. 

# Calculate the royle-nichols stats per species, study duration, 
# and time interval:

royle_nichols_stats = function(i){
  library(unmarked)
  
  # i is the species number:
  spec <- read.table('./data/processed/FSC and nonFSC data/species.txt')$V1[i]
  
  filename <- paste('./data/processed/FSC and nonFSC data/pres_abs', 
                    spec, '.txt')
  dets <- as.matrix(read.table(filename))
  cov <- read.table('./data/raw/FSC and nonFSC data/EffortPerCam.csv', sep=';', header=T)
  cov$Cluster <- factor(cov$Cluster)
  cov$Visibility = 3 - cov$Visibility
  
  # the output file:
  outputfile <- paste('./results/output/FSC and nonFSC data/royle_nichols_stats', 
                      spec, '.txt')
  outfile <- read.table(outputfile, header=T)
  
  # list the row numbers of all NAs so we can adjust this data:
  ix <- 1:length(outfile$species)
  ix <- ix[is.na(outfile$z)]
  
  for (i2 in ix){
    time_interval <- outfile$time_interval[i2]
    n_intervals <- floor(200/time_interval)
    stu_dur <- n_intervals*time_interval
    dets2 <- dets[,1:stu_dur]
    dets3 <- matrix(0, 455, n_intervals)
    x <- sort(rep(1:n_intervals, time_interval))
    for (i3 in 1:455) {
      dets3[i3,] <- tapply(dets2[i3,], x, max)
    }
    study_duration <- outfile$study_duration[i2]
    n_intervals <- (study_duration/time_interval)
    dets4 <- dets3[cov$Effort >= study_duration,1:n_intervals]
    cov2 <- cov[cov$Effort >= study_duration,]
    
    p_presence <- mean(rowSums(dets4) > 0)
    umf <- unmarkedFrameOccu(y=dets4, siteCovs=cov2)

    mod <- occuRN(~Visibility - 1 ~SiteType, umf)
    s1 <- summary(mod)
    
    data <- data.frame(species = spec, 
                       time_interval = time_interval, 
                       study_duration = study_duration, 
                       n_cams = length(cov2$Cam), 
                       z = s1$state[length(s1$state$Estimate),3], 
                       P = s1$state[length(s1$state$Estimate),4], 
                       p_presence=p_presence)

    outfile[i2,] <- data
  }
  
  write.table(outfile, outputfile, 
                  append=FALSE, col.names = TRUE, row.names = FALSE)
}

i = 1:35 # there are 35 species that we want to analyse

# run in parallel:
library(parallel)
cl   = makeCluster(10)
results = parSapply(cl, i, royle_nichols_stats)
stopCluster(cl)
