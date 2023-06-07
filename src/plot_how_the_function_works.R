# create an illustration of how the interval size is estimated 
# from the presence/absence matrix: 

pres_abs_matrix <- read.table('./data/processed/PresenceAbsence/PresAbs 9 10 .txt')
pres_abs_matrix <- pres_abs_matrix[,2:366]
pres_abs_matrix <- as.matrix(pres_abs_matrix)

pres_abs_matrix2 <- read.table('./data/processed/PresenceAbsence/PresAbs 11 10 .txt')
pres_abs_matrix2 <- pres_abs_matrix2[,2:366]
pres_abs_matrix2 <- as.matrix(pres_abs_matrix2)

pres_abs_matrix <- rbind(pres_abs_matrix, pres_abs_matrix2)

# add NA's randomly to the end of the matrix to make it look like real data:
#NA_start <- round(rnorm(50, 300, 30))
#for (i in 1:50){
#  if (NA_start[i] <= 365){
#    pres_abs_matrix[i,NA_start[i]:365] <- NA 
#  }
#}

library(ggplot2)
library(ggpubr)

# What is the survey effort per camera?
dets <- pres_abs_matrix
dets2 <- pres_abs_matrix
dets2[dets2 == 0] <- 1
dets2[is.na(dets2)] <- 0
survey_effort_per_cam <- rowSums(dets2)

# depending on the chosen maximum interval size,
# the number of cameras used in the analysis changes,
# and thus the survey effort and Ppresence change as well. 
# Thus, for a range of interval sizes, 
# we can estimate the z-value, and find out which 
# one we should use. 

x <- vector(length=0)
df <- data.frame(survey_effort = x,
                 Ppresence = x,
                 interval_size = x,
                 n_cams = x,
                 min_number_of_intervals = x,
                 z = x)

for (min_number_of_intervals in 2:30){
  max_interval_size <- 
    quantile(survey_effort_per_cam, 0.95)/min_number_of_intervals
  max_interval_size <- min(c(max_interval_size, 30))
  
  for (j in 1:max_interval_size){
    dets3 <- dets[survey_effort_per_cam >= (j*min_number_of_intervals),]
    survey_effort <- length(dets3[!is.na(dets3)])
    dets3[is.na(dets3)] <- 0
    Ppresence <- sum(rowSums(dets3) > 0) / length(dets3[,1])
    n_cams <- length(dets3[,1])
    
    df2 <- data.frame(survey_effort = survey_effort,
                      Ppresence = Ppresence,
                      interval_size = j,
                      n_cams = n_cams,
                      min_number_of_intervals = min_number_of_intervals,
                      z = NA)
    df <- rbind(df, df2)
  }
}

df$Ppresence[df$Ppresence == 1] <- 0.999
df <- df[df$Ppresence > 0,]

# 1. the sampling effort (= total camera days)
# 2. the proportion of cameras with detections of the species
# 3. the maximum interval size that can be used

x1 <- df$interval_size 
x2 <- 1/df$survey_effort
x3 <- log(df$Ppresence/(1 - df$Ppresence))
x4 <- 1/df$n_cams

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

df$z <- z
dT <- df$interval_size[df$z == max(df$z)]
nT <- df$min_number_of_intervals[df$z == max(df$z)]
#df$min_number_of_intervals <- factor(df$min_number_of_intervals)

windows(height=5, width=5)
ggplot(df, aes(x=interval_size, y=min_number_of_intervals, fill=1/z)) + 
  geom_raster() + 
  scale_fill_continuous(type='viridis', name='Model fit', trans='reverse', 
                        breaks=c(min(1/df$z),max(1/df$z)),labels=c("Best fit","Worst fit")) + 
  xlab('Interval size (days, dT)') + 
  ylab('Minimum number of intervals (nT)') + 
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
dT
nT

df$percentile_Fit <- rank(2/df$z) 

p2 <- ggplot(df, aes(x=interval_size, y=min_number_of_intervals, fill=percentile_Fit)) + 
  geom_raster() + 
  scale_fill_continuous(type='viridis', name='Model fit', trans='reverse',
                        breaks=c(0,0.5,1),labels=c("Minimum",0.5,"Maximum")) + 
  xlab('Interval size (days, dT)') + 
  ylab('Minimum number of intervals (nT)') + 
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
