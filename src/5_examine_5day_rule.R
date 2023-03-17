# How far off are the RN model results when using dT = 5 days 
# instead of optimal interval size?

library(ggplot2)

df <- read.table(
  './results/output/simulations/RoyleNicholsStats20230307 1 .txt', 
  header=T)
df$original_SD <- read.table(
  './results/output/simulations/originalStudyDuration.txt')$V1
df$dT <- read.table('./results/output/simulations/timeInterval.txt')$V1
df$Ppres <- round(df$Ppresence, 1)
df$cols <- 0 + 1*(df$P < 0.05) + 1*((df$P < 0.05)&(df$z < 0))
df$cols2 <- factor(df$cols, labels = c("No significant effect", 
                                       "Correct result", 
                                       "Incorrect result"))
df$sim_nr <- 1
df <- df[df$dT == 5,]

for (i in 2:10){
  df2 <- read.table(paste(
    './results/output/simulations/RoyleNicholsStats20230307', i, '.txt'), 
    header=T)
  
  df2$original_SD <- read.table(
    './results/output/simulations/originalStudyDuration.txt')$V1
  df2$dT <- read.table('./results/output/simulations/timeInterval.txt')$V1
  df2$Ppres <- round(df2$Ppresence, 1)
  
  df2$cols <- 0 + 1*(df2$P < 0.05) + 1*((df2$P < 0.05)&(df2$z < 0))
  df2$cols2 <- factor(df2$cols, labels = c("No significant effect", 
                                           "Correct result", 
                                           "Incorrect result"))
  df2$sim_nr <- i
  df2 <- df2[df2$dT == 5,]
  df <- rbind(df, df2)
}

# per combination of survey effort and Ppresence, we calculate the 
# percentage of RN models that resulted in a significant and correct
# answer: 

df$sampling_effort <- df$StudyDuration * df$nCams
df$samp_eff2 <- round(log(df$sampling_effort), 1)
df$Ppres <- ceiling(df$Ppresence*20)/20
df$z[is.na(df$z)] <- 0

group <- paste(df$samp_eff2, df$Ppres, sep='-')
n_tot <- tapply(df$z, group, length)
n_correct <- tapply(df$z > 1.96, group, sum)
p_correct <- n_correct/n_tot *100
sampling_effort <- tapply(df$samp_eff2, group, mean)
Ppresence <- tapply(df$Ppres, group, mean)
df2 <- data.frame(p_correct = p_correct,
                  sampling_effort = sampling_effort,
                  Ppresence = Ppresence)

labs <- c(50, 150, 400, 1000, 3000, 10000)
windows(height=5, width=5)
p1 <- ggplot(df2, aes(x=sampling_effort, y=Ppresence, fill=p_correct)) +
  geom_raster() + 
  scale_x_continuous(breaks = log(labs), labels=labs) +
  scale_fill_continuous(type='viridis', name='% correct') + 
  xlab('Sampling effort (total # camera days)') + 
  ylab(expression(
    paste('Proportion of cameras with detections (', P[presence], ')'))) +
  theme(legend.position = "top") 

df3 <- read.table('./results/output/simulations/highest_z_values.txt', 
                  header=TRUE)
df3$sampling_effort <- df3$StudyDuration * df3$nCams
df3$samp_eff2 <- round(log(df3$sampling_effort), 1)
df3$Ppres <- ceiling(df3$Ppresence*20)/20
df3$z[is.na(df3$z)] <- 0

group <- paste(df3$samp_eff2, df3$Ppres, sep='-')
n_tot <- tapply(df3$z, group, length)
n_correct <- tapply(df3$z > 1.96, group, sum)
p_correct <- n_correct/n_tot *100
sampling_effort <- tapply(df3$samp_eff2, group, mean)
Ppresence <- tapply(df3$Ppres, group, mean)
df4 <- data.frame(p_correct = p_correct,
                  sampling_effort = sampling_effort,
                  Ppresence = Ppresence)

p2 <- ggplot(df4, aes(x=sampling_effort, y=Ppresence, fill=p_correct)) +
  geom_raster() + 
  scale_x_continuous(breaks = log(labs), labels=labs) +
  scale_fill_continuous(type='viridis', name='% correct') + 
  xlab('Sampling effort (total # camera days)') + 
  ylab(expression(
    paste('Proportion of cameras with detections (', P[presence], ')'))) +
  theme(legend.position = "top") 

df2$type <- 'Using 5 day interval'
df4$type <- 'Using optimal interval size'
df5 <- rbind(df2, df4)

ggplot(df5, aes(x=sampling_effort, y=Ppresence, fill=p_correct)) +
  geom_raster() + 
  scale_x_continuous(breaks = log(labs), labels=labs) +
  scale_fill_continuous(type='viridis', name='% correct') + 
  xlab('Sampling effort (total # camera days)') + 
  ylab(expression(
    paste('Proportion of cameras with detections (', P[presence], ')'))) +
  facet_wrap(vars(type))

# for each value in df2, find the corresponding value
# in df4 and calculate the difference:
group1 <- paste(df2$sampling_effort, df2$Ppresence, sep='-')
group2 <- paste(df4$sampling_effort, df4$Ppresence, sep='-')

for (i in unique(group1)){
  p1 <- df2$p_correct[group1 == i]
  p2 <- df4$p_correct[group2 == i]
  if (length(p2) == 0){
    p2 <- 0
    x <- data.frame(p_correct =0,
                    sampling_effort = df2$sampling_effort[group1 == i],
                    Ppresence = df2$Ppresence[group1 == i],
                    type = 'Using optimal interval size')
    df5 <- rbind(df5, x)
  }
  x <- data.frame(p_correct = p2 - p1,
                  sampling_effort = df2$sampling_effort[group1 == i],
                  Ppresence = df2$Ppresence[group1 == i],
                  type = 'Difference')
  df5 <- rbind(df5, x)
}

sel <- df5$type != 'Difference'
ggplot(df5[sel,], aes(x=sampling_effort, y=Ppresence, fill=p_correct)) +
  geom_raster() + 
  scale_x_continuous(breaks = log(labs), labels=labs) +
  scale_fill_continuous(type='viridis', name='% correct') + 
  xlab('Sampling effort (total # camera days)') + 
  ylab(expression(
    paste('Proportion of cameras with detections (', P[presence], ')'))) +
  facet_wrap(vars(type))

sel <- df5$type == 'Difference'
ggplot(df5, aes(x=sampling_effort, y=Ppresence, fill=100 - p_correct)) +
  geom_raster() + 
  scale_x_continuous(breaks = log(labs), labels=labs) +
  scale_fill_continuous(type='viridis', name='Relative performance') + 
  xlab('Sampling effort (total # camera days)') + 
  ylab(expression(
    paste('Proportion of cameras with detections (', P[presence], ')'))) 
