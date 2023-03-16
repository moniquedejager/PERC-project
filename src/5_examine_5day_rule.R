# How far off are the RN model results when using dT = 5 days 
# instead of optimal interval size?

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
  df <- rbind(df, df2)
}

df <- df[df$dT == 5,]
# too be continued... 

