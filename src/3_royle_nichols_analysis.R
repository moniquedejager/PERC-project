library(ggplot2)
library(ggpubr)

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

df$logitPpres <- log(df$Ppresence/(1 - df$Ppresence))
df2 <- df[!is.infinite(df$logitPpres),]

# relation between average density and Ppresence:
df2$SD2 <- c('sampling one week', 
             'sampling one month', 
             'sampling three months',
             'sampling one year',
             rep('sampling one week', length(df2$dens1) - 4))
df2$SD2[df2$original_SD == 7] <- 'sampling one week'
df2$SD2[df2$original_SD == 28] <- 'sampling one month'
df2$SD2[df2$original_SD == 92] <- 'sampling three months'
df2$SD2[df2$original_SD == 365] <- 'sampling one year'
df2$SD2 <- factor(df2$SD2, labels = c('sampling one week', 
                                      'sampling one month', 
                                      'sampling three months',
                                      'sampling one year'))
sel <- (df2$nCams == 25)&(df2$dT == 1)&(df2$original_SD %in% c(7, 28, 92, 365))
df2$mDens <- (df2$dens1 + df2$dens2)/2
windows(height=6, width=6)

tiff(filename='./results/figures/simulations/Ppresence_per_average_density.tiff', 
     height=6, width=6, units='in', res=300)
ggplot(df2[sel,], aes(x=mDens/100, y=logitPpres)) + 
  geom_point() + 
  scale_x_continuous(trans='log10') + 
  facet_wrap(vars(SD2)) + 
  xlab(expression(paste('Average density (individuals ', km^-2, ' )'))) + 
  ylab('Proportion of cameras with detections (logit transformed)')
dev.off()

# Fit a model to the data that calculates the z-value from dT, survey effort, 
# and Ppresence (with dT transformed to fT (fraction of total study duration),
# which is logit-transformed as this ranges between 0 and 1. Ppresence is also
# logit transformed. 
# survey effort should be positively log-linearly related to the z-value 

y <- df2$z
x1 <- df2$dT 
x2 <- 1/(df2$StudyDuration * df2$nCams * 2)
x3 <- df2$logitPpres
x4 <- 1/(df2$nCams*2)

mod <- lm(y~x1*x2*x3*x4+ 
            I(x1^2)*x2*I(x3^2)*x4 + 
            I(x1^3)*x2*I(x3^3)*x4)
summary(mod)
# R2 = 0.583

s = summary(mod)

write.table(s$coefficients, 
            './results/output/simulations/to_estimate_z_values.txt', 
            append=FALSE, row.names = FALSE, col.names=TRUE)

# what is the best time interval size (with the highest z-score) per study duration, proportion of cameras with detections, and number of cameras?
df2$Ppres = ceiling(df2$Ppresence*20)/20
group = paste(df2$original_SD, df2$nCams, df2$Ppres, df2$sim_nr, sep='-')
df3 = df2[is.na(df2$dT),]
for (i in unique(group)) {
  df4 <- df2[group == i,]
  df4 <- df4[!is.na(df4$z),]
  df4 <- df4[(df4$z == max(df4$z)),]
  df3 <- rbind(df3,df4)
}
write.table(df3, './results/output/simulations/highest_z_values.txt', append=FALSE, col.names=TRUE, row.names = FALSE)

df3 <- read.table('./results/output/simulations/highest_z_values.txt', header=TRUE)

y <- df3$z
x1 <- 1/(df3$StudyDuration*df3$nCams*2)
x2 <- df3$logitPpres

mod <- lm(y~x1)
summary(mod)
mod1 <- lm(y~x1+x2)
summary(mod1)
mod2 <- lm(y~x1*x2 + I(x2^2)*x1 + I(x2^3)*x1 + I(x2^4))
summary(mod2)
mod3 <- lm(y~x1*x2 + I(x2^2)*x1 + I(x2^3)*x1 + I(x2^4)*x1)
summary(mod3)
anova(mod2, mod3)
# R2 = 0.76

df3$sampling_effort <- df3$StudyDuration * df3$nCams * 2
df3$est_z <- mod3$fitted.values

sel = (df3$est_z > 1.96)
labs = c(50, 150, 400, 1000, 3000, 10000)
p1 <- ggplot(df3[sel,], aes(x=round(log(sampling_effort), 1), y=Ppres, fill=est_z)) +
  geom_raster() + 
  ylim(0,1) + 
  scale_x_continuous(breaks = log(labs), labels=labs, limits=c(4,10)) +
  scale_fill_continuous(type='viridis', name='Estimated z-value') + 
  xlab('Sampling effort (total # camera days)') + 
  ylab(expression(paste('Proportion of cameras with detections (', P[presence], ')')))+
  theme(legend.position = "top") 

y <- df3$dT
mod <- glm(y~x1, family='poisson')
summary(mod)
mod1 <- glm(y~x1+x2, family='poisson')
summary(mod1)
mod2 <- glm(y~x1*x2 + I(x2^2)*x1 + I(x2^3)*x1 + I(x2^4), family='poisson')
summary(mod2)
mod3 <- glm(y~x1*x2 + I(x2^2)*x1 + I(x2^3)*x1 + I(x2^4)*x1, family='poisson')
summary(mod3)
s = summary(mod3)

with(summary(mod3), 1 - deviance/null.deviance)
# R2 = 0.22

df3$est_dT = ceiling(mod3$fitted.values)

p2 <- ggplot(df3[sel,], aes(x=round(log(sampling_effort), 1), y=Ppres, fill=est_dT)) +
  geom_raster() + 
  ylim(0, 1) + 
  scale_x_continuous(breaks = log(labs), labels=labs, limits=c(4, 10)) +
  scale_fill_continuous(type='viridis', trans='log10', name='Optimal interval size (days)') + 
  xlab('Sampling effort (total # camera days)') + 
  ylab('')+
  theme(legend.position = "top") 

# what is the difference between dT = 5 and optimal time interval? 
df4 <- df[df$dT == 5,]
# per combination of survey effort and Ppresence, we calculate the 
# percentage of RN models that resulted in a significant and correct
# answer: 

df4$sampling_effort <- df4$StudyDuration * df4$nCams * 2
df4$samp_eff2 <- round(log(df4$sampling_effort), 1)
df4$Ppres <- ceiling(df4$Ppresence*20)/20
df4$z[is.na(df4$z)] <- 0

group <- paste(df4$samp_eff2, df4$Ppres, sep='-')
n_tot <- tapply(df4$z, group, length)
n_correct <- tapply(df4$z > 1.96, group, sum)
p_correct <- n_correct/n_tot *100
sampling_effort <- tapply(df4$samp_eff2, group, mean)
Ppresence <- tapply(df4$Ppres, group, mean)
df5 <- data.frame(p_correct = p_correct,
                  sampling_effort = sampling_effort,
                  Ppresence = Ppresence)

df3$samp_eff2 <- round(log(df3$sampling_effort), 1)
df3$z[is.na(df3$z)] <- 0

group <- paste(df3$samp_eff2, df3$Ppres, sep='-')
n_tot <- tapply(df3$z, group, length)
n_correct <- tapply(df3$z > 1.96, group, sum)
p_correct <- n_correct/n_tot *100
sampling_effort <- tapply(df3$samp_eff2, group, mean)
Ppresence <- tapply(df3$Ppres, group, mean)
df6 <- data.frame(p_correct = p_correct,
                  sampling_effort = sampling_effort,
                  Ppresence = Ppresence)

df5$type <- 'Using 5 day interval'
df6$type <- 'Using optimal interval size'
df7 <- rbind(df5, df6)

# for each value in df2, find the corresponding value
# in df4 and calculate the difference:
group1 <- paste(df5$sampling_effort, df5$Ppresence, sep='-')
group2 <- paste(df6$sampling_effort, df6$Ppresence, sep='-')

for (i in unique(group1)){
  p1 <- df5$p_correct[group1 == i]
  p2 <- df6$p_correct[group2 == i]
  if (length(p2) == 0){
    p2 <- 0
    x <- data.frame(p_correct =0,
                    sampling_effort = df5$sampling_effort[group1 == i],
                    Ppresence = df5$Ppresence[group1 == i],
                    type = 'Using optimal interval size')
    df7 <- rbind(df7, x)
  } 
  performance <- p2 / (p1 + p2)
  if (is.na(performance)) { performance <- 0.5}
  if (performance < 0.5) { performance <- 0.5}
  x <- data.frame(p_correct = performance,
                  sampling_effort = df5$sampling_effort[group1 == i],
                  Ppresence = df5$Ppresence[group1 == i],
                  type = 'Difference')
  df7 <- rbind(df7, x)

}

sel <- df7$type == 'Difference'
p3 <- ggplot(df7[sel,], aes(x=sampling_effort, y=Ppresence, fill=p_correct)) +
  geom_raster() + 
  ylim(0,1)+
  scale_x_continuous(breaks = log(labs), labels=labs, limits=c(4, 10)) +
  scale_fill_continuous(type='viridis', name='Relative performance') + 
  xlab('Sampling effort (total # camera days)') + 
  ylab(expression(
    paste('Proportion of cameras with detections (', P[presence], ')'))) +
  theme(legend.position = "top") 

windows(height=5, width=12)
tiff(filename='./results/figures/simulations/estimated_z_values_optimal_interval_sizes_and_relative_performance.tiff', 
     height=5, width=12, units='in', res=300)
ggarrange(p1, p2 + rremove("ylab"), p3 + rremove("ylab"), 
          labels=c('A', 'B', 'C'), nrow=1)
dev.off()

sel <- df7$type != 'Difference'
p4 <- ggplot(df7[sel,], aes(x=sampling_effort, y=Ppresence, fill=p_correct)) +
  geom_raster() + 
  scale_x_continuous(breaks = log(labs), labels=labs) +
  scale_fill_continuous(type='viridis', name='% correct') + 
  xlab('Sampling effort (total # camera days)') + 
  ylab(expression(
    paste('Proportion of cameras with detections (', P[presence], ')'))) +
  facet_wrap(vars(type))
windows(height=5, width=9)
tiff(filename='./results/figures/simulations/p_correct.tiff', 
     height=5, width=9, units='in', res=300)
p4
dev.off()


# finding the good, the bad, and the ugly: 
# how often do combinations of study duration, time interval size, 
# ncams, and Ppresence result in correct results, non-significant results, 
# and incorrect results? 
df2 = df2[!is.na(df2$z),]
group <- paste(df2$dT, df2$original_SD, 
               df2$nCams, df2$dens1, df2$dens2, sep='-')
df5 <- data.frame(StudyDuration = tapply(df2$original_SD, group, mean),
                  nCams = tapply(df2$nCams, group, mean),
                  Ppresence = tapply(df2$Ppresence, group, mean),
                  dT = tapply(df2$dT, group, mean),
                  nCorrect = tapply((df2$z > 1.96), group, sum),
                  dens1 = tapply(df2$dens1, group, mean),
                  dens2 = tapply(df2$dens2, group, mean))
df5$Ppres <- round(df5$Ppresence, 1)

sel <- (df5$nCams == 25)
windows(height=10, width=10)
df5$dT2 <- exp(round(log(df5$dT)*2)/2)
df5$SD2 <- round(log(df5$StudyDuration)/log(365),1)

ggplot(df5[sel,], aes(x=SD2, y=dT2, fill=nCorrect)) + 
  geom_raster() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8, 1), labels=c(10, 30, 100, 365)) +
  facet_wrap(vars(Ppres)) + 
  scale_fill_continuous(type='viridis') + 
  xlab('Sampling period (days)') + 
  ylab('Time interval size (days)')+
  theme(legend.position = "top")

sel <- (df5$nCams == 25)&
  (df5$dens1 %in% sort(unique(df5$dens1))[c(1, 3, 6, 9, 12, 15,18)])&
  (df5$dens2 %in% sort(unique(df5$dens2))[c(1, 3, 6, 9, 12, 15,18)])

ggplot(df5[sel,], aes(x=SD2, y=dT2, fill=nCorrect*10)) + 
  geom_raster() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8, 1), labels=c(10, 30, 100, 365)) +
  facet_grid(cols=vars(round(dens1)), rows=vars(round(dens2))) + 
  scale_fill_continuous(type='viridis', name='% Correct') + 
  xlab('Sampling period (days)') + 
  ylab('Time interval size (days)')+
  theme(legend.position = "top")

mod = lm(nCorrect~SD2*dT2*dens1*dens2*nCams ,data=df5)
summary(mod)

mod = lm(nCorrect~SD2+dT2+dens1+dens2+nCams ,data=df5)
summary(mod)