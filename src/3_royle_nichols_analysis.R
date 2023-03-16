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

sel <- (df$nCams == 25)&(df$sim_nr == 1)
windows(height=10, width=10)
ggplot(df[sel,], aes(x=StudyDuration, y=dT, color=factor(cols2))) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') +
  facet_wrap(vars(Ppres)) + 
  scale_color_manual(values=c('grey70', 'darkolivegreen4','indianred4'), 
                     na.value='antiquewhite', name='') + 
  xlab('Sampling period (days)') + 
  ylab('Time interval size (days)')+
  theme(legend.position = "top") + 
  guides(colour = guide_legend(override.aes = list(size=10)))

sel <- (df$nCams == 5)&(df$StudyDuration == 365)&(df$dT == 1)
ggplot(df[sel,], aes(x=dens1, y=dens2, color=P)) +
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') +
  scale_color_continuous(type='viridis', trans='log10')

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

df$logitPpres = log(df$Ppresence/(1 - df$Ppresence))
df2 = df[!is.infinite(df$logitPpres),]

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

y = df2$z
x1 = df2$dT
x2 = df2$StudyDuration
x3 = df2$nCams
x4 = df2$logitPpres

mod = lm(y~x1 + x2 + x3 + x4 +
           x1:x2 + x1:x3 + x1:x4 + 
           x2:x3 + x2:x4 + x3:x4 + 
           x1:x2:x3 + x1:x3:x4 + x2:x3:x4 + x1:x2:x3:x4 +  
           
           I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) + 
           I(x1^2):I(x2^2) + I(x1^2):I(x3^2) + I(x1^2):I(x4^2) + 
           I(x2^2):I(x3^2) + I(x2^2):I(x4^2) + I(x3^2):I(x4^2) + 
           I(x1^2):I(x2^2):I(x3^2) + I(x1^2):I(x2^2):I(x4^2) + 
           I(x1^2):I(x3^2):I(x4^2) + I(x2^2):I(x3^2):I(x4^2) + 
           I(x1^2):I(x2^2):I(x3^2):I(x4^2) + 
           
           I(x1^3) + I(x2^3) + I(x3^3) + I(x4^3) + 
           I(x1^3):I(x2^3) + I(x1^3):I(x3^3) + I(x1^3):I(x4^3) + 
           I(x2^3):I(x3^3) + I(x2^3):I(x4^3) + I(x3^3):I(x4^3) + 
           I(x1^3):I(x2^3):I(x3^3) + I(x1^3):I(x2^3):I(x4^3) + 
           I(x1^3):I(x3^3):I(x4^3) + I(x2^3):I(x3^3):I(x4^3) + 
           I(x1^3):I(x2^3):I(x3^3):I(x4^3) + 
           
           I(x1^4) + I(x2^4) + I(x3^4) + I(x4^4) + 
           I(x1^4):I(x2^4) + I(x1^4):I(x3^4) + I(x1^4):I(x4^4) + 
           I(x2^4):I(x3^4) + I(x2^4):I(x4^4) + I(x3^4):I(x4^4) + 
           I(x1^4):I(x2^4):I(x3^4) + I(x1^4):I(x2^4):I(x4^4) + 
           I(x1^4):I(x3^4):I(x4^4) + I(x2^4):I(x3^4):I(x4^4) + 
           I(x1^4):I(x2^4):I(x3^4):I(x4^4))
summary(mod)

mod1 = lm(y~x1 + x2 + x3 + x4 +
           x1:x2 + x1:x3 + x1:x4 + 
           x2:x3 + x2:x4 + x3:x4 + 
           x1:x2:x3 + x1:x3:x4 + x2:x3:x4 + x1:x2:x3:x4 +  
           
           I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) + 
           I(x1^2):I(x2^2) + I(x1^2):I(x3^2) + I(x1^2):I(x4^2) + 
           I(x2^2):I(x3^2) + I(x2^2):I(x4^2) + I(x3^2):I(x4^2) + 
           I(x1^2):I(x2^2):I(x4^2) + 
           I(x1^2):I(x3^2):I(x4^2) + I(x2^2):I(x3^2):I(x4^2) + 
           I(x1^2):I(x2^2):I(x3^2):I(x4^2) + 
           
           I(x1^3) + I(x2^3) + I(x3^3) + I(x4^3) + 
           I(x1^3):I(x3^3) + I(x1^3):I(x4^3) + 
           I(x2^3):I(x3^3) + I(x2^3):I(x4^3) + I(x3^3):I(x4^3) + 
           I(x1^3):I(x2^3):I(x3^3) + I(x1^3):I(x2^3):I(x4^3) + 
           I(x1^3):I(x3^3):I(x4^3) + I(x2^3):I(x3^3):I(x4^3) + 
           I(x1^3):I(x2^3):I(x3^3):I(x4^3) + 
           
           I(x1^4) + I(x2^4) + I(x3^4) + I(x4^4) + 
           I(x1^4):I(x2^4) + I(x1^4):I(x3^4) + I(x1^4):I(x4^4) + 
           I(x2^4):I(x3^4) + I(x2^4):I(x4^4) + I(x3^4):I(x4^4) + 
           I(x1^4):I(x2^4):I(x3^4) + I(x1^4):I(x2^4):I(x4^4) + 
           I(x1^4):I(x3^4):I(x4^4) + I(x2^4):I(x3^4):I(x4^4) + 
           I(x1^4):I(x2^4):I(x3^4):I(x4^4))
s = summary(mod1)

# adjusted R squared = 0.604
data <-s$coefficients
write.table(data, './results/output/simulations/estimated_z_value.txt', append=F, col.names = TRUE, row.names=TRUE)


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

windows(height=12, width=4)
sel <- (df3$nCams == 25) & (df3$sim_nr %in% c(1, 2, 3, 4, 5))
ggplot(df3[sel,], aes(x=original_SD, y=Ppres, fill=P)) +
  geom_raster() + 
  scale_x_continuous(trans='log10') +
  facet_wrap(vars(sim_nr), ncol=1) + 
  scale_fill_continuous(type='viridis', trans='log10') + 
  xlab('Sampling period (days)') + 
  ylab('Ppresence')+
  theme(legend.position = "top") 

ggplot(df3[sel,], aes(x=original_SD, y=Ppres, fill=z)) +
  geom_raster() + 
  scale_x_continuous(trans='log10') +
  facet_wrap(vars(sim_nr), ncol=1) + 
  scale_fill_continuous(type='viridis') + 
  xlab('Sampling period (days)') + 
  ylab('Ppresence')+
  theme(legend.position = "top") 

ggplot(df3[sel,], aes(x=original_SD, y=Ppres, fill=dT)) +
  geom_raster() + 
  #scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') +
  facet_wrap(vars(sim_nr), ncol=1) + 
  scale_fill_continuous(type='viridis', trans='log10') + 
  xlab('Sampling period (days)') + 
  ylab('Ppresence')+
  theme(legend.position = "top") 

y <- df3$z
x1 <- df3$StudyDuration
x2 <- df3$nCams
x3 <- df3$logitPpres

mod0 <- lm(y~x1 + x2 + x3 + 
           x1:x2 + x1:x3 +  
           x2:x3 +
           x1:x2:x3 + 
           
           I(x1^2) + I(x2^2) + I(x3^2) + 
           I(x1^2):I(x2^2) + I(x1^2):I(x3^2) + 
           I(x2^2):I(x3^2) +
           I(x1^2):I(x2^2):I(x3^2) +  
           
           I(x1^3) + I(x2^3) + I(x3^3) + 
           I(x1^3):I(x2^3) + I(x1^3):I(x3^3) + 
           I(x2^3):I(x3^3) + 
           I(x1^3):I(x2^3):I(x3^3) + 
           
           I(x1^4) + I(x2^4) + I(x3^4) + 
           I(x1^4):I(x2^4) + I(x1^4):I(x3^4) + 
           I(x2^4):I(x3^4) + 
           I(x1^4):I(x2^4):I(x3^4))
summary(mod0)

mod1 <- lm(y~x1 + x2 + x3 + 
             x1:x2 + x1:x3 +  
             x2:x3 +
             x1:x2:x3 + 
             
             I(x1^2) + I(x2^2) + I(x3^2) + 
             I(x1^2):I(x2^2) + I(x1^2):I(x3^2) + 
             I(x2^2):I(x3^2) +
             
             I(x1^3) + I(x2^3) + I(x3^3) + 
             I(x1^3):I(x2^3) + I(x1^3):I(x3^3) + 
             I(x1^3):I(x2^3):I(x3^3) + 
             
             I(x1^4) + I(x2^4) + I(x3^4) + 
             I(x1^4):I(x2^4) + I(x1^4):I(x3^4) + 
             I(x2^4):I(x3^4) + 
             I(x1^4):I(x2^4):I(x3^4))
summary(mod1)

df3$est_z = mod1$fitted.values

sel = (df3$nCams %in% c(5, 10, 15, 20, 25)) & (df3$est_z > 1.96)
df3$SD2 <- round(log(df3$original_SD)/log(365),1)

p1 <- ggplot(df3[sel,], aes(x=SD2, y=Ppres, fill=est_z)) +
  geom_raster() + 
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8, 1), labels=c(10, 30, 100, 365)) +
  facet_wrap(vars(nCams), ncol=1) + 
  scale_fill_continuous(type='viridis', name='Estimated z-value') + 
  xlab('Sampling period (days)') + 
  ylab('Ppresence')+
  theme(legend.position = "top") 

y <- df3$dT
x1 <- df3$StudyDuration
x2 <- df3$nCams
x3 <- df3$logitPpres

mod0 <- glm(y~x1 + x2 + x3 + 
             x1:x2 + x1:x3 +  
             x2:x3 +
             x1:x2:x3 + 
             
             I(x1^2) + I(x2^2) + I(x3^2) + 
             I(x1^2):I(x2^2) + I(x1^2):I(x3^2) + 
             I(x2^2):I(x3^2) +
             I(x1^2):I(x2^2):I(x3^2) +  
             
             I(x1^3) + I(x2^3) + I(x3^3) + 
             I(x1^3):I(x2^3) + I(x1^3):I(x3^3) + 
             I(x2^3):I(x3^3) + 
             I(x1^3):I(x2^3):I(x3^3) + 
             
             I(x1^4) + I(x2^4) + I(x3^4) + 
             I(x1^4):I(x2^4) + I(x1^4):I(x3^4) + 
             I(x2^4):I(x3^4) + 
             I(x1^4):I(x2^4):I(x3^4), family='poisson')
summary(mod0)

mod2 <- glm(y~x1 + x2 + x3 + 
              x1:x2 + x1:x3 +  
              x2:x3 +
              x1:x2:x3 + 
              
              I(x1^2) + I(x2^2) + I(x3^2) + 
              I(x1^2):I(x2^2) + I(x1^2):I(x3^2) + 
              
              I(x1^3) + I(x2^3) + I(x3^3) + 
              I(x1^3):I(x2^3) + I(x1^3):I(x3^3) + 
              I(x2^3):I(x3^3) + 
              I(x1^3):I(x2^3):I(x3^3) + 
              
              I(x1^4) + I(x2^4) + I(x3^4) + 
              I(x1^4):I(x2^4) + I(x1^4):I(x3^4) + 
              I(x1^4):I(x2^4):I(x3^4), family='poisson')
summary(mod2)

df3$est_dT = ceiling(mod2$fitted.values)
p2 <- ggplot(df3[sel,], aes(x=SD2, y=Ppres, fill=est_dT)) +
  geom_raster() + 
  #scale_y_continuous(trans='log10') + 
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8, 1), labels=c(10, 30, 100, 365)) +
  facet_wrap(vars(nCams), ncol=1) + 
  scale_fill_continuous(type='viridis', trans='log10', name='Optimal interval size (days)') + 
  xlab('Sampling period (days)') + 
  ylab('Ppresence')+
  theme(legend.position = "top") 

ggarrange(p1, p2)

df3$dT_dif = abs(df3$dT - df3$est_dT)
ggplot(df3[sel,], aes(x=original_SD, y=Ppres, fill=dT_dif)) +
  geom_raster() + 
  #scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') +
  facet_wrap(vars(nCams), ncol=1) + 
  scale_fill_continuous(type='viridis', trans='log10',name='Difference interval size (days)') + 
  xlab('Sampling period (days)') + 
  ylab('Ppresence')+
  theme(legend.position = "top") 
