library(ggplot2)
library(ggpattern)
library(ggpubr)

spec <- read.table('./data/processed/FSC and nonFSC data/species.txt')$V1

filename <- paste('./results/output/FSC and nonFSC data/all_royle_nichols_stats', 
                  spec[1], '.txt')
df <- read.table(filename, header=TRUE)

for (i in 2:length(spec)){
  filename <- paste('./results/output/FSC and nonFSC data/all_royle_nichols_stats', 
                    spec[i], '.txt')
  df2 <- read.table(filename, header=TRUE)
  df <- rbind(df, df2)
}

df$cols = 0 + 1*(df$P < 0.05) + 1*((df$P < 0.05)&(df$z < 0))
df$cols2 = factor(df$cols, labels = c("No significant effect", "More abundant in non-FSC forest", "More abundant in FSC forest"))

# add the z-score: 
b <- read.table(
  './results/output/simulations/to_estimate_z_values.txt', 
  header=TRUE)$Estimate

x1 <- df$time_interval
x2 <- 1/df$survey_effort
x3 <- log(df$p_presence/(1 - df$p_presence))
x4 <- 1/df$n_cams[i]

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

df$z_score <- z

df$z2 <- df$z * -1
df$z2[df$z2 < -5] <- -5
df$z2[df$z2 > 5] <- 5
df$nT <- factor(df$min_number_of_intervals)

windows(height=10, width=13)

tiff(filename='./results/figures/FSC and nonFSC/results_per_species_and_dt_and_nT.tiff', 
     height=10, width=13, units='in', res=300)
ggplot(df, aes(x=time_interval, y=min_number_of_intervals, color=cols2)) + 
  geom_point(size=1.5) + 
  scale_color_manual(values = c("grey", "indianred4", "darkolivegreen4"), name='') + 
  scale_x_continuous(trans='log2') + 
  scale_y_continuous(trans='log2') + 
  facet_wrap(vars(species)) +
  xlab('Interval size (days)') + 
  ylab('Minimum number of intervals') +
  theme(legend.position = "top") 
dev.off()

windows(height=12, width=8)
tiff(filename='./results/figures/FSC and nonFSC/z-values_per_species_and_dt_and_nT.tiff', 
     height=12, width=8, units='in', res=300)
ggplot(df, aes(x=time_interval, y=z2, color=nT)) + 
  geom_line() + 
  scale_x_continuous(trans='log2') + 
  facet_wrap(vars(species), ncol=5) +
  xlab('Interval size (days)') + 
  ylab('z-value') +
  theme(legend.position = "top")+
  geom_hline(yintercept=1.96, linetype='dashed') + 
  geom_hline(yintercept=-1.96, linetype='dotted')
dev.off()

sel <- (df$species %in% c('Cricetomys emini', 
                          'Loxodonta cyclotis', 
                          'Mandrillus sphinx', 
                          'Gorilla gorilla'))&(df$min_number_of_intervals %in% 
                                                 c(2, 4, 8, 14, 24))
p1 <- ggplot(df[sel,], aes(x=time_interval, y=p_presence, color=nT)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(trans='log2') + 
  facet_wrap(vars(species), nrow=4) +
  xlab('Interval size (days)') + 
  ylab('Proportion of cameras with detections') +
  theme(legend.position = "top")

p2 <- ggplot(df[sel,], aes(x=time_interval, y=z*-1, color=nT)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(trans='log2') + 
  scale_y_continuous(position ='right') + 
  facet_wrap(vars(species), nrow=2) +
  xlab('Interval size (days)') + 
  ylab('z-value') +
  theme(legend.position = "top")+
  geom_hline(yintercept=1.96, linetype='dashed') + 
  geom_hline(yintercept=-1.96, linetype='dotted')

p5 <- ggplot(df[sel,], aes(x=time_interval, y=1/exp(z_score), color=nT)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(trans='log2') + 
  scale_y_continuous(trans='log10', position ='right') + 
  facet_wrap(vars(species), nrow=4) +
  xlab('Interval size (days)') + 
  ylab('estimated fit of Royle-Nichols occupancy model') +
  theme(legend.position = "top")+
  geom_hline(yintercept=1/exp(1.96), linetype='dashed') 
p5

p3 <- ggplot(df[sel,], aes(x=time_interval*min_number_of_intervals, y=z*-1, color=nT)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(trans='log2') + 
  facet_wrap(vars(species), nrow=4) +
  xlab('Interval size (days) * nT') + 
  ylab('z-value') +
  theme(legend.position = "top")+
  geom_hline(yintercept=1.96, linetype='dashed') + 
  geom_hline(yintercept=-1.96, linetype='dotted')

p4 <- ggplot(df[sel,], aes(x=survey_effort, y=z*-1, color=nT)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(trans='log2') + 
  facet_wrap(vars(species), nrow=4) +
  xlab('Sampling effort (days)') + 
  ylab('z-value') +
  theme(legend.position = "top")+
  geom_hline(yintercept=1.96, linetype='dashed') + 
  geom_hline(yintercept=-1.96, linetype='dotted')

tiff(filename='./results/figures/FSC and nonFSC/z-values_and_Ppresence_of_four_species_per_dt_and_nT.tiff', 
     height=10, width=6, units='in', res=300)
ggarrange(p1, p2, nrow=1, ncol=2, common.legend=TRUE)
dev.off()

tiff(filename='./results/figures/FSC and nonFSC/z-values_of_four_species_per_dt_and_nT.tiff', 
     height=5, width=5, units='in', res=300)
p2
dev.off()


tiff(filename='./results/figures/FSC and nonFSC/RNmodel_fit_and_Ppresence_of_four_species_per_dt_and_nT.tiff', 
     height=10, width=6, units='in', res=300)
ggarrange(p1, p5, nrow=1, ncol=2, common.legend=TRUE)
dev.off()

ggplot(df[sel,], aes(x=time_interval*min_number_of_intervals, y=survey_effort, color=nT)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(trans='log2') + 
  xlab('dT * nT') + 
  ylab('Sampling effort (days)') +
  theme(legend.position = "top")

p1 <- ggplot(df[sel,], aes(x=time_interval, y=survey_effort/1000, color=nT)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(trans='log2') + 
  scale_color_discrete(name='Minimum number of samples') + 
  xlab('Interval size (days)') + 
  ylab('Sampling effort (1,000 days)') +
  theme(legend.position = "top") 

p2 <- ggplot(df[sel,], aes(y=survey_effort/1000, x=n_cams, color=nT)) + 
  geom_line() + 
  geom_point() + 
  ylab('Sampling effort (1,000 days)') + 
  xlab('Number of cameras') +
  theme(legend.position = "top") 

tiff(filename='./results/figures/FSC and nonFSC/sampling_effort_per_dt_and_nT_and_ncams.tiff', 
     height=5, width=9, units='in', res=300)
ggarrange(p1, p2 + rremove("ylab"), common.legend = TRUE, labels=c('A', 'B'))
dev.off()

tiff(filename='./results/figures/FSC and nonFSC/sampling_effort_per_dt_and_nT.tiff', 
     height=5, width=5, units='in', res=300)
p1
dev.off()
