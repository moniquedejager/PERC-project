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

tiff(filename='./results/figures/FSC and nonFSC/z-values_per_species_and_dt_and_nT.tiff', 
     height=10, width=13, units='in', res=300)
df$z2 <- df$z * -1
df$z2[df$z2 < -5] <- -5
df$z2[df$z2 > 5] <- 5
df$nT <- factor(df$min_number_of_intervals)
ggplot(df, aes(x=time_interval, y=z2, color=nT)) + 
  geom_line() + 
  scale_x_continuous(trans='log2') + 
  facet_wrap(vars(species)) +
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
  facet_wrap(vars(species), nrow=4) +
  xlab('Interval size (days)') + 
  ylab('z-value') +
  theme(legend.position = "top")+
  geom_hline(yintercept=1.96, linetype='dashed') + 
  geom_hline(yintercept=-1.96, linetype='dotted')

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

