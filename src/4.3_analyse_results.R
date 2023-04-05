library(ggplot2)
library(ggpattern)
library(ggpubr)

spec <- read.table('./data/processed/FSC and nonFSC data/species.txt')$V1

filename <- paste('./results/output/FSC and nonFSC data/royle_nichols_stats', 
                  spec[1], '.txt')
df <- read.table(filename, header=TRUE)

for (i in 2:length(spec)){
  filename <- paste('./results/output/FSC and nonFSC data/royle_nichols_stats', 
                    spec[i], '.txt')
  df2 <- read.table(filename, header=TRUE)
  df <- rbind(df, df2)
}

df$cols = 0 + 1*(df$P < 0.05) + 1*((df$P < 0.05)&(df$z < 0))
df$cols2 = factor(df$cols, labels = c("No significant effect", "More abundant in non-FSC forest", "More abundant in FSC forest"))

summary(lm(time_interval~p_presence+min_number_of_intervals+survey_effort,data=df))

sel <- (df$type == 'optimal dT')
df2 <- df[sel,]
i <- order(df2$z)
df2 <- df2[i, ]
df2 <- rbind(df2, df[df$type == '5 day interval',])
df2$species = factor(df2$species, levels = unique(df2$species)) 

windows(height=10, width=7)
tiff(filename='./results/figures/FSC and nonFSC/z_per_species.tiff', 
          height=10, width=7, units='in', res=300)
ggplot(df2, aes(y=species, x=-1*z, pattern=type, fill=cols2)) + 
  geom_bar_pattern(stat = 'identity', position='dodge',
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  xlab('Z-value') + 
  ylab('') + 
  scale_fill_manual(values = c("grey", "indianred4", "darkolivegreen4")) + 
  guides(pattern = guide_legend(override.aes = list(fill = "white"), title=''),
         fill = guide_legend(override.aes = list(pattern = "none"), title=''))+
  theme(axis.text.x = element_text(angle=90), 
        legend.position = c(0.7, 0.9),
        legend.background = element_blank())
dev.off()

p1 <- ggplot(df2[df2$type=='optimal dT',], aes(x=survey_effort, y=p_presence, 
                                               color=time_interval)) + 
  geom_point(size=3) + 
  xlim(0, 40000) + 
  ylim(0, 1) + 
  scale_color_continuous(type='viridis', name='Optimal time interval (days)') + 
  ylab('Proportion of cameras with detections') + 
  xlab('Sampling effort (total # camera days)') + 
  theme(legend.position = c(0.25, 0.8),
        legend.background = element_blank())

p2 <- ggplot(df2[df2$type=='optimal dT',], 
             aes(x=survey_effort, y=p_presence, color=cols2)) + 
  geom_point(size=3) + 
  xlim(0, 40000) + 
  ylim(0, 1) + 
  scale_color_manual(values = c("grey", "indianred4", "darkolivegreen4"), 
                     name='') + 
  ylab('Proportion of cameras with detections') + 
  xlab('Sampling effort (total # camera days)') + 
  theme(legend.position = c(0.25, 0.9),
        legend.background = element_blank())

windows(height=5, width=10)
tiff(filename='./results/figures/FSC and nonFSC/time_interval_and_effect_per_survey_effort_and_p_presence.tiff', 
     height=5, width=10, units='in', res=300)
ggarrange(p1, p2 + rremove("ylab"), labels=c('A', 'B'))
dev.off()

sel <- (df$type == 'optimal dT')
df2 <- df[sel,]

windows(height=5, width=5)
tiff(filename='./results/figures/FSC and nonFSC/Ppresence_per_optimal_interval_size_and_no_intervals.tiff', 
     height=5, width=5, units='in', res=300)
ggplot(df2, aes(x=time_interval, y=min_number_of_intervals, color=p_presence)) + 
  geom_point(size=3) + 
  scale_color_continuous(type='viridis', name='Proportion of cameras with detections',
                         limits = c(0, 1)) +
  xlab('Optimal interval size (days)') + 
  ylab('Minimum number of intervals') + 
  scale_y_continuous(trans='log10') + 
  theme(legend.position = "top") 
dev.off()

# create a table with all useful information per species:
write.csv(df, 'results/tables/FSC and nonFSC/results_per_species.csv')
