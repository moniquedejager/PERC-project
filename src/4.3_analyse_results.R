library(ggplot2)

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

windows(height=9, width=9)
ggplot(df[df$p_presence > 0,], aes(x=study_duration, y=time_interval, color=factor(cols2))) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') +
  facet_wrap(vars(species)) + 
  scale_color_manual(values=c('grey70', 'indianred4', 'darkolivegreen4'), na.value='antiquewhite', name='') + 
  xlab('Sampling period (days)') + 
  ylab('Time interval size (days)')+
  theme(legend.position = "top") + 
  guides(colour = guide_legend(override.aes = list(size=10)))

ggplot(df[df$p_presence > 0,], aes(x=study_duration, y=time_interval, color=p_presence)) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') +
  facet_wrap(vars(species)) + 
  scale_color_continuous(type='viridis', name='') + 
  xlab('Sampling period (days)') + 
  ylab('Time interval size (days)')+
  theme(legend.position = "top") 

df$cols3 = df$cols
df$cols3[is.na(df$cols)] = 4
df$cols3 = factor(df$cols3, labels = c("No significant effect", "More abundant in non-FSC forests", "More abundant in FSC forests", "NA"))

m = tapply((df$cols3 == "More abundant in FSC forests"), df$species, sum)
m = m - tapply((df$cols3 == "More abundant in non-FSC forests"), df$species, sum)
i = order(m)

plot(m, tapply(df$p_presence, df$species, max))
i = i[tapply(df$p_presence, df$species, max) > 0.1]

df2 = df[df$species %in% sort(unique(df$species))[i],]
df2$species = factor(df2$species, levels = sort(unique(df$species))[i]) 
df2$cols3 = factor(df2$cols3, levels = c("More abundant in non-FSC forests", "No significant effect", "NA","More abundant in FSC forests"))

windows(height=10, width=8)
ggplot(df2, aes(y=species, fill=cols3)) + 
  geom_bar() +
  scale_fill_manual(values=c('indianred4', 'grey70', 'antiquewhite','darkolivegreen4'), name='') +
  ylab('') +
  theme(legend.position = "top")

mod = lm(df$p_presence~df$species*df$time_interval*df$study_duration)
summary(mod)

windows(height=9, width=9)
ggplot(df2[df2$p_presence > 0,], aes(x=study_duration, y=time_interval, color=factor(cols2))) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') +
  facet_wrap(vars(species)) + 
  scale_color_manual(values=c('grey70', 'indianred4', 'darkolivegreen4'), na.value='antiquewhite', name='') + 
  xlab('Sampling period (days)') + 
  ylab('Time interval size (days)')+
  theme(legend.position = "top") + 
  guides(colour = guide_legend(override.aes = list(size=10)))

windows(height=9, width=9)
ggplot(df2, aes(x=study_duration, y=time_interval, color=p_presence)) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') +
  facet_wrap(vars(species)) + 
  scale_color_continuous(type='viridis', name='') + 
  xlab('Sampling period (days)') + 
  ylab('Time interval size (days)')+
  theme(legend.position = "top") 

windows(height=9, width=9)
ggplot(df2, aes(x=study_duration, y=time_interval, color=n_cams)) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') +
  facet_wrap(vars(species)) + 
  scale_color_continuous(type='viridis', name='') + 
  xlab('Sampling period (days)') + 
  ylab('Time interval size (days)')+
  theme(legend.position = "top") 
