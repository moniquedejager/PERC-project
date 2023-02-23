# plot the results

simNr = 1
outputfile = paste('./data/processed/PresenceAbsence/RoyleNicholsStats', simNr, '.txt')
df = read.table(outputfile, header=T)
df = df[df$nCams == 25,]
df$effect = (0 - 1*(df$z < -1.959964) + 1*(df$z > 1.959964))
df$effect[df$effect == 0] = 'No effect'  
df$effect[df$effect == -1] = 'Incorrect relation'
df$effect[df$effect == 1] = 'Correct relation'
df$effect = factor(df$effect, levels=c('No effect', 'Incorrect relation', 'Correct relation'))

library(ggplot2)
ggplot(df, aes(x=dens1, y=dens2, fill=effect)) + geom_raster() +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10')
