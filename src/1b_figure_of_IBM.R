# make a figure that explains the IBM:

densities <- round(exp(seq(1.5, 6, 0.25))) # vector with 19 different densities,
# in number of individuals per 100km2

# two examples, one with a high density, and one with a low density:

camTrapX <- rep(seq(300, 1900, 400), 5)
camTrapY <- sort(camTrapX)
camTrapPos <- paste(camTrapX, camTrapY, sep='-')
df <- data.frame(site = 'Low density site', 
                 type = 'Camera trap', 
                 x = camTrapX,
                 y = camTrapY)
n <- densities[1]*12.25
x <- runif(n, 0, 7000)
y <- runif(n, 0, 7000)
df2 <- data.frame(site = 'Low density site', 
                  type = 'Individual', 
                  x = x,
                  y = y)
df <- rbind(df, df2)

n <- densities[19]
x <- runif(n, 0, 2000)
y <- runif(n, 0, 2000)
df2 <- data.frame(site = 'High density site', 
                  type = 'Individual', 
                  x = x,
                  y = y)
df <- rbind(df, df2)
df2 <- data.frame(site = 'High density site', 
                 type = 'Camera trap', 
                 x = camTrapX,
                 y = camTrapY)
df <- rbind(df, df2)
df2 <- data.frame(site = 'High density site', 
                  type = 'max', 
                  x = 2000,
                  y = 2000)
df <- rbind(df, df2)
df2 <- data.frame(site = 'Low density site', 
                  type = 'max', 
                  x = 7000,
                  y = 7000)
df <- rbind(df, df2)

library(ggplot2)

windows(height=4, width=7)
ggplot(df, aes(x=x*5/1000, y=y*5/1000, color=type)) + 
  geom_point(aes(shape=type), size=0.8) + 
  facet_wrap(vars(site), scales='free') + 
  xlab('X position (km)') + 
  ylab('Y position (km)') + 
  scale_color_manual(values=c('black', 'lightsalmon', 'white')) +
  scale_shape_manual(values=c(15, 8, 3)) + 
  theme(legend.position="none")
  
