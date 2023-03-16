library(ggplot2)

sim_nr <- 1
outputfile <- paste('./results/output/simulations/royle_nichols_stats_with_NAs', 
                    sim_nr, '.txt')
df <- read.table(outputfile, header=T)
df$sim_nr <- sim_nr

for (sim_nr in 2:10){
  outputfile <- paste('./results/output/simulations/royle_nichols_stats_with_NAs', 
                      sim_nr, '.txt')
  df2 <- read.table(outputfile, header=T)
  df2$sim_nr <- sim_nr
  df <- rbind(df,  df2)
}
df$sim_nr <- factor(df$sim_nr)

windows(height=5, width=10)
ggplot(df, aes(x=study_duration, y=z, color=factor(sim_nr))) +
  geom_point() + 
  geom_hline(yintercept = 1.96, linetype='dashed') + 
  facet_wrap(vars(type))

tapply(df$study_duration[df$type == 'Without NAs'], df$sim_nr[df$type == 'Without NAs'], max)

mod <- lm(z~study_duration*factor(sim_nr)*type,data=df)
summary(mod)

library(lme4)
mod <- lmer(z~study_duration + type + (1|sim_nr), data=df)
summary(mod)

df <- df[(!is.na(df$z))&(df$study_duration <= 140),]
tapply(df$z > 1.96, df$type, sum) / tapply(df$z, df$type, length)
