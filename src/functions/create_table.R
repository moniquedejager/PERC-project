# 1. Create a file that contains all the parameter combinations:
i1 <- 1:19
i2 <- 1:19

df <- expand.grid(i1 = i1, i2 = i2, time_interval = time_interval, study_duration = study_duration, n_cams = n_cams)
# exclude all occurrences where i1 >= i2:
df <- df[df$i1 >= df$i2,]
# exclude all occurrences where time_interval > study_duration / 2
df <- df[df$time_interval < (df$study_duration / 2),]
# write the parameter combinations to a file:
write.table(df, './data/temp/parameter_combinations_royle_nichols_stats.txt', append = FALSE, row.names = FALSE, col.names = TRUE)