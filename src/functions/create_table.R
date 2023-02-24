# We create a file that contains all the parameter combinations,
# which we can read in when running the royle_nichols_stats model 
# in parallel. 

# There are 19 simulations with different densities for which we
# have presence/absence matrices of 365 days in 1 day intervals.
i1 <- 1:19
i2 <- 1:19

# To get all possible combinations of the parameter values, 
# we create a dataframe containing these
df <- expand.grid(i1 = i1, i2 = i2, time_interval = time_interval, study_duration = study_duration, n_cams = n_cams)

# Some parameter value combinations are impossible and are excluded
df <- df[df$i1 >= df$i2,]
df <- df[df$time_interval < (df$study_duration / 2),]

write.table(df, './data/temp/parameter_combinations_royle_nichols_stats.txt', append = FALSE, row.names = FALSE, col.names = TRUE)