# with the presence/absence matrices per density, 
# we can use the Royle-Nichols occupany model to 
# examine whether it provides a significant difference between
# pairs of simulations with different densities.

# We can do this with different parameter settings:
# - density-combinations
# - time intervals
# - study durations
# - number of camera traps
# do this for 10 simulation runs, using 10 parallel sessions...

# 1. Create a file that contains all the parameter combinations:
# for the time interval,
# we generally use the range unique(round(1.1^(4:54)))
# for study duration,
# we generally use the range round(exp(log(365)*(10:30/30)))
# for the number of camera traps per location,
# we generally use 25:5.
time_interval <- 1
study_duration <- 365
n_cams <- 20
source('./src/functions/create_table.R')

# 2. 
source('./src/functions/royle_nichols_stats.R')

# 3. run in parallel:
library(parallel)
cl <- makeCluster(10)
results <- parSapply(cl, 1:10, royle_nichols_stats)
stopCluster(cl)
