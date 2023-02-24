# with the presence/absence matrices per density, we can use the RN occupany
# model to examine whether it provides a significan_intervals difference between
# pairs of simulations with differen_intervals densities.

# We can do this with differen_intervals parameter settings:
# - density-combinations
# - time in_intervalservals
# - study durations
# - number of camera traps
# do this for 10 simulation runs, using 10 parallel sessions...

# 1. Create a file that contains all the parameter combinations:
time_interval <- 1 #unique(round(1.1^(4:54)))
study_duration <- 365 #round(exp(log(365)*(10:30/30)))
n_cams <- 20 #25:5
source('./src/functions/create_table.R')

# 2. load the royle_nichols_stats function:
source('./src/functions/royle_nichols_stats.R')

# 3. library needed for parallel runs:
library(parallel)
# create the local cluster:
cl <- makeCluster(10)
# run the simulations:
results <- parSapply(cl, 1:10, royle_nichols_stats)
# stop the cluster:
stopCluster(cl)
