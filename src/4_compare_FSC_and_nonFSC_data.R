# To compare the camera trapping events between FSC and non-FSC forests,
# we need to do the following:

# 4.1: Create presence/absence matrices per species, per day, and
# per camera trap:
source('./src/4.1_create_pres_abs_matrices.R')

# 4.2: Calculate the royle-nichols stats per species, study duration, 
# and time interval:
source('/src/4.2_royle_nichols_stats.R')