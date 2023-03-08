# This R-script contains all the elements needed for this project: 
# 1. The individual-based model that generates camera-trapping events #####
source('./src/1_IBM.R') 
# 2. The royle-nichols model that is run on many combinations of: #####
#    - study duration
#    - number of camera traps
#    - time interval size
#    - density at location 1
#    - density at location 2 #####
source('/src/2_royle_nichols_stats.R')
# 3. The analysis of the royle-nichols results #####
source('./src/3_royle_nichols_analysis.R')
# 4. The application of the analysis results on real data #####





