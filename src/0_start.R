# This R-script contains all the elements needed for this project: 
# 1. The individual-based model that generates camera-trapping events
# 2. The script to turn these data into presence/absence matrices
# 3. The royle-nichols model that is run on many combinations of:
#    - study duration
#    - number of camera traps
#    - time interval size
#    - density at location 1
#    - density at location 2
# 4. The analysis of the royle-nichols results
# 5. The application of the analysis results on real data

# 1. The individual-based model that generates camera-trapping events #####
source('./src/1_IBM.R')

