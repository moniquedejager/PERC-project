# 1. create an outputfile:
# source('./src/functions/create_outputfile.R')

outputfile <- paste("./data/processed/PresenceAbsence/RoyleNicholsStats",
                    sim_nr,
                    ".txt")

# if the outputfile does not yet exist, create it
# and insert the column names: 
if (!file.exists(outputfile)) {
  x <- vector(length = 0)
  data <- data.frame(
    dens1 = x,
    dens2 = x,
    time_interval = x,
    study_duration = x,
    n_cams = x,
    z = x,
    P = x,
    p_presence = x)
  write.table(data, outputfile, col.names = TRUE, row.names = FALSE)
}