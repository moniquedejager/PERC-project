create_outputfile <- function(sim_nr) {
  # Create an outputfile. If the outputfile does not yet exist, create it
  # and insert the column names.  
  
  outputfile <- paste("./data/processed/PresenceAbsence/RoyleNicholsStats",
                      sim_nr,
                      ".txt")
  
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
}
