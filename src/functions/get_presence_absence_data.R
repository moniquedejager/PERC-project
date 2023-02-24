get_pres_abs_data <- function(i1, i2, sim_nr) {
  #' Get the presence/absence data from the files.
  df1 <- as.matrix(
    read.table(paste("./data/processed/PresenceAbsence/PresAbs",
                     i1, sim_nr, ".txt")))
  dens1 <- df1[1, 1]
  df1 <- df1[, 2:366]
  
  df2 <- as.matrix(
    read.table(paste("./data/processed/PresenceAbsence/PresAbs",
                     i2, sim_nr, ".txt")))
  dens2 <- df2[1, 1]
  df2 <- df2[, 2:366]
  
  cov <- data.frame(
    location = factor(c(rep("A", 25), rep("B", 25))),
    density = c(rep(dens1, 25), rep(dens2, 25))
  )
  dets <- rbind(df1, df2)
  return(cbind(cov, dets))
}
