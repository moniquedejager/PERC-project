adjust_matrices <- function(study_duration, time_interval, dets, cov) {
  # We need to adjust the detection and covariates matrices to the
  # parameters used in the i'th row of the values table
  
  n_intervals <- floor(study_duration / time_interval)
  study_duration <- n_intervals * time_interval
  
  dets2 <- dets[, 1:study_duration]
  dets3 <- matrix(0, 50, n_intervals)
  
  # voor iedere rij in dets2, presence/absence samenvoegen per
  # tijdsintervalserval...
  x <- sort(rep(1:n_intervals, time_interval))
  for (i3 in 1:50) {
    dets3[i3, ] <- tapply(dets2[i3, ], x, max)
  }
  
  # randomly select cameras to use...
  cams_used <- sort(rank(runif(25))[1:n_cams])
  # using n_cams random camera traps.
  dets3 <- dets3[c(cams_used, cams_used + 25), ]
  cov2 <- cov[c(cams_used, cams_used + 25), ]
  return(cbind(cov2, dets3))
}
