# with the presence/absence matrices per density, we can use the RN occupany
# model to examine whether it provides a significan_intervals difference between
# pairs of simulations with differen_intervals densities.

# We can do this with differen_intervals parameter settings:
# - density-combinations
# - time in_intervalservals
# - study durations
# - number of camera traps
# do this for 10 simulation runs, using 10 parallel sessions...

# TODO iteraties eruit halen en parameter combinaties inlezen uit file!

royle_nichols_stats <- function(sim_nr) {
  library(unmarked)

  x <- vector(length = 0)
  data <- data.frame(
    dens1 = x,
    dens2 = x,
    time_in_intervalserval = x,
    study_duration = x,
    n_cams = x,
    z = x,
    P = x,
    p_presence = x)

  outputfile <- paste("./data/processed/PresenceAbsence/RoyleNicholsStats",
                      sim_nr,
                      ".txt")

  if (!file.exists(outputfile)) {
    write.table(data, outputfile, col.names = TRUE, row.names = FALSE)
  }

  for (i1 in 1:18){
    for (i2 in (i1 + 1):19){
      study_duration <- 365
      n_cams <- 25

      print(paste(i1, i2))
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

      for (n_cams in 25) { # 25:5)
        for (study_duration in 50) { # round(exp(log(365)*(10:30/30))))
          for (time_in_intervalserval in 1) { # unique(round(1.1^(4:54))))
            if (time_in_intervalserval <= study_duration / 2) {
              sd <- study_duration
              n_intervals <- floor(sd / time_in_intervalserval)
              sd <- n_intervals * time_in_intervalserval

              dets2 <- dets[, 1:sd]
              dets3 <- matrix(0, 50, n_intervals)

              # voor iedere rij in dets2, presence/absence samenvoegen per
              # tijdsintervalserval...
              x <- sort(rep(1:n_intervals, time_in_intervalserval))
              for (i3 in 1:(n_cams * 2)) {
                dets3[i3, ] <- tapply(dets2[i3, ], x, max)
              }

              # randomly select cameras to use...
              cams_used <- sort(rank(runif(25))[1:n_cams])
              # using n_cams random camera traps.
              dets3 <- dets3[c(cams_used, cams_used + 25), ]
              cov2 <- cov[c(cams_used, cams_used + 25), ]

              p_presence <- mean(rowSums(dets3) > 0)

              umf <- unmarkedFrameOccu(y = dets3, siteCovs = cov2)
              m1 <- occuRN(~1 ~ location, umf)
              s1 <- summary(m1)

              data <- data.frame(
                dens1 = dens1,
                dens2 = dens2,
                time_in_intervalserval = 1,
                study_duration = sd,
                n_cams = n_cams,
                z = s1$state$z[2],
                P = s1$state$`P(>|z|)`[2],
                p_presence = p_presence)

              write.table(
                data, outputfile,
                append = TRUE,
                col.names = FALSE,
                row.names = FALSE)
            }
          }
        }
      }
    }
  }
}

# library needed for parallel runs:
library(parallel)

# create the local cluster:
cl <- makeCluster(10)

# run the simulations:
results <- parSapply(cl, 1:10, royle_nichols_stats)
# stop the cluster:
stopCluster(cl)

library(lin_intervalsr)
