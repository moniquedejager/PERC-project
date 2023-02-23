# with the presence/absence matrices per density, we can use the RN occupany
# model to examine whether it provides a significant difference between
# pairs of simulations with different densities.

# We can do this with different parameter settings:
# - density-combinations
# (also same densties against each other to examine false significant effects)
# - time intervals
# - study durations
# - number of camera traps
# do this for 10 simulation runs, using 10 parallel sessions...

royle_nichols_stats <- function(simNr) {
  library(unmarked)

  x <- vector(length = 0)
  DATA <- data.frame(dens1 = x, dens2 = x, dT = x, StudyDuration = x, nCams = x, z = x, P = x, Ppresence = x)
  outputfile <- paste("./data/processed/PresenceAbsence/RoyleNicholsStats", simNr, ".txt")
  if (!file.exists(outputfile)) {
    write.table(DATA, outputfile, col.names = T, row.names = F)
  }

  for (i1 in 1:18)
  {
    for (i2 in (i1 + 1):19)
    {
      studyDuration <- 365
      nCams <- 25

      print(paste(i1, i2))
      df1 <- as.matrix(read.table(paste("./data/processed/PresenceAbsence/PresAbs", i1, simNr, ".txt")))
      dens1 <- df1[1, 1]
      df1 <- df1[, 2:366]

      df2 <- as.matrix(read.table(paste("./data/processed/PresenceAbsence/PresAbs", i2, simNr, ".txt")))
      dens2 <- df2[1, 1]
      df2 <- df2[, 2:366]

      cov <- data.frame(
        location = factor(c(rep("A", 25), rep("B", 25))),
        density = c(rep(dens1, 25), rep(dens2, 25))
      )
      Dets <- rbind(df1, df2)

      for (nCams in 25) # 25:5)
      {
        for (studyDuration in 50) # round(exp(log(365)*(10:30/30))))
        {
          for (dT in 1) # unique(round(1.1^(4:54))))
          {
            if (dT <= studyDuration / 2) {
              SD <- studyDuration
              nT <- floor(SD / dT)
              SD <- nT * dT

              Dets2 <- Dets[, 1:SD]
              Dets3 <- matrix(0, 50, nT)
              # voor iedere rij in Dets2, presence/absence samenvoegen per tijdsinterval...
              x <- sort(rep(1:nT, dT))
              for (i3 in 1:(nCams * 2)) {
                Dets3[i3, ] <- tapply(Dets2[i3, ], x, max)
              }

              # randomly select cameras to use...
              camsUsed <- sort(rank(runif(25))[1:nCams]) # using nCams random camera traps.
              Dets3 <- Dets3[c(camsUsed, camsUsed + 25), ]
              cov2 <- cov[c(camsUsed, camsUsed + 25), ]

              Ppresence <- mean(rowSums(Dets3) > 0)

              umf <- unmarkedFrameOccu(y = Dets3, siteCovs = cov2)
              m1 <- occuRN(~1 ~ location, umf)
              s1 <- summary(m1)

              DATA <- data.frame(dens1 = dens1, dens2 = dens2, dT = 1, StudyDuration = SD, nCams = nCams, z = s1$state$z[2], P = s1$state$`P(>|z|)`[2], Ppresence = Ppresence)
              write.table(DATA, outputfile, append = T, col.names = F, row.names = F)
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

# blabla
library(lintr)
