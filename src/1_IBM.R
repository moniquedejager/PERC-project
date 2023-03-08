cam_trap = function(sim_nr)
{
  # with this script, we model camera trapping events of a species
  
  # parameters that can be changed between simulations:
  speed <- 1   # in 5 m/timestep (1 timestep = 1/3 minute)
  densities <- round(exp(seq(1.5, 6, 0.25))) # vector with 19 different densities,
  # in # individuals per 100km2
  
  # per density, we simulate camera trapping events:
  for (i_dens in length(densities)) {
    # We consider periodic boundary conditions.  
    # In case of low densities (< 10 individuals/100km2), 
    # Individuals move randomly within a space of 7,000 x 7,000 patches, 
    # representing an area of 35 x 35 km 
    # (= 1225 km2, one cell representing a 5 x 5 m area).
    # In case of high densities (> 10 individuals/100km2), 
    # Individuals move randomly within a space of 2,000 x 2,000 patches, 
    # representing an area of 10 x 10 km 
    # (= 100 km2, one cell representing a 5 x 5 m area).)
    dens <- densities[i_dens]
    xmin <- 0
    ymin <- 0
    if (dens < 10) { 
      xmax <- 7000
      ymax <- 7000
      n <- dens * 12.25
    } else {
      xmax <- 2000
      ymax <- 2000
      n <- dens
    }
    
    # Individuals are randomly distributed across the simulation area: 
    x <- runif(n, xmin, xmax)
    y <- runif(n, ymin, ymax)
    
    # set up the 25 camera traps:
    camTrapX <- rep(seq(300, 1900, 400), 5)
    camTrapY <- sort(camTrapX)
    camTrapPos <- paste(camTrapX, camTrapY, sep='-')
    
    # Each simulation represents a period of 365 days. Each time step 
    # represents a 1/3 minute; within a minute, each individual moves 
    # in a random direction. Default speed is 15 m/min. This 
    # can be changed. During nighttime, the individuals stop moving 
    # (= 11 hours per day). 
    # Per day and camera trap, we record whether a camera trapping 
    # event occurred.
    pres_abs <- matrix(0, 25, 365)
    for (i in 1:(365*11*60*3))
    {
      day <- ceiling(i/3/60/11)
      a <- runif(n, 0, 2*pi)
      dx <- speed*cos(a)
      dy <- speed*sin(a)
      x <- x + dx
      y <- y + dy
      
      # if the new location is outside the area, the individual will 
      # appear on the other side of the area:
      x <- x + xmax*(x <= xmin) - xmax*(x > xmax)
      y <- y + ymax*(y <= ymin) - ymax*(y > ymax)
      
      # every timestep, a photo can be made at the camera trap locations:
      pos <- paste(round(x), round(y), sep='-')
      sighting <- (1:25)[camTrapPos %in% pos]
      for (i2 in sighting) {
        pres_abs[i2,day] <- 1
      }
    }
    
    # write data to a file:
    filename <- paste('./data/raw/IBM results/presence_absence', 
                      i_dens, 
                      sim_nr, 
                      '.txt')
    DATA <- cbind(dens, pres_abs)
    write.table(DATA, filename, 
                append=FALSE, col.names = FALSE, row.names = FALSE)
  }
}

# use parallel runs to increase computational speed:
library(parallel)
n_cores <- detectCores()
cl <- makeCluster(min(20, n_cores))
results <- parSapply(cl, 1:100, cam_trap)
stopCluster(cl)


