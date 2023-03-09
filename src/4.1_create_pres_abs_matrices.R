library(lubridate)
m1 <- read.table('./data/raw/FSC and nonFSC data/Congo_data_2021_09_10.csv', 
                 sep=',', header=TRUE)
m2 <- read.table('./data/raw/FSC and nonFSC data/Gabon_data_2021_09_09.csv', 
                 sep=',', header=TRUE)
m <- rbind(m1, m2)

# how many cameras? And what's their survey period?
n_cams <- length(unique(m$Unique_CT_Name))
max_effort <- max(m$Effort..days..New)
max_effort_per_cam <- tapply(m$Effort..days..New, m$Unique_CT_Name, max)

# we create a template, which can be used for each species:
pres_abs_template <- matrix(NA, n_cams, max_effort)
for (i in 1:n_cams){
  pres_abs_template[i, 1:max_effort_per_cam[i]] <- 0
}
cov <- data.frame(cam = sort(unique(m$Unique_CT_Name)), 
                           cluster = tapply(m$Cluster, m$Unique_CT_Name, min),
                           site_type = NA)
for (i in cov$cam){
  cov$site_type[cov$cam == i] <- m$Type.of.Site[m$Unique_CT_Name == i][1]
}
write.table(cov, './data/processed/FSC and nonFSC data/covariates.txt', 
            append=FALSE, col.names = TRUE, row.names = FALSE)

# we will only focus on mammals in this study:
m <- m[m$Class == "MAMMALIA",]
m <- m[m$Species != "",]
m <- m[m$Species != "unknown",]
m <- m[m$Species != "unknown ",]

genus <- gsub(' ','',m$Genus)
species <- gsub(' ', '', m$Species)
species <- paste(genus, species, sep=' ')
n <- tapply(species, species, length)

# write data to new input files, as summarized per species per camera per day:
uSpec <- sort(unique(species))
uSpec <- uSpec[n > 150]
write.table(uSpec, './data/processed/FSC and nonFSC data/species.txt',
            append=FALSE, row.names = FALSE, col.names = FALSE)
for (spec in uSpec)
{
  m2 <- m[species == spec,]
  
  date1 <- as.Date(m2$Camera.Start.Date.New)
  date2 <- as.Date(m2$Photo.Date.New)
  x <- interval(ymd(date1),ymd(date2))
  m2$x <- x %/% days(1)
  
  group <- paste(m2$Unique_CT_Name, m2$x, sep='-')
  i <- 1:length(group)
  i <- tapply(i, group, min)
  
  m3 <- m2[i,]
  pres_abs <- pres_abs_template

  for (i in 1:length(m3$x)) {
    pres_abs[m3$Unique_CT_Name[i] == cov$cam, m3$x[i]] <- 1
  }
  
  filename <- paste('./data/processed/FSC and nonFSC data/pres_abs', 
                    spec, '.txt')
  write.table(pres_abs, filename, 
              row.names=FALSE, col.names=FALSE, append=FALSE)
}


