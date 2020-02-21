setwd('E:\\Certhia\\certhia_MS\\UR_Total')
r <-read.csv('E:\\Certhia\\certhia_MS\\UR_Total\\UR_total.csv')
colnames(r) <- c('species', 'Longitude', 'Latitude')
u <- unique(r$species)
head(r)

i=1

cols = 1:3 # columns extracted from dataset - species, long, lat
for (i in 1:length(u)) {
  sptable <- r[r$species == u[i], cols]
  write.csv(sptable, paste(paste(u[i], collapse = '_'), ".csv", sep = ""), row.names = FALSE)
  }
  

#Load variables 
#cropping layers in batch------------------
setwd('E:\\Certhia\\bio_2-5m_bil')
require(rgdal)
library(raster)
shape <- readOGR(dsn = "E:\\Certhia\\certhia_MS\\",layer = "certhia") #SHAPEFILE YOU CREATED FOR PROJECTION AREA
path <- "E:\\Certhia\\bio_2-5m_bil\\"
varaibles_list <-list.files(path = path, pattern = ".bil", # vector of variables
                            full.names = TRUE)
variables <- stack(varaibles_list) #create a stack

#variables <- getData("worldclim", var = "bio", res = "2.5")[[c(1, 12)]]

# masking variables
var_mask <- mask(crop(variables, shape), shape)

## names for layers
rnames <- paste0("E:/Certhia/certhia_MS/Certhia_G/", names(variables), ".asc") # users select the format

## saving layers in new folder
sav <- lapply(1:nlayers(var_mask), function(x) {
  writeRaster(var_mask[[x]], filename = rnames[x], format = "ascii", overwrite=T) # change format accordingly
})

# Principal components naalysis--------------

library(kuenm)
suppressWarnings({
  if(!require(raster)){
    install.packages("raster")
    library(raster)
  }
}) 

# simple raster PCA
## functions help
help(kuenm_rpca)

## preparing function's arguments
var_folder <- "E:\\Certhia\\certhia_MS\\Certhia_G" # name of folder with variables to be combined in distinct sets
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil 
scalev <- TRUE # scale variables
writer <- TRUE # save results
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
out_folder <- "E:\\Certhia\\certhia_MS\\certhia_pca" # name of folder that will contain the sets 
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters

## runing PCA
kuenm_rpca(variables = var_folder, in.format = in_format, var.scale = scalev, 
           write.result = writer, out.format = out_format, out.dir = out_folder,
           n.pcs = n_pcs)

#ENMTOOLSS------------
install.packages("devtools")
library(devtools)
install_github("danlwarren/ENMTools")
library(ENMTools)
env.files <- list.files(path = "E:\\Certhia\\certhia_MS\\ENMtools\\env", pattern = "pc", full.names = TRUE)
env <- stack(env.files)
names(env) <- c("pc1", "pc2", "pc3", "pc4")
env <- setMinMax(env)
#americana-------------
america <- enmtools.species(species.name = "america", 
                            presence.points = read.csv("E:\\Certhia\\certhia_MS\\UR_Total\\Certhia americana.csv")[,2:3])
america$range <- background.raster.buffer(america$presence.points, 50000, mask = env)
crs(america$range) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
america$background.points <- background.points.buffer(points = america$presence.points,
                                                   radius = 20000, n = 1000, mask = env[[1]])
america.glm <- enmtools.glm(species = america, env = env, test.prop = 0.2)
visualize.enm(america.glm, env, layers = c("pc1", "pc2"), plot.test.data = TRUE)
raster.breadth(america.glm)

#brachydactyla-------------
brachy <- enmtools.species(species.name = "brachy", 
                            presence.points = read.csv("E:\\Certhia\\certhia_MS\\UR_Total\\Certhia brachydactyla.csv")[,2:3])
brachy$range <- background.raster.buffer(brachy$presence.points, 50000, mask = env)
crs(brachy$range) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
brachy$background.points <- background.points.buffer(points = brachy$presence.points,
                                                      radius = 20000, n = 1000, mask = env[[1]])
brachy.glm <- enmtools.glm(species = brachy, env = env, test.prop = 0.2)
visualize.enm(brachy.glm, env, layers = c("pc1", "pc2"), plot.test.data = TRUE)

#discolor-----------
discolor <- enmtools.species(species.name = "discolor", 
                           presence.points = read.csv("E:\\Certhia\\certhia_MS\\UR_Total\\Certhia discolor.csv")[,2:3])
discolor$range <- background.raster.buffer(discolor$presence.points, 50000, mask = env)
discolor$background.points <- background.points.buffer(points = discolor$presence.points,
                                                     radius = 20000, n = 1000, mask = env[[1]])
discolor.glm <- enmtools.glm(species = discolor, env = env, test.prop = 0.2)
visualize.enm(discolor.glm, env, layers = c("pc1", "pc2"), plot.test.data = TRUE)

#familiaris----------
familiaris <- enmtools.species(species.name = "familiaris", 
                             presence.points = read.csv("E:\\Certhia\\certhia_MS\\UR_Total\\Certhia familiaris.csv")[,2:3])
familiaris$range <- background.raster.buffer(familiaris$presence.points, 50000, mask = env)
familiaris$background.points <- background.points.buffer(points = familiaris$presence.points,
                                                       radius = 20000, n = 1000, mask = env[[1]])
familiaris.glm <- enmtools.glm(species = familiaris, env = env, test.prop = 0.2)

#himalayana-----------
himalayana <- enmtools.species(species.name = "himalayana", 
                               presence.points = read.csv("E:\\Certhia\\certhia_MS\\UR_Total\\Certhia himalayana.csv")[,2:3])
himalayana$range <- background.raster.buffer(himalayana$presence.points, 50000, mask = env)
himalayana$background.points <- background.points.buffer(points = himalayana$presence.points,
                                                         radius = 20000, n = 1000, mask = env[[1]])
himalayana.glm <- enmtools.glm(species = himalayana, env = env, test.prop = 0.2)


#hodgsoni-----------
hodgsoni <- enmtools.species(species.name = "hodgsoni", 
                               presence.points = read.csv("E:\\Certhia\\certhia_MS\\UR_Total\\Certhia hodgsoni.csv")[,2:3])
hodgsoni$range <- background.raster.buffer(hodgsoni$presence.points, 50000, mask = env)
hodgsoni$background.points <- background.points.buffer(points = hodgsoni$presence.points,
                                                         radius = 20000, n = 1000, mask = env[[1]])
hodgsoni.glm <- enmtools.glm(species = hodgsoni, env = env, test.prop = 0.2)

#manipurensis-----------
mani <- enmtools.species(species.name = "mani", 
                             presence.points = read.csv("E:\\Certhia\\certhia_MS\\UR_Total\\Certhia manipurensis.csv")[,2:3])
mani$range <- background.raster.buffer(mani$presence.points, 50000, mask = env)
mani$background.points <- background.points.buffer(points = mani$presence.points,
                                                       radius = 20000, n = 1000, mask = env[[1]])
mani.glm <- enmtools.glm(species = mani, env = env, test.prop = 0.2)

#manipurensis---------
mani <- enmtools.species(species.name = "mani", 
                         presence.points = read.csv("E:\\Certhia\\certhia_MS\\UR_Total\\Certhia manipurensis.csv")[,2:3])
mani$range <- background.raster.buffer(mani$presence.points, 50000, mask = env)
mani$background.points <- background.points.buffer(points = mani$presence.points,
                                                   radius = 20000, n = 1000, mask = env[[1]])
mani.glm <- enmtools.glm(species = mani, env = env, test.prop = 0.2)

#manipurensis---------
mani <- enmtools.species(species.name = "mani", 
                         presence.points = read.csv("E:\\Certhia\\certhia_MS\\UR_Total\\Certhia manipurensis.csv")[,2:3])
mani$range <- background.raster.buffer(mani$presence.points, 50000, mask = env)
mani$background.points <- background.points.buffer(points = mani$presence.points,
                                                   radius = 20000, n = 1000, mask = env[[1]])
mani.glm <- enmtools.glm(species = mani, env = env, test.prop = 0.2)

#nipalensis----------
nipa <- enmtools.species(species.name = "nipa", 
                         presence.points = read.csv("E:\\Certhia\\certhia_MS\\UR_Total\\Certhia nipalensis.csv")[,2:3])
nipa$range <- background.raster.buffer(nipa$presence.points, 50000, mask = env)
crs(nipa$range) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
nipa$background.points <- background.points.buffer(points = nipa$presence.points,
                                                   radius = 5000, n = 500, mask = env[[1]])
nipa.glm <- enmtools.glm(species = nipa, env = env, test.prop = 0.2)
plot(nipa.glm)
#tianquanensis-----------
tia <- enmtools.species(species.name = "tia", 
                         presence.points = read.csv("E:\\Certhia\\certhia_MS\\UR_Total\\Certhia tianquanensis.csv")[,2:3])
tia$range <- background.raster.buffer(tia$presence.points, 50000, mask = env)
crs(tia$range) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(tia$range)
tia$background.points <- background.points.buffer(points = tia$presence.points,
                                                   radius = 5000, n = 100, mask = env[[1]])
tia.glm <- enmtools.glm(species = tia, env = env, test.prop = 0.2)
visualize.enm(tia.glm, env, layers = c("pc1", "pc2"), plot.test.data = TRUE)
#niche breath---------