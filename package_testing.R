# To be tested: gen3sis, FAVITES, SLiM, MESS (TESS?)

###################
##### gen3sis #####
###################

# O Hagen, B Flueck, F Fopp, JS Cabral, F Hartig, M Pontarp, TF Rangel, L Pellissier (2021) gen3sis: A general engine for eco-evolutionary simulations of the processes that shape Earth’s biodiversity. PLOS Biology. doi:10.1371/journal.pbio.3001340

library(gen3sis)
library(raster)
??gen3sis

# Available data:
  # South America <40MYA with temperature, precipitation, and area (not sure what area represents)
  # Worldwide temperature data? The documentation for this is not clear
  # Can use 'raster' package to acquire bricks (.grd) of existing data

# Data handling:
  # Seems relatively easy to feed in paleogeographic rasters and format into landscapes (a gen3sis object) with create_input_landscape()
  # Seems relatively easy to pull information out of landscapes for use in our models (landscapes are just arrays)

# Simulation:
  # Need to specify cost function defining how dispersal between sites is penalized per grid cell *in addition* to distance, will then create distance matrices based on landscape objects and cost functions

## Creating raster bricks out of the data from the gen3sis GitHub
datapath <- "./genesis/rasters"
temperature_brick <- brick(paste(datapath,"/temperature_rasters.grd",sep=""))
aridity_brick <- brick(paste(datapath,"/aridity_rasters.grd",sep=""))
area_brick <- brick(paste(datapath,"/area_rasters.grd",sep=""))

## Putting the raster bricks into list form
landscapes_list <- list(temp = NULL, arid = NULL, area = NULL)
for (i in 1:nlayers(temperature_brick)) {
  landscapes_list$temp <- c(landscapes_list$temp, temperature_brick[[i]])
  landscapes_list$arid <- c(landscapes_list$arid, aridity_brick[[i]])
  landscapes_list$area <- c(landscapes_list$area, area_brick[[i]])
}

## Establishing a simple cost function
cost_function_null <- function(source, habitable_src, dest, habitable_dest) {
  return(1/1000)
}

## A more reasonable cost function
cost_function_water <- function(source, habitable_src, dest, habitable_dest) {
  if (!all(habitable_src, habitable_dest)) {
    return(2/1000)
  } else {
    return(1/1000)
  }
}

## Creating a landscape object out of the raster list and cost function
create_input_landscape(landscapes = landscapes_list, cost_function = cost_function_water, directions = 8, output_directory =  "./genesis/simulation/landscape", timesteps = paste0(seq(65, 0, by = -1), "Ma"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", calculate_full_distance_matrices = F)

## Creating a config object from the example config file
config_object <- create_input_config("./genesis/simulation/config_southamerica.R")

## Running the example simulation
sim <- run_simulation(config = config_object,
                          landscape = "./genesis/simulation/landscape",
                          output_directory="./genesis/simulation/output",
                          call_observer = 4,
                          verbose=0)

#################
##### slimr #####
#################

# https://github.com/rdinnager/slimr

library(slimr)

# Have to install SLiM first --> there is an install tool for this on Mac

# Example slim script from slimr github
slim_script(
  slim_block(initialize(),
             {
               ## set the overall mutation rate
               initializeMutationRate(1e-7);
               ## m1 mutation type: neutral
               initializeMutationType("m1", 0.5, "f", 0.0);
               ## g1 genomic element type: uses m1 for all mutations
               initializeGenomicElementType("g1", m1, 1.0);
               ## uniform chromosome of length 100 kb
               initializeGenomicElement(g1, 0, 99999);
               ## uniform recombination along the chromosome
               initializeRecombinationRate(1e-8);
             }),
  slim_block(1,
             {
               sim.addSubpop("p1", 500);
             }),
  slim_block(10000,
             {
               sim.simulationFinished();
             })
) -> script_1
slim_run(script_1)


