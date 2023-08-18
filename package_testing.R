# To be tested: gen3sis, FAVITES, SLiM, MESS (TESS?)

###################
##### gen3sis #####
###################

# O Hagen, B Flueck, F Fopp, JS Cabral, F Hartig, M Pontarp, TF Rangel, L Pellissier (2021) gen3sis: A general engine for eco-evolutionary simulations of the processes that shape Earth’s biodiversity. PLOS Biology. doi:10.1371/journal.pbio.3001340

library(gen3sis)
library(raster)
??gen3sis
setwd("~/Projects/simulation_tests/package_testing-master")
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

gen3sis::save_phylogeny()

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
create_input_landscape(landscapes = landscapes_list,  
                       cost_function = cost_function_water, 
                       directions = 8, output_directory =  "./genesis/simulation/landscape", 
                       timesteps = paste0(seq(65, 0, by = -1), "Ma"), 
                       crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                       calculate_full_distance_matrices = F)

## Creating a config object from the example config file
config_object <- create_input_config("./genesis/simulation/config_southamerica.R")


#modify observer function to save phylogeny and character matrices at end of each time step.
config_object$gen3sis$general$end_of_timestep_observer=function(data, vars, config){
  save_species()
  save_traits()
  save_phylogeny()
  plot_richness(data$all_species, data$landscape)
  # example 1 plot over simulation
  # par(mfrow=c(2,3))
  # plot_raster_single(data$landscape$environment[,"temp"], data$landscape, "temp", NA)
  # plot_raster_single(data$landscape$environment[,"arid"], data$landscape, "arid", NA)
  # plot_raster_single(data$landscape$environment[,"area"], data$landscape, "area", NA)
  # plot_richness(data$all_species, data$landscape)
  # plot_species_presence(data$all_species[[1]], data$landscape)
  # plot(0,type='n',axes=FALSE,ann=FALSE)
  # mtext("STATUS",1)
  # example 2 plot over simulations saving plots
  # plot_richness(data$all_species, data$landscape)
  # plot_landscape(data$landscape)
  
}





## Running the example simulation
sim <- run_simulation(config = config_object,
                          landscape = "./genesis/simulation/landscape",
                          output_directory="./genesis/simulation/output",
                          call_observer = 4,
                          verbose=0)


#get phylogeny from specific time step
phy=read.nexus("genesis/simulation/output/default_config/202308180408-06045/phylogeny/phylogeny_t_0.nex")

#get traits of all populations of each species at each timestep
traits=readRDS("genesis/simulation/output/default_config/202308180408-06045/traits/traits_t_0.rds")


sgen3sis=readRDS("genesis/simulation/output/default_config/202308180408-06045/sgen3sis.rds")
sgen3sis$summary$phylo_summary

#################
##### slimr #####
#################

# https://github.com/rdinnager/slimr
library(devtools)
install_github("https://github.com/rdinnager/slimr")

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

out=slim_run(script_1, keep_all_output = T, capture_output = T, throw_error = T)


# can specify either a pedigree of individuals or tree sequences (tree for each loci)

#// The following constants need to be defined:
#  // - outfile: path to save the tree sequence of the simulation.
#// - popsize: population size.
#// - num_gens: run simulation for num_gens ticks.
#// - infile: a string with path to tree sequence to start from; optional.

slim_script(
  slim_block(
initialize(), {
  initializeSLiMModelType("WF");
  initializeTreeSeq();
  initializeMutationRate(1e-8);
  initializeMutationType("m1", 0.5, "f", -0.01);
  initializeGenomicElementType("g1", m1, 0.1);
  initializeGenomicElement(g1, 0, 1e6-1);
  initializeRecombinationRate(1e-9);
}),


slim_block(1, late(), {
#  // if no input tree sequence is provided, then start a subpopulation
  if (infile == "") {
    p = sim.addSubpop("p1", popsize);
  } else {
#    // reloading must happen in late()
    sim.readFromPopulationFile(infile);
    parent = sim.subpopulations[0];
    p = sim.addSubpopSplit(max(sim.subpopulations.id) + 1, popsize, parent);
    parent.setSubpopulationSize(0);
  }
  p.name = popname;
}
),

slim_block(
#// schedule the end of the simulation
1, late(), {
  finaltick = num_gens + community.tick;
  community.rescheduleScriptBlock(s0, ticks=finaltick);
}
),
#// event that saves the tree sequence

slim_block(
s0, 1000, late(), {
  sim.treeSeqRememberIndividuals(sim.subpopulations.individuals);
  sim.treeSeqOutput(outfile);
}
))-> script_1
