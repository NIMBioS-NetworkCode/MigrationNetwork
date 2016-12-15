#################################
##       NETWORK CODE          ##
#################################
## This code first sets up the network variables to be used in the Network Simulation.
## It then calls NetworkSetup.R to set up the network 
## and calls NetworkSimulation.R to run simulations until an equilibrium population is reached.
## Finally, it calls NetworkOutputs.R to display some basic results.

library(rJava)
library(XLConnect)
library(lattice)
library(R.matlab)
library(xtable)
library(data.table)


### USERS MUST CHANGE THIS PART OF THE CODE TO SET UP SPECIES-SPECIFIC NETWORK INFO ###

#########################################
### SET SPECIES SPECIFIC NETWORK INFO ###
#########################################

### VERY BASIC NETWORK EXAMPLE ###

# IMPORTANT NOTE #
# Users should specify network model functions f_(i,t), p_(ij,t), and s_(ij,t) in the file:
# SpeciesFunctions.R, even if these are constants they must be set as such.

### User specified information specific to the simulation ###

SIMNAME <- "Baseline1" # Specifies the name of the subfolder where 
# the network data, outputs, and .RData file are saved
# This is also the subfolder where the network input files are stored: 
#   ./SIMNAME/network_inputs_NETNAME.xlsx

NETNAME <- c("example") # Give a distinct name for each class as used in input files
# Order is important here we would index [[1]] = class 1 and [[2]] = class 2 in alpha and beta.

seasons <- 2 # Number of seasons or steps in one annul cycle. 
# This must match number of spreadsheets in input files

num_nodes <- 3 # Number of nodes in the network
# This must match the number of initial conditions given in input files

ERR <- .01 # Error tolerance for convergence. 
# To test convergence, we compare total population of all classes in the current season
# to the matching season in the previous year.

tmax <- 200 # Maximum number of steps to take - assume non convergence if t=tmax

OUTPUTS <- TRUE # TRUE = Process final outputs, FALSE = Do not process just run the sumulation.

## For debugging your model equations ##
SILENT <- TRUE # TRUE = Do not print data to console - silence outputs.
# FALSE = Print population data and network function data to the Console for debugging.






### Users should not need to interact with the code below ###

################
## SIMULATION ##
################

# Set the working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Set location of source code
netcode <- c("../NetworkCode1.1/")
# Clear the workspace reserving needed network input variables
base_variables <- c("seasons", "num_nodes", "NETNAME", "tmax", "SIMNAME", "ERR", "OUTPUTS", "SILENT","netcode","base_variables")
rm(list=setdiff(ls(), base_variables))

### SET UP THE NETWORK(S) ###
source(paste(netcode,"NetworkSetup.R",sep=""))

### RUN THE BASELINE SIMULATION ###
print(paste("Running", SIMNAME, sep=" "))
source(paste(netcode,"NetworkSimulation.R",sep=""))

########################  
###  PROCESS OUTPUTS ###
########################

if (OUTPUTS == T){
  source(paste(netcode,"NetworkOutputs.R",sep=""))
}


######################
### SAVE THE DATA ####
######################

save.image(file=paste(SIMNAME,"/",SIMNAME, ".RData", sep = ""))

