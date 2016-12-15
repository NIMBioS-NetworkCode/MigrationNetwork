#################################
##    MODEL FUNCTIONS          ##
#################################

############## CALCULATE f -> NODE DYNAMICS UPDATE #################

f_function <- function() {
  # TYPICAL VARIABLES USED (N,alpha,type,ind,t)
  # N - contains the population arriving at the node at all previous time steps - N[["time_step"]]$"class_name"
  # NOTE if you specify "class_name" the it will look to that specific class's population each time
  # NOTE if you set "class_name"=[[type]] it will get the numbers for the current calculation's class
  # alpha - contains all information about node dynamics - alpha[["class_number"]][["season"]]$"variable_name"
  # type - gives the class number as an integer - the class types are defined in NETNAME 
  # NOTE if NETNAME <- c("seeds","plants") then type 1 = seeds and type 2 = plants
  # ind - gives the season number as an integer
  # t - gives the time step
  
  # Initialized Parts of f_update
  trans <- 0*N[[t]][[type]]
  repro <- 0*N[[t]][[type]]
  
  
  # Density dependent adult survival rate:
  SI_A <- alpha[[2]][[ind]]$s0 # Baseline adult survival rate type=2
  SI_C <- alpha[[1]][[ind]]$s0 # Baseline calf survival rate type=1
  K <- alpha[[1]][[ind]]$K+alpha[[2]][[ind]]$K # Total carying capacity
  Na <- N[[t]][[2]] # Number of adults pre
  Nc <- N[[t]][[1]] # Number of calves pre
  R <- alpha[[2]][[ind]]$r # Adult reproductive rate
  
  SI_A <- sqrt(SI_A*exp(-0.219*((Na+Nc)/(K))^(3.77))) # Adult node survival rate - density dependent
  
  
  # Node survival properties:
  if(type==1){ # calves
    survi <- SI_C*Nc
    if(ind == 2){ 
      survi <- 0*Nc # All one year calves transition to adults in the summer.
    }
  }
  
  if(type==2){ # adult
    survi <- SI_A*Na
  }
    
  # Node class transitions (reproduction) (transition)
  if(ind == 2){ # summer
    if(type == 1){ # calves
      repro <- R*SI_A*Na # Number of calves born
      }
    if(type == 2){ # adults
      calves_survival <- SI_C*exp(1-(Na+Nc)/(K)) # Density dependent survival
      trans <-calves_survival*Nc*0.8 # Number of calves transitioning to adult
    }
  }
  
  f_new <- survi+repro+trans
  # Units of f_new should be total number of migrants leaving the node
  
  return(f_new)
}
########################################################################## 

############## CALCULATE p -> PATH TRANSITION RATES #################
p_function <- function(){
  # TYPICAL VARIABLES USED (N,f_update,alpha,beta,type,ind,t)
  # N - contains the population arriving at the nodes at all previous time steps - N[["time_step"]]$"class_name"
  # NOTE this population count is before accounting for node dynamics
  # f_update - contains the population after node dyanamics have been applied for the current time step - f_update[["class_number"]]
  # alpha - contains all information about node dynamics - alpha[["class_number"]][["season"]]$"variable_name"
  # beta - contains all information about path dynamics - beta[["class_number"]][["season]]$"variable_name"
  # NOTE edge transition rates are found in beta[["class_number"]][["season]]$p_ij
  # type - gives the class number as an integer - the class types are defined in NETNAME 
  # for example: if NETNAME <- c("seeds","plants") then type 1 = seeds and type 2 = plants
  # ind - gives the season number as an integer
  # t - gives the time step
  # NOTE if it is required to change any variable globally then use "old_value" <<- "new_value"
  
  p <- beta[[type]][[ind]]$p_ij
  n <- length(p)
  resident <- matrix(0,n,n)
  
  if(ind == 1){ # winter
    if(t==1){
      Npre3 <- N[[1]][[type]][3]
      Npre2 <- N[[1]][[type]][2]
    
      resident[3,3] <- Npre3-0.13/0.87*Npre2
    }
    else{
      resident[3,3] <- M[[t-1]][[type]][3,3]
    }
    
    Npre3 <- N[[t]][[type]][3]
    if(Npre3 > 0){
      p[3,3] <- resident[3,3]/Npre3
      p[3,1] <- 1 - resident[3,3]/Npre3
    }
    else{
      p[3,3] <- 0
      p[3,1] <- 0
    }
  }
  
  p_new <- p
  
  return(p_new)
  # p_new should be unitless and represent the edge transition probabilities for the given step
  
}
##########################################################################

############## CALCULATE s -> PATH SURVIVAL RATES #################
s_function <- function(){
  # TYPICAL VARIABLES USED (N,f_update,alpha,beta,type,ind,t)
  # N - contains the population arriving at the nodes at all previous time steps - N[["time_step"]]$"class_name"
  # NOTE this population count is before accounting for node dynamics
  # f_update - contains the population after node dyanamics have been applied for the current time step - f_update[["class_number"]]
  # alpha - contains all information about node dynamics - alpha[["class_number"]][["season"]]$"variable_name"
  # beta - contains all information about path dynamics - beta[["class_number"]][["season]]$"variable_name
  # NOTE edge survival rates are found in beta[["class_number"]][["season]]$s_ij
  # type - gives the class number as an integer - the class types are defined in NETNAME 
  # NOTE if NETNAME <- c("seeds","plants") then type 1 = seeds and type 2 = plants
  # ind - gives the season number as an integer
  # t - gives the time step
  # NOTE if it is required to change any variable globally then use "old_value" <<- "new_value"
  
  s <- beta[[type]][[ind]]$s_ij
  
  s_new <- s
  # s_new should be unitless and represent the edge survival for the given step
  
  return(s_new)
}
##########################################################################
