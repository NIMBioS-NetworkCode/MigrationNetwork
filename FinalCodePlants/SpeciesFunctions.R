#################################
##    MODEL FUNCTIONS          ##
#################################

### PLANTS ###

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
  
  t_node <- alpha[[type]][[ind]]$t_node # Transition between classes
  r_node <- alpha[[type]][[ind]]$r_node # Reproductive rate
  psi_node <- alpha[[type]][[ind]]$psi_node
  K_node <- alpha[[type]][[ind]]$K_node
  
  plant_pop <- N[[t]]$plants  
  plant_survival <- alpha$plants[[ind]]$s_node
  
  seed_pop <- N[[t]]$seeds
  seed_survival <- alpha$seeds[[ind]]$s_node
  
  # SEED UPDATE - type=1 implies seeds
  if( type==1){
    trans <- t_node*seed_survival*seed_pop #Transitioning seeds - attempt germination
    repro <- r_node*plant_pop #Seeds produced by fruiting plants
    survi <- seed_survival*seed_pop #Surviving seeds - seed bank
  }
  
  # PLANT UPDATE - type=2 implies plants
  if( type==2){
    
    seed_trans <- t_node*seed_survival*seed_pop #Total number of transitioning seeds
    Psi <- psi_node*exp(-(plant_pop+seed_trans)/K_node) #Density dependent successful germination rate
    trans <- Psi*seed_trans #Total number of new plants
    repro <- 0 #No plants produced by reproduction
    survi <- plant_survival*plant_pop #Surviving plants
  }
  
  f_new <- survi + repro + trans  
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
  
  p <- beta[[type]][[ind]]$p_ij # Constant transition probabilities
  
  p_new <- p
  # p_new should be unitless and represent the edge transition probabilities for the given step
  
  return(p_new)
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
  
  s <- beta[[type]][[ind]]$s_ij # Constant edge survival rates
  
  s_new <- s
  # s_new should be unitless and represent the edge survival for the given step
  
  return(s_new)
}
##########################################################################
