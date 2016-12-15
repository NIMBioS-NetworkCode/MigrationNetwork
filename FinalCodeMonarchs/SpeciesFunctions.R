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
  
  eggs <- 0.75*358
  mw_param <- alpha[[type]][[ind]]$m
  n <- length(mw_param)
  POP <- N[[t]]
  sL <- matrix(0,1,n)
  
  for (i in 1:n){
    if(mw_param[i]==0) sL[i] <- 0 else sL[i] <- (0.0327*(1+exp(-1.0175)))/(1+exp(-1.0175+0.1972*eggs*POP[[1]][i]/mw_param[i]))
  }
  
  Survi <- alpha[[type]][[ind]]$sA
  sP <- alpha[[type]][[ind]]$sP
  
  f_new <- Survi*(1 + eggs*sP*sL)*N[[t]][[type]]
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
  
  s <- beta[[type]][[ind]]$s_ij
  
  s_new <- s
  # s_new should be unitless and represent the edge survival for the given step
  
  return(s_new)
}
##########################################################################
