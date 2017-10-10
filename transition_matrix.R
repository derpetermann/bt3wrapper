transition_matrix <- function(mcmc_res, feature, type = "median", prob = T, years = 100) 
  {
  mat_names = list(c('N','Y'),c('N','Y'))
  
  # Define the rate matrix
  x <- matrix(c(-mcmc_res[[feature]][["ARD"]][["ARD_rate"]][[type]]$qNY,  
                mcmc_res[[feature]][["ARD"]][["ARD_rate"]][[type]]$qNY,
                mcmc_res[[feature]][["ARD"]][["ARD_rate"]][[type]]$qYN, 
                -mcmc_res[[feature]][["ARD"]][["ARD_rate"]][[type]]$qYN),
         byrow = T, ncol = 2, dimnames = mat_names)
  
  # Compute the state probabilities after xx years
  if (prob) {
    m <- matrix(matexpo(years*x), byrow=F, ncol=2, dimnames=mat_names)
  } 
  else {m <- x}
  
  return(m)}
