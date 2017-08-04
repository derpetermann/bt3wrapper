# TO DO: Description

# IN:
# sum_tree
# mcmc_res
# feature
# output

# OUT
# treemaps  the actual output from the function
# OR
# plot      the visualization of the treemaps

simmap_plot <- function (sum_tree, 
                         mcmc_res, 
                         feature, 
                         output = c("plot", "treemaps")){
  
  require(phytools)
  
  # Create transition matrix
  q_mat <- matrix(c(-median,
                    median,
                    median,
                    -median),
                  byrow = T, ncol = 2, dimnames = list(c('N','Y'), c('N','Y')))
  
  
  states <- sum_tree$data[, feature]
  names(states) <- sum_tree$data$Language
  
  sum_tree$tree$tip.label <- as.character(sum_tree$data$Language[
    match(sum_tree$tree$tip.label, sum_tree$data$taxa)])
  
  # Compute simmap
  treemaps <- make.simmap(sum_tree$tree, states, Q = q_mat, pi = 'estimated',
                          nsim = 100, message = F)
  map <- densityMap(treemaps, plot=F, res=300)
  map$cols[1:length(map$cols)] <- rev(map$cols[1:length(map$cols)])
  
  if(missing(output) || toupper(output) == "TREEMAPS"){
    return(treemaps)
  }
  else if(toupper(output) == "PLOT"){
    # Plot simmap
    plot(map, outline=T, lwd=c(3,3),
         fsize=c(0.6,0.6), ftype='reg',
         legend=2000, leg.txt=c("0",paste0("Posterior Probability (Y)"),"1"))
  }
  
}