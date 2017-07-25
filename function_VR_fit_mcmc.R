# The function runs an MCMC for Bayesian inference of the evolution of a discrete feature in BT3 
# using the equal rates (ER),  all rates different (ARD) or Variable rates (VR) model with various 
# repetitions.
# The function structures the user input according to the needs of BT3, sends it to BT3, 
# initializes the discrete character evolution and reads its output to R

# IN:
# tree:       a phylogenetic tree of class phylo or multiPhylo
# data:       a data.frame consisting of two columns:
#              - one containing the feature
#              - one named taxa, refering to the tip.label in the phylogenetic tree, ensuring the match between tree and data
# model:      either ER, ARD or VR:
#              - ER assumes the rate of gaining (qNY) or losing a feature (qNY) is the same
#              - ARD assumes the rates to differ
# it:         number of MCMC iterations
# bi:         number of iterations discarded as burn in
# sa:         numbers of iterations to skip between taking samples
# pa:         prior distribution of the parameters (for ARD assumed the same for qNY and qYN)
# rm:         logical, indicates whether to clean up ad remove all temporary files
# reuse_tree  logical, indicates whether to reuse the tree or remove it after each run
# rep_mcmc    number of replications of the Discrete_BT3_VR_simple function

# OUT: 
# VR_test     all ouputs from each individual model
# VR_part_no  the number of partitions in every model
# VR_summary  summary statistics of the replicated model
# VR_PIE    
# VR_best     the model with the highest likelihood

VR_fit_mcmc <- function(feature, tree, prior, iterations_mcmc, 
                        burn_in, sampling_rate, reuse_tree=TRUE,
                        rep_mcmc) {
  
  
  test_VR <- replicate(n = rep_mcmc, 
                       expr = Discrete_BT3_VR_simple(data = tree$data[, c("taxa", feature)], 
                                                     tree = tree$tree, model = "VR",
                                                     it = iterations_mcmc, bi = burn_in, sa = sampling_rate, pa = prior, 
                                                     reuse_tree = reuse_tree),
                       simplify = F)
  
  
  PIE_VR <- apply(do.call(rbind,lapply(test_VR, function(x) 
    tail(x$log_file, 1)[, 8:9])), 2, median)
  
  names(PIE_VR)  <- c("Prob.N","Prob.Y")
  
  #number of partitions
  VR_part_no<-as.vector(do.call(rbind,lapply(test_VR, function(x) 
    nrow(x$var_rates$tree_partitions))))
  
  #list of partitions
  VR_part_list<-do.call(rbind,lapply(test_VR, function(x) 
    list(x$var_rates$tree_partitions)))
  
  names(VR_part_list)<-paste("partitions_rep_","",1:rep_mcmc, sep = "")
  
  
  summary_VR <- list(
    median = data.frame(
      qNY = median(do.call(rbind,lapply(test_VR, function(x) tail(x$log_file$qNY, 1)))),
      qYN = median(do.call(rbind,lapply(test_VR, function(x) tail(x$log_file$qYN, 1))))),
    mean =  data.frame(
      qNY = mean(do.call(rbind,lapply(test_VR, function(x) tail(x$log_file$qNY, 1)))),
      qYN = mean(do.call(rbind,lapply(test_VR, function(x) tail(x$log_file$qYN, 1))))),
    sd = data.frame(
      qNY = sd(do.call(rbind,lapply(test_VR, function(x) tail(x$log_file$qNY, 1)))),
      qYN = sd(do.call(rbind,lapply(test_VR, function(x) tail(x$log_file$qYN, 1))))),
    mad = data.frame(
      qNY = mad(do.call(rbind,lapply(test_VR, function(x) tail(x$log_file$qNY, 1)))),
      qYN = mad(do.call(rbind,lapply(test_VR, function(x) tail(x$log_file$qYN, 1))))))
  
  #compute best 
  #first all the likelihoods
  all_Lh<-do.call(rbind,lapply(test_VR, function(x) tail(x$log_file$Lh, 1)))
  #position of best
  pos_best<-which(all_Lh == max(all_Lh, na.rm = TRUE))
  Lh_best<-max(all_Lh, na.rm = TRUE)
  best_info<-list(Lh_best, pos_best)
  names(best_info)<-c("Highest_Lh", "pos_highest_Lh")
  
  
  best_VR<-list(best_info, test_VR[[pos_best]]$var_rates)
  names(best_VR)<-c("Lh_Info", "Var_Rates_best")
  
  
  return(list(VR_test=test_VR,
              VR_part_no = VR_part_no,
              VR_part_list = VR_part_list,
              VR_summary=summary_VR,
              VR_PIE=PIE_VR,
              VR_best = best_VR
  ))
}