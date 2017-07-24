# The function runs an MCMC for Bayesian inference of the evolution of a discrete feature in BT3
# The function structures the user input according to the needs of BT3, sends it to BT3, 
# initializes the discrete character evolution and reads its output to R

# IN:
# tree:       a phylogenetic tree of class phylo or multiPhylo
# data:       a data.frame consisting of two columns:
#              - one containing the feature
#              - one named taxa, refering to the tip.label in the phylogenetic tree, ensuring the match between tree and data
# model:      either ER or ARD:
#              - ER assumes the rate of gaining (qNY) or losing a feature (qNY) is the same
#              - ARD assumes the rates to differ
# it:         number of MCMC iterations
# bi:         number of itarations discarded as burn in
# sa:         numbers of iterations to skip between taking samples
# pa:         prior distribution of the parameters (for ARD assumed the same for qNY and qYN)
# rm:         logical, indicates whether to clean up ad remove all temporary files
# reuse_tree  logical, indicates whether to reuse the tree or remove it after each run

# OUT: 
# log_file    summary statistics of the MCMC results
# lh          the marginal likelihood of the model


Discrete_BT3_bayesian <- function (tree, data, model, it = 1e+05, bi = 50000, sa = 100, 
                                   pa = NULL, rm = T, reuse_tree = F, st = 100, st_it = 10000) {
  
  # Create a temporary folder where the output from BT will be stored
  
  bayes_traits_folder <- file.path(getwd(), "bt3_temp")
  if (!file.exists(bayes_traits_folder)) {dir.create(bayes_traits_folder)}
  
  # Structure the input information
  
  if (class(tree) == "phylo") {tree$node.label <- NULL}
  
  if (class(tree) == "phylo") {treelabs = tree$tip.label}
  else if (class(tree) == "multiPhylo") {treelabs = attributes(tree)$TipLabel}
  else {stop("Tree must be of class phylo or multiPhylo")}
  
  if (!(class(data[, 1]) %in% c("character", "factor"))) {stop("First column of data should contain species names.")}
  
  if (!(class(data[, 2]) %in% c("character", "factor")) ||  !all(unique (data[,2]) %in% c("N","Y"))) {
    stop("Second column of data should contain presence (Y) and absence (N) of a feature. 
         NAs cannot be handled!")}
  
  if (length(setdiff(treelabs, data[, 1])) > 0) {
    stop(paste("No match found in the data:", 
               paste(setdiff(tree$tip.label,data[, 1]), collapse = ", ")))}
  
  if (length(setdiff(data[, 1], treelabs)) > 0) {
    stop(paste("No match found in the phylogeny:", 
               paste(setdiff(data[,1], tree$tip.label), collapse = ", ")))}
  
  if (length(setdiff(treelabs, data[, 1])) > 0 | 
      length(setdiff(data[, 1], treelabs)) > 0) {
    stop("Species in your phylogeny and data must match up exactly.")}
  
  if (ncol(data) > 2) {stop("Please provide only one trait.")}
  
  if (!exists(".BayesTraitsPath") | !file.exists(.BayesTraitsPath)) {
    stop("Must define '.BayesTraitsPath' to be the path to BayesTraitsV3 on your computer. 
         For example: .BayesTraitsPath <- User/Desktop/BayesTraitsV3")
  }
  
  # Write input file
  # 1 makes BT3 use a Bayesian MCMC and 2 selects the Discrete model
  input <- c (1, 2)
  
  # Choose a model (ER or ARD)
  if (model == 'ER') {
    # Restrict all parameters to be the same as "qYN"
    input <-  c(input, "resall qYN")
    print(input)} 
  else if (model == 'ARD') {}
  else {stop("Please choose between one of the two models: ER, ARD")}
  
  input = c(input, paste("it", format(it, scientific = F)))
  input = c(input, paste("bi", format(bi, scientific = F)))
  input = c(input, paste("sa", format(sa, scientific = F)))
  
  if (!is.null(pa)) {input = c(input, paste("pa", pa))}
  
  input = c(input, paste0("lf ", bayes_traits_folder, "/BTout"))
  input = c(input, paste0("Stones ", st, " ", st_it))
  input = c(input, "run")
  input_filepath <- file.path(bayes_traits_folder, "/inputfile.txt")
  
  write(input, file = input_filepath)
  
  # Write tree
  tree_filepath <- file.path(bayes_traits_folder, "/BT.current.tree.nex")
  
  if (!reuse_tree | !file.exists(tree_filepath)) {ape::write.nexus(tree, file = tree_filepath, translate = T)}
  
  # Write data
  data_filepath <- file.path(bayes_traits_folder, "/BT.current.data.txt")
  write.table(data, data_filepath, quote = F, col.names = F, row.names = F)
  
  # Execute BT 
  # For Windows
  if (Sys.info()['sysname'] == "Windows") {
    print ("Windows OS detected")  
    shell(cmd = paste(.BayesTraitsPath, 
                      tree_filepath, data_filepath, "<", input_filepath), 
          flag="/c", intern=FALSE, wait=TRUE,
          translate=TRUE, mustWork=FALSE)}
  
  # For Mac
  else if (Sys.info()['sysname'] == "Darwin") {
    print("Mac OS detected")
    system(paste(.BayesTraitsPath, 
                 paste0(bayes_traits_folder, "/BT.current.tree.nex"), 
                 paste0(bayes_traits_folder, "/BT.current.data.txt"),
                 paste0("< ", bayes_traits_folder, "/inputfile.txt")))}
  
  else {stop("Your OS is not supported")} 
  
  # Structure output from BT
  # log file
  log_filepath <- file.path(bayes_traits_folder, "/BTout.Log.txt")
  Skip = grep("Tree No", scan(file = log_filepath, what = "c", 
                              quiet = T, sep = "\n", blank.lines.skip = FALSE)) - 1
  
  log_file = read.table(log_filepath, skip = Skip, sep = "\t", 
                        quote = "\"", header = TRUE)
  log_file = log_file[, -ncol(log_file)]
  
  # Stepping stone sampler (marginal likelkihood)
  stones_filepath <- file.path(bayes_traits_folder, "BTout.Stones.txt")
  
  if (file.exists(stones_filepath)) {
    stones <- read.table(stones_filepath, sep='\t')
    stones <- as.character(stones[nrow(stones), 1])
    stones <- as.numeric(stones)}
  
  # Schedule file
  schedule_filepath <- file.path(bayes_traits_folder, "BTout.Schedule.txt")
  
  skip_schedule <- grep("Accepted", scan(file = schedule_filepath, 
                                         what = "c", quiet = T, sep = "\n", 
                                         blank.lines.skip = FALSE)) - 1
  
  schedule <- read.table(schedule_filepath, skip = skip_schedule, 
                         sep = "\t", quote = "\"", header = TRUE)
  
  if (mean(schedule$X..Accepted < 0.2) > 0.5 & mean(schedule$X..Accepted > 0.4) > 0.5) {
    prop_below <- 100 * round(mean(schedule$X..Accepted < 0.2), 2)
    prop_above <- 100 * round(mean(schedule$X..Accepted > 0.4), 2)
    warning(paste0("The acceptance rate was below .20 in ", 
                   prop_below, "% and above .40 in ", prop_above, 
                   "% of the iterations!"), call. = F)}
  
  # remove all files from temporary folder and delete folder
  if (rm) {
    file.remove(log_filepath)
    file.remove(input_filepath)
    file.remove(data_filepath)
    
    if (file.exists(stones_filepath)) {file.remove(stones_filepath)}
    if (file.exists(schedule_filepath)) {file.remove(schedule_filepath)}
    if (!reuse_tree) {file.remove(tree_filepath) unlink(bayes_traits_folder, recursive = T)}}
  
  return(list(log_file = log_file, lh = stones))}
