# The function runs an MCMC for Bayesian inference of the evolution of a discrete feature in BT3 
# using the equal rates (ER),  all rates different (ARD) or Variable rates (VR) model
# The function structures the user input according to the needs of BT3, sends it to BT3, 
# initializes the discrete character evolution and reads its output to R

# IN:
# tree       a phylogenetic tree of class phylo or multiPhylo
# data        a data.frame consisting of two columns:
#              - one containing the feature
#              - one named taxa, refering to the tip.label in the phylogenetic tree, ensuring the match between tree and data
# model       either ER or ARD:
#              - ER assumes the rate of gaining (qNY) or losing a feature (qNY) is the same
#              - ARD assumes qNY and qYN to differ
#              - VR assumes qNY and qYN differ and to vary across the phylogeny
# it          number of MCMC iterations
# bi          number of iterations discarded as burn in
# sa          numbers of iterations to skip between taking samples
# pa          prior distribution of the parameters (for ARD assumed the same for qNY and qYN)
# rm          logical, indicates whether to clean up ad remove all temporary files
# reuse_tree  logical, indicates whether to reuse the tree or remove it after each run

# OUT: 
# log_file    summary statistics of the MCMC results
# lh          the marginal likelihood of the model
# var_rates   summary statistics for the variable rates model (if model = VR)
# node_info   nodes and branches where a rate change occurs (if model = VR)

Discrete_BT3_VR_bayesian <- function (tree, data, model, it = 1e+05, bi = 50000, 
                                      sa = 100, pa = NULL, rm = T, reuse_tree = F, 
                                      st = 100, st_it = 10000) {
  
  # Create a temporary folder where the output from BT will be stored
  
  bayes_traits_folder <- file.path(getwd(), "bt3_temp")
  
  if (!file.exists(bayes_traits_folder)) {dir.create(bayes_traits_folder)}
  
  # Structure the input information
  
  if (class(tree) == "phylo") {tree$node.label <- NULL}
  
  if (class(tree) == "phylo") {treelabs = tree$tip.label}
  else if (class(tree) == "multiPhylo") {treelabs = attributes(tree)$TipLabel}
  else {stop("Tree must be of class phylo or multiPhylo")}
  
  if (!(class(data[, 1]) %in% c("character", "factor"))) {stop("First column of data should contain species names.")}
  if (!(class(data[, 2]) %in% c("character", "factor")) || !all(unique (data[,2]) %in% c("N","Y"))) {
    stop("Second column of data should contain presence (Y) and absence (N) of a feature. 
         NAs cannot be handled!")}
  
  if (length(setdiff(treelabs, data[, 1])) > 0) {
    stop(paste("No match found in the data:", 
               paste(setdiff(tree$tip.label,
                             data[, 1]), collapse = ", ")))}
  
  if (length(setdiff(data[, 1], treelabs)) > 0) {
    stop(paste("No match found in the phylogeny:", 
               paste(setdiff(data[,1], tree$tip.label), collapse = ", ")))}
  
  if (length(setdiff(treelabs, data[, 1])) > 0 | 
      length(setdiff(data[, 1], treelabs)) > 0) {
    stop("Species in your phylogeny and data must match up exactly.")}
  
  if (ncol(data) > 2) {
    stop("Please provide only one trait.")}
  
  if (!exists(".BayesTraitsPath") | !file.exists(.BayesTraitsPath)) {
    stop("Must define '.BayesTraitsPath' to be the path 
         to BayesTraitsV3 on your computer. For example: 
         .BayesTraitsPath <- User/Desktop/BayesTraitsV3")}
  
  # Write input file
  # 1 makes BT3 use a Bayesian MCMC and 2 selects the Discrete model
  input <- c (1, 2)
  
  # Choose a model (ER, ARD, VR)
  if (model == 'ER') {
    # Restrict all parameters to be the same as "qYN"
    input <-  c(input, "resall qYN")}
  
  else if (model == 'ARD') {}
  
  else if (model == 'VR') {input = c(input, "VarRates")}
  
  else {stop("Please choose between one of the three models: ER, ARD, VR")}
  
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
  
  if (!reuse_tree | !file.exists(tree_filepath)) {
    ape::write.nexus(tree, file = tree_filepath, 
                     translate = T)}
  
  # Write data
  data_filepath <- file.path(bayes_traits_folder, "/BT.current.data.txt")
  write.table(data, data_filepath, quote = F, col.names = F, row.names = F)
  
   # Execute BT 
  # For Windows
  if (Sys.info()['sysname'] == "Windows") {
    print ("Windows OS detected")  
    invisible(shell(cmd = paste(.BayesTraitsPath, 
                                tree_filepath, data_filepath, "<", input_filepath), 
                                flag="/c", intern=TRUE, wait=TRUE,
                                translate=TRUE, mustWork=FALSE))}
  
  # For Mac
  else if (Sys.info()['sysname'] == "Darwin") {
    print("Mac OS detected")
    invisible(system(paste(.BayesTraitsPath, 
                 paste0(bayes_traits_folder, "/BT.current.tree.nex"), 
                 paste0(bayes_traits_folder, "/BT.current.data.txt"),
                 paste0("< ", bayes_traits_folder, "/inputfile.txt"))))}
  
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
  
  if (model == "VR"){
    var_rates_filepath <- file.path(bayes_traits_folder, "BTout.VarRates.txt")
    
    if (file.exists(var_rates_filepath)) {
      #the first line shows the number of languages
      num_lang<-as.numeric(readLines(var_rates_filepath, n = 1))
      
      #read all the languages
      lang_names <- read.table(var_rates_filepath, sep='\t', fill=T, colClasses = "character",
                               skip = 1, nrows = num_lang)
      #add colnames
      colnames(lang_names)<-c("tip_id", "lang_name")
      
      #read the number of nodes
      ancest_n<-as.numeric(read.table(var_rates_filepath, sep='\t', fill=T, colClasses = "character",
                                      skip = (num_lang+1), nrows = 1))
      
      #read all nodes
      ancest_info<-read.table(var_rates_filepath, sep='\t', fill=T, colClasses = "character",
                              skip = (num_lang+2), nrows = ancest_n)
      
      #all tips
      tip_info<-subset(ancest_info, V3 <= 1)
      tip_info<-tip_info[,c(1,4)]
      colnames(tip_info)<-c("node_id", "tip_id")
      lang_names<-merge(lang_names, tip_info, all.x = TRUE, by = "tip_id")
      #remove all tips from the nodes
      ancest_info<-subset(ancest_info, V3 > 1)
      ancest_info<-ancest_info[,c(1,3:(3+num_lang))]
      #add colnames
      colnames(ancest_info)[1:2]<-c("node_id", "num_descendants")
      colnames(ancest_info)[3:ncol(ancest_info)]<-paste("descendant_","",1:(ncol(ancest_info)-2),
                                                        sep="")
      #what to give back
      node_info<-list(tip_names = lang_names, node_structure = ancest_info)
      
      
      nc <- max(count.fields(var_rates_filepath, sep = '\t'))
      var_rates_log <- read.table(var_rates_filepath, sep='\t', 
                                  col.names = paste("v", 1:nc,sep="."), fill=T, 
                                  colClasses = "character")
      
      var_rates_log_ll <- var_rates_log[nrow(var_rates_log),]
      
      no_pram <- as.integer(var_rates_log_ll[,4])
      n_inf <- var_rates_log_ll[, 8:(4 * no_pram + 7)]
      
      var_rates <- vector("list", 2)
      names(var_rates) <- c("meta", "tree_partitions")
      
      colnames(var_rates_log_ll)[1:7]<- var_rates_log[which(var_rates_log == "Sigma^2", 
                                                            arr.ind = TRUE)[1],][1:7]
      
      var_rates[[1]] <- var_rates_log_ll[, 1:7]
      
      if(no_pram != 0){
        m1 <- matrix(1:ncol(n_inf), ncol = no_pram)
        lst <- split(m1, row(m1))
        m1 <- reshape(n_inf, varying=lst, direction='long')
        
        colnames(m1)[2:5]<- var_rates_log[which(var_rates_log == "Sigma^2",
                                                arr.ind = TRUE)[1],][8:11]
        var_rates[[2]] <- m1[,2:5]}
      
      else{var_rates[[2]] <- NA}}
    
    var_trees_filepath <- file.path(bayes_traits_folder, "BTout.Output.trees")}
  
  # remove all files from temporary folder and delete folder
  
  if (rm) {
    
    file.remove(log_filepath)
    file.remove(input_filepath)
    file.remove(data_filepath)
    
    if (file.exists(stones_filepath)) {file.remove(stones_filepath)}
    
    if (file.exists(schedule_filepath)) {file.remove(schedule_filepath)}
    
    if (model == "VR"){
      if (file.exists(var_rates_filepath)) {file.remove(var_rates_filepath)}
      if (file.exists(var_trees_filepath)) {file.remove(var_trees_filepath)}}
    
    if (!reuse_tree) {
      file.remove(tree_filepath)
      unlink(bayes_traits_folder, recursive = T)}}
  
  
  if (model == 'VR') {
    return(list(log_file = log_file, lh = stones, var_rates = var_rates, 
                node_info = node_info))}
  else {return(list(log_file = log_file, lh = stones))}}
