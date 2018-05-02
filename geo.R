# The function runs the Geo model in Bayes Traits 3 (BT3)
# It structures the user input according to the needs of BT3, sends it to BT3, 
# initializes geo model and reads its output to R

# IN:
# tree       a phylogenetic tree of class phylo
# data       a data.frame containing three columns
#            taxa: the tip.labels, ensuring the match between tree and data
#            long: the longitude of the tips
#            lat: the latitude of the tips
# it         number of MCMC iterations
# bi         number of iterations discarded as burn in
# sa         numbers of iterations to skip between taking samples
# pr         prior distribution of the parameters (alpha and scale)
# rm         logical, indicates whether to clean up ad remove all temporary files
# reuse_tree logical, indicates whether to reuse the tree or remove it after each run

# OUT: 
# log_file          summary statistics of the MCMC results
# lh                the marginal likelihood of the model
# ancestral_states  summary statistics of the geo  model
# ancestral_nodes   coordinates of the ancestral nodes

Geo_BT3 <- function (tree, data, it = 1e+05, 
                     bi = 5000, sa = 100, pr = NULL,
                     silent = TRUE, rm = T, 
                     reuse_tree = F, st = 100, st_it = 10000) {
  
  # Create a temporary folder where the output from BT will be stored
  
  dir.create(file.path(getwd(), "bt3_temp"))
  bayes_traits_folder <- file.path(getwd(), "bt3_temp")
  
  # Structure the input information
  
  if (class(tree) == "phylo") {
    tree$node.label <- NULL
  }
  if (class(tree) == "phylo") {
    treelabs = tree$tip.label
  }
  else if (class(tree) == "multiPhylo") {
    stop("Tree must be of class phylo. multiPhylo is currently not supported by the Geo model in BT3")
  }
  else {
    stop("Tree must be of class phylo")
  }
  if (!(class(data$taxa) %in% c("character", "factor"))) {
    stop("Names should contain species names.")
  }
  if (length(setdiff(treelabs, data$taxa)) > 0) {
    stop(paste("No match found in the data:", 
               paste(setdiff(tree$tip.label,
                             data$taxa), collapse = ", ")))
  }
  if (length(setdiff(data$taxa, treelabs)) > 0) {
    stop(paste("No match found in the phylogeny:", 
               paste(setdiff(taxa, tree$tip.label), collapse = ", ")))
  }
  if (length(setdiff(treelabs, data$taxa)) > 0 | 
      length(setdiff(data$taxa, treelabs)) > 0) {
    stop("Species in your phylogeny and data must match up exactly.")
  }
  if (!is.numeric(data$long)) {
    stop("Numeric coordinates expected (long).")
  }
  if (!is.numeric(data$lat)) {
    stop("Numeric coordinates expected (lat).")
  }
  if (!exists(".BayesTraitsPath") | !file.exists(.BayesTraitsPath)) {
    stop("Must define '.BayesTraitsPath' to be the path 
         to BayesTraitsV3 on your computer. For example: 
         .BayesTraitsPath <- User/Desktop/BayesTraitsV3")
  }
  
  # Write input file
  
  model = 13
  input = c(model)
  input = c(input, paste("it", format(it, scientific = F)))
  input = c(input, paste("bi", format(bi, scientific = F)))
  input = c(input, paste("sa", format(sa, scientific = F)))
  
  if (!is.null(pr)) {
    for (i in 1:length(pr)) {
      input = c(input, paste("prior", pr[i]))
    }
  }
  
  input = c(input, paste0("lf ", bayes_traits_folder, "/BTout"))
  input = c(input, paste("Stones ", st, " ", st_it))
  input = c(input, "run")
  input_filepath <- file.path(bayes_traits_folder, "/inputfile.txt")
  
  write(input, file = input_filepath)
  
  # Write tree
  
  tree_filepath <- file.path(bayes_traits_folder, "/BT.current.tree.nex")
  
  if (!reuse_tree | !file.exists(tree_filepath)) {
    ape::write.nexus(tree, file = tree_filepath, 
                     translate = T)
  }
  
  # Write data
  
  data_filepath <- file.path(bayes_traits_folder, "/BT.current.data.txt")
  data <- data.frame(data$taxa, data$long, data$lat)
  write.table(data, data_filepath, quote = F, col.names = F, row.names = F)
  
  # Execute BT 
  # For Windows
  
  if (Sys.info()['sysname'] == "Windows") {
    print ("Windows OS detected")  
    shell(cmd = paste(.BayesTraitsPath, tree_filepath, data_filepath, "<", input_filepath), 
          flag="/c", intern=FALSE, wait=TRUE,
          translate=TRUE, mustWork=FALSE)
  }
  
  # For Mac
  
  else if (Sys.info()['sysname'] == "Darwin") {
    print("Mac OS detected")
    system(paste(.BayesTraitsPath, paste0(bayes_traits_folder, "/BT.current.tree.nex"), 
                 paste0(bayes_traits_folder, "/BT.current.data.txt"),
                 paste0("< ", bayes_traits_folder, "/inputfile.txt")), ignore.stdout = silent)
  }
  else {
    stop("Your OS is not supported")
  } 
  
  # Structure output from BT
  # log file
  
  log_filepath <- file.path(bayes_traits_folder, "/BTout.Log.txt")
  Skip <- grep("Tree No", scan(file = log_filepath, what = "c", 
                               quiet = T, sep = "\n", blank.lines.skip = FALSE)) - 1
  
  log_file <- read.table(log_filepath, skip = Skip, sep = "\t", 
                         quote = "\"", header = TRUE)
  log_file <- log_file[, -ncol(log_file)]
  
  # stepping stone sampler (marginal likelkihood)
  
  stones_filepath <- file.path(bayes_traits_folder, "BTout.Stones.txt")
  
  if (file.exists(stones_filepath)) {
    stones <- read.table(stones_filepath, sep='\t')
    stones <- as.character(stones[nrow(stones), 1])
    stones <- as.numeric(stones)
  }
  
  # schedule file
  
  schedule_filepath <- file.path(bayes_traits_folder, "BTout.Schedule.txt")
  
  
  skip_schedule <- grep("Accepted", 
                        scan(file = schedule_filepath, 
                             what = "c", quiet = T, sep = "\n", blank.lines.skip = FALSE)) - 1
  schedule <- read.table(schedule_filepath, skip = skip_schedule, 
                         sep = "\t", quote = "\"", header = TRUE)
  
  if (mean(schedule$X..Accepted < 0.2) > 0.5 & mean(schedule$X..Accepted > 0.4) > 0.5) {
    prop_below <- 100 * round(mean(schedule$X..Accepted < 0.2), 2)
    prop_above <- 100 * round(mean(schedule$X..Accepted > 0.4), 2)
    warning(paste0("The acceptance rate was below .20 in ", 
                   prop_below, "% and above .40 in ", prop_above, 
                   "% of the iterations!"), call. = F)
  }
  
  # Ancestral states file
  
  anc_states_filepath <- file.path(bayes_traits_folder, 
                                   "BTout.AncStates.txt")
  
  if (file.exists(anc_states_filepath)) {
    doc_break = grep("Itter", scan(file = anc_states_filepath, 
                                   what = "c", quiet = T, sep = "\n", 
                                   blank.lines.skip = F)) - 1
    anc_states <- read.table(anc_states_filepath, 
                             skip = doc_break, sep = '\t', header = T)
    
    anc_col_names <- read.table(anc_states_filepath, 
                                sep = '\t', nrow = 1, skip = doc_break, 
                                stringsAsFactors = F)
    
    anc_nodes_df <- read.table(anc_states_filepath, sep = '\t',
                               header = F, nrow = doc_break, 
                               fill = T, na.strings=c("","NA"))
    
    row.names(anc_nodes_df) <- anc_nodes_df[, 1]
    anc_nodes_df <- anc_nodes_df[, 2:ncol(anc_nodes_df)]
    anc_nodes <-  setNames(split(anc_nodes_df,
                                 seq(nrow(anc_nodes_df))), 
                           rownames(anc_nodes_df))
    
    anc_nodes <- lapply(anc_nodes, function(x) x[!is.na(x)])
    summary_anc_states <- tail(anc_states[, 1:ncol(anc_states)], 1) 
    colnames(summary_anc_states) <- anc_col_names
  }
  
  # remove all files from temporary folder and delete folder 
  
  if (rm) {
    file.remove(log_filepath)
    file.remove(input_filepath)
    file.remove(data_filepath)
    file.remove(anc_states_filepath)
    if (file.exists(stones_filepath)) {
      file.remove(stones_filepath)
    }
    if (file.exists(schedule_filepath)) {
      file.remove(schedule_filepath)
    }
    if (!reuse_tree) {
      file.remove(tree_filepath)
      unlink(bayes_traits_folder, recursive = T) 
    }
  }
  
  return(list(log_file = log_file, lh = stones, ancestral_states = summary_anc_states, 
              ancestral_nodes = anc_nodes))
  
  }