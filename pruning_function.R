# The function prunes a phylogenetic tree according to a set of language data, 
# such that it keeps only those parts in the tree and data that match by glottocode and for which there are data. 
# The tree must have the glottocode as a tip.label, the data must have a column containing the glottocode.

# IN:
# data:     a data.frame with linguistic features for a set of languages
# tree:     a phylogenetic tree of class phylo or multiPhylo with tip.labels containing the glottocode
# prune_by: names of the features in the data used for pruning 
#           (i.e. if one! of these features of a language is NA, the corresponding tip is removed from the tree)
# glott:    name of the column in the data containing the glottocode  

# OUT: 
# a list containing the pruned data and the pruned tree

pruning <- function(data, tree, pruned_by, glott) {
  
  # check if the class is a phylogeny 
  if(class(tree) %in% 'phylo') { 
    tree.glottocodes <- tree$tip.label 
  } 
  # .. or multiphylogeny
  else if (class(tree) %in% 'multiPhylo') {
    
    attributes(tree)$TipLabel <- tree[[1]]$tip.label	  
    tree.glottocodes <- tree$TipLabel$tip.label
  }
  else {stop("Tree must be of class phylo or multiPhylo")}
  
  # keep only those data that are in the tree AND are not NA in all data columns
  # the glottocode is the glue for pruning
  vg.data <- data[rowSums(is.na(data[c(pruned_by)])) == 0, ]
  vg.data <- subset(vg.data, vg.data[[glott]] %in% tree.glottocodes)
  vg.data$taxa <- vg.data[[glott]]
  
  if(class(tree) %in% 'phylo') {
    tree.data <- drop.tip(tree, setdiff(tree$tip.label, vg.data[[glott]]))
  }
  
  else {
    tree.data <- lapply(tree, function(t) {
      drop.tip(t, setdiff(t$tip.label, vg.data[[glott]]))
    })
    class(tree.data) <- 'multiPhylo'	
    attributes(tree.data)$TipLabel <- tree.data[[1]]$tip.label
  }
  
  return(list(data = vg.data[, names(vg.data) != glott], tree=tree.data))
  
}