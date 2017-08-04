# This function computes the original NEXUS ID's for the partitions found by
# the Discrete_BT3_VR function.

# IN
# tree:           a phylogenetic tree of class phylo or multiPhylo, same that was also
#                 used for the previous function
# vr_result       the result of the Discrete_BT3_VR function

# OUT
# partitions      the individual vr partitions as found by Discrete_BT3_VR with
#                 additional NEXUS IDs
# matching_table  a matching table of the NEXUS IDs and the IDs produced by the
#                 mcmc function
node_id_match<-function(tree, vr_result){
  
  require(phangorn)
  
  ancest_info<-vr_result$node_info$node_structure
  
  lang_names<-vr_result$node_info$tip_names
  
  partitions<-vr_result$var_rates$tree_partitions
  
  #create vector of node labels
  node_label<-ancest_info[,1]
  
  n_length<-ancest_info[,2]
  
  ancest_info<-ancest_info[,3:ncol(ancest_info)]
  
  #create new df
  new <- ancest_info
  
  #test all entries against the lookup table lang_names
  new[] <- lang_names$lang_name[match(unlist(ancest_info), lang_names$tip_id)]
  
  node_list <- vector("list", nrow(new))
  
  names(node_list) <- node_label
  
  for(i in 1: nrow(new)){
    
    node_list[[i]]<-as.vector(t(new[i,(1:n_length[i])]))
    
  }
  
  #order all nodes and create df
  ordered_nodes <- unique(tree$tree$edge[, 1])
  ordered_nodes <- data.frame(ordered_nodes)
  row.names(ordered_nodes) <- ordered_nodes[, 1]
  
  # apply descendants (package phangorn) to ordered nodes
  descendants <- sapply(rownames(ordered_nodes),
                        function(o) {as.vector(Descendants(tree$tree,
                                                           node = ordered_nodes[o, ],
                                                           type = "tips")[[1]])},
                        USE.NAMES = T, simplify = F)
  
  #check tip labels
  descendants_tip_labels <- sapply(descendants, function (d)
    tree$tree$tip.label[d], USE.NAMES = T, simplify = F)
  
  #match the tables
  match_table <- as.data.frame(sapply(node_list, 
                                      function(k) sapply(descendants_tip_labels,
                                                         function (d) identical(sort(d), sort(k)),
                                                         USE.NAMES = T), USE.NAMES = T))
  #simplify the table
  simple_match_table <- as.data.frame(apply(match_table, 2,
                                            function(u) paste(names(which(u)), collapse = "," )))
  
  #add names    
  colnames(simple_match_table) <- c("Node ID Nexus")
  simple_match_table$`Node ID` <- rownames(simple_match_table)
  
  #add matching tips
  tip_match<-lang_names[,2:3]
  
  #now check the correct match of tips in the nexus tree
  is_tip <- test_tree$tree$edge[,2] <= length(test_tree$tree$tip.label)
  ordered_tips <- data.frame(node_id_nexus=as.factor(test_tree$tree$edge[is_tip, 2]))
  ordered_tips$lang_name<-test_tree$tree$tip.label[ordered_tips$node_id_nexus]
  #match them
  tip_match<-merge(tip_match, ordered_tips, by = "lang_name", all.x = TRUE)
  #change order
  tip_match<-tip_match[,c(3,2)]
  colnames(tip_match)<-c("Node ID Nexus", "Node ID")
  
  #merge them
  simple_match_table<-rbind(simple_match_table, tip_match)
  simple_match_table<-simple_match_table[order(simple_match_table$`Node ID`),] 
  
  
  
  partitions<-merge(x = partitions, y = simple_match_table, by = "Node ID", all.x = TRUE)
  partitions$`Node ID Nexus`<- as.numeric(levels(partitions$`Node ID Nexus`))[partitions$`Node ID Nexus`]
  
  return(list(partitions = partitions, 
              matching_table = simple_match_table))
  
}