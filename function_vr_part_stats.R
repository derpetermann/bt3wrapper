# The function computes how often a partition is occuring at each NEXUS ID branch. 
# It therefore takes the result of the vr mcmc as well as the original tree used for  
# the vr mcmc as input. It further matches the IDs given by the Discrete_VR_Simple
# function with the NEXUS IDs found in the original tree.

# IN:
# mcmc_res        the result of the VR_fit_mcmc function
# tree       a phylogenetic tree of class phylo or multiPhylo that was also used for
#                 VR_fit_mcmc function

# OUT
# part_stats      summary of the most frequent found partitions (at which Branch ID)
# matching_table  a matching table of the NEXUS IDs and the IDs produced by the
#                 mcmc function

vr_part_stats<-function(mcmc_res, tree){
  require(phangorn)
  
  #number of repetitions
  num_rep<-length(mcmc_res$VR_test)
  
  #list of all found partitions
  all_partitions<-mcmc_res$VR_part_list
  #remove all cases where no partitions where found
  all_partitions<-all_partitions[!is.na(all_partitions)]
  #take the information from the node structure
  ancest_info<-mcmc_res$VR_test[[1]]$node_info$node_structure
  #store all the language names
  lang_names<-mcmc_res$VR_test[[1]]$node_info$tip_names
  
  #add all partitions to a dataframe
  part_list<-data.frame(do.call(rbind, Map(data.frame, all_partitions)))
  colnames(part_list)<-colnames(all_partitions$partitions_rep_1)
  
  #create vector of node labels
  node_label<-ancest_info[,1]
  
  n_length<-ancest_info[,2]
  
  #remove the now unnecessary information
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
  
  
  #add the information to the part_list
  part_list<-merge(x = part_list, y = simple_match_table, by = "Node ID", all.x = TRUE)
  part_list$`Node ID Nexus`<- as.numeric(levels(
    part_list$`Node ID Nexus`))[part_list$`Node ID Nexus`]
  
  #now we set the nodes to the branches that lead to the nodes
  #we first create dataframes for only the branches/nodes
  part_list_nodes<-subset(part_list, `Node / Branch` == "Node")
  part_list_branch<-subset(part_list, `Node / Branch` == "Branch")[,5]
  
  #create a lookup dataframe to check for parental nodes
  parent_df<-data.frame(tree$tree$edge)
  colnames(parent_df)<-c("parent", "Node ID Nexus")
  
  #add the parental nodes back to the partitions
  part_list_nodes<-merge(part_list_nodes, parent_df, by = "Node ID Nexus", all.x = TRUE)
  part_list_nodes<-part_list_nodes[,6]
  
  #combine the two
  part_list_table<-data.frame(table(c(part_list_branch, part_list_nodes)))
  
  #calculate the percental occurence of each partition
  part_list_table$Freq<-part_list_table$Freq/num_rep
  #change the column names
  colnames(part_list_table)<-c("Branch ID Nexus", "Perc")
  
  #order them so that the most occuring partition is at the top
  part_list_table<-part_list_table[order(part_list_table$Perc, decreasing = TRUE),]
  
  
  return(list(part_stats = part_list_table,
              matching_table = simple_match_table))
  
}