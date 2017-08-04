# The function visualizes the findings of the Discrete_BT3_VR function. 
# All found partitions are colored differently.

# IN:
# tree              a phylogenetic tree of class phylo or multiPhylo
# vr_partitions     the results of the node_id_match function

# OUT
# the visualization
tree_visu<-function(tree, 
                    vr_partitions){
  require(ggtree)
  require(phangorn)
  require(RColorBrewer)
  require(ggplot2)
  
  #only need the two last rows that contain the branch/node information(4) and the nexus id
  match2<-vr_partitions[[1]][,4:5]
  
  #add number of descendant
  match2$no_desc<-lengths(Descendants(tree$tree, match2$`Node ID Nexus`, type = "all"))
  
  #sort by number of descendants
  match2<-match2[order(-match2$no_desc),]
  
  #create a colorlist
  if(nrow(vr_partitions[[1]])<= 9){
    color_list<-brewer.pal(9, "Set1")
  } else if (nrow(vr_partitions[[1]])> 9 & nrow(vr_partitions)<= 12){
    color_list<-brewer.pal(12, "Set3")
  } else {
    color_list<-c(brewer.pal(8, "Dark2"),brewer.pal(8,"Accent"))
  }
  
  #create a df to store the nodes and their respective colors
  color_df<-data.frame(node=1:(Nnode(tree$tree)+length(tree$tree$tip.label)),
                       color = "black",stringsAsFactors = FALSE)
  
  #assign all vr_partitions a color
  for(i in 1:nrow(match2)){
    desc<-unlist(Descendants(tree$tree, match2$`Node ID Nexus`[i], type = "all"))
    color_df[color_df$node %in% desc,]$color<-rep(color_list[i], length(desc))
  }
  
  #now, all vr_partitions are colored up to the node where they change,
  #however, we also need to color certain branches (when the change is at a branch,
  #and not at a node)
  
  #first, we identify all vr_partitions at branches and add the nodes to a new dataframe
  branch_ancest<-data.frame(node = match2[match2$`Node / Branch` == "Branch",]$`Node ID Nexus`)
  #when there is a change at a branch, we need to color the branch that leads to the node
  #stored in the new dataframe
  #we therefore first identify the color at the node
  branch_ancest<-merge(x = branch_ancest, y = color_df, all.x = TRUE)
  #then we identify the ancestor
  branch_ancest$ancest<-Ancestors(tree$tree,branch_ancest$node, type = "parent")
  #to add it back correctly, we have to order it
  branch_ancest<-branch_ancest[order(branch_ancest$ancest),]
  
  #we add the color information back
  color_df[color_df$node %in% branch_ancest$ancest,]$color<-branch_ancest$color
  
  visu<-ggtree(tree$tree) %<+% color_df + aes(color=I(color))+ 
    geom_text2(aes(subset=!isTip, label=node, color = "black"), hjust=-.3)+ 
    geom_tiplab(color = "black")+
    theme(legend.position="right")+xlim(0,6000)
  
  return(visu)
  
}