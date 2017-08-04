# The function visualizes the findings of the vr_part_stats function. 
# Every branch is colored according the frequency it was found in the 
# mcmc function.

# IN:
# tree              a phylogenetic tree of class phylo or multiPhylo
# part_stat_result  the results of the vr_part_stats function
# color             the color that should be used for the visualization
#                   (either blue, red or green, standard is blue)
# OUT
# the visualization

vr_part_stat_visu<-function(tree, 
                            part_stat_result, 
                            color = c("blue", "green", "red")){
  require(ggtree)
  require(phangorn)
  require(RColorBrewer)
  require(ggplot2)
  
  #create a df to store the nodes and their respective colors
  color_df<-data.frame(node=1:(Nnode(tree$tree)+length(tree$tree$tip.label)))
  
  color_df<-merge(x = color_df, y = part_stat_result, by.x = "node", by.y = "Branch ID Nexus",
                  all.x = TRUE)
  #set NA to 0
  color_df[is.na(color_df)] <- 0
  
  if(missing(color) || toupper(color) == "BLUE"){
    pal <- "Blues"
  }
  else if(toupper(color) == "RED"){
    pal <- "Reds"
  }
  else if(toupper(color) == "GREEN"){
    pal <- "Greens"
  }
  
  m<-ggtree(tree$tree,aes(color=Perc), size = 1.5) %<+% color_df +
    scale_color_distiller(palette = pal, guide = "colourbar", direction = 1,
                          name = "Change Prob. [%]")+
    geom_tiplab(size=4, color="black")+
    theme(legend.position="left",
          legend.title = element_text(size=10, face="bold")) + 
    xlim(NA, 6200)+
    geom_text2(aes(subset=!isTip, label=node), hjust=-.3)
  
  return(m)
}