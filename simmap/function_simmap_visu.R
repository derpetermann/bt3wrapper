# The function visualizes the results of the simmap and the edge_change_rates function.

# IN
# sim_map_res   the results of the simmap_plot function
# edge_changes  the results of the edge_change_rates function
# type          either gain, loss, presence, absence, all or highest, 
#               standard is "all"

# OUT
# The visualizations
simmap_visu<-function(sim_map_res, 
                      edge_changes, 
                      type = c("gain", "loss", "presence", "absence", "all", "highest")){
  require(ggtree)
  require(RColorBrewer)
  #take first tree as output tree
  tree<-sim_map_res[[1]]
  
  #create new df
  color_df<-data.frame(node=1:(Nnode(tree)+length(tree$tip.label)))
  
  color_df<-merge(x = color_df, y = edge_changes[,c(3,5:8,12)], by.x = "node", by.y = "to",
                  all.x = TRUE)
  
  #remove na's
  color_df<-subset(color_df,!is.na(max))
  #change factor order
  color_df$max<-factor(color_df$max, levels = c("Y", "Y_N", "N_Y", "N"), exclude = NA)
  
  if(missing(type) || toupper(type) == "ALL"){
    v1<-ggtree(tree,aes(color=N_Y), size = 1.5) %<+% color_df +
      scale_color_distiller(palette = "Blues", guide = "colourbar", direction = 1,
                            name = "Gain Prob. [%]")+
      geom_tiplab(size=4, color="black")+
      theme(legend.position="left",
            legend.title = element_text(size=10, face="bold")) + 
      xlim(NA, 8900)
    
    v2<-ggtree(tree,aes(color=Y_N), size = 1.5) %<+% color_df +
      scale_color_distiller(palette = "Reds", guide = "colourbar", direction = 1,
                            name = "Loss Prob. [%]")+
      geom_tiplab(size=4, color="black")+
      theme(legend.position="left",
            legend.title = element_text(size=10, face="bold")) + 
      xlim(NA, 8900)  
    v3<-ggtree(tree,aes(color=N), size = 1.5) %<+% color_df +
      scale_color_distiller(palette = "Greys", guide = "colourbar", direction = 1,
                            name = "Absence Prob. [%]")+
      geom_tiplab(size=4, color="black")+
      theme(legend.position="left",
            legend.title = element_text(size=10, face="bold")) + 
      xlim(NA, 9100)
    v4<-ggtree(tree,aes(color=Y), size = 1.5) %<+% color_df +
      scale_color_distiller(palette = "Greys", guide = "colourbar", direction = 1,
                            name = "Presence Prob. [%]")+
      geom_tiplab(size=4, color="black")+
      theme(legend.position="left",
            legend.title = element_text(size=10, face="bold")) + 
      xlim(NA, 9100)  
    v5<-ggtree(tree,aes(color=max), size = 1.5) %<+% color_df +
      scale_colour_brewer(palette = "RdYlBu", name = "Type", breaks = c("Y", "Y_N", "N_Y", "N"),
                          labels = c("Presence", "Loss", "Gain", "Absence"),
                          direction = -1)+ 
      geom_tiplab(size=4, color="black")+
      theme(legend.position="left",
            legend.title = element_text(size=14, face="bold")) + 
      xlim(NA, 8800) 
    return(list(v1, v2, v3, v4, v5))
  }
  
  else if(toupper(type) == "GAIN"){
    v<-ggtree(tree,aes(color=N_Y), size = 1.5) %<+% color_df +
      scale_color_distiller(palette = "Blues", guide = "colourbar", direction = 1,
                            name = "Gain Prob. [%]")+
      geom_tiplab(size=4, color="black")+
      theme(legend.position="left",
            legend.title = element_text(size=10, face="bold")) + 
      xlim(NA, 8900)
    return(v)
  }
  
  else if(toupper(type) == "LOSS"){
    v<-ggtree(tree,aes(color=Y_N), size = 1.5) %<+% color_df +
      scale_color_distiller(palette = "Reds", guide = "colourbar", direction = 1,
                            name = "Loss Prob. [%]")+
      geom_tiplab(size=4, color="black")+
      theme(legend.position="left",
            legend.title = element_text(size=10, face="bold")) + 
      xlim(NA, 8900)
    return(v)
  }
  
  else if(toupper(type) == "ABSENCE"){
    v<-ggtree(tree,aes(color=N), size = 1.5) %<+% color_df +
      scale_color_distiller(palette = "Greys", guide = "colourbar", direction = 1,
                            name = "Absence Prob. [%]")+
      geom_tiplab(size=4, color="black")+
      theme(legend.position="left",
            legend.title = element_text(size=10, face="bold")) + 
      xlim(NA, 9100)
    return(v)
  }
  
  else if(toupper(type) == "PRESENCE"){
    v<-ggtree(tree,aes(color=Y), size = 1.5) %<+% color_df +
      scale_color_distiller(palette = "Greys", guide = "colourbar", direction = 1,
                            name = "Presence Prob. [%]")+
      geom_tiplab(size=4, color="black")+
      theme(legend.position="left",
            legend.title = element_text(size=10, face="bold")) + 
      xlim(NA, 9100)
    return(v)
  }
  else if(toupper(type) == "HIGHEST"){
    v<-ggtree(tree,aes(color=max), size = 1.5) %<+% color_df +
      scale_colour_brewer(palette = "RdYlBu", name = "Type", breaks = c("Y", "Y_N", "N_Y", "N"),
                          labels = c("Presence", "Loss", "Gain", "Absence"),
                          direction = -1)+ 
      geom_tiplab(size=4, color="black")+
      theme(legend.position="left",
            legend.title = element_text(size=14, face="bold")) + 
      xlim(NA, 8800)
    return(v)
  }
  
}
