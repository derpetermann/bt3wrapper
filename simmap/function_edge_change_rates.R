# This funciton takes the output from the sim_map function and computes
# statistics about the changes at each individual edge, i.e. branch. 

# IN
# sim_map_res   the results of the sim_map function

# OUT
# edge_change   statistics about the changes at each edge 

edge_change_rates<-function(sim_map_res){
  #take first element to have a look at tree
  tree<-sim_map_res[[1]]
  
  #create dataframe to store all the information
  edge_change<-data.frame(tree$edge)
  #add edge length
  edge_change$length<-tree$edge.length
  
  colnames(edge_change)[1:2]<-c("from", "to")
  #add edge ids
  edge_change$edge_id<-rownames(edge_change)
  
  #create list so that we can store all segment information
  segments<-vector("list", nrow(sim_map_res[[1]]$edge))
  
  for(i in 1:length(segments)){
    segments[[i]]<-vector("list",length(sim_map_res))
    for(j in 1:length(sim_map_res)){
      segments[[i]][[j]]<-names(sim_map_res[[j]][[5]][[i]])
    }
    names(segments)<-rownames(edge_change)
  }
  
  for(k in 1:nrow(edge_change)){
    #for one
    edge_change$Y[k]<-length(segments[[k]][segments[[k]]=="Y"])/100
    edge_change$N[k]<-length(segments[[k]][segments[[k]]=="N"])/100
    #for 2
    cond1 <- sapply(segments[[k]], function(x) length(x) == 2)
    seg1a<-segments[[k]][cond1]
    seg1b<-lapply(seg1a, `[[`, 1)
    edge_change$Y_N[k]<-length(seg1b[seg1b=="Y"])/100
    edge_change$N_Y[k]<-length(seg1b[seg1b=="N"])/100
    
    #for 3
    cond2 <- sapply(segments[[k]], function(x) length(x) == 3)
    seg2a<-segments[[k]][cond2]
    seg2b<-lapply(seg2a, `[[`, 1)
    edge_change$Y_N_Y[k]<-length(seg2b[seg2b=="Y"])/100
    edge_change$N_Y_N[k]<-length(seg2b[seg2b=="N"])/100
  }
  #the rest
  edge_change$rest<-1-rowSums(edge_change[,5:10])
  edge_change<-edge_change[,c(4,1:3,5:11)]
  tmp<-edge_change[,5:11]
  edge_change$max<-colnames(tmp)[apply(tmp,1,which.max)]
  
  return(edge_change)
}