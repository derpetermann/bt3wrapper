boxp_root<-function(results){

  res_df<-data.frame(type = character(),
                     Prob_Y = numeric(),
                     Prob_N = numeric(),
                     stringsAsFactors = FALSE)
  
  row.length <- as.numeric(nrow(results[[1]]$Root_dist))
  t_names <- names(results)
  
  for(i in 1:length(results)){

      df_add<-data.frame(type = t_names[i],
                       Prob_Y = results[[i]]$Root_dist[[1]],
                       Prob_N = results[[i]]$Root_dist[[2]])
      
      res_df<-rbind(res_df, df_add, stringsAsFactors=FALSE)
      
      
  }
  
  df_melt <- melt(res_df, id.vars = "type")
  
  df_melt <- subset(df_melt, !is.na(value))
  
  vis <- ggplot(df_melt, aes(x=type, y=value, fill = variable)) +
                geom_boxplot()+
                labs(x=NULL, y = "Probability at the root")+
                scale_fill_manual(name = NULL, values=c("blue", "red", "darkgrey"),
                                  labels = c("present (Y)", "absent (N)"))+
                theme(axis.text=element_text(size=10),
                      axis.title=element_text(size=14),
                      legend.text = element_text(size = 10))
  return(vis)}
