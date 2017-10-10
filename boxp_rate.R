boxp_rate<-function(results){

  res_df<-data.frame(type = character(),
                     qNY = numeric(),
                     qYN = numeric(),
                     rate = numeric(),
                     stringsAsFactors = FALSE)
  
  row.length<-as.numeric(nrow(results[[1]]$Rate_dist))
  t_names<-names(results)
  
  for(i in 1:length(results)){
    if(length(results[[i]]$Rate_dist) == 1){

      df_add<-data.frame(type = t_names[i],
                         qNY = NA,
                         qYN = NA,
                         rate = results[[i]]$Rate_dist[[1]])
      
      res_df <- rbind(res_df, df_add, stringsAsFactors=FALSE)
      
      
    }
    else if(length(results[[i]]$Rate_dist) == 2){
      
      df_add<-data.frame(type = t_names[i],
                         qNY = results[[i]]$Rate_dist[[1]],
                         qYN = results[[i]]$Rate_dist[[2]],
                         rate = NA)
      
      res_df <- rbind(res_df, df_add, stringsAsFactors = FALSE)
    }
  }
  
  
  df_melt<-melt(res_df, id.vars = "type")
  
  df_melt<-subset(df_melt, !is.na(value))
  
  vis <- ggplot(df_melt, aes(x=type, y=value, fill = variable)) +
                geom_boxplot()+
                labs(x=NULL, y = "rate estimate")+
                scale_fill_manual(name = NULL, values=c("blue", "red", "darkgrey"))+
                theme(axis.text=element_text(size=10),
                      axis.title=element_text(size=12),
                      legend.text = element_text(size = 10))
        
  return(vis)
}
