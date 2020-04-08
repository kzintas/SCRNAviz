#function for checking alternate core edge

check_alternate<-function(core_df, all_df)
{
  assign_df<- data.frame("Sink_node"=character(),
                         "Assign_from"=character(),
                         "Assign_to"=character(),
                         stringsAsFactors=FALSE)
  for(i in seq(1,nrow(core_df))){
    
    for(j in seq(1,nrow(all_df))){
      if(core_df$to_node[i]==all_df$to_node[j] && core_df$from_node[i]!=all_df$from_node[j] && all_df$is_core[j]==TRUE){
        assign_df[nrow(assign_df)+1,]<- c(core_df$to_node[i],core_df$from_node[i],all_df$from_node[j])
        
      }
    }
  }
  return (assign_df)
}

#Apply function
assign_df<-check_alternate(core_df,graph_df)


#split the DF into resolution and cluster name

Names<-colnames(assign_df)
for(names in Names){
  assign_df<-assign_df %>%
    separate(names, c(paste(names,"_res"),paste( names,"_clust")), "C")
}

#c(paste(Names,"_res")


#for each row in the core edge false dataset, find out corresponding source and sink in pbmc, and change accordingly
for(i in seq(1:nrow(assign_df))){
  sink_res=match(assign_df$`Sink_node _res`[i],colnames(pbmc@meta.data))
  source_res=match(assign_df$`Assign_from _res`[i],colnames(pbmc@meta.data))
  for(j in 1:nrow(pbmc@meta.data)){
    if(pbmc@meta.data[j,sink_res] ==assign_df$`Sink_node _clust`[i]  && pbmc@meta.data[j,source_res] ==assign_df$`Assign_from _clust`[i]){
      print(pbmc@meta.data[j,source_res])
      pbmc@meta.data[j,source_res]=assign_df$`Assign_to _clust`[i]
      print(pbmc@meta.data[j,source_res])
    }
    
  }
}