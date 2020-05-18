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
    separate(names, c(paste0(names,"_res"),paste0( names,"_clust")), "C")
}

#c(paste(Names,"_res")


#for each row in the core edge false dataset, find out corresponding source and sink in pbmc, and change accordingly

for(i in seq(1:nrow(assign_df))){
  sink_res=match(assign_df$`Sink_node_res`[i],colnames(pbmc@meta.data))
  source_res=match(assign_df$`Assign_from_res`[i],colnames(pbmc@meta.data))
  for(j in 1:nrow(pbmc@meta.data)){
    if(pbmc@meta.data[j,sink_res]==assign_df$Sink_node_clust[i]  && pbmc@meta.data[j,source_res] ==assign_df$Assign_from_clust[i]){
      print(pbmc@meta.data[j,source_res])
      pbmc@meta.data[j,source_res]=assign_df$`Assign_to_clust`[i]
      print(pbmc@meta.data[j,source_res])
    }
    
  }
}



graph_modified<-clustree(pbmc ,prop_filter=0.2, return="graph")

plot(graph_modified, layout=layout_as_tree, vertex.size=4,
     vertex.label.dist=1,  edge.arrow.size=0.5)

core_df <- subset(graph_df[c("from_node","to_node", "is_core")], is_core == FALSE)
core_df <-core_df[order(core_df$to_node), ]
modified_graph_df<-as_long_data_frame(graph_modified)

View(modified_graph_df)


n(neighbors(graph_modified, 1, mode="out"))

for (nodes in V(graph_modified)){
  a<-neighbors(graph_modified, nodes, mode="out")
  b<-neighbors(graph_modified, nodes, mode="in")
  if(length(a)==length(b))
    print(nodes)
}
V(graph_modified)
vertex.attributes(graph_modified, 4)
edge.attributes(graph_modified, 4)
degree(graph_modified, mode="out")
degree(graph_modified, mode="in")

x<-ego(graph_modified, nodes=1, mode="out", mindist=1)

print(x)
length(x)

dfs(graph_modified,1)

distance_table(graph_modified)
max(distances(graph_modified, to=1))
x<-distances(graph_modified, to=1)
a<-table(x)
type(x)
df<-as.data.frame(x)
df$observation <- 1:nrow(df) 
graph2<-graph_modified
p<- unique(x)
for(val in unique(x)){
  #m<- df %>% filter(V1==val)#select(df,df$V1==val)
  m<-which(x==val)
  print(m)
  
  
}


for(i in seq(2,length(p))){
  #m<- df %>% filter(V1==val)#select(df,df$V1==val)
  m<-which(x==p[i-1])
  print(p[i])
  print(m)
  
}
x[14]

unique(x)



for(i in seq(2,length(p))){
  #m<- df %>% filter(V1==val)#select(df,df$V1==val)
  m<-which(x==p[i-1])
  print(p[i])
  print(m)
  
}