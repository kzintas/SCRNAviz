#function for checking alternate core edge

check_alternate<-function(sub_df, all_df)
{
  assign_df<- data.frame("Sink_node_res"=character(),
                         "Assign_from"=character(),
                         "Assign_to"=character(),
                         stringsAsFactors=FALSE)
  for(i in seq(1,nrow(sub_df))){
    
    for(j in seq(1,nrow(all_df))){
      if(sub_df$to_node[i]==all_df$to_node[j] && sub_df$from_node[i]!=all_df$from_node[j] && all_df$is_core[j]==TRUE){
        assign_df[nrow(assign_df)+1,]<- c(sub_df$to_node[i],sub_df$from_node[i],all_df$from_node[j])
        
      }
    }
  }
  assign_df
}

#split the DF into resolution and cluster name
manip_names<-function(Df)
{
  Names<-colnames(Df)
  print(Names)
  for(names in Names){
    Df<-Df %>%
      separate(names, c(paste0(names,"_res"),paste0( names,"_clust")), "C")
  }
  Df
}

#for each row in the core edge false dataset, find out corresponding source and sink in pbmc, and change accordingly
change_assignment<- function(Df, Seurat_obj){
  for(i in seq(1:nrow(Df))){
    sink_res=match(Df$`Sink_node_res_res`[i],colnames(Seurat_obj@meta.data))
    source_res=match(Df$`Assign_from_res`[i],colnames(Seurat_obj@meta.data))
    
    #lapply(Seurat_obj@meta.data, change_src_sink, df_row=Df[1,])
    for(j in 1:nrow(Seurat_obj@meta.data)){
      #change_src_sink(Seurat_obj@meta.data[j,],df_row=Df[1,])
      if(Seurat_obj@meta.data[j,sink_res]==Df$Sink_node_res_clust[i]  && Seurat_obj@meta.data[j,source_res] ==Df$Assign_from_clust[i]){
        #print(Seurat_obj@meta.data[j,source_res])
        Seurat_obj@meta.data[j,][source_res]<-Df$`Assign_to_clust`[i]
        
      }
      
    }
  }
  return (Seurat_obj)
  #TODO-> separate assignment into another function
}

#Prune the tree so only core edges remain
prune_tree<-function(graph_Df){
  repeat{
    graph_Df <-graph_Df[ !duplicated(names(graph_Df)) ]
    k=nrow(graph_Df[graph_Df$is_core==FALSE,])
    print(k)
    if(nrow(graph_Df[graph_Df$is_core==FALSE,])==0) break
    #select rows where CORE is False
    not_core_df<-  graph_Df%>%
      filter(is_core==FALSE)%>%  
      select(from_node,to_node,is_core) 
    
    #Apply function
    assign_df<-check_alternate(not_core_df,graph_Df)
    assign_df<- manip_names(assign_df)
    pbmc<- change_assignment(assign_df, pbmc)
    Graph<-clustree(pbmc , return="graph")
    graph_Df<-as_long_data_frame(Graph)
    
    
  }
  
  return (Graph)
}

#collapse TREE

collapse_tree<- function(Original_graph){
  node_dists<-distances(Original_graph, to=1)
  layered_dist<- unique(distances(Original_graph, to=1))
  delete_set<-vector('numeric')
  for(i in seq(length(layered_dist),2)){  
    #Get distances for each layer
    
    prev_layer<-which(node_dists==layered_dist[i-1])
    current_layer<-which(node_dists==layered_dist[i])
    
    if(length(prev_layer)==length(current_layer)){
      for(j in seq(1,length(current_layer))){
        #find parent
        nei_in<-neighbors(Original_graph,prev_layer[j],mode="in")
        #find_Original_Child
        nei_out<-neighbors(Original_graph,prev_layer[j],mode="out")
        for(vertices in nei_out){
         
          Original_graph<-Original_graph+edge(nei_in,vertices)#%>%
          from_edge<-get.edge.ids(Original_graph,c(as.integer(nei_in),prev_layer[j]))
          
          to_edge<-get.edge.ids(Original_graph, c(prev_layer[j],vertices))
          
          new_edge<-get.edge.ids(Original_graph, c(as.integer(nei_in),vertices))
          print(E(graph_modified)[[from_edge]])
          E(graph_modified)[[from_edge]]
        }
      }
      delete_set<-c(delete_set,prev_layer)

    }
  }
  Original_graph<-delete_vertices(Original_graph, delete_set)
  
  return(Original_graph)
}


#Plot tree to check if implementation is correct
plotting_collapsed_tree<-function(Graph_original){
  new_graph <- graph.empty() + vertices(V(Graph_original))
  Edges_original<-get.edgelist(Graph_original)
  
  #RE ORDER Edges
  Edges_original[order(Edges_original[,1]),]
  new_graph <- new_graph + edges(as.vector(t(Edges_original)))
  plot(new_graph, layout=layout_as_tree, vertex.size=4,
       vertex.label.dist=1,  edge.arrow.size=0.5)
}


#Get pruned tree with only true core edges
modified_graph<-prune_tree(graph_df)
#Data Frame of modified tree
modified_graph_df<-as_long_data_frame(modified_graph)


collapsed_graph <-collapse_tree(modified_graph)
plot(collapsed_graph, layout=layout_as_tree, vertex.size=4,
     vertex.label.dist=1,  edge.arrow.size=0.5)

collapsed_graph_df<-as_long_data_frame(collapsed_graph)



## REST of the code is not required to run the analysis
plotting_collapsed_tree(collapsed_graph)

plot(modified_graph, layout=layout_as_tree, vertex.size=4,
     vertex.label.dist=1,  edge.arrow.size=0.5)



###Tried some stuff

#c(paste(Names,"_res")
plot(graph, layout=layout_as_tree, vertex.size=4,
     vertex.label.dist=1,  edge.arrow.size=0.5)



for(i in seq(1:nrow(assign_df))){
  sink_res=match(assign_df$`Sink_node_res_res`[i],colnames(pbmc@meta.data))
  source_res=match(assign_df$`Assign_from_res`[i],colnames(pbmc@meta.data))
  for(j in 1:nrow(pbmc@meta.data)){
    if(pbmc@meta.data[j,sink_res]==assign_df$Sink_node_res_clust[i]  && pbmc@meta.data[j,source_res] ==assign_df$Assign_from_clust[i]){
      print(pbmc@meta.data[j,source_res])
      pbmc@meta.data[j,source_res]=assign_df$`Assign_to_clust`[i]
      print(pbmc@meta.data[j,source_res])
    }
    
  }
}


#See information of data frame
str(graph_df)
colnames(graph_df)

not_core_df <- subset(graph_df[c("from_clust","to_clust","from_RNA_snn_res.","to_RNA_snn_res.", "is_core")], is_core == FALSE)
not_core_df[order(not_core_df$to_node, decreasing=TRUE), ]




graph_modified<-clustree(pbmc ,prop_filter=0.2, return="graph")

# Drop duplicated columns
graph_df <-graph_df[ !duplicated(names(graph_df)) ]
colnames(graph_df)

#select rows where CORE is False
not_core_df<-  graph_df%>%
              filter(is_core==FALSE)%>%  
              select(from_node,to_node,is_core) 
  
modified_graph_df<-as_long_data_frame(graph_modified)

View(modified_graph_df)


##TRies Tree collapsing using data frame


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

#newDF = df[ core_df$from_clust == graph_df$from_clust , ]

assign_df[1,]$Sink_node_res_res
pbmc@meta.data[2,][9]<- 3
x<-pbmc@meta.data[2,]
x[9]

pbmc@meta.data[[2]]

###NOT WORKING
#for each row in the core edge false dataset, find out corresponding source and sink in pbmc, and change accordingly
change_src_sink<- function(seurat_row, df_row){
  print(seurat_row)
  print(df_row)
  sink_res=match(df_row$`Sink_node_res_res`,colnames(seurat_row))
  source_res=match(df_row$`Assign_from_res`,colnames(seurat_row))
  
  print(sink_res)
  print(source_res)
  if(seurat_row[sink_res]==df_row$Sink_node_res_clust  && seurat_row[source_res] ==df_row$Assign_from_clust[i]){
    seurat_row[source_res]=df_row$`Assign_to_clust`
  }
  seurat_row[source_res]
}