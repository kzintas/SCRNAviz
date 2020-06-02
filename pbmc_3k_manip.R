#function for checking alternate core edge

check_alternate <- function(sub_df, all_df)
{
  assign_df <- data.frame(
    "Sink_node_res" = character(),
    "Assign_from" = character(),
    "Assign_to" = character(),
    stringsAsFactors = FALSE
  )
  for (i in seq(1, nrow(sub_df))) {
    for (j in seq(1, nrow(all_df))) {
      if (sub_df$to_node[i] == all_df$to_node[j] &&
          sub_df$from_node[i] != all_df$from_node[j] &&
          all_df$is_core[j] == TRUE) {
        assign_df[nrow(assign_df) + 1, ] <-
          c(sub_df$to_node[i],
            sub_df$from_node[i],
            all_df$from_node[j])
        
      }
    }
  }
  assign_df
}

#split the DF into resolution and cluster name
manip_names <- function(Df)
{
  Names <- colnames(Df)
  #print(Names)
  for (names in Names) {
    Df <- Df %>%
      separate(names, c(paste0(names, "_res"), paste0(names, "_clust")), "C")
  }
  Df
}

#for each row in the core edge false dataset, find out corresponding source and sink in Seurat, and change accordingly
change_assignment <- function(Df, Seurat_obj) {
  for (i in seq(1:nrow(Df))) {
    sink_res = match(Df$`Sink_node_res_res`[i], colnames(Seurat_obj@meta.data))
    source_res = match(Df$`Assign_from_res`[i], colnames(Seurat_obj@meta.data))
    
    #lapply(Seurat_obj@meta.data, change_src_sink, df_row=Df[1,])
    for (j in 1:nrow(Seurat_obj@meta.data)) {
      #change_src_sink(Seurat_obj@meta.data[j,],df_row=Df[1,])
      if (Seurat_obj@meta.data[j, sink_res] == Df$Sink_node_res_clust[i]  &&
          Seurat_obj@meta.data[j, source_res] == Df$Assign_from_clust[i]) {
        #print(Seurat_obj@meta.data[j,source_res])
        Seurat_obj@meta.data[j, ][source_res] <-
          Df$`Assign_to_clust`[i]
        
      }
      
    }
  }
  return (Seurat_obj)
  #TODO-> separate assignment into another function
}

#Prune the tree so only core edges remain
prune_tree <- function(graph_Df, Seurat_obj) {
  repeat {
    graph_Df <- graph_Df[!duplicated(names(graph_Df))]
    k = nrow(graph_Df[graph_Df$is_core == FALSE, ])
    #print(k)
    if (nrow(graph_Df[graph_Df$is_core == FALSE, ]) == 0)
      break
    #select rows where CORE is False
    not_core_df <-  graph_Df %>%
      filter(is_core == FALSE) %>%
      select(from_node, to_node, is_core)
    
    #Apply function
    assign_df <- check_alternate(not_core_df, graph_Df)
    assign_df <- manip_names(assign_df)
    Seurat_obj <- change_assignment(assign_df, Seurat_obj)
    Graph <- clustree(Seurat_obj , return = "graph")
    graph_Df <- as_long_data_frame(Graph)
    
    
  }
  
  return (Graph)
}


#collapse TREE


collapse_tree <- function(Original_graph) {
  #Distances of nodes and layers
  node_dists <- distances(Original_graph, to = 1)
  layered_dist <- unique(distances(Original_graph, to = 1))
  
  #Which vertices and edges we want to delete
  delete_set_edges <- vector('numeric')
  delete_set_vertices <- vector('numeric')
  
  #Vertex and edge lists of the graph from which we will construct the collapsed graph
  ver_list <- as_data_frame(Original_graph, what = "vertices")
  edge_list <- as_data_frame(Original_graph, what = "edges")
  
  i <- length(layered_dist)
  while (i >= 2) {
    prev_layer <- which(node_dists == layered_dist[i - 1])
    current_layer <- which(node_dists == layered_dist[i])
    
    if (length(prev_layer) == length(current_layer)) {
      while (length(prev_layer) == length(current_layer)) {
        delete_set_vertices <- c(delete_set_vertices, prev_layer)
        i <- i - 1
        
        prev_layer <- which(node_dists == layered_dist[i - 1])
        
      }
      
      for (j in seq(1, length(prev_layer))) {
        #find parent
        nei_in <- prev_layer[j]
        
        #find_Original_Child
        temp_out <-
          neighbors(Original_graph, prev_layer[j], mode = "out")
        
        nei_out <- c()
        
        for (nodes in temp_out) {
          temp_node <- neighbors(Original_graph, nodes, mode = "out")
          repeat {
            if (temp_node %in% current_layer) {
              nei_out <- c(nei_out, temp_node)
              break
            }
            
            temp_node <-
              neighbors(Original_graph, temp_node, mode = "out")
          }
        }
        
        for (k in seq (1, length(nei_out))) {
          from_edge <-
            get.edge.ids(Original_graph, c(prev_layer[j], as.integer(temp_out[k])))
          
          to_edge <-
            incident(Original_graph, as.integer(nei_out[k]), mode = "in")
          
          get_src_attrs <-
            as.list(sapply(edge_attr(Original_graph), '[[', from_edge))
          
          get_sink_attrs <-
            as.list(sapply(edge_attr(Original_graph), '[[', to_edge))
          
          get_src_attrs["to_clust"] <- get_sink_attrs["to_clust"]
          
          get_src_attrs["to_RNA_snn_res."] <-
            get_sink_attrs["to_RNA_snn_res."]
          
          get_src_attrs <-
            append(get_src_attrs,
                   list(from = prev_layer[j], to = as.integer(nei_out[k])),
                   0)
          
          edge_list <- rbind(edge_list, get_src_attrs)
          
          
        }
      }
      
    }
    i <- i - 1
  }
  
  for (nodes in delete_set_vertices) {
    adj_edge <- incident(Original_graph, nodes, mode = "all")
    delete_set_edges <- c(delete_set_edges, adj_edge)
    
  }
  
  ver_list <- ver_list[-delete_set_vertices, ]
  id <- rownames(ver_list)
  ver_list <- cbind(name = id, ver_list)
  
  edge_list <- edge_list[-delete_set_edges,]
  
  dummy_graph <-
    graph_from_data_frame(edge_list, directed = TRUE, vertices = ver_list)
  dummy_graph <- as_tbl_graph(dummy_graph)
  return(dummy_graph)
}




#Plot tree to check if implementation is correct
plotting_collapsed_tree <- function(Graph_original) {
  new_graph <- graph.empty() + vertices(V(Graph_original))
  Edges_original <- get.edgelist(Graph_original)
  
  #RE ORDER Edges
  Edges_original[order(Edges_original[, 1]), ]
  new_graph <- new_graph + edges(as.vector(t(Edges_original)))
  plot(
    new_graph,
    layout = layout_as_tree,
    vertex.size = 4,
    vertex.label.dist = 1,
    edge.arrow.size = 0.5
  )
}

#take a clustree as graph and Seurat object as input

reassign_and_collapse <- function(graph, Seurat_obj) {
  graph_df <- as_long_data_frame(graph)
  
  #Get pruned tree with only true core edges
  modified_graph <- prune_tree(graph_df, Seurat_obj)
  #Data Frame of modified tree
  modified_graph_df <- as_long_data_frame(modified_graph)
  
  
  collapsed_graph <- collapse_tree(modified_graph)
  plot(
    collapsed_graph,
    layout = layout_as_tree,
    vertex.size = 4,
    vertex.label.dist = 1,
    edge.arrow.size = 0.5
  )
  
  collapsed_graph_df <- as_long_data_frame(collapsed_graph)
}