library(clustree)
library(SingleCellExperiment)
library(TreeSummarizedExperiment)
library(dplyr)
library(tidygraph)
library(ggraph)
library(igraph)
library(Seurat)
library(tidyr)
library(stringr)
library(metagenomeSeq)
library(msd16s)
library(S4Vectors)
library(metavizr)
# install scater https://bioconductor.org/packages/release/bioc/html/scater.html
library(scater)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
library(loomR)
library(patchwork)
library(scRNAseq)
library(SC3)
library(scran)
library(scater)


#function for checking alternate core edge

check_alternate <- function(sub_df, all_df)
{
  assign_df <- data.frame(
    "Sink_node_res" = character(),
    "Sink_node_clust" = character(),
    "Assign_from_res" = character(),
    "Assign_from_clust" = character(),
    "Assign_to_res" = character(),
    "Assign_to_clust" = character(),
    stringsAsFactors = FALSE
  )
  for (i in seq(1, nrow(sub_df))) {
    for (j in seq(1, nrow(all_df))) {
      if (sub_df$to_node[i] == all_df$to_node[j]) {
        assign_df[nrow(assign_df) + 1, ] <-
          c(
            sub_df$to_cluster[i],
            sub_df$to_clust[i],
            sub_df$from_cluster[i],
            sub_df$from_clust[i],
            all_df$from_cluster[j],
            all_df$from_clust[j]
          )
        
      }
    }
  }
  #print(assign_df)
  assign_df
}


#for each row in the core edge false dataset, find out corresponding source and sink in Seurat, and change accordingly
change_assignment <- function(Df, Cluster_obj) {
  for (i in seq(1:nrow(Df))) {
    #print(Df$`Sink_node_res`[i])
    #print(colnames(Cluster_obj))
    match_str <- paste0('cluster', Df$`Sink_node_res`[i])
    sink_res = match(paste0('cluster', Df$`Sink_node_res`[i]),
                     colnames(Cluster_obj))
    source_res = match(paste0('cluster', Df$`Assign_from_res`[i]),
                       colnames(Cluster_obj))
    #print(match_str)
    #print(source_res)
    #print(Df$`Sink_node_res`[i])
    #print(colnames(Cluster_obj))
    
    
    #lapply(Seurat_obj@meta.data, change_src_sink, df_row=Df[1,])
    for (j in 1:nrow(Cluster_obj)) {
      #change_src_sink(Seurat_obj@meta.data[j,],df_row=Df[1,])
      #print(Cluster_obj[[j, sink_res]])
      #print(Df$Sink_node_clust[i])
      if (as.numeric(as.character(Cluster_obj[[j, sink_res]])) == as.numeric(Df$Sink_node_clust[i])  &&
          as.numeric(as.character(Cluster_obj[[j, source_res]])) == as.numeric(Df$Assign_from_clust[i])) {
        #print(Seurat_obj@meta.data[j,source_res])
        Cluster_obj[[j, source_res]] <-
          as.factor(Df$`Assign_to_clust`[i])
      }
      
    }
  }
  return (Cluster_obj)
  #TODO-> separate assignment into another function
}

check_cycle<- function(pruned_graph){
  
  delete_set_vertices <- vector('numeric')
  ver_list <- V(pruned_graph)
  #print(ver_list)
  for (nodes in ver_list) {
    adj_edge <- incident(pruned_graph, nodes, mode = "in")
    #print(adj_edge)
    #delete_set_edges <- c(delete_set_edges, adj_edge)
      if (length(adj_edge)>1){
          #print(as_ids(adj_edge[[1]]))
          #str(adj_edge)
          #print(nodes)
          adj_ver<- adjacent_vertices(pruned_graph,nodes, mode="in")
          adj_ver<- as_ids(adj_ver[[1]])
          str(adj_ver)
          remove_node<- sample(adj_ver,1)
          print(remove_node)
          
          remove_edge <- incident(pruned_graph, remove_node, mode = "all")
          print(remove_edge)
          pruned_graph<- delete_edges(pruned_graph,remove_edge)
          delete_set_vertices <- c(delete_set_vertices, remove_node)
          #pruned_graph<-delete.vertices(pruned_graph,remove_node)
          #print(as_tbl_graph(pruned_graph))
          #  delete_set_edges <- c(delete_set_edges, adj_edge)
          
          #print(adj_ver[[1]][1])
      }
  }
  pruned_graph<-delete.vertices(pruned_graph,delete_set_vertices)
  #as_tbl_graph(pruned_graph)
  
}

#Prune the tree so only core edges remain
prune_tree <- function(graph_Df, cluster_df) {
  repeat {
    graph_Df <- graph_Df[!duplicated(names(graph_Df))]
    #print(graph_Df)
    print(nrow(graph_Df[graph_Df$is_core == FALSE,]))
    #No False edges acyclic tree
    if (nrow(graph_Df[graph_Df$is_core == FALSE,]) == 0)
      break
    
    #Apply function
    assign_df <-
      check_alternate(graph_Df[graph_Df$is_core == FALSE,], graph_Df[graph_Df$is_core == TRUE,])
    
    cluster_df <- change_assignment(assign_df, cluster_df)
    #print(cluster_df)
    Graph <-
      clustree(
        cluster_df ,
        prefix = "cluster",
        prop_filter = 0,
        return = "graph"
      )
    
    graph_Df <- as_long_data_frame(Graph)
    
    
    
  }
  Graph<-check_cycle(Graph)
  
  objList <-
    list("Cluster_obj" = cluster_df, "Clustree_obj" = Graph)
  #return (Graph)
}



#collapse TREE


collapse_tree <- function(Original_graph) {
  #Distances of nodes and layers
  node_dists <- distances(Original_graph, to = 1)
  layered_dist <- unique(distances(Original_graph, to = 1))
  
  #Which vertices and edges we want to delete
  #delete_set_edges <- vector('numeric')
  delete_set_vertices <- vector('numeric')
  
  #Vertex and edge lists of the graph from which we will construct the collapsed graph
  ver_list <- as_data_frame(Original_graph, what = "vertices")
  #edge_list <- as_data_frame(Original_graph, what = "edges")
  
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
      
      
    }
    i <- i - 1
  }
  
  #for (nodes in delete_set_vertices) {
  #  adj_edge <- incident(Original_graph, nodes, mode = "all")
  #  delete_set_edges <- c(delete_set_edges, adj_edge)
  
  #}
  
  ver_list <- ver_list[-delete_set_vertices,]
  #id <- rownames(ver_list)
  #ver_list <- cbind(name = id, ver_list)
  
  #edge_list <- edge_list[-delete_set_edges, ]
  
  #dummy_graph <-
  #  graph_from_data_frame(edge_list, directed = TRUE, vertices = ver_list)
  #dummy_graph <- as_tbl_graph(dummy_graph)
  return(ver_list)
}



checkIfNotTree <- function(cluster_df) {
  cols <- colnames(cluster_df)
  if (length(unique(cluster_df[[1]])) > 1) {
    cluster_df$root <- "AllClusters"
    cols <- c("root", cols)
  }
  print(unique(cluster_df[[1]]))
  cluster_df[cols]
}


rename_clusters <- function(clusdata) {
  clusnames <- colnames(clusdata)
  clusnames <-
    as.numeric(unlist(regmatches(
      clusnames,
      gregexpr(
        "[[:digit:]]+\\.*[[:digit:]]*",
        clusnames
      )
    )))
  clusnames<- seq(length(colnames(clusdata)))
  #as.numeric(gsub("[^\\d]+\\.*[^\\d]", "", clusnames, perl = TRUE))
  clusnames <- paste0("cluster", clusnames)
  print(clusnames)
  #colnames(clusdata) <- unique(clusnames[(order(clusnames))])
  colnames(clusdata) <- clusnames
  clusdata
}

reassign_and_collapse <-
  function(cluster_df, count_matrix) {
    cluster_df <- rename_clusters(cluster_df)
    clustree_graph <-
      clustree(
        cluster_df ,
        prefix = "cluster",
        prop_filter = 0,
        return = "graph"
      )
    graph_df <- as_long_data_frame(clustree_graph)
    
    #Get pruned tree with only true core edges
    
    #cluster_df <- rename_clusters(cluster_df)
    
    modified_obj <- prune_tree(graph_df, cluster_df)
    #Data Frame of modified tree
    #modified_graph_df <- as_long_data_frame(modified_graph)
    
    
    #modified graph and seurat object
    modified_graph = modified_obj$Clustree_obj
    clusters <-  modified_obj$Cluster_obj
    
    
    collapsed_graph <- collapse_tree(modified_graph)
    
    #collapsed_graph_df <- as_long_data_frame(collapsed_graph)
    cluster_names <-
      unique(sapply(strsplit(collapsed_graph$node, "C"), '[', 1))
    
    #clusters <- modified_Seurat@meta.data
    clusters <- clusters[, cluster_names]
    
    #[[:digit:]]+\\.*[[:digit:]]
    #clusnames <-
    #  as.numeric(gsub("[^\\d]+\\.*[^\\d]", "", names(clusters), perl = TRUE))
    #clusnames <- paste0("clust", names(clusters))
    #print(clusnames)
    #clusnames <-str_replace(names(clusters),pattern = "RNA_snn_res.",replacement = "Clust")
    #names(clusters) <- clusnames
    
    for (clusnames in names(clusters)) {
      clusters[[clusnames]] <-
        paste(clusnames, clusters[[clusnames]], sep = 'C')
    }
    
    #print(clusnames)
    samples <- rownames(clusters)
    clusters <- cbind(clusters, samples)
    
    clusters <- checkIfNotTree(clusters)
    
    #print(clusters)
    
    tree <- TreeIndex(clusters)
    rownames(tree) <- rownames(clusters)
    
    TreeSE_obj <-
      TreeSummarizedExperiment(SimpleList(counts = count_matrix), colData = tree)
    
  }



#setwd("./Documents/RScripts/TreeSE/")
load("pbmc_clustree.Rdata")
##pbmc example
clusterdata <- pbmc@meta.data
clusterdata <- clusterdata %>%
  select(starts_with("RNA_snn"))

str(clusterdata)
#Create PBMC object and graph
#graph <- clustree(pbmc , prop_filter = 0, return = "graph")

#Call func from here
pbmc_TreeSE <-
  reassign_and_collapse(clusterdata, GetAssayData(pbmc))
#str(GetAssayData(pbmc))

#Save Object
save(pbmc_TreeSE, file = "pbmc_TreeSE.Rdata")


#Vusaualization
app <- startMetaviz()

icicle_plot <-
  app$plot(pbmc_TreeSE, datasource_name = "SCRNA", tree = "col")

#unscaled heatmap
#heatmap <- app$chart_mgr$revisualize(chart_type = "HeatmapPlot", chart = icicle_plot)


# Identify the 10 most highly variable genes
pbmc <-
  FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#Generate Heatmap
facetZoom <- app$get_ms_object(chart_id_or_object = icicle_plot)

ms_list <- facetZoom$get_measurements()
top100 <- head(VariableFeatures(pbmc), 100)
subset_ms_list <- Filter(function(ms)
  ms@id %in% top100, ms_list)

app$chart_mgr$visualize(chart_type = "HeatmapPlot", measurements = subset_ms_list)

#Stop app and check if stopped
app$stop_app()
app$is_server_closed()


#check if Aggregate works
aggr <- aggregateTree(pbmc_TreeSE, selectedLevel = 3, by = "col")

icicle_plot <-
  app$plot(aggr, datasource_name = "SCRNA", tree = "col")
