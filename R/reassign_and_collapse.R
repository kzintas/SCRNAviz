check_alternate <- function(sub_df, all_df)
{
  #check where clusters from not-core edges should actually belong
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

  assign_df
}



change_assignment <- function(Df, Cluster_obj) {
  #for each row in the core edge false dataset, find out corresponding
  #source and sink in Seurat, and change accordingly
  for (i in seq(1:nrow(Df))) {
    sink_res = match(paste0('cluster', Df$`Sink_node_res`[i]),
                     colnames(Cluster_obj))
    source_res = match(paste0('cluster', Df$`Assign_from_res`[i]),
                       colnames(Cluster_obj))


    for (j in 1:nrow(Cluster_obj)) {
      if (as.numeric(as.character(Cluster_obj[[j, sink_res]])) == as.numeric(Df$Sink_node_clust[i])  &&
          as.numeric(as.character(Cluster_obj[[j, source_res]])) == as.numeric(Df$Assign_from_clust[i])) {
        Cluster_obj[[j, source_res]] <-
          as.factor(Df$`Assign_to_clust`[i])
      }

    }
  }
  return (Cluster_obj)

}

check_cycle <- function(pruned_graph) {
  #Check if a diamond or circle still exists in tree
  delete_set_vertices <- vector('numeric')
  ver_list <- V(pruned_graph)

  for (nodes in ver_list) {
    adj_edge <- incident(pruned_graph, nodes, mode = "in")

    if (length(adj_edge) > 1) {
      #Not a Tree yet
      adj_ver <-
        adjacent_vertices(pruned_graph, nodes, mode = "in")
      adj_ver <- as_ids(adj_ver[[1]])
      str(adj_ver)
      remove_node <- sample(adj_ver, 1)

      remove_edge <-
        incident(pruned_graph, remove_node, mode = "all")

      pruned_graph <- delete_edges(pruned_graph, remove_edge)
      delete_set_vertices <-
        c(delete_set_vertices, remove_node)

    }
  }
  pruned_graph <- delete.vertices(pruned_graph, delete_set_vertices)
  #as_tbl_graph(pruned_graph)

}


prune_tree <- function(graph_Df, cluster_df) {
  #Prune the tree so only core edges remain
  print("Remove not core edges")
  repeat {
    #drop duplicated columns (Throws error otherwise)

    graph_Df <- graph_Df[!duplicated(names(graph_Df))]

    #See number of core edges at each ites
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
  Graph <- check_cycle(Graph)

  objList <-
    list("Cluster_obj" = cluster_df, "Clustree_obj" = Graph)
  #return (Graph)
}





collapse_tree <- function(original_graph) {
  #collapse TREE
  #Distances of nodes and layers
  node_dists <- distances(original_graph, to = 1)
  layered_dist <- unique(distances(original_graph, to = 1))

  #Which vertices and edges we want to delete

  delete_set_vertices <- vector('numeric')

  #Vertex and edge lists of the graph from which we will construct the collapsed graph
  ver_list <-
    igraph::as_data_frame(original_graph, what = "vertices")


  i <- length(layered_dist)
  while (i >= 2) {
    prev_layer <- which(node_dists == layered_dist[i - 1])
    current_layer <- which(node_dists == layered_dist[i])

    if (length(prev_layer) == length(current_layer)) {
      while (length(prev_layer) == length(current_layer)) {
        delete_set_vertices <- c(delete_set_vertices, prev_layer)
        i <- i - 1

        prev_layer <-
          which(node_dists == layered_dist[i - 1])

      }
    }
    i <- i - 1
  }



  ver_list <- ver_list[-delete_set_vertices,]
  #print(ver_list)
  return(ver_list)
}




plot_tree <-
  function(hierarchydf) {
    print("abc")
    hierarchydf <-
      hierarchydf[,!colnames(hierarchydf) %in% c("samples", "otu_index")]

    df <- data.frame(from = numeric(), to = numeric())
    for (i in seq(ncol(hierarchydf) - 1)) {
      edges <- hierarchydf %>%
        dplyr::rename(from = colnames(hierarchydf)[i],
               to = colnames(hierarchydf)[i + 1]) %>%
        select(i, i + 1) %>%
        unique

      df <- rbind(df, edges)
    }

    mygraph <- graph_from_data_frame(df)
    plot(mygraph, layout=layout_as_tree, vertex.size=4,
         vertex.label.dist=1,  edge.arrow.size=0.5)
    #print(mygraph)
  #  ggraph::ggraph(mygraph, layout = 'dendrogram', circular = FALSE) +
  #    ggraph::geom_edge_diagonal() +
  #    ggraph::geom_node_point(show.legend=TRUE)
  }

checkIfNotTree <- function(cluster_df) {
  #Handle Forests
  cols <- colnames(cluster_df)
  if (length(unique(cluster_df[[1]])) > 1) {
    cluster_df$root <- "AllClusters"
    cols <- c("root", cols)
  }

  cluster_df[cols]
}


check_unique_parent <- function(clusterdata) {
  #check if user provided data has unique parents at each level
  for (i in seq(2, ncol(clusterdata))) {
    childs <- unique(clusterdata[[i]])
    for (values in childs) {
      subsetted_list <- clusterdata %>%
        filter(.data[[colnames(clusterdata)[[i]]]] == values)
      parent <- length(unique(subsetted_list[[i - 1]]))
      if (parent > 1) {
        stop("Not a tree, some nodes with multiple parents in level", i)
      }
      #cat(nrow(subsetted_list)," ",parent, "\n")
    }
  }
}



rename_clusters <- function(clusdata) {
  #convert cluster names to 'cluster1', 'cluster2'...
  clusnames <- seq(length(colnames(clusdata)))

  clusnames <- paste0("cluster", clusnames)
  #print(clusnames)
  #colnames(clusdata) <- unique(clusnames[(order(clusnames))])
  colnames(clusdata) <- clusnames
  clusdata
}


#' Creates the `TreeSE` object from provided hierarchy dataframe and count matrix
#'
#' This module returns the Tree Summarized Experiment if the dataframe is tree, throws error otherwise
#'
#' @param Dataframe containing cluster information at different resolutions
#' @param matrix Dense or sparse matrix containing the count matrix
#' @return TreeSummarizedExperiment Object that can be visualized with metavizr
#' @export
#' @examples
#'
simplified_treese <-
  function(cluster_df, count_matrix) {
    cluster_df <- rename_clusters(cluster_df)
    clusters <- cluster_df

    for (clusnames in names(clusters)) {
      clusters[[clusnames]] <-
        paste(clusnames, clusters[[clusnames]], sep = 'C')
    }

    check_unique_parent(clusters)
    samples <- rownames(clusters)
    clusters <- cbind(clusters, samples)

    clusters <- checkIfNotTree(clusters)

    #print(clusters)

    tree <- TreeIndex(clusters)
    #str(tree)
    #View(tree)
    rownames(tree) <- rownames(clusters)

    TreeSE_obj <-
      TreeSummarizedExperiment(SimpleList(counts = count_matrix), colData = tree)
    plot_tree(TreeSE_obj@colData@hierarchy_tree)
    TreeSE_obj


  }

#' Creates the `TreeSE` object from provided dataframe and count matrix
#'
#' This module Reassigns Clusters in Dataframe Object
#'
#' @param Dataframe containing cluster information at different resolutions
#' @param matrix Dense or sparse matrix containing the count matrix
#' @return TreeSummarizedExperiment Object that can be visualized with metavizr
#' @export
#' @examples
#'

reassign_and_collapse <-
  function(cluster_df, count_matrix) {
    cluster_df <- rename_clusters(cluster_df)

    #Create clustree object as graph and get the dataframe
    clustree_graph <-
      clustree(
        cluster_df ,
        prefix = "cluster",
        prop_filter = 0,
        return = "graph"
      )
    graph_df <- as_long_data_frame(clustree_graph)

    #Get pruned tree with only true core edges

    modified_obj <- prune_tree(graph_df, cluster_df)

    #modified graph and seurat object
    modified_graph = modified_obj$Clustree_obj
    clusters <-  modified_obj$Cluster_obj


    collapsed_graph <- collapse_tree(modified_graph)

    #collapsed_graph_df <- as_long_data_frame(collapsed_graph)
    cluster_names <-
      unique(sapply(strsplit(collapsed_graph$node, "C"), '[', 1))

    clusters <- clusters[, cluster_names]

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
    plot_tree(TreeSE_obj@colData@hierarchy_tree)
    TreeSE_obj

  }


#' Wrapper for `Seurat`
#'
#' This module Reassigns Clusters in Dataframe Object
#'
#' @param Seurat_object S4 object of class Seurat containing cluster information at different resolutions
#' @return TreeSummarizedExperiment Object that can be visualized with metavizr
#' @export
#' @examples
#'

visualizeSeurat <-
  function(Seurat_object) {
    clusterdata <- Seurat_object@meta.data
    #print(colnames(clusterdata))
    clusterdata <-
      clusterdata[, grep("*snn*", colnames(clusterdata))]
    #print(colnames(clusterdata))
    pbmc_TreeSE <-
      reassign_and_collapse(clusterdata, GetAssayData(Seurat_object))

  }


#' Wrapper for `SingleCellExperiment``
#'
#' This module Reassigns Clusters in Dataframe Object
#'
#' @param `SingleCellExperiment` S4 object of class Single Cell Experiment containing cluster information at different resolutions
#' @return TreeSummarizedExperiment Object that can be visualized with metavizr
#' @export
#' @examples
#'
visualizeSingleCellExperiment <-
  function(Sce_object) {
    clusterdata <- colData(Sce_object)
    #colnames(clusterdata)
    #clusterdata <- clusterdata[ , grep("sc3_", colnames(clusterdata))]
    clusterdata <-
      clusterdata[, grep("cluster", colnames(clusterdata))]
    count <- counts(Sce_object)
    rownames(count) <- rownames(counts(Sce_object))

    TreeSE <-
      reassign_and_collapse(as.data.frame(clusterdata@listData), count)


  }


#' Finds top variable genes `TreeSE`
#'
#' Finds the n top variable genes within the genes present in Summarized Experiment object
#'
#' @param `TreeSE` S4 object of class TreeSummarizedExperiment
#' @param `number` Number of top genes to be calculated
#' @return TreeSummarizedExperiment Object with added top_variable_gene information in metadata
#' @export
#' @examples
#'
#'
find_top_variable_genes <-
  function(treeseobject, number) {
    dec.Tree_SE <- modelGeneVar(assays(treeseobject)$counts)
    top_n <- getTopHVGs(dec.Tree_SE, n = number)
    treeseobject@metadata[['top_variable']] <- top_n
    treeseobject
  }
