library(TreeSummarizedExperiment)
library(metagenomeSeq)
library(msd16s)
library(S4Vectors)
library(clustree)
library(Seurat)
library(SingleCellExperiment)
library(metavizr)

pbmc_TreeSE <- ImportFromSeurat(pbmc, graph_modified)
pbmc_TreeSE
#pbmc@assays$RNA@
  
ImportFromSeurat <- function(seurat, clustree, cluster_names = NULL) {
    # get clusters from seurat
    clusters <- seurat@meta.data
    if(is.null(cluster_names)) {
      cluster_names <- colnames(clusters)
      cluster_names <- cluster_names[grepl("res", cluster_names)]
    }
    
    clusters <- clusters[, cluster_names]
    
    # get cluster nodes from clustree
    clusnodes <- as.data.frame(with_graph(clustree, .N()))
    clusnames <- clusnodes$cluster
    names(clusnames) <- clusnodes$node
    
    # parse and create a dataframish structure
    clusters <- lapply(colnames(clusters), function(cn) {
      col <- clusters[, cn]
      col <- paste0(cn, "C", col)
      for (i in names(clusnames)) {
        col <- replace(col, col==i, paste0(i, "-cluster", clusnames[[i]]))
      }
      col
    })
    
    names(clusters) <- cluster_names
    clusters <- as.data.frame(clusters)
    clusters$samples <- colnames(GetAssayData(seurat))
    cluster_names <- c(cluster_names, "samples")
    
    tree <- TreeIndex(clusters, cluster_names)
    rownames(tree) <- colnames(GetAssayData(seurat))
    
    TreeSummarizedExperiment(SimpleList(counts = GetAssayData(seurat)), colData = tree)
}

aggregateTree(pbmc_TreeSE, selectedLevel=3, by="col")

pbmc_sce<-as.SingleCellExperiment(pbmc)

? colData

"""

sce <- SingleCellExperiment(assays = list(counts = sc_example$counts,
                                          logcounts = sc_example$logcounts),
                            colData = sc_example$sc3_clusters,
                              reducedDims = SimpleList(TSNE = sc_example$tsne))
"""  
  counts <- t(assays(pbmc_sce)$counts)
  tree <- as.data.frame(colData(pbmc_TreeSE))
  rownames(tree) <- rownames(counts)
  genes <- colnames(counts)
  genes <- as.data.frame(genes)
  rownames(genes) <- colnames(counts)
  
  sExp <- newMRexperiment(counts, featureData = AnnotatedDataFrame( as.data.frame(tree)), phenoData = AnnotatedDataFrame(genes))

#sExp <- newMRexperiment(pbmc_TreeSE)
# ignore host if using the server, else need to run metaviz locally
app <- startMetaviz()

app$plot(sExp, datasource_name = "single_cell", 
         feature_order = colnames(fData(sExp)),
         type="LeafCounts")