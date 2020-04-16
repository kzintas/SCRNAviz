library(TreeSummarizedExperiment)
library(metagenomeSeq)
library(msd16s)
library(S4Vectors)
library(clustree)
library(Seurat)
library(SingleCellExperiment)
library(metavizr)



#Seurat function to check validity

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

#aggregateTree(pbmc_TreeSE, selectedLevel=3, by="col")

#convert seurat to sce
pbmc_sce<-as.SingleCellExperiment(pbmc)

#call function
pbmc_TreeSE <- ImportFromSeurat(pbmc, graph_modified)
pbmc_TreeSE
#pbmc@assays$RNA@
#sExp<-  SeuratToMR(pbmc, graph_modified)

#? colData
'''
SeuratToMR <- function(seurat, clustree, cluster_names = NULL) {
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
  
  #tree <- TreeIndex(clusters, cluster_names)
  #rownames(tree) <- colnames(GetAssayData(seurat))
  
  Seurat_sce<-as.SingleCellExperiment(seurat)
  
  #TreeSummarizedExperiment(SimpleList(counts = GetAssayData(seurat)), colData = tree)
  counts <- t(assays(Seurat_sce)$counts)
  tree <- as.data.frame(clusters, cluster_names)
  rownames(tree) <- rownames(counts)
  genes <- colnames(counts)
  genes <- as.data.frame(genes)
  rownames(genes) <- colnames(counts)
  
  newMRexperiment(counts, featureData = AnnotatedDataFrame( as.data.frame(tree)), phenoData = AnnotatedDataFrame(genes))
  
}




sce <- SingleCellExperiment(assays = list(counts = sc_example$counts,
                                          logcounts = sc_example$logcounts),
                            colData = sc_example$sc3_clusters,
                              reducedDims = SimpleList(TSNE = sc_example$tsne))
'''  

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


data(mouseData)
counts <- MRcounts(mouseData)
hierarchy <- fData(mouseData)
tree <- TreeIndex(hierarchy)
mbiome <- TreeSummarizedExperiment(SimpleList(counts=counts), rowData=tree, colData=pData(mouseData))
app <- startMetaviz()
icicle_plot <- app$plot(mbiome, datasource_name="mmssdd", tree = "row")

icicle_plot <- app$plot(pbmc_TreeSE, datasource_name="SCRNA", tree = "col")

app <- startMetaviz(host = "http://localhost:7782")
app$stop_app()
