# install scater https://bioconductor.org/packages/release/bioc/html/scater.html
library(scater)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
library(loomR)
library(Seurat)
library(patchwork)
library(scRNAseq)


cell416B <- LunSpikeInData('416b')

cell416B <- CreateSeuratObject(counts = cell416B@assays@data@listData[["counts"]],
                          project = "sce")

cell416B<- NormalizeData(cell416B)
all.genes <- rownames(cell416B)
cell416B <- ScaleData(cell416B)
cell416B <- FindVariableFeatures(object = cell416B)
cell416B <- RunPCA(cell416B, features = VariableFeatures(object = cell416B))
cell416B <- FindNeighbors(cell416B, dims = 1:10)
cell416B <- FindClusters(cell416B, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), print.output = 0, save.SNN = TRUE)

#Create PBMC object and graph
graph_416B<-clustree(cell416B , prop_filter=0, return="graph")

clustree(cell416B)


cell416B_TreeSE<- reassign_and_collapse(graph, cell416B)


save(cell416B_TreeSE, file="cell416B_TreeSE.Rdata")


app <- startMetaviz()

icicle_plot <-
  app$plot(cell416B_TreeSE, datasource_name = "SCRNA", tree = "col")
heatmap <- app$chart_mgr$revisualize(chart_type = "HeatmapPlot", chart = icicle_plot)

app$stop_app()
app$is_server_closed()

#aggr<-aggregateTree(pbmc_TreeSE, selectedLevel=3, by="col")
#icicle_plot <-  app$plot(aggr, datasource_name = "SCRNA", tree = "col")

#str(clusters)


#remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
