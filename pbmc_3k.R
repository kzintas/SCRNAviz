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

#setwd("./Documents/RScripts/TreeSE/")
# data set link - https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

#PBMC Seuratcode
pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200,
                           project = "pbmc3k")
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), print.output = 0, save.SNN = TRUE)

#create Clustree
clustree(pbmc)


graph<-clustree(pbmc , return="graph")

graph


#Pruning
subgraph<-graph%>% filter(is_core==TRUE)
subgraph


#Plotting
plot(subgraph, layout=layout_as_tree, vertex.size=4,
     vertex.label.dist=1,  edge.arrow.size=0.5)

plot(graph, layout=layout_as_tree, vertex.size=4,
     vertex.label.dist=1,  edge.arrow.size=0.5)


#Viewing
names(graph)
names(pbmc)
str(pbmc)
levels(pbmc$RNA_snn_res.0.3)
class(pbmc$RNA_snn_res.0.3)


levels<-data.frame(pbmc$RNA_snn_res.0.4,pbmc$RNA_snn_res.0.5)

#convert graph to dataframe
graph_df<-as_long_data_frame(graph)
nrow(graph_df)

#core_df<- data.frame()
#extract rows with false core
core_df <- subset(graph_df[c("from_node","to_node", "is_core")], is_core == FALSE)
core_df <-core_df[order(core_df$to_node, decreasing=TRUE), ]
pbmc$seurat_clusters

#core_df<-core_df[order(to)]
#for (i in 1:nrow(graph_df)){
  
#  if(graph_df[i,"is_core"]== FALSE)
#    rbind(core_df, graph_df[i,"from_node"], graph_df[i,"to_node"])
    #print(graph_df[i, "from_node"])
#}



#Saving
save(pbmc,graph,graph_df,pbmc_sce,pbmc_TreeSE, file="pbmc_clustree2.Rdata")
save(pbmc,graph, file="pbmc_clustree.Rdata")

unlink("pbmc_clustree.Rdata")

#loading
load("pbmc_clustree.Rdata")
load("pbmc_clustree2.Rdata")
