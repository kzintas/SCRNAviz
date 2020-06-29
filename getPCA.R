compute_tsne<- function(treeseobject){
  tsne <- Rtsne(t(as.matrix(assays(treeseobject)$counts)))
  treeseobject@metadata[['tsne']]<- tsne
  treeseobject
}


getPCA=function(treeseobject) {
  " Compute PCA over all features for given samples
    
    \\describe{
    \\item{measurements}{Samples to compute PCA over}
    \\item{start}{Start of feature range to query }
    \\item{end}{End of feature range to query}
    }
    "
  if(is.null(metadata(treeseobject)$tsne)){
    treeseobject<- compute_tsne(treeseobject)
  }

  #SE_tsne<-Rtsne(t(as.matrix(assays(treeseobject)$counts)))
  #print(SE_tsne)
  measurements<-  metadata(treeseobject)$tsne
  data <- list()
  #tsne_plot <- data.frame(x = SE_tsne$Y[,1], y = SE_tsne$Y[,2])
  #ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))
  for (col in seq(colnames(treeseobject))) {
    #print(col)
    temp
    <- list(sample_id = colnames(treeseobject)[col], PC1 = measurements$Y[col,1], PC2 = measurements$Y[col,2])
    data[[col]] <- temp
  }
  print(data)
  
  #result <- list(data = unname(data), pca_variance_explained = ord$sdev[1:2])
  #return(result)
  return(data)
}

json <- toJSON( data )

