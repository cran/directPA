#' Perturbation Plot 3D
#'
#' This function takes in a matrix of test statistics with two columns (3-dimensional space) and the 
#' annotation list such as pathway annotation or kinase-substrate annotation, and visualize the enrichment
#' of pathways or kinases in direction specific manner.
#' 
#' @usage perturbPlot3d(Tc, annotation, minSize=5, ...)
#' @param Tc a numeric matrix. The columns are genes or phosphorylation sites and the columns are treatments 
#' vs control statistics.
#' @param annotation a list with names correspond to pathways or kinases and elements correspond to genes or
#' substrates belong to each pathway or kinase, respectively.
#' @param minSize the size of annotation groups to be considered for calculating enrichment. Groups 
#' that are smaller than the minSize will be removed from the analysis.
#' @param ... parameters for controling the plot.
#' @return a list of coordinates for pathways or kinases
#' @export
#' 
perturbPlot3d <- function(Tc, annotation, minSize=5, ...) {
  
  # step 1. convert statistics into z-scores
  Tc.zscores <- apply(Tc, 2, function(x){qnorm(rank(x)/(nrow(Tc)+1))})
  
  # step 2. filter the groups that are smaller than the minimun cutoff
  DE = lapply(annotation, function(x){
    if(sum(rownames(Tc.zscores) %in% x) >= minSize) {
      X <- Tc.zscores[rownames(Tc.zscores)%in%x,]
      n = nrow(X)
      Z1 = sum(X[,1])/sqrt(n)
      Z2 = sum(X[,2])/sqrt(n)
      Z3 = sum(X[,3])/sqrt(n)
      list(Z1=Z1, Z2=Z2, Z3=Z3)
    }
  })
  
  # step3. filter DE that has 0 element
  DE <- DE[which(sapply(DE, length) != 0)]
  Z1 <- unlist(sapply(DE, function(x){x[1]}))
  Z2 <- unlist(sapply(DE, function(x){x[2]}))
  Z3 <- unlist(sapply(DE, function(x){x[3]}))
  
  # visualization
  plot3d(Z1, Z2, Z3, col="darkblue", pch=16, xlab=colnames(Tc)[1], ylab=colnames(Tc)[2], ...)
  text3d(Z1, Z2, Z3, names(DE), col="black", cex=1)
  abclines3d(x=0, y=0, z=0, a=diag(3), lwd=3, col="gold")
  
  ## return the results
  result <- list()
  result$Z1 <- Z1
  result$Z2 <- Z2
  result$Z3 <- Z3
  return(result)
}