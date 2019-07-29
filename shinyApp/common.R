source("importData.R")

downloadNotInstalled <- function(x)
  {
  for(i in x)
    {
    if(!require(i, character.only = TRUE))
      {
      source("https://bioconductor.org/biocLite.R")
      biocLite(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}

requiredPackages = c("shiny", "vegan", "RColorBrewer", "plotly", "ggplot2",
                    "gplots", "metagenomeSeq", "biomformat", "matrixStats", "rafalib")

downloadNotInstalled(requiredPackages)

MRobj = hmp_metag

calcPCs <- function(obj, tran = TRUE, comp = 1:2, norm = TRUE, log = TRUE, 
  usePCA = TRUE, useDist = FALSE, distfun = stats::dist,
  dist.method = "euclidian", n = NULL, ...)
  {
    mat = returnAppropriateObj(obj, norm, log)
    
    if (useDist == FALSE & usePCA == FALSE) 
        stop("Classical MDS requires distances")
    
    if (is.null(n)) 
        n = min(nrow(mat), 1000)
    
    if (length(comp) > 2) 
        stop("Can't display more than two components")
    
    otusToKeep <- which(rowSums(mat) > 0)
    otuVars <- rowSds(mat[otusToKeep, ])
    otuIndices <- otusToKeep[order(otuVars, decreasing = TRUE)[seq_len(n)]]
    mat <- mat[otuIndices, ]
    
    if (tran == TRUE) 
    {
        mat = t(mat)
    }
    
    if (useDist == TRUE)
    {
        d <- distfun(mat, method = dist.method)
    }
    
    else 
    {
        d = mat
    }
    
    if (usePCA == FALSE) 
    {
        ord = cmdscale(d, k = max(comp))
        xl = paste("MDS component:", comp[1])
        yl = paste("MDS component:", comp[2])
    }
    
    else 
    {
        pcaRes <- prcomp(d)
        ord <- pcaRes$x
        vars <- pcaRes$sdev^2
        vars <- round(vars/sum(vars), 5) * 100
        xl <- sprintf("PCA %s: %.2f%% variance", colnames(ord)[comp[1]], vars[comp[1]])
        yl <- sprintf("PCA %s: %.2f%% variance", colnames(ord)[comp[2]], vars[comp[2]])
    }
    
    list(ord[,comp], xl, yl)
}
