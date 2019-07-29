#METHOD 1

DESeq2_call <- function(group1, group2, countMatrix, phenotypes, phenoChar) {
  
  suppressMessages(library(DESeq2))
  suppressMessages(library(qvalue))
  
  saveprefix <- "DEseq2"
  plimit <- 0.05 
  setupinfo <- data.frame("sample_ID" = phenotypes$External.ID, "Condition" = phenotypes[[phenoChar]])
  setupinfo$Condition <- factor(setupinfo$Condition, levels <- c(group1, group2))
  
  
  dds <- DESeqDataSetFromMatrix(countData = round(countMatrix), colData = setupinfo, design =~Condition)
  dds <- dds[, dds$Condition %in% c(group1, group2)]
  dds <- DESeq(dds, parallel = FALSE)
  
  res <- results(dds)
  res <- DESeq2::results(dds, alpha = plimit)
  results <- as.data.frame(res)
  results$foldChange <- 2 ^ results$log2FoldChange
  results <- results[c("foldChange", "log2FoldChange", "pvalue", "padj")]
  colnames(results) <- c("foldChange", "logFC", "PValue", "FDR")
  
  normalized.counts <- as.data.frame(counts(dds, normalized = TRUE))
  
  results <- merge(normalized.counts, results, by.x = 0, by.y = 0, all = FALSE)
  results <- results[order(results$FDR, results$PValue), ]
  
  name = paste0(saveprefix, ".Results.csv")
  
  qv <- qvalue(results$PValue)
  
  results$qValue <- qv$qvalues
  results$Local.FDR <- qv$lfdr
  
  return(results)
}



#METHOD 2

limma_voom_call <- function(group1, group2, countMatrix, phenotypes, phenoChar) {
  
  suppressMessages(library(limma))
  suppressMessages(library(qvalue))
  suppressMessages(library(edgeR))
  
  saveprefix <- "limma.voom"
  plimit <- 0.05 
  setupinfo <- data.frame("sample_ID" = phenotypes$External.ID, "Condition" = phenotypes[[phenoChar]])
  setupinfo$Condition <- factor(setupinfo$Condition, levels <- c(group1, group2))
  #setupinfo <- setupinfo[complete.cases(setupinfo), ]
  
  limma.design <- model.matrix(~factor(setupinfo$Condition, levels = c(group1, group2)))
  
  y <- DGEList(counts = countMatrix, group = factor(setupinfo$Condition, levels = c(group1, group2)))
  y <- edgeR::calcNormFactors(y)
  remove <- is.na(y$samples$group)
  y <- y[, !remove]
  v <- voom(y, limma.design, plot = FALSE)
  
  fit <- lmFit(v, limma.design)
  fit <- eBayes(fit) #this is the hierarchical test procedure
  
  results <- topTable(fit, coef = 2, number = nrow(fit), adjust.method = "BH", sort.by = "none")
  results$logFC <- results$logFC
  results$foldChange <- 2 ^ results$logFC
  results <- results[, c("foldChange", "logFC", "P.Value", "adj.P.Val")]
  
  colnames(results) <- c("foldChange", "logFC", "PValue", "FDR")
  setupinfo <- setupinfo[complete.cases(setupinfo), ]
  group.exp <- round(2 ^ (as.data.frame.EList(v, row.names = rownames(1:nrow(countMatrix)))), digits = 0)
  
  Group1.Mean <- rowMeans(group.exp[, which(setupinfo$Condition == group1)])
  Group2.Mean <- rowMeans(group.exp[, which(setupinfo$Condition == group2)])
  All.Mean <- (Group1.Mean + Group2.Mean) / 2
  Means <- cbind(Group1.Mean, Group2.Mean, All.Mean)
  
  colnames(Means) <- c(paste0("G1.", group1), paste0("G2.", group2), "All.Mean")
  
  results <- cbind(rownames(pd), group.exp, Means, results)
  
  #name = paste0(saveprefix, ".", group1, ".vs.", group2, ".", saveprefix, "PlotMA.jpg")
  
  results <- results[order(results$FDR, results$PValue), ]
  
  name = paste0(saveprefix, ".", group1, ".vs.", group2, ".", saveprefix, ".Results.csv")
  
  qv <- qvalue(results$PValue)
  
  results$qValue <- qv$qvalues
  results$Local.FDR <- qv$lfdr
  
  return(results)
}




#METHOD 3

limma_voom_weights_call <- function(group1, group2, countMatrix, phenotypes, phenoChar) {
  
  suppressMessages(library(limma))
  suppressMessages(library(qvalue))
  suppressMessages(library(edgeR))
  
  saveprefix <- "limma.voom.Weights"
  plimit <- 0.05 
  setupinfo <- data.frame("sample_ID" = phenotypes$External.ID, "Condition" = phenotypes[[phenoChar]])
  setupinfo$Condition <- factor(setupinfo$Condition, levels <- c(group1, group2))
  
  limma.design <- model.matrix(~factor(setupinfo$Condition, levels = c(group1, group2)))
  
  y <- DGEList(counts = countMatrix, group = factor(setupinfo$Condition, levels = c(group1, group2)))
  y <- edgeR::calcNormFactors(y)
  remove <- is.na(y$samples$group)
  y <- y[, !remove]
  v <- voomWithQualityWeights(y, limma.design, plot = FALSE)
  fit <- lmFit(v, limma.design)
  fit <- eBayes(fit)
  
  results <- topTable(fit, coef = 2, number = nrow(fit), adjust.method = "BH", sort.by = "none")
  results$logFC <- results$logFC
  results$foldChange <- 2 ^ results$logFC
  results <- results[, c("foldChange", "logFC", "P.Value", "adj.P.Val")]
  colnames(results) <- c("foldChange", "logFC", "PValue", "FDR")
  
  setupinfo <- setupinfo[complete.cases(setupinfo), ]
  group.exp <- round(2 ^ (as.data.frame.EList(v)), digits = 0)
  Group1.Mean <- rowMeans(group.exp[, which(setupinfo$Condition == group1)])
  Group2.Mean <- rowMeans(group.exp[, which(setupinfo$Condition == group2)])
  All.Mean <- (Group1.Mean + Group2.Mean) / 2
  Means <- cbind(Group1.Mean, Group2.Mean, All.Mean)
  colnames(Means) <- c(paste0("G1.", group1), paste0("G2.", group2), "All.Mean")
  
  results <- cbind(rownames(pd), group.exp, Means, results )
  
  name = paste0(saveprefix, ".", group1, ".vs.", group2, ".", saveprefix, "PlotMA.jpg")
  
  results <- results[order(results$FDR, results$PValue), ]
  
  name = paste0(saveprefix, ".", group1, ".vs.", group2, ".", saveprefix, ".Results.csv")
  qv <- qvalue(results$PValue)
  
  results$qValue <- qv$qvalues
  results$Local.FDR <- qv$lfdr
  
  return(results)
}




#METHOD 4

edgeR_RLRTRobust_call <- function(group1, group2, countMatrix, phenotypes, phenoChar) {
  
  suppressMessages(library(edgeR))
  suppressMessages(library(qvalue))
  
  saveprefix <- "RLRTRobust"
  
  plimit <- 0.05 
  setupinfo <- data.frame("sample_ID" = phenotypes$External.ID, "Condition" = phenotypes[[phenoChar]])
  setupinfo$Condition <- factor(setupinfo$Condition, levels <- c(group1, group2))
  
  design <- model.matrix(~factor(setupinfo$Condition, levels = c(group1, group2)))
  
  dge <- DGEList(countMatrix, group = factor(setupinfo$Condition, levels = c(group1, group2)))
  remove <- is.na(dge$samples$group)
  dge <- dge[, !remove]
  dge <- edgeR::calcNormFactors(dge)
  dge <- estimateGLMRobustDisp(dge, design = design)
  
  fit <- glmFit(dge, design = design)
  lrt <- glmLRT(fit)
  
  results <- topTags(lrt, n = Inf)
  results <- results$table
  results$foldChange <- 2 ^ results$logFC
  results <- results[, c("logFC", "foldChange", "PValue", "FDR")]
  results <- cbind(rownames(pd), results)
  results <- results[order(results$FDR, results$PValue), ]
  
  name = paste0(saveprefix, ".", group1, ".vs.", group2, ".", saveprefix, ".Results.csv")
  
  qv <- qvalue(results$PValue)
  
  results$qValue <- qv$qvalues
  results$Local.FDR <- qv$lfd
  
  return(results)
}




#METHOD 5

edgeR_RQLF_call <- function(group1, group2, countMatrix, phenotypes, phenoChar) {
  
  suppressMessages(library(edgeR))
  suppressMessages(library(qvalue))
  
  saveprefix<-"RQLF"
  
  plimit <- 0.05 
  setupinfo <- data.frame("sample_ID" = phenotypes$External.ID, "Condition" = phenotypes[[phenoChar]])
  setupinfo$Condition <- factor(setupinfo$Condition, levels <- c(group1, group2))
  
  design <- model.matrix(~factor(setupinfo$Condition, levels = c(group1, group2)))
  
  dge <- DGEList(countMatrix, group = factor(setupinfo$Condition, levels = c(group1, group2)))
  remove <- is.na(dge$samples$group)
  dge <- dge[, !remove]
  dge <- edgeR::calcNormFactors(dge)
  dge <- estimateDisp(dge, design = design)
  
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit)
  
  results <- topTags(qlf, n = Inf)
  results <- results$table
  results$foldChange <- 2 ^ results$logFC
  results <- results[, c("logFC", "foldChange", "PValue", "FDR")]
  results <- cbind(rownames(pd), results)
  results <- results[order(results$FDR, results$PValue), ]
  
  name = paste0(saveprefix, ".", group1, ".vs.", group2, ".", saveprefix, ".Results.csv")
  
  qv <- qvalue(results$PValue)
  
  results$qValue <- qv$qvalues
  results$Local.FDR <- qv$lfd
  
  return(results)
}




#METHOD 6

metagenomeSeq_call <- function(group1, group2, countMatrix, phenotypes, phenoChar) {
  #rm(group1, group2,  phenoChar)
  suppressMessages(library(metagenomeSeq))
  suppressMessages(library(qvalue))
  # countMatrix <- mats
  saveprefix <- "MSeq"
  plimit <- 0.05 
  # countMatrix<- mats
  # group1 <- "CD"
  # group2 <- "UC"
  # phenoChar <- "diagnosis"
  setupinfo2 <- data.frame("Condition" = phenotypes[[phenoChar]])
  rownames(setupinfo2) <- rownames(phenotypes)
  setupinfo2 <- subset(setupinfo2, Condition == group1 | Condition == group2)
  setupinfo2$Condition <- factor(setupinfo2$Condition, levels <- c(group1, group2))
  countMatrix <- countMatrix[, colnames(countMatrix) %in% rownames(setupinfo2)]
  countMatrix <- countMatrix[, order(colnames(countMatrix))]
  # rownames(countMatrix) <- x$TaxID
  obj <- newMRexperiment(countMatrix, phenoData = new("AnnotatedDataFrame", data = setupinfo2))
  
  # phenoData(obj) = AnnotatedDataFrame(setupinfo2) #new("AnnotatedDataFrame", data = setupinfo2))
  p <- cumNormStatFast(obj)
  obj <- cumNorm(obj, p = p)
  
  mod <- model.matrix(~factor(setupinfo2$Condition, levels = c(group1, group2)), data = pData(obj))
  res <- fitFeatureModel(obj = obj, mod = mod, coef = 2)
  
  results <- MRtable(obj = res, number = Inf, by = 2, adjustMethod = "BH")
  # results$logFC <- results$logFC
  results$foldChange <- 2 ^ results$logFC
  results <- data.frame(TaxID = rownames(results), results)
  colnames(results) <- c("TaxID", "present.in.N.samples.in.group.0", "present.in.N.samples.in.group.1", 
                         "Counts.in.group.0", "Counts.in.group.1", "logFC", "se", "PValue", "FDR", "foldChange")
  
  results <- results[order(results$FDR, results$PValue), ]
  
  name = paste0(saveprefix, ".", group1, ".vs.", group2, ".", saveprefix, ".Results.csv")
  
  qv <- qvalue(results$PValue)
  
  results$qValue <- qv$qvalues
  results$Local.FDR <- qv$lfdr
  return(results)
}



returnMatrix <- function(x, input)
{
    inputFeature = input[[input$level]]
    if(input$norm == TRUE)
    {
      mat = log2(aggTax(x, lvl = input$level, out = 'matrix', norm = input$norm) + 1)
    } else
      {
      mat = aggTax(x, lvl = input$level, out = 'matrix', norm = input$norm)
      }
    dat = list(mat = mat, inputFeature = inputFeature)
    return(dat)
}

shinyServer(function(input, output) 
{
  output$plotly1 <- renderPlotly({
    if(input$pd != "None") coll = factor(pData(MRobj)[, input$pd])
    else coll = 1

      plist = list(level = input$level, norm = input$norm,
            Scientific_Name = input$Scientific_Name, Species = input$Species, Genus = input$Genus, Family = input$Family,
            Order = input$Order, Class = input$Class, Phylum = input$Phylum)
      vals = returnMatrix(MRobj, input = plist)
      mat  = vals$mat
      main = inputFeature = vals$inputFeature

    df = data.frame(tax_abundance = mat[inputFeature,], sample = colnames(mat), phenotype = coll)
    if(input$box)
    {
      p = ggplot(df, aes(x = phenotype, y = tax_abundance, fill = phenotype)) + geom_boxplot() + ggtitle(main) + theme_bw() + theme(axis.text.x = element_blank()) #+ ylab(ylabel)
    } else 
      {
        p = ggplot(df, aes(x = sample, y = tax_abundance, color = phenotype)) + geom_point() + ggtitle(main) + theme_bw() + theme(axis.text.x = element_blank()) #+ ylab(ylabel)
      }
    ggplotly(p)    
  })
  
  output$pcaPlot <- renderPlotly({
    pd = factor(pData(MRobj)[, input$ppd])
    x = calcPCs(MRobj, n = nrow(MRobj), pch = 21, bg = pd, usePCA = input$pcaOrMds,
      comp = c(input$dimensionx, input$dimensiony),
      useDist = input$useDist, distfun = vegan::vegdist, dist.method = input$distance)
    df = data.frame(xpc = round(x[[1]][, 1], 3), ypc = round(x[[1]][, 2], 3), samples = rownames(x[[1]]), pd)
    a <- list(title = x[[2]]); b <- list(title = x[[3]])
    plot_ly(df, x = ~xpc, y = ~ypc, color = ~pd, text = ~samples, type = "scatter") %>% layout(xaxis = a, yaxis = b)
  })
  
  output$diversity <- renderPlotly({
    pdata = pData(MRobj)[, input$dpd]
    mat = t(MRcounts(MRobj, norm = FALSE, log = FALSE))
    H = vegan::diversity(mat, index = input$diversity)
    ggplot(data.frame(H, pdata), aes( x = pdata, y = H, fill = pdata)) + geom_boxplot() + ggtitle("Diversity index") + theme_bw()
  })
  
  output$ntaxa <- renderPlotly({
    pdata = factor(pData(MRobj)[, input$spd])
    ntaxa = colSums(MRcounts(MRobj) > 0)
    p = ggplot(data.frame(ntaxa, samples = colnames(MRobj), pdata), aes(x = samples, y = ntaxa, color = pdata)) + geom_point() + theme_bw() + theme(axis.text.x = element_blank())
    ggplotly(p)
  })
  
  output$sparsity <- renderPlotly({
    pdata = factor(pData(MRobj)[, input$spd])
    ntaxa = colSums(MRcounts(MRobj) > 0)
    p = ggplot(data.frame(ntaxa, pdata), aes(x = pdata, y = ntaxa, fill = pdata)) + geom_boxplot() + ggtitle("Number of taxa detected") + theme_bw()
    ggplotly(p)
  })

  output$sparsityPCA <- renderPlot({
    keep = 1:ncol(MRobj)
    dat = MRobj[, keep]
    useDist = input$useDist

    pd = factor(pData(dat)[, input$spd])
    res = plotOrd(dat, n = 200, pch = 21, bg = pd, usePCA = input$pcaOrMds,
      comp = c(1, 2), useDist = useDist, distfun = vegan::vegdist,
      dist.method = input$distance)
  
    ntaxa = colMeans(MRcounts(dat) > 0)
    plot(y = res[, 1], x = ntaxa, ylab = "PC1", xlab = "Detection rate",
      bg = pd, pch = 21)
  })
  
  
  output$diversityTable <- renderTable({
    pdata = pData(MRobj)[, input$dpd]
    mat = t(MRcounts(MRobj, norm = FALSE, log = FALSE))
    H = vegan::diversity(mat, index = input$diversity)

    muH =as.vector(by(H, pdata, mean))
    sdH =by(H, pdata, sd)

    divs = rbind(muH, sdH)
    rownames(divs) = c("Mean", "SD")
    colnames(divs) = levels(pdata)
    divs
  })
  
  output$table <- renderDataTable({
          fd = fData(MRobj)[, -c(1, 2, 9, 10)]
          Index = 1:nrow(fd)
          Rownames = rownames(fd)
          fd = cbind(Index, Rownames, fd)
          return(fd)
    })

  
  output$plotHeatmap <- renderPlot({
    keep = 1:ncol(MRobj)
    dat = MRobj[, keep]
    trials = factor(pData(dat)[, input$hpd])
    heatmapColColors = brewer.pal(12, "Set3")[as.integer(trials)];
    heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50);
    if(input$heatmap_level == 'OTU')
    {
      subset = MRcounts(dat, norm = TRUE, log = TRUE)
      rownames(subset) = paste(rownames(subset), fData(dat)[, "Genus"], sep = ":")
    } else 
      {
      subset = log2(aggTax(dat, lvl = input$heatmap_level, norm = TRUE, out = 'matrix') + 1)
      }
    if (input$heatmap_level == "Phylum") {validate(need(input$heatNumber <= as.integer(length(unique(x$Phylum))), "Maximum number of OTUs at the Phylum level 16"))}
    else if (input$heatmap_level == "Class") {validate(need(input$heatNumber <= as.integer(length(unique(x$Class))), "Maximum number of OTUs at the Class level 30"))}
    else if (input$heatmap_level == "Order") {validate(need(input$heatNumber <= as.integer(length(unique(x$Order))), "Maximum number of OTUs at the Order level 58"))}
    else if (input$heatmap_level == "Family") {validate(need(input$heatNumber <= as.integer(length(unique(x$Family))), "Maximum number of OTUs at the Family level 112"))}
    else if (input$heatmap_level == "Genus") {validate(need(input$heatNumber <= as.integer(length(unique(x$Genus))), "Maximum number of OTUs at the Genus level 230"))}
    else if (input$heatmap_level == "Species") {validate(need(input$heatNumber <= as.integer(length(unique(x$Species))), "Maximum number of OTUs at the Species level 750"))}
    else if (input$heatmap_level == "Scientific_Name") {validate(need(input$heatNumber <= as.integer(length(unique(x$Scientific_Name))), "Maximum number of OTUs at the Strain | sub-Species level 1401"))}
    plotMRheatmap(subset, n = input$heatNumber, fun = input$heat,
              main = "Microbial Abundances Heatmap",
              cexRow = .95, cexCol = 0.4, trace = "none",
              dendrogram = "column", key = TRUE,
              lwid = c(1, 4), lhei = c(1, 4),
              margins = c(2, 10),
              col = heatmapCols, ColSideColors = heatmapColColors)
    legend("left", fill = unique(heatmapColColors), legend = unique(trials), title = "Label: ")
  })

  output$group1 <- renderUI({
        selectInput("g1", "Group 1:", choices = as.character(na.omit(unique(pData(MRobj)[, input$phen]))))
  })
  
  output$group2 <- renderUI({
        selectInput("g2", "Group 2:", choices = as.character(na.omit(unique(pData(MRobj)[, input$phen]))))
  })
  
  output$demethod <- renderUI({
        selectInput("dem", "DE Method:", choices = as.character(c("DESeq2", "limma voom", "limma voom weights", "edgeR RLRTRobust", "edgeR RQLF", "metagenomeSeq")))
  })
  
  
  rt <- reactive({
        if (input$dem == "DESeq2") {
          DEresults <- DESeq2_call(input$g1, input$g2, mats, phenotypes, input$phen)
        } else if (input$dem == "limma voom") {
          DEresults <- limma_voom_call(input$g1, input$g2, mats, phenotypes, input$phen)
        } else if (input$dem == "limma voom weights") {
          DEresults <- limma_voom_weights_call(input$g1, input$g2, mats, phenotypes, input$phen)
        } else if (input$dem == "edgeR RLRTRobust") {
          DEresults <- edgeR_RLRTRobust_call(input$g1, input$g2, mats, phenotypes, input$phen)
        } else if (input$dem == "edgeR RQLF") {
          DEresults <- edgeR_RQLF_call(input$g1, input$g2, mats, phenotypes, input$phen)
        } else if (input$dem == "metagenomeSeq") {
          DEresults <- metagenomeSeq_call(input$g1, input$g2, mats, phenotypes, input$phen)
        }
  })
  
  output$deanalysis <- renderDataTable({
      rt()
  })
  
  
  output$downloadData <- downloadHandler(
        filename = function() {"DE_results.csv"},
        content = function(fname){
        write.table(rt(), fname, sep = "\t", row.names = FALSE)
        }
  )
  
  output$phenolist <- renderDataTable({
        pd = as.matrix(pData(MRobj))
        cbind(colnames(MRobj), pd)
  })

  output$featurelist <- renderDataTable({
        fd = as.matrix(fData(MRobj)[, -10])
        cbind(rownames(MRobj), fd)
  })

  output$clusterSequences <- renderText({
    plist = list(level = input$level, norm = input$norm, level = input$level,
            species = input$species, genus = input$genus, family = input$family,
            order = input$order, class = input$class, phylum = input$phylum)
    inputFeature = plist[[plist$level]]
    kk = which(rowSums(MRcounts(MRobj) > 0) > 20)
    k = which(fData(MRobj)[kk, input$level] == inputFeature)
    otuids = paste(sprintf(" > OTU_%s", rownames(MRobj)[kk[k]]), "\n", sep = "")
    seqs = as.character(fData(MRobj)[kk[k], 10])
    paste(otuids, seqs, collapse = "\n", sep = "")
  })
  
})
