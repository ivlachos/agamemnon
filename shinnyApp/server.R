returnMatrix<-function(x,input){
    inputFeature=input[[input$level]]
    if(input$norm == TRUE){
      mat = log2(aggTax(x,lvl=input$level,out='matrix',norm=input$norm)+1)
    } else{
      mat = aggTax(x,lvl=input$level,out='matrix',norm=input$norm)
    }
    dat = list(mat=mat,inputFeature = inputFeature)
    return(dat)
}

shinyServer(function(input, output) {
  output$plotly1 <- renderPlotly({
    if(input$pd!="None") coll = factor(pData(MRobj)[,input$pd])
    else coll = 1

    if(input$level!="OTU"){
      plist = list(level=input$level,norm=input$norm,
            Species=input$Species,Genus=input$Genus,Family=input$Family,
            Order=input$Order,Class=input$Class,Phylum=input$Phylum)
      vals = returnMatrix(MRobj,input=plist)
      mat  = vals$mat
      main = inputFeature = vals$inputFeature
    } else {
      if(input$norm == TRUE){
        mat = MRcounts(MRobj,norm=TRUE,log=TRUE)
      } else {
        mat = MRcounts(MRobj,norm=FALSE,log=FALSE)
      }
      main = inputFeature = input$feature
    }
    # ylabel = ifelse(input$norm,yes=expression(bold("Log"[2]*" Abundance")),no="No. raw reads") 
    df= data.frame(tax_abundance = mat[inputFeature,],sample = colnames(mat),phenotype = coll)
    if(input$box){
      p = ggplot(df,aes(x=phenotype,y=tax_abundance,fill=phenotype)) + geom_boxplot() + ggtitle(main) + theme_bw()+ theme(axis.text.x=element_blank()) #+ ylab(ylabel)
    } else {
      p = ggplot(df,aes(x=sample,y=tax_abundance,color=phenotype)) + geom_point() + ggtitle(main) + theme_bw()+ theme(axis.text.x=element_blank()) #+ ylab(ylabel)
    }
    ggplotly(p)    
  })
  output$pcaPlot <- renderPlotly({
    pd = factor(pData(MRobj)[,input$ppd])
    x = calcPCs(MRobj,n=nrow(MRobj),pch=21,bg=pd,usePCA=input$pcaOrMds,
      comp=c(input$dimensionx,input$dimensiony),
      useDist=input$useDist,distfun=vegan::vegdist,dist.method=input$distance)
    df = data.frame(xpc = round(x[[1]][,1],3),ypc = round(x[[1]][,2],3),samples = rownames(x[[1]]),pd)
    a <- list(title = x[[2]]); b <- list(title = x[[3]])
    plot_ly(df,x = ~xpc,y =~ypc,color = ~pd, text=~samples,type="scatter") %>% layout(xaxis=a,yaxis=b)
  })
  output$diversity <- renderPlotly({
    pdata = pData(MRobj)[,input$dpd]
    mat = t(MRcounts(MRobj,norm=FALSE,log=FALSE))
    H = vegan::diversity(mat,index=input$diversity)
    ggplot(data.frame(H,pdata),aes(x=pdata,y=H,fill=pdata)) + geom_boxplot() + ggtitle("Diversity index") + theme_bw()
  })
  
  output$ntaxa <- renderPlotly({
    pdata = factor(pData(MRobj)[,input$spd])
    ntaxa = colSums(MRcounts(MRobj)>0)
    p = ggplot(data.frame(ntaxa,samples=colnames(MRobj),pdata),aes(x=samples,y=ntaxa,color=pdata)) + geom_point()+ theme_bw()+ theme(axis.text.x=element_blank())
    ggplotly(p)
  })
  
  output$sparsity <- renderPlotly({
    pdata = factor(pData(MRobj)[,input$spd])
    ntaxa = colSums(MRcounts(MRobj)>0)
    p = ggplot(data.frame(ntaxa,pdata),aes(x=pdata,y=ntaxa,fill=pdata)) + geom_boxplot() + ggtitle("Number of taxa detected") + theme_bw()
    ggplotly(p)
  })

  output$sparsityPCA <- renderPlot({
    keep = 1:ncol(MRobj)
    dat = MRobj[,keep]
    useDist = input$useDist

    pd = factor(pData(dat)[,input$spd])
    res = plotOrd(dat,n=200,pch=21,bg=pd,usePCA=input$pcaOrMds,
      comp=c(1,2),useDist=useDist,distfun=vegan::vegdist,
      dist.method=input$distance)
  
    ntaxa = colMeans(MRcounts(dat)>0)
    plot(y=res[,1],x=ntaxa,ylab="PC1",xlab="Detection rate",
      bg=pd,pch=21)
  })
  output$diversityTable <- renderTable({
    pdata = pData(MRobj)[,input$dpd]
    mat = t(MRcounts(MRobj,norm=FALSE,log=FALSE))
    H = vegan::diversity(mat,index=input$diversity)

    muH =as.vector(by(H,pdata,mean))
    sdH =by(H,pdata,sd)

    divs = rbind(muH,sdH)
    rownames(divs) = c("Mean","SD")
    colnames(divs) = levels(pdata)
    divs
  })
  
  output$table <- renderDataTable({
          fd = fData(MRobj)[,-c(1,2,9,10)]
          Index = 1:nrow(fd)
          Rownames = rownames(fd)
          fd = cbind(Index,Rownames,fd)
          return(fd)
    })

  output$plotRare<-renderPlot({

    cl = factor(pData(dat)[,input$rpd])
    numFeatures = colSums(MRcounts(dat)>0)
    totalCounts = libSize(dat)
    plot_ly(data.frame(numFeatures,totalCounts,cl),x=~totalCounts,y=~numFeatures,color=~cl)    
  })

  output$plotHeatmap<-renderPlot({
    keep = 1:ncol(MRobj)
    dat = MRobj[,keep]
    trials = factor(pData(dat)[,input$hpd])
    heatmapColColors=brewer.pal(12,"Set3")[as.integer(trials)];
    heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50);
    if(input$heatmap_level=='OTU'){
      subset = MRcounts(dat,norm=TRUE,log=TRUE)
      rownames(subset) = paste(rownames(subset),fData(dat)[,"Genus"],sep=":")
    } else {
      subset = log2(aggTax(dat,lvl=input$heatmap_level,norm=TRUE,out='matrix')+1)
    }
    plotMRheatmap(subset,n=input$heatNumber,fun=input$heat,
              main="Bacterial Abundance Heatmap",
              cexRow=.95,cexCol=0.4,trace="none",
              dendrogram="column", key=TRUE,
              lwid=c(1,4), lhei=c(1,4),
              margins=c(2,10),
              col = heatmapCols,ColSideColors = heatmapColColors)
    legend("left",fill=unique(heatmapColColors),legend=unique(trials),title="Column labels")
  })

  output$phenolist<- renderDataTable({
        pd = as.matrix(pData(MRobj))
        cbind(colnames(MRobj),pd)
  })

  output$featurelist<- renderDataTable({
        fd = as.matrix(fData(MRobj)[,-10])
        cbind(rownames(MRobj),fd)
  })

  output$clusterSequences<-renderText({
    plist = list(level=input$level,norm=input$norm,level=input$level,
            species=input$species,genus=input$genus,family=input$family,
            order=input$order,class=input$class,phylum=input$phylum)
    inputFeature=plist[[plist$level]]
    kk = which(rowSums(MRcounts(MRobj)>0)>20)
    k = which(fData(MRobj)[kk,input$level]==inputFeature)
    otuids = paste(sprintf(">OTU_%s",rownames(MRobj)[kk[k]]),"\n",sep="")
    seqs = as.character(fData(MRobj)[kk[k],10])
    paste(otuids,seqs,collapse="\n",sep="")
  })
  
})
