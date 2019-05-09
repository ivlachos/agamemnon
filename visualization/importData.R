setwd("../results/final/all")
library("Biobase")
library("metagenomeSeq")

i = list.files()[1]
x = read.csv(i,
	sep="\t",
	stringsAsFactors=FALSE)
pd = x[,-ncol(x)]

pd = data.frame(pd)
pd = AnnotatedDataFrame(pd)

mats = sapply(list.files(),function(i){
			x = read.csv(i,
			sep="\t",
			stringsAsFactors=FALSE)[,10]})

rownames(mats) = rownames(pd)
colnames(mats) = gsub("\\.tab","",colnames(mats))

hmp_metag = newMRexperiment(mats,featureData = pd)
phenoData(hmp_metag) = AnnotatedDataFrame(data.frame(id = colnames(hmp_metag)))
hmp_metag = filterData(hmp_metag,present=1,depth=1)

