
```

library(GEOquery)
library(Biobase)
library(AnnotationDbi)
library(BiocGenerics)
library(dplyr)
library(hgu133a.db)
library(readr)
library(limma)
library(knitr)
library(pvca)
library(stringr)
library(topGO)
library(ReactomePA)
library(ggnewscale)
library(clusterProfiler)
library(genefilter)
library(arrayQualityMetrics)
library(oligo)
library(ggplot2)
#memory.limit(size= 20000)
#Sys.setenv("VROOM_CONNECTION_SIZE" = 520000 * 2) 
url.stats4bioinfo <- "http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics"
source(file.path(url.stats4bioinfo, 'R-files/config.R'))
dir.base <- 'http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics'
dir.base.ASG1 <- 'http://pedagogix-tagc.univ-mrs.fr/courses/ASG1'
source(file.path(dir.base, 'R-files', 'config.R'))
setwd(dir.results)
print(paste("Result directory", dir.results))
```

### Download the data


``` 

gset <- getGEO("GSE68468", GSEMatrix = TRUE, getGPL = TRUE, AnnotGPL = TRUE)
length(gset)
gse <- gset[[1]]

```

### Transform

Log transformation will be needed throughout every step, thus this will occur first.

```
exprs(gse) <- log2(exprs(gse))
```

### Filter

```
pheno <- pData(gse)
head(pheno) ## check the columns of interest again.
if ((rownames(pheno)) == colnames(exprs(gse)))  
{
  print("TRUE")
} 

### 8 Groups, to reduce redundancy the only the main one is included
```

```
kable(table(pheno$characteristics_ch1.4), format = "markdown", caption = "Summary of samples in the experiment", col.names=c("Histology", "Sample Count"), align="lr")
```

```
filter1 <- colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: colon cancer"]

if 
(filter1 == colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: colon cancer"]){
   print("TRUE")
}

coloncancer <- gse[, filter1]

table(coloncancer@phenoData$characteristics_ch1.4)

filter2 <- colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal colon mucosa"]

if 
(filter2 == colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal colon mucosa"]){
  print("TRUE")
}

normalcolon <- gse[, filter2]

vs1 <- colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer" | gse@phenoData@data$"characteristics_ch1.4"=="histology: normal colon mucosa"]

if (vs1 == colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer" | gse@phenoData@data$"characteristics_ch1.4"=="histology: normal colon mucosa"])
{
  print("TRUE")
}

cancerVSnormal <- gse[, vs1]

table(cancerVSnormal@phenoData$characteristics_ch1.4)
```

### Quality control

```
plotDensities(cancerVSnormal, legend = FALSE, main = "Colon cancer versus normal colon density plot")

#aarrayQualityMetrics(cancerVSnormal, outdit = tempdir(), do.logtransform = FALSE)

oligo::boxplot(cancerVSnormal,  main = "Colon cancer versus normal colon boxplot")
```

### Quantile normalise

```
QcancerVSnormal <- cancerVSnormal
exprs(QcancerVSnormal) <- normalizeBetweenArrays(exprs(QcancerVSnormal))
oligo::boxplot(QcancerVSnormal)
```

```
par(mfrow=c(2,1))
oligo::boxplot(cancerVSnormal, main = "Before normalisation")
oligo::boxplot(QcancerVSnormal, main = "After normalisation")
```

```
par(mfrow=c(1,2))
plotDensities(cancerVSnormal, legend=FALSE, "Density plot before normalisation", main = "Before normalisation")
plotDensities(QcancerVSnormal, legend=FALSE, "Density plot after normalisation", main = "After normalisation")
```

### Dimensionality reduction

```
EX <- Biobase::exprs(QcancerVSnormal)
PCA1 <- prcomp(t(EX), scale = FALSE)
percentVar <- round(100*PCA1$sdev^2/sum(PCA1$sdev^2),1)
sd_ratio <- sqrt(percentVar[2]/percentVar[1])
dataGG <- data.frame(PC1 = PCA1$x[,1], PC2 = PCA1$x[,2],
            Histology = Biobase::pData(QcancerVSnormal)$characteristics_ch1.4)
ggplot(dataGG, aes(PC1, PC2)) + 
          geom_point(aes(colour = Histology)) + 
          ggtitle("PCA plot of normalised data") +
          xlab(paste0("PC1, VarExp:", percentVar[1], "%")) +
          ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          coord_fixed(ratio = sd_ratio) +
          scale_color_manual(values = c("darkorange2", "dodgerblue4"))

```

### Inspection of batch effects

```
pct_threshold <- 0.6
colnames(pData(QcancerVSnormal))
batch.factors <-c("source_name_ch1", "characteristics_ch1.4")
pvca0bj <- pvcaBatchAssess(QcancerVSnormal, batch.factors, pct_threshold)

sds <- apply(exprs(QcancerVSnormal), 1, sd)
sds0 <- sort(sds)
plot(1:length(sds0), sds0, main = "Distribution of variability for all genes", 
      sub="Vertical lines represent 90% and 95% percentiles",
      xlab="Gene index (from least to most variable)",
      ylab="Standard deviation")
 abline(v=length(sds)*c(0.9,0.95))
 
```

### More filtering and Annotation

```
AQcancervsnormal <- QcancerVSnormal
annotation(AQcancervsnormal) <- "hgu133a.db"
filtered <- nsFilter(AQcancervsnormal, 
                      require.entrez = TRUE, 
                      remove.dupEntrez = TRUE,
                      var.filter = TRUE,
                      var.func = IQR,
                      var.cutoff = 0.25,
                      filterByQuantile = TRUE,
                      feature.exclude = "^AFFX")
print(filtered$filter.log)
filtered_cancervsnorm <- filtered$eset
nrow(filtered_cancervsnorm)
head(filtered_cancervsnorm@featureData$"Gene Symbol")

```

### Linear modelling.

```
metadata <- pData(filtered_cancervsnorm)
design <- model.matrix(~0+metadata$characteristics_ch1.4)
colnames(design) <- c("Cancer", "Normal")
fit <- lmFit(filtered_cancervsnorm, design)
contrasts <- makeContrasts(Cancer - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
class(fit2)
results <- decideTests(fit2)
summary(results)
topTab <- topTable(fit2, number = nrow(fit2), adjust = "fdr")
cancerVSnormTopTable <- topTab[, ! names(topTab) %in% c("GB_ACC", "SPOT_ID", "Species.Scientific.Name", "Annotation.Date", "Sequence.Source", "Target.Description", "Representative.Public.ID", "Refeq.Transcript.ID", "Gene.Ontology.Biological.Process", "Gene.Ontology.Cellular.Component", "Gene.Ontology.Molecular.Function")]
cancerVSnormTopTable <- cancerVSnormTopTable[, ! names(cancerVSnormTopTable) %in% c("Sequence.Type", "Gene.Title", "RefSeq.Transcript.ID")]

results1 <- decideTests(fit2, method = "separate", adjust.method = "fdr", p.value = 0.1, lfc = 1)

sum.res.rows <- apply(abs(results1), 1, sum)

res.selected <- results1[sum.res.rows!=0,]
```

```

colnames(cancerVSnormTopTable)
names(cancerVSnormTopTable)[names(cancerVSnormTopTable) == "ENTREZ_GENE_ID"] <- "ENTREZID"

listOfTables <- list(ColonCancerVSNormalColon = cancerVSnormTopTable)
listofSelected <- list()

for (i in 1:length(listOfTables)) 
  {
  topTable <- listOfTables[[i]]
  whichGenes <- topTable["adj.P.Val"] < 0.15
  selectedIDs <- rownames(topTable)[whichGenes]
  EntrezIDs <- AnnotationDbi::select(hgu133a.db, selectedIDs, c("ENTREZID"))
  EntrezIDs <- EntrezIDs$ENTREZID
  listofSelected[[i]] <- EntrezIDs
   names(listofSelected)[i] <- names(listOfTables)[i]
}
sapply(listofSelected, length)
listOfSelected <- listofSelected


```

### Pathway Enrichment
#### Method 1

```

mapped_genes2GO <- mappedkeys(hgu133aGO)
mapped_genes2KEGG <- mappedkeys(hgu133aPATH)
mapped_genes <- union(mapped_genes2GO, mapped_genes2KEGG)
library(ReactomePA)

listOfData <- listOfSelected[1]
comparisonsNames <- names(listOfData)
universe <- mapped_genes

#for (i in 1:length(listOfData)){
#genesIn <- listOfData[[i]]
#comparison <- comparisonsNames[i]
#enrich.result <- enrichPathway(gene = genesIn,
                                #pvalueCutoff = 0.05,
                               # qvalueCutoff = 0.9,
                                # readable = T,
                                # pAdjustMethod = "BH",
                                # organism = "human",
                                # universe = universe)
  
#  cat("##################################")
 # cat("\nComparison: ", comparison,"\n")
 # print(head(enrich.result))
 
  # if (length(rownames(enrich.result@result)) != 0) {
 #  write.csv(as.data.frame(enrich.result), 
# file =paste0("./results/","ReactomePA.Results.",comparison,".csv"), 
            # row.names = FALSE)
   
#pdf(file=paste0("./results/","ReactomePABarplot.",comparison,".pdf"))
# print(barplot(enrich.result, showCategory = 15, font.size = 4, 
# title = paste0("Reactome Pathway Analysis for ", comparison,". Barplot")))
 # dev.off()
  
#  pdf(file = paste0("./results/","ReactomePAcnetplot.",comparison,".pdf"))
 #    print(cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
    #      vertex.label.cex = 0.75))
  # dev.off()
  # }
# }
```


#### Method 2

```
back <- subset(cancerVSnormTopTable, adj.P.Val < 0.1)$ID
backidx <- genefilter::genefinder(filtered_cancervsnorm, as.character(back), method = "manhattan", scale = "none")
backidx <- sapply(backidx, function(x)x$indices)
backs <- featureNames(filtered_cancervsnorm )[backidx]
backs <- setdiff(backs, back)
intersect(backs, back)

IDs <- rownames(cancerVSnormTopTable)
un <- IDs %in% c(backs, back)
sl <- IDs %in% back

genes <- sl[un]
genes <- factor(as.integer(sl[un]))
names(genes) <- IDs[un]

entrez1 <- mapIds(hgu133a.db, 
                 keys = rownames(results),
                 keytype = "PROBEID",
                 colum = "ENTREZID")

reactome <- enrichPathway(gene = entrez1[back], 
                          universe = entrez1[c(back, 
                                              backs)],
                          organism = "human",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.9, 
                          readable = TRUE)
reactome @result$Description <- paste0(str_sub(
    reactome @result$Description, 1, 20),
    "...")



top_Go_data <- new("topGOdata", ontology = "BP", allGenes = genes, nodeSize = 10, annot = annFUN.db, affyLib = "hgu133a.db")

topGoElim <- runTest(top_Go_data, algorithm = "elim", statistic = "Fisher")

topGoClassic <- runTest(top_Go_data, algorithm = "classic", statistic = "Fisher")


resTopGo <- GenTable(top_Go_data, Fisher.elim = topGoElim, Fisher.classic = topGoClassic,
                     orderBy = "Fisher.elim", topNodes = 100)

genesTopGo <- printGenes(top_Go_data, whichTerms = resTopGo$GO.ID, chip = "hgu133a.db", geneCutOff = 1000)

resTopGo$sig_genes <- sapply(genesTopGo, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"], ";"),
        collapse = "")
})

```

#### Method 3

```
##KEGG

entrez <- entrez1
enrich_kegg <- kegga(fit2, geneid = entrez, species = "Hs")
topKEGG(enrich_kegg)

##GO
enrich_go <- goana(fit2, geneid = entrez, species = "Hs")
topGO(enrich_go, ontology = "BP")

## Filtering for genes involved in the biological processes
bpGO <- cancerVSnormTopTable
x <- org.Hs.egGO2ALLEGS
Rkeys(x) <- "GO:0008150"
ID <- mappedLkeys(x)
i <- bpGO$ENTREZID %in% ID
bpGO <- bpGO[i,] 

## Filtering for genes involved in xenobiotic metabolic process
xmpGO <- cancerVSnormTopTable
x1 <- org.Hs.egGO2ALLEGS
Rkeys(x1) <- "GO:0006805"
ID1 <- mappedLkeys(x1)
i1 <- xmpGO$ENTREZID %in% ID1
i1_d <- as.numeric(i1)
xmpGO <- xmpGO[i1,]


```

#### Method 4

```
library(clusterProfiler)
## CC - cellular  component, MF - molecular function, BP - biological process.

geneList <- cancerVSnormTopTable[,4]
names(geneList) <- as.character(cancerVSnormTopTable[,3])
geneList <- sort(geneList, decreasing = TRUE)

ggoCC <- groupGO(gene = entrez,
               OrgDb = org.Hs.eg.db,
               ont = "CC",
               level = 3,
               readable = TRUE)

ggoMF <- groupGO(gene = entrez,
               OrgDb = org.Hs.eg.db,
               ont = "MF",
               level = 3,
               readable = TRUE)

ggoBP <- groupGO(gene = entrez,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               level = 3,
               readable = TRUE)
```


### Over-representation analysis

```
##egoCC <- enrichGO(gene = entrez,
                ##  OrgDb = org.Hs.eg.db,
                 ## universe = names(geneList), 
                 ## ont = "CC",
                 ## pAdjustMethod = "BH",
                 ## pvalueCutoff = 0.01,
                 ## qvalueCutoff = 0.05,
                 ## readable = TRUE)

##egoMF <- enrichGO(gene = entrez,
                 ## OrgDb = org.Hs.eg.db,
                 ## ont = "MF",
                 ## pAdjustMethod = "BH",
                 ## pvalueCutoff = 0.01,
                 ## qvalueCutoff = 0.05,
                 ## readable = TRUE)

##egoBP <- enrichGO(gene = entrez,
                ##  OrgDb = org.Hs.eg.db,
                 ## ont = "BP",
                 ## pAdjustMethod = "BH",
                 ## pvalueCutoff = 0.01,
                 ## qvalueCutoff = 0.05,
                 ## readable = TRUE)

##gene.df <- bitr(entrez, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)

##egoCC2 <- enrichGO(gene = gene.df$ENSEMBL,
                  ## OrgDb = org.Hs.eg.db,
                 ##  keyType = 'ENSEMBL',
                 ##  ont = "CC",
                 ##  pAdjustMethod = "BH",
                 ##  pvalueCutoff = 0.01,
                  ## qvalueCutoff = 0.05)

##egoMF2 <- enrichGO(gene = gene.df$ENSEMBL,
                 ##  OrgDb = org.Hs.eg.db,
                 ##  keyType = 'ENSEMBL',
                ##   ont = "MF",
                 ##  pAdjustMethod = "BH",
                 ##  pvalueCutoff = 0.01,
                 ##  qvalueCutoff = 0.05)

##egoBP2 <- enrichGO(gene = gene.df$ENSEMBL,
                 ##  OrgDb = org.Hs.eg.db,
                  ## keyType = 'ENSEMBL',
                  ## ont = "BP",
                  ## pAdjustMethod = "BH",
                  ## pvalueCutoff = 0.01,
                  ## qvalueCutoff = 0.05)

##GOaCC <- gseGO(geneList = geneList,
              ## OrgDb = org.Hs.eg.db,
              ## ont = "CC",
             ##  minGSSize = 100,
             ##  maxGSSize = 500,
             ##  pvalueCutoff = 0.05,
             ##  verbose = FALSE)
```

```{r}
##goplot(egoCC)
##cnetplot(egoCC, showCategory = 10)
```

### Differential expression with t testing

```
packages <- c("knitr")
for (pkg in packages) {
  if (!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

pkg <- "qvalue"
if (!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
  source("http://bioconductor.org/biocLite.R")
  biocLite();
  biocLite(pkg)
}
library(pkg, character.only = TRUE)

url.stats4bioinfo <- "http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics"
source(file.path(url.stats4bioinfo, 'R-files/config.R'))
dir.base <- 'http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics'
dir.base.ASG1 <- 'http://pedagogix-tagc.univ-mrs.fr/courses/ASG1'
source(file.path(dir.base, 'R-files', 'config.R'))
setwd(dir.results)
print(paste("Result directory", dir.results))
dir.base <- 'http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics'
dir.base.ASG1 <- 'http://pedagogix-tagc.univ-mrs.fr/courses/ASG1'
source(file.path(dir.base, 'R-files', 'config.R'))
setwd(dir.results)
print(paste("Result directory", dir.results))
source(file.path(dir.util, "util_student_test_multi.R"))
group.of.interest <- "histology: colon cancer"
```

### New filter - Different cut off

```
LDRfilter <- nsFilter(AQcancervsnormal, 
                      require.entrez = TRUE, 
                      remove.dupEntrez = TRUE,
                      var.filter = TRUE,
                      var.func = IQR,
                      var.cutoff = 0.75,
                      filterByQuantile = TRUE,
                      feature.exclude = "^AFFX")
print(LDRfilter$filter.log)
LDR_canVSnorm <- LDRfilter$eset
nrow(LDR_canVSnorm)
head(LDR_canVSnorm@featureData$"Gene Symbol")


table(LDR_canVSnorm$characteristics_ch1.4)
      
canVSnormal <- exprs(LDR_canVSnorm)

canVSnorm <- as.data.frame(canVSnormal)

vv <- LDR_canVSnorm$characteristics_ch1.4

colnames(canVSnorm) <- paste(colnames(canVSnorm), vv, sep = "_")
```

### Single Student Test. Random gene to test hypothesis

```

n.Cancer <- 186
n.Normal <- 55

g <- sample(1:nrow(canVSnorm), 1)

g.profile <-as.vector(as.matrix(canVSnorm[g,]))
plot.col <- c('histology: colon cancer' = '#4444BB', 'histology: normal colon mucosa' = '#FFFF88')
```

```
par(mar = c(5.1, 4.1, 4.1, 11))
barplot(g.profile, col=plot.col[vv], main = "One random gene expression in the data", ylab = "Relative expression", xlab = "Individual samples n=65" , xlim = c(1, 65))
legend('topright', c("histology: colon cancer", "histology: normal colon mucosa"), col=plot.col[c("histology: colon cancer", "histology: normal colon mucosa")], pch=15, bty="o", bg='white', cex = 0.65)
```

```
sample.Cancer <- g.profile[vv == "histology: colon cancer"]
sample.Normal <- g.profile[vv == "histology: normal colon mucosa"]

mean.est.Cancer <- mean(sample.Cancer)
mean.est.Normal <- mean(sample.Normal)

sample.sd.Cancer <- sd(sample.Cancer) * sqrt((n.Cancer-1)/n.Cancer)
sample.sd.Normal <- sd(sample.Normal) * sqrt((n.Normal-1)/n.Normal)

sd.est.Cancer <- sd(sample.Cancer)
sd.est.Normal <- sd(sample.Normal)


sd.err.est.Cancer <- sd(sample.Cancer) / sqrt(n.Cancer)
sd.err.est.Normal <- sd(sample.Normal) / sqrt(n.Normal)

diff.sd.est <- sqrt((n.Cancer * sample.sd.Cancer^2 + n.Normal * sample.sd.Normal^2) * (1/n.Cancer + 1/n.Normal) / (n.Cancer+n.Normal-2))

d <- abs(mean.est.Cancer - mean.est.Normal)
t.obs.Student <- d/diff.sd.est
P.val.Student <- 2*pt(q=t.obs.Student,
                      df=n.Cancer+n.Normal-2,
                      lower.tail = F)
t.student <- t.test(sample.Cancer, sample.Normal, var.equal = TRUE)

print(t.student)

### Welch's

t.welch <- t.test(sample.Cancer, sample.Normal, var.equal = FALSE)
 
print(t.welch)

t.statistics <- vector()
P.values <- vector()

for (g in 1:nrow(canVSnorm)) {
  print(paste("Random gene", g))
  g.profile <- as.vector(canVSnorm[g,])
  sample.Cancer <- g.profile[vv == "histology: colon cancer"]
  sample.Normal <- g.profile[vv == "histology: normal colon mucosa"]
  t <- t.test(sample.Cancer, sample.Normal)
  t.statistics <- append(t.statistics, t$statistic)
  P.values <- append(P.values, t$p.value)
}
print(P.values)

### Multiple variance test

canvsnorm.t.result <- t.test.multi(canVSnorm, vv)
dim(canvsnorm.t.result)
names(canvsnorm.t.result)
sum(canvsnorm.t.result$E.value <=1)
canvsnorm.E <- canVSnorm[canvsnorm.t.result$E.value <=1,]


```

#### T test multivariance function
```
t.test.multi  function (x, ## a matrix or data.frame
                          cl, ##  vector of classes
                          P.threshold=NA, ## select objects below a given threshold of P.value
                          E.threshold=NA, ## select objects below a given threshold of E.value
                          FDR.threshold=NA, ## select objects below a given threshold of FDR
			  robust.est=F, ## use robust estimators for central tendency and dispersion
			  verbosity=1, ## Level of verbosity
                          volcano.plot = T, ## Draw a volcano plot
                          alternative = "two.sided", ## Supported: "two.sided", "less", "greater"
                          ... ## additional parameters are passed to the function volcano.plot()
                          ) {
  ## Report starting time
  if (verbosity >= 1) { print (paste(date(), "Multiple t-test started", sep= " - ")) }
  
  ## Dimensions of the data set
  n <- nrow(x)
  p <- ncol(x)

  ## check the dimensions of the class vector
  if (length(cl) != p) {
    stop (paste ('The number of columns of the dataset (',
                 p,
                 ') should be the same as the length of the class vector (',
                 length(cl),
                 ')'))
  }

  classes <- unique(cl)
  if (robust.est) {
    means.1 <- apply(x[,cl==classes[1]],1,median,na.rm=T)
    means.2 <- apply(x[,cl==classes[2]],1,median,na.rm=T)
    iqr.1 <- apply (x[,cl==classes[1]], 1, IQR, na.rm=T)
    iqr.2 <- apply (x[,cl==classes[2]], 1, IQR, na.rm=T)
    sd.est.1 <- iqr.1/(qnorm(0.75) - qnorm(0.25))
    sd.est.2 <- iqr.2/(qnorm(0.75) - qnorm(0.25))
    var.est.1 <- sd.est.1^2
    var.est.2 <- sd.est.2^2
  } else {
    means.1 <- apply(x[,cl==classes[1]],1,mean,na.rm=T)
    means.2 <- apply(x[,cl==classes[2]],1,mean,na.rm=T)
    var.est.1 <- apply(x[,cl==classes[1]],1,var,na.rm=T)
    var.est.2 <- apply(x[,cl==classes[2]],1,var,na.rm=T)
    sd.est.1 <- sqrt(var.est.1)
    sd.est.2 <- sqrt(var.est.2)
  }


  ## 2012-05-02: restored means.diff = means.1 - means.2 for the sake
  ## of consistency with R funciton t.test()
  means.diff <- means.1 - means.2
  ##  means.diff <- means.2 - means.1
  ##  hist(means.diff,breaks=50)
  
  n.1 <- sum(cl == classes[1])
  n.2 <- sum(cl == classes[2])

  ## Calculate observed t value
  st.err.diff <- sqrt(var.est.1/n.1 + var.est.2/n.2)
  t.obs <- means.diff/st.err.diff
  ## Calculate degrees of freedom with Welch's formula  
  df.welch <- (var.est.1/n.1 + var.est.2/n.2)^2 / ((var.est.1/n.1)^2/(n.1-1) + (var.est.2/n.2)^2/(n.2-1))

  ## Calculate P-value and E-value
#  P.value.normal.approx <- 2*pnorm(abs(t.obs),lower.tail=F)
  if (alternative == "greater") {
    P.value <- pt(t.obs,df.welch,lower.tail=FALSE)
  } else if (alternative == "less") {
    P.value <- pt(t.obs,df.welch,lower.tail=TRUE)
  } else if (alternative == "two.sided") {
    P.value <- 2*pt(abs(t.obs),df.welch,lower.tail=F)
  } else {
    stop('Invalid alternative option for t.test.multi(). Supported: "two.sided", "less", or "greater"')
  }
  E.value <- P.value*nrow(x)
  sig <- -log(E.value)/log(10)

  multi.corr <- multitest.corrections(P.value, plots=FALSE)

  ## Collect all statistics in a data frame
  result <- data.frame(
                       means.1,
                       means.2,
                       means.diff,
                       var.est.1,
                       var.est.2,
                       sd.est.1,
                       sd.est.2,
                       st.err.diff,
                       t.obs,
                       df.welch,
#                       P.value.normal.approx,
                       P.value,
                       E.value,
		       sig)
  result$fwer <- multi.corr$multitest.table$fwer
  result$q.value <- multi.corr$multitest.table$qval.Storey
  result$fdr <- multi.corr$multitest.table$fdr
  result$rank <- multi.corr$multitest.table$rank

  if (robust.est) {
    names(result)[1:7] <- c(
                            paste("median.", classes[1], sep=""),
                            paste("median.", classes[2], sep=""),
                            "medians.diff",
                            paste("var.est.", classes[1], sep=""),
                            paste("var.est.", classes[2], sep=""),
                            paste("sd.est.", classes[1], sep=""),
                            paste("sd.est.", classes[2], sep="")
                            )
  } else {
    names(result)[1:7] <- c(
                            paste("mean.", classes[1], sep=""),
                            paste("mean.", classes[2], sep=""),
                            "means.diff",
                            paste("var.est.", classes[1], sep=""),
                            paste("var.est.", classes[2], sep=""),
                            paste("sd.est.", classes[1], sep=""),
                            paste("sd.est.", classes[2], sep="")
                            )
  }

  
  ## Filtering on P-value and E-value thresholds
  if (!is.na(P.threshold)) {
    result <- result[result$P.value < P.threshold,]
  }
  if (!is.na(E.threshold)) {
    result <- result[result$E.value < E.threshold,]
  }
  if (!is.na(FDR.threshold)) {
    result <- result[result$fdr < FDR.threshold,]
  }

  ## Report done time
  if (verbosity >= 1) { print (paste(date(), "Multiple t-test done", sep= " - ")) }

  ##  plot(P.value.normal.approx,P.value,panel.first=grid(col='#0000ff'),log="xy")
  return (result)
}

```

### Prep for LDA

```
CVSNMpheno <- phenoData(LDR_canVSnorm)
colnames(CVSNMpheno)
sample.labels <- as.vector(CVSNMpheno@data$characteristics_ch1.4)
sample.colors <- plot.col
group.descriptions <- rbind(sample.labels, sample.colors[vv])
group.descriptions <- t(group.descriptions)
group.descriptions <- as.data.frame(group.descriptions)
group.labels <- group.descriptions$sample.labels
names(group.labels) <- row.names(group.descriptions)
group.colors <- group.descriptions$V2
names(group.colors) <- row.names(group.descriptions)
group.legend <- paste(sep ="", group.labels, row.names(group.descriptions))

## Variance per gene

var.per.gene <- apply(canVSnorm, 1, var)
sd.per.gene <- apply(canVSnorm, 1, sd)

```

```

par(mar = c(2, 2, 2, 2))
par(mfrow=c(1,2))

hist(var.per.gene, breaks=100, col="#BBFFDD", main="Variance per gene", xlab="Variance", ylab="Number of genes")
hist(sd.per.gene, breaks=100, col="#BBFFDD", main="Standard dev. per gene", xlab="Standard deviation", ylab="Number of genes")

```

```

## Select 20 top-ranking gens for training

genes.by.decr.var <- sort(var.per.gene,decreasing=TRUE)
top.nb <- 20
genes.selected.by.var <- names(genes.by.decr.var[1:top.nb])

## Rank by cross-sample variance

gene.ranks <- data.frame(var=var.per.gene)
gene.ranks$var.rank <- rank(-gene.ranks$var, ties.method='random')
head(gene.ranks, n=10)
```

```
g1 <- 236
g2 <- 1213

maxvar.g1 <- names(genes.by.decr.var[1])
maxvar.g2 <- names(genes.by.decr.var[2])

x <- as.vector(as.matrix(canVSnorm[maxvar.g1,]))
y <- as.vector(as.matrix(canVSnorm[maxvar.g2,]))

plot(x,y,
      col=sample.colors,
      type='n',
      panel.first=grid(col='black'), 
      main="2 genes with the highest variance", 
      xlab=paste('gene', maxvar.g1), 
      ylab=paste('gene', maxvar.g2))

text(x, y,labels=sample.labels,col=sample.colors,pch=0.5)
legend('topright',col=group.colors, 
         legend=group.legend,pch=1,cex=0.6,bg='white',bty='o')

var.per.gene <- apply(canVSnorm, 1, var)
sd.per.gene <- apply(canVSnorm, 1, sd)


par(mfrow=c(1,2))
boxplot(x ~ sample.labels, las=1,
        horizontal=TRUE,
        main=maxvar.g1, 
        col="#BBBBBB")
data.frame(group.colors)

boxplot(y ~ sample.labels, las=1,
        horizontal=TRUE,
        main=maxvar.g2, col="#BBBBBB")
par(mfrow=c(1,1))

```
### Load Welch's
```


source(file.path(dir.util, "util_student_test_multi.R"))
group.of.interest <- "histology: colon cancer"
one.vs.others <- sample.labels
## one.vs.others[sample.labels != group.of.interest] <- "histology: normal colon mucosa"


welch.one.vs.others <- t.test.multi(canVSnorm, 
                                    one.vs.others,
                                    volcano.plot = FALSE)

kable(head(welch.one.vs.others), caption = "Head of the Welch result table. Each row corrresponds to one probeset (gene), each column to one statistics used for the Welch test.")

test.name <- paste(group.of.interest, '.vs.others.sig', sep='')
gene.ranks[,test.name] <- welch.one.vs.others$sig
gene.ranks[,paste(test.name, ".rank", sep="")] <- 
    rank(-welch.one.vs.others$sig, ties.method='random')

group.of.interest1 <- as.data.frame(group.of.interest)
ample.labels1 <- as.data.frame(sample.labels)
one.vs.others2 <- as.vector(one.vs.others)

## Apply the Welch test for the 3 other majority groups
## Not necessary for the cancer vs normal dataset, as there are only two groups that we are looking at.

for (group.of.interest1 in c("histology: colon cancer", "histology: normal colon mucosa")) {
    print(paste("Selecting differentially expressed genes for", group.of.interest1, "versus others"))
    one.vs.others2 <- sample.labels
    one.vs.others2[sample.labels != group.of.interest1] <- "histology: normal colon mucosa"
    
    welch.one.vs.others <- t.test.multi(canVSnorm,                           one.vs.others2,
                         volcano.plot = FALSE)

    test.name <- paste(group.of.interest1, '.vs.others.sig', sep='')
    gene.ranks[,test.name] <- welch.one.vs.others$sig
    gene.ranks[,paste(test.name, ".rank", sep="")] <- 
      rank(-welch.one.vs.others$sig, ties.method='random')
}


write.table(gene.ranks, file=file.path(dir.results, 'gene_ranks.tab'), sep='\t', quote=F, col.names=NA)

```

### ANOVA ordering

```

g <- 1234 ## random gene

g.expr <- unlist(canVSnorm[g,])

g.for.anova <- data.frame("expr"=g.expr, "group"=sample.labels)

g.aov.result <- aov(formula = expr ~ group, data = g.for.anova)


g.anova.result <- anova(lm(formula = expr ~ group, data = g.for.anova))



pval <- as.numeric(unlist(g.anova.result)["Pr(>F)1"])


eval <- pval * nrow(canVSnorm)



g.anova.summary <- data.frame("g"=g, 
                     "name"=row.names(canVSnorm[g,]),
                     "pval"=pval,
                     "eval"=eval,
                     "sig"=-log(eval, base=10))

```

### PCA
```

expr.prcomp <- prcomp(t(canVSnorm))

attributes(expr.prcomp) 
names(expr.prcomp) 

sd.per.pc <- expr.prcomp$sdev

var.per.pc <- sd.per.pc^2

sd.per.pc.percent <- sd.per.pc/sum(sd.per.pc)

var.per.pc.percent <- var.per.pc/sum(var.per.pc)

```

### Linear Discriminant Analysis (LDA)

```
library(MASS)

one.vs.others.lda.allvars <- lda(t(canVSnorm), one.vs.others, CV=FALSE)

## Incorrect usage - as we used a training set to train on itself.

predict.lda.allvars <- predict(object = one.vs.others.lda.allvars, newdata = t(canVSnorm))
hits <- sum(one.vs.others == predict.lda.allvars$class)
errors <- sum(one.vs.others != predict.lda.allvars$class)
total <- hits + errors
(hit.rate.internal <- hits / total)
(error.rate.internal <- errors / total)
```

### Leave-one-out. Each element of training set will be left out to evaluate a classifier trained with all other elements.

```

one.vs.others.lda.allvars.loo <- lda(t(canVSnorm),one.vs.others,CV=TRUE)

(hit.rate.loo <- sum(one.vs.others == one.vs.others.lda.allvars.loo$class) / total)
(error.rate.loo <- 1 - hit.rate.loo) ## classification rate falls with extra testing
```


### Random expectation for the hit rate. Estimator
```

random.hit.rates <- vector()
for (rep in 1:10000) {
  random.hit.rates <- append(random.hit.rates, sum(one.vs.others == sample(one.vs.others)) / total)
}
(random.hit.rates.mean <- mean(random.hit.rates))

prior <- as.vector(table(one.vs.others))/length(one.vs.others)
(hit.rate.expect <- sum(prior^2))

hist(random.hit.rates, breaks=(0:total)/total, col="lightgrey", 
     freq=TRUE,
     main="Hit rate analysis",
     xlab="Hit rate",
     ylab="Frequency")
arrows(x0=hit.rate.loo, y0 = 1000, x1=hit.rate.loo, y1=100, 
       col="darkgreen", lwd=2, code=2, , length=0.2, angle=20)
arrows(x0=hit.rate.internal, y0 = 1000, x1=hit.rate.internal, y1=100, 
       col="red", lwd=2, code=2, , length=0.2, angle=20)
arrows(x0=hit.rate.expect, y0 = 1000, x1=hit.rate.expect, y1=100, 
       col="darkblue", lwd=2, code=2, , length=0.2, angle=20)

legend("topleft", legend=c("random", "random expectation", "LOO", "internal"), 
       col=c("grey", "darkblue", "darkgreen", "red"), 
       lwd=c(4,2,2,2))

```

### Selecting a subset of the variables (feature selection)
```

welch.Bo.vs.others <- t.test.multi(canVSnorm, 
                                   one.vs.others,
                                   volcano.plot = FALSE)

welch.Bo.vs.others.sorted <- welch.Bo.vs.others[order(welch.Bo.vs.others$sig, decreasing=TRUE),]

sorted.names <- rownames(welch.Bo.vs.others.sorted)

welch.Bo.vs.others[sorted.names[0:20], c("E.value","sig")]

g1 <- sorted.names[1]
g2 <- sorted.names[2]
x <- as.vector(as.matrix(canVSnorm[g1,]))
y <- as.vector(as.matrix(canVSnorm[g2,]))

top.variables <- 20
selected.genes <- sorted.names[1:top.variables]

one.vs.others.lda.classifier <- lda(t(canVSnorm[selected.genes,]),one.vs.others,CV=FALSE) 

one.vs.others.lda.loo <- lda(t(canVSnorm[selected.genes,]),one.vs.others,CV=TRUE) 

one.vs.others.loo.predicted.class <- as.vector(one.vs.others.lda.loo$class)

(one.vs.others.lda.loo.xtab <- table(one.vs.others, one.vs.others.loo.predicted.class))

(one.vs.others.lda.loo.xtab2 <- table(sample.labels, one.vs.others.loo.predicted.class))


(nb.hits <- sum(na.omit(hits)))

(nb.pred <- length(na.omit(hits)))

(hit.rate <- nb.hits / nb.pred )

```

### Multi-group classification

```

multigroup.lda.loo <- lda(t(canVSnorm[selected.genes,]),sample.labels,CV=TRUE) 

multigroup.loo.predicted.class <- as.vector(multigroup.lda.loo$class)

table(multigroup.loo.predicted.class)

(multigroup.lda.loo.xtab <- table(sample.labels, multigroup.loo.predicted.class))

library(lattice)
levelplot(multigroup.lda.loo.xtab)

hits <- sample.labels == multigroup.loo.predicted.class

errors <- sample.labels != multigroup.loo.predicted.class

(nb.hits <- sum(na.omit(hits)))

(nb.pred <- length(na.omit(hits)))

(hit.rate <- nb.hits / nb.pred )


## Random expectation of a hit rate

n <- length(sample.labels)
n.goi<- sum(sample.labels == group.of.interest1)
n.others <- sum(sample.labels != group.of.interest1)
prior <- c("goi"=n.goi/n,
           "others" = n.others/n)

exp.hits <- c("goi"=(n.goi*prior["goi"]),
              "other"=(n.others*prior["others"]))

print(exp.hits)

exp.hit.rate <- sum(prior^2)

print(exp.hit.rate)

multi.samples.per.class <- unlist(table(sample.labels))

multi.prior <- (multi.samples.per.class)/sum(multi.samples.per.class)
(multi.expect.hit.rate <- sum(multi.prior^2))

## Training a classifier with permuted labels

sample.labels.perm <- as.vector(sample(sample.labels))

table(sample.labels, sample.labels.perm)

lda.loo.labels.perm <- lda(t(canVSnorm[selected.genes,]),sample.labels.perm,CV=TRUE) 

loo.predicted.class.labels.perm <- as.vector(lda.loo.labels.perm$class)
lda.loo.labels.perm.xtab <- table(sample.labels.perm, loo.predicted.class.labels.perm)

hits.label.perm <- sample.labels.perm == loo.predicted.class.labels.perm
(nb.hits.label.perm <- sum(na.omit(hits.label.perm)))

(nb.pred.label.perm <- length(na.omit(hits.label.perm)))

(hit.rate.label.perm <- nb.hits.label.perm / nb.pred.label.perm )


## Label permutation test for two-group classification

sample.labels.perm.2gr <- as.vector(sample(one.vs.others))
table(one.vs.others, sample.labels.perm.2gr)

permuted.equal <- sum(diag(table(one.vs.others, sample.labels.perm.2gr)))
(permuted.equal.rate <- permuted.equal/length(one.vs.others))

lda.loo.labels.perm.2gr <- lda(t(canVSnorm[selected.genes,]),sample.labels.perm.2gr,CV=TRUE) 

loo.predicted.class.labels.perm.2gr <- as.vector(lda.loo.labels.perm.2gr$class)
lda.loo.labels.perm.2gr.xtab <- table(sample.labels.perm.2gr, loo.predicted.class.labels.perm.2gr)
print(lda.loo.labels.perm.2gr.xtab)


(hit.rate.label.perm.2gr <- sum(diag(lda.loo.labels.perm.2gr.xtab)) / sum(lda.loo.labels.perm.2gr.xtab))
```

### Quadratic discriminant analys (QDA)

```

(n.variables <- length(selected.genes))

one.vs.others.qda.loo <- qda(t(canVSnorm[selected.genes,]),one.vs.others,CV=TRUE)

table(one.vs.others, one.vs.others.qda.loo$class)

(one.vs.others.qda.loo.hit.rate <- sum(one.vs.others == one.vs.others.qda.loo$class)/n)

label.perm.one.vs.others.qda.loo <- qda(t(canVSnorm[selected.genes,]),
                                        sample(one.vs.others, replace=FALSE),  CV=TRUE) 
  

table(one.vs.others, label.perm.one.vs.others.qda.loo$class)

(label.perm.one.vs.others.qda.loo.hit.rate <- 
   sum(one.vs.others == label.perm.one.vs.others.qda.loo$class)/n)


one.vs.others.lda.loo <- lda(t(canVSnorm[selected.genes,]),one.vs.others,CV=TRUE) 

(one.vs.others.lda.loo.hit.rate <- sum(one.vs.others == one.vs.others.lda.loo$class)/n)

```

### Connect to Python

Done via 'reticulate'. Download anaconda to local computer, and in the R terminal create:
```
conda create -n py3.8 python=3.8 scikit-learn pandas numpy matplotlib
conda env list
```
Back in console:
```
library(reticulate)
reticulate::conda_list()
use_condaenv("PATH*", required = TRUE)
```
*PATH = local computer path that matches the location of virtual enivronment created on the conda.
```
py_config()

repl_python() >>> Interactive interface

```
```
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn
import GEOparse

gse = GEOparse.get_GEO(geo="GSE68468", destdir="./")

for gsm_name, gsm in gse.gsms.items():
    print("Name: ", gsm_name)
    print("Metadata:",)
    for key, value in gsm.metadata.items():
        print(" - %s : %s" % (key, ", ".join(value)))
    print ("Table data:",)
    print (gsm.table.head())
    break
  
for gpl_name, gpl in gse.gpls.items():
    print("Name: ", gpl_name)
    print("Metadata:",)
    for key, value in gpl.metadata.items():
        print(" - %s : %s" % (key, ", ".join(value)))
    print("Table data:",)
    print(gpl.table.head())
    break
  
  

df = pd.read_excel(r'C:/Users/greta/Documents/LDR/data.xlsx')
vv = pd.read_excel(r'C:/Users/greta/Documents/LDR/types.xlsx') 
vv = vv.iloc[:,1:]
vv.head()
 
## SVM classifier

array = df.values
array1 = vv.values

X = array[:,1:241]
Y = array[:,241]

from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression

model = LogisticRegression()
rfe = RFE(model, 3)
fit = rfe.fit(X, Y)

### To be continued
```

### Connection to Spark and MLib

```
install.packages("spkarlyr", repos = "http://cran.us.r-project.org")
library(sparklyr)
spark_available_versions() ## check for the latest versions available at the time.
spark_install("3.2")
Sys.setenv(JAVA_HOME="")
system("java -version") ## check the java version installed.

```

### MLib Decision trees and Random Forest
coloncancer <- gse[, filter1]

table(coloncancer@phenoData$characteristics_ch1.4)

filter2 <- colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal colon mucosa"]

if 
(filter2 == colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal colon mucosa"]){
  print("TRUE")
}

normalcolon <- gse[, filter2]

vs1 <- colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer" | gse@phenoData@data$"characteristics_ch1.4"=="histology: normal colon mucosa"]

if (vs1 == colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer" | gse@phenoData@data$"characteristics_ch1.4"=="histology: normal colon mucosa"])
{
  print("TRUE")
}

cancerVSnormal <- gse[, vs1]

table(cancerVSnormal@phenoData$characteristics_ch1.4)
```

### Quality control

```
plotDensities(cancerVSnormal, legend = FALSE, main = "Colon cancer versus normal colon density plot")

#aarrayQualityMetrics(cancerVSnormal, outdit = tempdir(), do.logtransform = FALSE)

oligo::boxplot(cancerVSnormal,  main = "Colon cancer versus normal colon boxplot")
```

### Quantile normalise

```
QcancerVSnormal <- cancerVSnormal
exprs(QcancerVSnormal) <- normalizeBetweenArrays(exprs(QcancerVSnormal))
oligo::boxplot(QcancerVSnormal)
```

```
par(mfrow=c(2,1))
oligo::boxplot(cancerVSnormal, main = "Before normalisation")
oligo::boxplot(QcancerVSnormal, main = "After normalisation")
```

```
par(mfrow=c(1,2))
plotDensities(cancerVSnormal, legend=FALSE, "Density plot before normalisation", main = "Before normalisation")
plotDensities(QcancerVSnormal, legend=FALSE, "Density plot after normalisation", main = "After normalisation")
```

### Dimensionality reduction

```
EX <- Biobase::exprs(QcancerVSnormal)
PCA1 <- prcomp(t(EX), scale = FALSE)
percentVar <- round(100*PCA1$sdev^2/sum(PCA1$sdev^2),1)
sd_ratio <- sqrt(percentVar[2]/percentVar[1])
dataGG <- data.frame(PC1 = PCA1$x[,1], PC2 = PCA1$x[,2],
            Histology = Biobase::pData(QcancerVSnormal)$characteristics_ch1.4)
ggplot(dataGG, aes(PC1, PC2)) + 
          geom_point(aes(colour = Histology)) + 
          ggtitle("PCA plot of normalised data") +
          xlab(paste0("PC1, VarExp:", percentVar[1], "%")) +
          ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          coord_fixed(ratio = sd_ratio) +
          scale_color_manual(values = c("darkorange2", "dodgerblue4"))

```

### Inspection of batch effects

```
pct_threshold <- 0.6
colnames(pData(QcancerVSnormal))
batch.factors <-c("source_name_ch1", "characteristics_ch1.4")
pvca0bj <- pvcaBatchAssess(QcancerVSnormal, batch.factors, pct_threshold)

sds <- apply(exprs(QcancerVSnormal), 1, sd)
sds0 <- sort(sds)
plot(1:length(sds0), sds0, main = "Distribution of variability for all genes", 
      sub="Vertical lines represent 90% and 95% percentiles",
      xlab="Gene index (from least to most variable)",
      ylab="Standard deviation")
 abline(v=length(sds)*c(0.9,0.95))
 
```

### More filtering and Annotation

```
AQcancervsnormal <- QcancerVSnormal
annotation(AQcancervsnormal) <- "hgu133a.db"
filtered <- nsFilter(AQcancervsnormal, 
                      require.entrez = TRUE, 
                      remove.dupEntrez = TRUE,
                      var.filter = TRUE,
                      var.func = IQR,
                      var.cutoff = 0.25,
                      filterByQuantile = TRUE,
                      feature.exclude = "^AFFX")
print(filtered$filter.log)
filtered_cancervsnorm <- filtered$eset
nrow(filtered_cancervsnorm)
head(filtered_cancervsnorm@featureData$"Gene Symbol")

```

### Linear modelling.

```
metadata <- pData(filtered_cancervsnorm)
design <- model.matrix(~0+metadata$characteristics_ch1.4)
colnames(design) <- c("Cancer", "Normal")
fit <- lmFit(filtered_cancervsnorm, design)
contrasts <- makeContrasts(Cancer - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
class(fit2)
results <- decideTests(fit2)
summary(results)
topTab <- topTable(fit2, number = nrow(fit2), adjust = "fdr")
cancerVSnormTopTable <- topTab[, ! names(topTab) %in% c("GB_ACC", "SPOT_ID", "Species.Scientific.Name", "Annotation.Date", "Sequence.Source", "Target.Description", "Representative.Public.ID", "Refeq.Transcript.ID", "Gene.Ontology.Biological.Process", "Gene.Ontology.Cellular.Component", "Gene.Ontology.Molecular.Function")]
cancerVSnormTopTable <- cancerVSnormTopTable[, ! names(cancerVSnormTopTable) %in% c("Sequence.Type", "Gene.Title", "RefSeq.Transcript.ID")]

results1 <- decideTests(fit2, method = "separate", adjust.method = "fdr", p.value = 0.1, lfc = 1)

sum.res.rows <- apply(abs(results1), 1, sum)

res.selected <- results1[sum.res.rows!=0,]
```

```

colnames(cancerVSnormTopTable)
names(cancerVSnormTopTable)[names(cancerVSnormTopTable) == "ENTREZ_GENE_ID"] <- "ENTREZID"

listOfTables <- list(ColonCancerVSNormalColon = cancerVSnormTopTable)
listofSelected <- list()

for (i in 1:length(listOfTables)) 
  {
  topTable <- listOfTables[[i]]
  whichGenes <- topTable["adj.P.Val"] < 0.15
  selectedIDs <- rownames(topTable)[whichGenes]
  EntrezIDs <- AnnotationDbi::select(hgu133a.db, selectedIDs, c("ENTREZID"))
  EntrezIDs <- EntrezIDs$ENTREZID
  listofSelected[[i]] <- EntrezIDs
   names(listofSelected)[i] <- names(listOfTables)[i]
}
sapply(listofSelected, length)
listOfSelected <- listofSelected


```

### Pathway Enrichment
#### Method 1

```

mapped_genes2GO <- mappedkeys(hgu133aGO)
mapped_genes2KEGG <- mappedkeys(hgu133aPATH)
mapped_genes <- union(mapped_genes2GO, mapped_genes2KEGG)
library(ReactomePA)

listOfData <- listOfSelected[1]
comparisonsNames <- names(listOfData)
universe <- mapped_genes

#for (i in 1:length(listOfData)){
#genesIn <- listOfData[[i]]
#comparison <- comparisonsNames[i]
#enrich.result <- enrichPathway(gene = genesIn,
                                #pvalueCutoff = 0.05,
                               # qvalueCutoff = 0.9,
                                # readable = T,
                                # pAdjustMethod = "BH",
                                # organism = "human",
                                # universe = universe)
  
#  cat("##################################")
 # cat("\nComparison: ", comparison,"\n")
 # print(head(enrich.result))
 
  # if (length(rownames(enrich.result@result)) != 0) {
 #  write.csv(as.data.frame(enrich.result), 
# file =paste0("./results/","ReactomePA.Results.",comparison,".csv"), 
            # row.names = FALSE)
   
#pdf(file=paste0("./results/","ReactomePABarplot.",comparison,".pdf"))
# print(barplot(enrich.result, showCategory = 15, font.size = 4, 
# title = paste0("Reactome Pathway Analysis for ", comparison,". Barplot")))
 # dev.off()
  
#  pdf(file = paste0("./results/","ReactomePAcnetplot.",comparison,".pdf"))
 #    print(cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
    #      vertex.label.cex = 0.75))
  # dev.off()
  # }
# }
```


#### Method 2

```
back <- subset(cancerVSnormTopTable, adj.P.Val < 0.1)$ID
backidx <- genefilter::genefinder(filtered_cancervsnorm, as.character(back), method = "manhattan", scale = "none")
backidx <- sapply(backidx, function(x)x$indices)
backs <- featureNames(filtered_cancervsnorm )[backidx]
backs <- setdiff(backs, back)
intersect(backs, back)

IDs <- rownames(cancerVSnormTopTable)
un <- IDs %in% c(backs, back)
sl <- IDs %in% back

genes <- sl[un]
genes <- factor(as.integer(sl[un]))
names(genes) <- IDs[un]

entrez1 <- mapIds(hgu133a.db, 
                 keys = rownames(results),
                 keytype = "PROBEID",
                 colum = "ENTREZID")

reactome <- enrichPathway(gene = entrez1[back], 
                          universe = entrez1[c(back, 
                                              backs)],
                          organism = "human",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.9, 
                          readable = TRUE)
reactome @result$Description <- paste0(str_sub(
    reactome @result$Description, 1, 20),
    "...")



top_Go_data <- new("topGOdata", ontology = "BP", allGenes = genes, nodeSize = 10, annot = annFUN.db, affyLib = "hgu133a.db")

topGoElim <- runTest(top_Go_data, algorithm = "elim", statistic = "Fisher")

topGoClassic <- runTest(top_Go_data, algorithm = "classic", statistic = "Fisher")


resTopGo <- GenTable(top_Go_data, Fisher.elim = topGoElim, Fisher.classic = topGoClassic,
                     orderBy = "Fisher.elim", topNodes = 100)

genesTopGo <- printGenes(top_Go_data, whichTerms = resTopGo$GO.ID, chip = "hgu133a.db", geneCutOff = 1000)

resTopGo$sig_genes <- sapply(genesTopGo, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"], ";"),
        collapse = "")
})

```

#### Method 3

```
##KEGG

entrez <- entrez1
enrich_kegg <- kegga(fit2, geneid = entrez, species = "Hs")
topKEGG(enrich_kegg)

##GO
enrich_go <- goana(fit2, geneid = entrez, species = "Hs")
topGO(enrich_go, ontology = "BP")

## Filtering for genes involved in the biological processes
bpGO <- cancerVSnormTopTable
x <- org.Hs.egGO2ALLEGS
Rkeys(x) <- "GO:0008150"
ID <- mappedLkeys(x)
i <- bpGO$ENTREZID %in% ID
bpGO <- bpGO[i,] 

## Filtering for genes involved in xenobiotic metabolic process
xmpGO <- cancerVSnormTopTable
x1 <- org.Hs.egGO2ALLEGS
Rkeys(x1) <- "GO:0006805"
ID1 <- mappedLkeys(x1)
i1 <- xmpGO$ENTREZID %in% ID1
i1_d <- as.numeric(i1)
xmpGO <- xmpGO[i1,]


```

#### Method 4

```
library(clusterProfiler)
## CC - cellular  component, MF - molecular function, BP - biological process.

geneList <- cancerVSnormTopTable[,4]
names(geneList) <- as.character(cancerVSnormTopTable[,3])
geneList <- sort(geneList, decreasing = TRUE)

ggoCC <- groupGO(gene = entrez,
               OrgDb = org.Hs.eg.db,
               ont = "CC",
               level = 3,
               readable = TRUE)

ggoMF <- groupGO(gene = entrez,
               OrgDb = org.Hs.eg.db,
               ont = "MF",
               level = 3,
               readable = TRUE)

ggoBP <- groupGO(gene = entrez,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               level = 3,
               readable = TRUE)
```


### Over-representation analysis

```
##egoCC <- enrichGO(gene = entrez,
                ##  OrgDb = org.Hs.eg.db,
                 ## universe = names(geneList), 
                 ## ont = "CC",
                 ## pAdjustMethod = "BH",
                 ## pvalueCutoff = 0.01,
                 ## qvalueCutoff = 0.05,
                 ## readable = TRUE)

##egoMF <- enrichGO(gene = entrez,
                 ## OrgDb = org.Hs.eg.db,
                 ## ont = "MF",
                 ## pAdjustMethod = "BH",
                 ## pvalueCutoff = 0.01,
                 ## qvalueCutoff = 0.05,
                 ## readable = TRUE)

##egoBP <- enrichGO(gene = entrez,
                ##  OrgDb = org.Hs.eg.db,
                 ## ont = "BP",
                 ## pAdjustMethod = "BH",
                 ## pvalueCutoff = 0.01,
                 ## qvalueCutoff = 0.05,
                 ## readable = TRUE)

##gene.df <- bitr(entrez, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)

##egoCC2 <- enrichGO(gene = gene.df$ENSEMBL,
                  ## OrgDb = org.Hs.eg.db,
                 ##  keyType = 'ENSEMBL',
                 ##  ont = "CC",
                 ##  pAdjustMethod = "BH",
                 ##  pvalueCutoff = 0.01,
                  ## qvalueCutoff = 0.05)

##egoMF2 <- enrichGO(gene = gene.df$ENSEMBL,
                 ##  OrgDb = org.Hs.eg.db,
                 ##  keyType = 'ENSEMBL',
                ##   ont = "MF",
                 ##  pAdjustMethod = "BH",
                 ##  pvalueCutoff = 0.01,
                 ##  qvalueCutoff = 0.05)

##egoBP2 <- enrichGO(gene = gene.df$ENSEMBL,
                 ##  OrgDb = org.Hs.eg.db,
                  ## keyType = 'ENSEMBL',
                  ## ont = "BP",
                  ## pAdjustMethod = "BH",
                  ## pvalueCutoff = 0.01,
                  ## qvalueCutoff = 0.05)

##GOaCC <- gseGO(geneList = geneList,
              ## OrgDb = org.Hs.eg.db,
              ## ont = "CC",
             ##  minGSSize = 100,
             ##  maxGSSize = 500,
             ##  pvalueCutoff = 0.05,
             ##  verbose = FALSE)
```

```{r}
##goplot(egoCC)
##cnetplot(egoCC, showCategory = 10)
```

### Differential expression with t testing

```
packages <- c("knitr")
for (pkg in packages) {
  if (!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

pkg <- "qvalue"
if (!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
  source("http://bioconductor.org/biocLite.R")
  biocLite();
  biocLite(pkg)
}
library(pkg, character.only = TRUE)

url.stats4bioinfo <- "http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics"
source(file.path(url.stats4bioinfo, 'R-files/config.R'))
dir.base <- 'http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics'
dir.base.ASG1 <- 'http://pedagogix-tagc.univ-mrs.fr/courses/ASG1'
source(file.path(dir.base, 'R-files', 'config.R'))
setwd(dir.results)
print(paste("Result directory", dir.results))
dir.base <- 'http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics'
dir.base.ASG1 <- 'http://pedagogix-tagc.univ-mrs.fr/courses/ASG1'
source(file.path(dir.base, 'R-files', 'config.R'))
setwd(dir.results)
print(paste("Result directory", dir.results))
source(file.path(dir.util, "util_student_test_multi.R"))
group.of.interest <- "histology: colon cancer"
```

### New filter - Different cut off

```
LDRfilter <- nsFilter(AQcancervsnormal, 
                      require.entrez = TRUE, 
                      remove.dupEntrez = TRUE,
                      var.filter = TRUE,
                      var.func = IQR,
                      var.cutoff = 0.75,
                      filterByQuantile = TRUE,
                      feature.exclude = "^AFFX")
print(LDRfilter$filter.log)
LDR_canVSnorm <- LDRfilter$eset
nrow(LDR_canVSnorm)
head(LDR_canVSnorm@featureData$"Gene Symbol")


table(LDR_canVSnorm$characteristics_ch1.4)
      
canVSnormal <- exprs(LDR_canVSnorm)

canVSnorm <- as.data.frame(canVSnormal)

vv <- LDR_canVSnorm$characteristics_ch1.4

colnames(canVSnorm) <- paste(colnames(canVSnorm), vv, sep = "_")
```

### Single Student Test. Random gene to test hypothesis

```

n.Cancer <- 186
n.Normal <- 55

g <- sample(1:nrow(canVSnorm), 1)

g.profile <-as.vector(as.matrix(canVSnorm[g,]))
plot.col <- c('histology: colon cancer' = '#4444BB', 'histology: normal colon mucosa' = '#FFFF88')
```

```
par(mar = c(5.1, 4.1, 4.1, 11))
barplot(g.profile, col=plot.col[vv], main = "One random gene expression in the data", ylab = "Relative expression", xlab = "Individual samples n=65" , xlim = c(1, 65))
legend('topright', c("histology: colon cancer", "histology: normal colon mucosa"), col=plot.col[c("histology: colon cancer", "histology: normal colon mucosa")], pch=15, bty="o", bg='white', cex = 0.65)
```

```
sample.Cancer <- g.profile[vv == "histology: colon cancer"]
sample.Normal <- g.profile[vv == "histology: normal colon mucosa"]

mean.est.Cancer <- mean(sample.Cancer)
mean.est.Normal <- mean(sample.Normal)

sample.sd.Cancer <- sd(sample.Cancer) * sqrt((n.Cancer-1)/n.Cancer)
sample.sd.Normal <- sd(sample.Normal) * sqrt((n.Normal-1)/n.Normal)

sd.est.Cancer <- sd(sample.Cancer)
sd.est.Normal <- sd(sample.Normal)


sd.err.est.Cancer <- sd(sample.Cancer) / sqrt(n.Cancer)
sd.err.est.Normal <- sd(sample.Normal) / sqrt(n.Normal)

diff.sd.est <- sqrt((n.Cancer * sample.sd.Cancer^2 + n.Normal * sample.sd.Normal^2) * (1/n.Cancer + 1/n.Normal) / (n.Cancer+n.Normal-2))

d <- abs(mean.est.Cancer - mean.est.Normal)
t.obs.Student <- d/diff.sd.est
P.val.Student <- 2*pt(q=t.obs.Student,
                      df=n.Cancer+n.Normal-2,
                      lower.tail = F)
t.student <- t.test(sample.Cancer, sample.Normal, var.equal = TRUE)

print(t.student)

### Welch's

t.welch <- t.test(sample.Cancer, sample.Normal, var.equal = FALSE)
 
print(t.welch)

t.statistics <- vector()
P.values <- vector()

for (g in 1:nrow(canVSnorm)) {
  print(paste("Random gene", g))
  g.profile <- as.vector(canVSnorm[g,])
  sample.Cancer <- g.profile[vv == "histology: colon cancer"]
  sample.Normal <- g.profile[vv == "histology: normal colon mucosa"]
  t <- t.test(sample.Cancer, sample.Normal)
  t.statistics <- append(t.statistics, t$statistic)
  P.values <- append(P.values, t$p.value)
}
print(P.values)

### Multiple variance test

canvsnorm.t.result <- t.test.multi(canVSnorm, vv)
dim(canvsnorm.t.result)
names(canvsnorm.t.result)
sum(canvsnorm.t.result$E.value <=1)
canvsnorm.E <- canVSnorm[canvsnorm.t.result$E.value <=1,]


```

#### T test multivariance function
```
t.test.multi  function (x, ## a matrix or data.frame
                          cl, ##  vector of classes
                          P.threshold=NA, ## select objects below a given threshold of P.value
                          E.threshold=NA, ## select objects below a given threshold of E.value
                          FDR.threshold=NA, ## select objects below a given threshold of FDR
			  robust.est=F, ## use robust estimators for central tendency and dispersion
			  verbosity=1, ## Level of verbosity
                          volcano.plot = T, ## Draw a volcano plot
                          alternative = "two.sided", ## Supported: "two.sided", "less", "greater"
                          ... ## additional parameters are passed to the function volcano.plot()
                          ) {
  ## Report starting time
  if (verbosity >= 1) { print (paste(date(), "Multiple t-test started", sep= " - ")) }
  
  ## Dimensions of the data set
  n <- nrow(x)
  p <- ncol(x)

  ## check the dimensions of the class vector
  if (length(cl) != p) {
    stop (paste ('The number of columns of the dataset (',
                 p,
                 ') should be the same as the length of the class vector (',
                 length(cl),
                 ')'))
  }

  classes <- unique(cl)
  if (robust.est) {
    means.1 <- apply(x[,cl==classes[1]],1,median,na.rm=T)
    means.2 <- apply(x[,cl==classes[2]],1,median,na.rm=T)
    iqr.1 <- apply (x[,cl==classes[1]], 1, IQR, na.rm=T)
    iqr.2 <- apply (x[,cl==classes[2]], 1, IQR, na.rm=T)
    sd.est.1 <- iqr.1/(qnorm(0.75) - qnorm(0.25))
    sd.est.2 <- iqr.2/(qnorm(0.75) - qnorm(0.25))
    var.est.1 <- sd.est.1^2
    var.est.2 <- sd.est.2^2
  } else {
    means.1 <- apply(x[,cl==classes[1]],1,mean,na.rm=T)
    means.2 <- apply(x[,cl==classes[2]],1,mean,na.rm=T)
    var.est.1 <- apply(x[,cl==classes[1]],1,var,na.rm=T)
    var.est.2 <- apply(x[,cl==classes[2]],1,var,na.rm=T)
    sd.est.1 <- sqrt(var.est.1)
    sd.est.2 <- sqrt(var.est.2)
  }


  ## 2012-05-02: restored means.diff = means.1 - means.2 for the sake
  ## of consistency with R funciton t.test()
  means.diff <- means.1 - means.2
  ##  means.diff <- means.2 - means.1
  ##  hist(means.diff,breaks=50)
  
  n.1 <- sum(cl == classes[1])
  n.2 <- sum(cl == classes[2])

  ## Calculate observed t value
  st.err.diff <- sqrt(var.est.1/n.1 + var.est.2/n.2)
  t.obs <- means.diff/st.err.diff
  ## Calculate degrees of freedom with Welch's formula  
  df.welch <- (var.est.1/n.1 + var.est.2/n.2)^2 / ((var.est.1/n.1)^2/(n.1-1) + (var.est.2/n.2)^2/(n.2-1))

  ## Calculate P-value and E-value
#  P.value.normal.approx <- 2*pnorm(abs(t.obs),lower.tail=F)
  if (alternative == "greater") {
    P.value <- pt(t.obs,df.welch,lower.tail=FALSE)
  } else if (alternative == "less") {
    P.value <- pt(t.obs,df.welch,lower.tail=TRUE)
  } else if (alternative == "two.sided") {
    P.value <- 2*pt(abs(t.obs),df.welch,lower.tail=F)
  } else {
    stop('Invalid alternative option for t.test.multi(). Supported: "two.sided", "less", or "greater"')
  }
  E.value <- P.value*nrow(x)
  sig <- -log(E.value)/log(10)

  multi.corr <- multitest.corrections(P.value, plots=FALSE)

  ## Collect all statistics in a data frame
  result <- data.frame(
                       means.1,
                       means.2,
                       means.diff,
                       var.est.1,
                       var.est.2,
                       sd.est.1,
                       sd.est.2,
                       st.err.diff,
                       t.obs,
                       df.welch,
#                       P.value.normal.approx,
                       P.value,
                       E.value,
		       sig)
  result$fwer <- multi.corr$multitest.table$fwer
  result$q.value <- multi.corr$multitest.table$qval.Storey
  result$fdr <- multi.corr$multitest.table$fdr
  result$rank <- multi.corr$multitest.table$rank

  if (robust.est) {
    names(result)[1:7] <- c(
                            paste("median.", classes[1], sep=""),
                            paste("median.", classes[2], sep=""),
                            "medians.diff",
                            paste("var.est.", classes[1], sep=""),
                            paste("var.est.", classes[2], sep=""),
                            paste("sd.est.", classes[1], sep=""),
                            paste("sd.est.", classes[2], sep="")
                            )
  } else {
    names(result)[1:7] <- c(
                            paste("mean.", classes[1], sep=""),
                            paste("mean.", classes[2], sep=""),
                            "means.diff",
                            paste("var.est.", classes[1], sep=""),
                            paste("var.est.", classes[2], sep=""),
                            paste("sd.est.", classes[1], sep=""),
                            paste("sd.est.", classes[2], sep="")
                            )
  }

  
  ## Filtering on P-value and E-value thresholds
  if (!is.na(P.threshold)) {
    result <- result[result$P.value < P.threshold,]
  }
  if (!is.na(E.threshold)) {
    result <- result[result$E.value < E.threshold,]
  }
  if (!is.na(FDR.threshold)) {
    result <- result[result$fdr < FDR.threshold,]
  }

  ## Report done time
  if (verbosity >= 1) { print (paste(date(), "Multiple t-test done", sep= " - ")) }

  ##  plot(P.value.normal.approx,P.value,panel.first=grid(col='#0000ff'),log="xy")
  return (result)
}

```

### Prep for LDA

```
CVSNMpheno <- phenoData(LDR_canVSnorm)
colnames(CVSNMpheno)
sample.labels <- as.vector(CVSNMpheno@data$characteristics_ch1.4)
sample.colors <- plot.col
group.descriptions <- rbind(sample.labels, sample.colors[vv])
group.descriptions <- t(group.descriptions)
group.descriptions <- as.data.frame(group.descriptions)
group.labels <- group.descriptions$sample.labels
names(group.labels) <- row.names(group.descriptions)
group.colors <- group.descriptions$V2
names(group.colors) <- row.names(group.descriptions)
group.legend <- paste(sep ="", group.labels, row.names(group.descriptions))

## Variance per gene

var.per.gene <- apply(canVSnorm, 1, var)
sd.per.gene <- apply(canVSnorm, 1, sd)

```

```

par(mar = c(2, 2, 2, 2))
par(mfrow=c(1,2))

hist(var.per.gene, breaks=100, col="#BBFFDD", main="Variance per gene", xlab="Variance", ylab="Number of genes")
hist(sd.per.gene, breaks=100, col="#BBFFDD", main="Standard dev. per gene", xlab="Standard deviation", ylab="Number of genes")

```

```

## Select 20 top-ranking gens for training

genes.by.decr.var <- sort(var.per.gene,decreasing=TRUE)
top.nb <- 20
genes.selected.by.var <- names(genes.by.decr.var[1:top.nb])

## Rank by cross-sample variance

gene.ranks <- data.frame(var=var.per.gene)
gene.ranks$var.rank <- rank(-gene.ranks$var, ties.method='random')
head(gene.ranks, n=10)
```

```
g1 <- 236
g2 <- 1213

maxvar.g1 <- names(genes.by.decr.var[1])
maxvar.g2 <- names(genes.by.decr.var[2])

x <- as.vector(as.matrix(canVSnorm[maxvar.g1,]))
y <- as.vector(as.matrix(canVSnorm[maxvar.g2,]))

plot(x,y,
      col=sample.colors,
      type='n',
      panel.first=grid(col='black'), 
      main="2 genes with the highest variance", 
      xlab=paste('gene', maxvar.g1), 
      ylab=paste('gene', maxvar.g2))

text(x, y,labels=sample.labels,col=sample.colors,pch=0.5)
legend('topright',col=group.colors, 
         legend=group.legend,pch=1,cex=0.6,bg='white',bty='o')

var.per.gene <- apply(canVSnorm, 1, var)
sd.per.gene <- apply(canVSnorm, 1, sd)


par(mfrow=c(1,2))
boxplot(x ~ sample.labels, las=1,
        horizontal=TRUE,
        main=maxvar.g1, 
        col="#BBBBBB")
data.frame(group.colors)

boxplot(y ~ sample.labels, las=1,
        horizontal=TRUE,
        main=maxvar.g2, col="#BBBBBB")
par(mfrow=c(1,1))

```
### Load Welch's
```


source(file.path(dir.util, "util_student_test_multi.R"))
group.of.interest <- "histology: colon cancer"
one.vs.others <- sample.labels
## one.vs.others[sample.labels != group.of.interest] <- "histology: normal colon mucosa"


welch.one.vs.others <- t.test.multi(canVSnorm, 
                                    one.vs.others,
                                    volcano.plot = FALSE)

kable(head(welch.one.vs.others), caption = "Head of the Welch result table. Each row corrresponds to one probeset (gene), each column to one statistics used for the Welch test.")

test.name <- paste(group.of.interest, '.vs.others.sig', sep='')
gene.ranks[,test.name] <- welch.one.vs.others$sig
gene.ranks[,paste(test.name, ".rank", sep="")] <- 
    rank(-welch.one.vs.others$sig, ties.method='random')

group.of.interest1 <- as.data.frame(group.of.interest)
ample.labels1 <- as.data.frame(sample.labels)
one.vs.others2 <- as.vector(one.vs.others)

## Apply the Welch test for the 3 other majority groups
## Not necessary for the cancer vs normal dataset, as there are only two groups that we are looking at.

for (group.of.interest1 in c("histology: colon cancer", "histology: normal colon mucosa")) {
    print(paste("Selecting differentially expressed genes for", group.of.interest1, "versus others"))
    one.vs.others2 <- sample.labels
    one.vs.others2[sample.labels != group.of.interest1] <- "histology: normal colon mucosa"
    
    welch.one.vs.others <- t.test.multi(canVSnorm,                           one.vs.others2,
                         volcano.plot = FALSE)

    test.name <- paste(group.of.interest1, '.vs.others.sig', sep='')
    gene.ranks[,test.name] <- welch.one.vs.others$sig
    gene.ranks[,paste(test.name, ".rank", sep="")] <- 
      rank(-welch.one.vs.others$sig, ties.method='random')
}


write.table(gene.ranks, file=file.path(dir.results, 'gene_ranks.tab'), sep='\t', quote=F, col.names=NA)

```

### ANOVA ordering

```

g <- 1234 ## random gene

g.expr <- unlist(canVSnorm[g,])

g.for.anova <- data.frame("expr"=g.expr, "group"=sample.labels)

g.aov.result <- aov(formula = expr ~ group, data = g.for.anova)


g.anova.result <- anova(lm(formula = expr ~ group, data = g.for.anova))



pval <- as.numeric(unlist(g.anova.result)["Pr(>F)1"])


eval <- pval * nrow(canVSnorm)



g.anova.summary <- data.frame("g"=g, 
                     "name"=row.names(canVSnorm[g,]),
                     "pval"=pval,
                     "eval"=eval,
                     "sig"=-log(eval, base=10))

```

### PCA
```

expr.prcomp <- prcomp(t(canVSnorm))

attributes(expr.prcomp) 
names(expr.prcomp) 

sd.per.pc <- expr.prcomp$sdev

var.per.pc <- sd.per.pc^2

sd.per.pc.percent <- sd.per.pc/sum(sd.per.pc)

var.per.pc.percent <- var.per.pc/sum(var.per.pc)

```

### Linear Discriminant Analysis (LDA)

```
library(MASS)

one.vs.others.lda.allvars <- lda(t(canVSnorm), one.vs.others, CV=FALSE)

## Incorrect usage - as we used a training set to train on itself.

predict.lda.allvars <- predict(object = one.vs.others.lda.allvars, newdata = t(canVSnorm))
hits <- sum(one.vs.others == predict.lda.allvars$class)
errors <- sum(one.vs.others != predict.lda.allvars$class)
total <- hits + errors
(hit.rate.internal <- hits / total)
(error.rate.internal <- errors / total)
```

### Leave-one-out. Each element of training set will be left out to evaluate a classifier trained with all other elements.

```

one.vs.others.lda.allvars.loo <- lda(t(canVSnorm),one.vs.others,CV=TRUE)

(hit.rate.loo <- sum(one.vs.others == one.vs.others.lda.allvars.loo$class) / total)
(error.rate.loo <- 1 - hit.rate.loo) ## classification rate falls with extra testing
```


### Random expectation for the hit rate. Estimator
```

random.hit.rates <- vector()
for (rep in 1:10000) {
  random.hit.rates <- append(random.hit.rates, sum(one.vs.others == sample(one.vs.others)) / total)
}
(random.hit.rates.mean <- mean(random.hit.rates))

prior <- as.vector(table(one.vs.others))/length(one.vs.others)
(hit.rate.expect <- sum(prior^2))

hist(random.hit.rates, breaks=(0:total)/total, col="lightgrey", 
     freq=TRUE,
     main="Hit rate analysis",
     xlab="Hit rate",
     ylab="Frequency")
arrows(x0=hit.rate.loo, y0 = 1000, x1=hit.rate.loo, y1=100, 
       col="darkgreen", lwd=2, code=2, , length=0.2, angle=20)
arrows(x0=hit.rate.internal, y0 = 1000, x1=hit.rate.internal, y1=100, 
       col="red", lwd=2, code=2, , length=0.2, angle=20)
arrows(x0=hit.rate.expect, y0 = 1000, x1=hit.rate.expect, y1=100, 
       col="darkblue", lwd=2, code=2, , length=0.2, angle=20)

legend("topleft", legend=c("random", "random expectation", "LOO", "internal"), 
       col=c("grey", "darkblue", "darkgreen", "red"), 
       lwd=c(4,2,2,2))

```

### Selecting a subset of the variables (feature selection)
```

welch.Bo.vs.others <- t.test.multi(canVSnorm, 
                                   one.vs.others,
                                   volcano.plot = FALSE)

welch.Bo.vs.others.sorted <- welch.Bo.vs.others[order(welch.Bo.vs.others$sig, decreasing=TRUE),]

sorted.names <- rownames(welch.Bo.vs.others.sorted)

welch.Bo.vs.others[sorted.names[0:20], c("E.value","sig")]

g1 <- sorted.names[1]
g2 <- sorted.names[2]
x <- as.vector(as.matrix(canVSnorm[g1,]))
y <- as.vector(as.matrix(canVSnorm[g2,]))

top.variables <- 20
selected.genes <- sorted.names[1:top.variables]

one.vs.others.lda.classifier <- lda(t(canVSnorm[selected.genes,]),one.vs.others,CV=FALSE) 

one.vs.others.lda.loo <- lda(t(canVSnorm[selected.genes,]),one.vs.others,CV=TRUE) 

one.vs.others.loo.predicted.class <- as.vector(one.vs.others.lda.loo$class)

(one.vs.others.lda.loo.xtab <- table(one.vs.others, one.vs.others.loo.predicted.class))

(one.vs.others.lda.loo.xtab2 <- table(sample.labels, one.vs.others.loo.predicted.class))


(nb.hits <- sum(na.omit(hits)))

(nb.pred <- length(na.omit(hits)))

(hit.rate <- nb.hits / nb.pred )

```

### Multi-group classification

```

multigroup.lda.loo <- lda(t(canVSnorm[selected.genes,]),sample.labels,CV=TRUE) 

multigroup.loo.predicted.class <- as.vector(multigroup.lda.loo$class)

table(multigroup.loo.predicted.class)

(multigroup.lda.loo.xtab <- table(sample.labels, multigroup.loo.predicted.class))

library(lattice)
levelplot(multigroup.lda.loo.xtab)

hits <- sample.labels == multigroup.loo.predicted.class

errors <- sample.labels != multigroup.loo.predicted.class

(nb.hits <- sum(na.omit(hits)))

(nb.pred <- length(na.omit(hits)))

(hit.rate <- nb.hits / nb.pred )


## Random expectation of a hit rate

n <- length(sample.labels)
n.goi<- sum(sample.labels == group.of.interest1)
n.others <- sum(sample.labels != group.of.interest1)
prior <- c("goi"=n.goi/n,
           "others" = n.others/n)

exp.hits <- c("goi"=(n.goi*prior["goi"]),
              "other"=(n.others*prior["others"]))

print(exp.hits)

exp.hit.rate <- sum(prior^2)

print(exp.hit.rate)

multi.samples.per.class <- unlist(table(sample.labels))

multi.prior <- (multi.samples.per.class)/sum(multi.samples.per.class)
(multi.expect.hit.rate <- sum(multi.prior^2))

## Training a classifier with permuted labels

sample.labels.perm <- as.vector(sample(sample.labels))

table(sample.labels, sample.labels.perm)

lda.loo.labels.perm <- lda(t(canVSnorm[selected.genes,]),sample.labels.perm,CV=TRUE) 

loo.predicted.class.labels.perm <- as.vector(lda.loo.labels.perm$class)
lda.loo.labels.perm.xtab <- table(sample.labels.perm, loo.predicted.class.labels.perm)

hits.label.perm <- sample.labels.perm == loo.predicted.class.labels.perm
(nb.hits.label.perm <- sum(na.omit(hits.label.perm)))

(nb.pred.label.perm <- length(na.omit(hits.label.perm)))

(hit.rate.label.perm <- nb.hits.label.perm / nb.pred.label.perm )


## Label permutation test for two-group classification

sample.labels.perm.2gr <- as.vector(sample(one.vs.others))
table(one.vs.others, sample.labels.perm.2gr)

permuted.equal <- sum(diag(table(one.vs.others, sample.labels.perm.2gr)))
(permuted.equal.rate <- permuted.equal/length(one.vs.others))

lda.loo.labels.perm.2gr <- lda(t(canVSnorm[selected.genes,]),sample.labels.perm.2gr,CV=TRUE) 

loo.predicted.class.labels.perm.2gr <- as.vector(lda.loo.labels.perm.2gr$class)
lda.loo.labels.perm.2gr.xtab <- table(sample.labels.perm.2gr, loo.predicted.class.labels.perm.2gr)
print(lda.loo.labels.perm.2gr.xtab)


(hit.rate.label.perm.2gr <- sum(diag(lda.loo.labels.perm.2gr.xtab)) / sum(lda.loo.labels.perm.2gr.xtab))
```

### Quadratic discriminant analys (QDA)

```

(n.variables <- length(selected.genes))

one.vs.others.qda.loo <- qda(t(canVSnorm[selected.genes,]),one.vs.others,CV=TRUE)

table(one.vs.others, one.vs.others.qda.loo$class)

(one.vs.others.qda.loo.hit.rate <- sum(one.vs.others == one.vs.others.qda.loo$class)/n)

label.perm.one.vs.others.qda.loo <- qda(t(canVSnorm[selected.genes,]),
                                        sample(one.vs.others, replace=FALSE),  CV=TRUE) 
  

table(one.vs.others, label.perm.one.vs.others.qda.loo$class)

(label.perm.one.vs.others.qda.loo.hit.rate <- 
   sum(one.vs.others == label.perm.one.vs.others.qda.loo$class)/n)


one.vs.others.lda.loo <- lda(t(canVSnorm[selected.genes,]),one.vs.others,CV=TRUE) 

(one.vs.others.lda.loo.hit.rate <- sum(one.vs.others == one.vs.others.lda.loo$class)/n)

```

### Connect to Python

Done via 'reticulate'. Download anaconda to local computer, and in the R terminal create:
```
conda create -n py3.8 python=3.8 scikit-learn pandas numpy matplotlib
conda env list
```
Back in console:
```
library(reticulate)
reticulate::conda_list()
use_condaenv("PATH*", required = TRUE)
```
*PATH = local computer path that matches the location of virtual enivronment created on the conda.
```
py_config()

repl_python() >>> Interactive interface

```
```
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn
import GEOparse

gse = GEOparse.get_GEO(geo="GSE68468", destdir="./")

for gsm_name, gsm in gse.gsms.items():
    print("Name: ", gsm_name)
    print("Metadata:",)
    for key, value in gsm.metadata.items():
        print(" - %s : %s" % (key, ", ".join(value)))
    print ("Table data:",)
    print (gsm.table.head())
    break
  
for gpl_name, gpl in gse.gpls.items():
    print("Name: ", gpl_name)
    print("Metadata:",)
    for key, value in gpl.metadata.items():
        print(" - %s : %s" % (key, ", ".join(value)))
    print("Table data:",)
    print(gpl.table.head())
    break
  
  

df = pd.read_excel(r'C:/Users/greta/Documents/LDR/data.xlsx')
vv = pd.read_excel(r'C:/Users/greta/Documents/LDR/types.xlsx') 
vv = vv.iloc[:,1:]
vv.head()
 
## SVM classifier

array = df.values
array1 = vv.values

X = array[:,1:241]
Y = array[:,241]

from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression

model = LogisticRegression()
rfe = RFE(model, 3)
fit = rfe.fit(X, Y)

### To be continued
```

### Connection to Spark and MLib

```
install.packages("spkarlyr", repos = "http://cran.us.r-project.org")
library(sparklyr)
spark_available_versions() ## check for the latest versions available at the time.
spark_install("3.2")
Sys.setenv(JAVA_HOME="")
system("java -version") ## check the java version installed.

```

