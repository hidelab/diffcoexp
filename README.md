diffcoexp
=========
Differential coexpression analysis

##### Wenbin Wei, Sandeep Amberkar, Winston Hide, March 2, 2018

## 1. Description

This package identifies differentially coexpressed links (DCLs) and differentially coexpressed genes (DCGs). DCLs are gene pairs with significantly different correlation coefficients under two conditions (de la Fuente 2010, Jiang et al., 2016). DCGs are genes with significantly more DCLs than by chance (Yu et al., 2011, Jiang et al., 2016). It takes two gene expression matrices or data frames under two conditions as input, calculates gene-gene correlations under two conditions and compares them with Fisher's Z transformation(Fisher 1915 and Fisher 1921). It filters gene pairs with the thresholds for correlation coefficients and their adjusted p value as well as the thresholds for the difference between the two correlation coefficients and its adjusted p value. It identifies DCGs using binomial probability model (Jiang et al., 2016).

The main steps are as follows:

a). Correlation coefficients and p values of all gene pairs under two conditions are calculated.

b). The difference between the correlation coefficients  under two conditions are calculated and the p value is calculated using Fisher's Z-transformation.

c). p values are adjusted.

d). Gene pairs (links) coexpressed in at least one condition are identified using the criteria that at least one of the correlation coefficients under two conditions having absolute value greater than the threshold *rth* and the adjusted p value less than the threshold *qth*. The links that meet the criteria are included in co-expressed links (CLs).

e). Differentially coexpressed links (gene pairs) are identified from CLs using the criteria that the absolute value of the difference between the two correlation coefficients is greater than the threshold *r.diffth* and the adjusted p value is less than the threshold *q.diffth*. The links that meet the criteria are included in DCLs.

f). The DCLs are classified into three categories: *same signed*, *diff signed*, or *switched opposites*. *same signed* indicates that the gene pair has same signed correlation coefficients under both conditions. *diff signed* indicates that the gene pair has oppositely signed correlation coefficients under two conditions and only one of them meets the criteria that the absolute correlation coefficient is greater than the threshold *rth* and adjusted p value less than the threshold *qth*. *switched opposites* indicates that the gene pair has oppositely signed correlation coefficients under two conditions and both of them meet the criteria that the absolute correlation coefficients are greater than the threshold *rth* and adjusted p value less than the threshold *qth*.

g). All the genes in DCLs are tested for their enrichment of DCLs, i.e, whether they have more DCLs than by chance using binomial probability model (Jiang et al., 2016). Those with adjusted p value less than the threshold *q.dcgth* are included in DCGs.

## 2. Installation and removal
From April 2018, this package will be available from Bioconductor and can be
installed within R as follows:
```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("dicoexp")
```
To install this package from GitHub, start R and enter:
```R
library(devtools)
devtools::install_git("git://github.com/hidelab/diffcoexp.git", branch = "master")
```
The above method does not build and install vignette. To install the package with vignette, enter the following from command line:
```
git clone https://github.com/hidelab/diffcoexp.git
R CMD build diffcoexp
R CMD check diffcoexp_0.99.1.tar.gz
R CMD INSTALL diffcoexp_0.99.1.tar.gz
```
To remove this package, start R and enter:
```R
remove.packages("diffcoexp")
```

## 3.Input and output of *diffcoexp* function
The main function of this package is *diffcoexp* function. The first two arguments, *exprs.1* and *exprs.2*, are normalized gene expression data under two conditions with rows as genes and columns as samples. They should be objects of classes *SummarizedExperiment*, *data.frame* or *matrix*. Both should have the same number of genes in the same order. The third argument *r.method* is passed to the *cor* function of the *WGCNA* package as argument *method*, details of which can be found by typing
```R
help(cor, WGCNA)
```
The fourth argument *q.method* is passed to the *p.adjust* function of the *stats* package as argument *method*, details of which can be found by typing
```R
help(p.adjust, stats)
```
Details of other arguments of *diffcoexp* function can be found by typing
```R
help(diffcoexp, diffcoexp)
```
The output of *diffcoexp* function is a list of two data frames, one for differentially co-expressed links (DCLs), the other for differentially co-expressed genes (DCGs). Further details of the output can be seen on the help page.

##4. Analysis and interpretation of DCGs and DCLs
DCGs are a list of genes and therefore can be further analysed using other tools such as FGNet (https://bioconductor.org/packages/release/bioc/html/FGNet.html), clusterProfiler (https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) and enrichr (http://amp.pharm.mssm.edu/Enrichr/). DCLs are a list of differentially co-expressed gene pairs and can be assembled into a differential coexpression network. The network is scale-free but not smallworld (Hsu et al., 2017). The network can be visualized and analyzed using igraph (https://cran.r-project.org/web/packages/igraph/index.html). DCLs can also be further analyzed to identify upstream causal regulators using other tools such as DCGL v2.0 (Yang et al., 2013).

## 5. Example

This example illustrates the workflow of downloading gene expression data from GEO and identifying differentially coexpressed links (DCLs) and differentially coexpressed genes (DCGs).

```R
library(GEOquery)
gse4158 <- getGEO("GSE4158")
exprs<-exprs(gse4158[[1]])
keep<-rowSums(is.na(exprs)) < ncol(exprs)/5
exprs<-exprs[keep,]
dim(exprs)
GPL3415<-getGEO("GPL3415")
exprs<-data.frame(ID=rownames(exprs), exprs)
exprs<-merge(GPL3415@dataTable@table, exprs, by.x="ID", by.y="ID")
colnames(exprs)
exprs<-exprs[, c(7, 11:36)]
exprs<-aggregate(exprs[, -1], by=list(Gene=exprs$ORF), FUN=mean, na.action = na.omit)
rownames(exprs)<-exprs$Gene
exprs<-exprs[, -1]
```
Analysis of all the genes (6104) will take about 20 minutes on a computer with 8 cores and 16GB RAM.
```R
exprs.1<-exprs[, c(1:14)]
exprs.2<-exprs[, c(15:26)]
library(diffcoexp)
allowWGCNAThreads()
res=diffcoexp(exprs.1 = exprs.1, exprs.2 = exprs.2, r.method = "spearman" )
```
The results are a list of two data frames, one for differentially co-expressed links (DCLs, gene pairs), one for differentially co-expressed genes (DCGs).
```R
str(res)
sessionInfo()
```
## References
de la Fuente A (2010). From “differential expression” to “differential networking” –
identification of dysfunctional regulatory networks in diseases. *Trends in Genetics*, 26(7):326-33.

Fisher, R. A. (1915). Frequency distribution of the values of the correlation coefficient in samples of an indefinitely large population. *Biometrika*, 10 (4): 507–521. 

Fisher, R. A. (1921). On the 'probable error' of a coefficient of correlation deduced from a small sample. *Metron*, 1: 3–32.

Hsu C-L, Juan H-F, Huang H-C (2015). Functional analysis and characterization of differential coexpression networks. *Scientific Reports*, 5: 13295

Jiang Z, Dong X, Li Z-G, He F, Zhang Z (2016). Differential coexpression analysis reveals extensive rewiring of Arabidopsis gene coexpression in response to Pseudomonas syringae infection. *Scientific Reports*, 6(1):35064.

Yang J, Yu H, Liu B-H, Zhao Z, Liu L, Ma L-X, et al. (2013) DCGL v2.0: An R package for unveiling differential regulation from differential co-expression. *PLoS ONE*, 8(11):e79729.

Yu H, Liu B-H, Ye Z-Q, Li C, Li Y-X, Li Y-Y (2011). Link-based quantitative methods to identify differentially coexpressed genes and gene pairs. *BMC bioinformatics*, 12(1):315.
