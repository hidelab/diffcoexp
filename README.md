diffcoexp
=========
Differential coexpression analysis

## 1. Description

This package identifies differentially coexpressed links (DCLs) and differentially coexpressed genes (DCGs). DCLs are gene pairs with significantly different correlation coefficients under two conditions (de la Fuente 2010, Jiang et al., 2016). DCGs are genes with significantly more DCLs than by chance (Yu et al., 2011, Jiang et al., 2016). It takes two gene expression matrices or data frames under two conditions as input, calculates gene-gene correlations under two conditions and compare them with Fisher's Z transformation. It filters gene pairs with the thresholds for correlation coefficients and their adjusted p value as well as the thresholds for the absolute value of the difference between the two correlation coefficients and its adjusted p value. It identifies DCGs using binomial probability model (Jiang et al., 2016).

The main steps are as follows:

a). In this step, gene pairs (links) are filtered using the criteria that at least one of the correlation coefficients under two conditions having absolute value greater than the threshold and the adjusted p value less than the threshold. The links that meet the criteria are included in All.links.

b). In this step, gene pairs (links) are furhter filtered using the criteria that the absolute value of the difference between the two correlation coefficients is greater the threshold and the adjusted p value is less than the threshold. The links that meet the criteria are included in DCLs and DC.links.

c). In this step, the DCLs are classified into three categories: "same signed", "diff signed", or "switched opposites". "same signed" indicates that the gene pair has same signed correlation coefficients under both conditions. "diff signed" indicates that the gene pair has opposite signed correlation coefficients under two conditions and only one of them meets the criteria that the absolute correlation coefficients greater than the threshold and adjusted p value less than the threshold. "switched opposites" indicates that the gene pair has opposite signed correlation coefficients under two conditions and both of them meet the criteria that the absolute correlation coefficients greater than the threshold and adjusted p value less than the threshold.

d). In this step, all the genes in DCLs are tested for their enrichment, i.e, whether they have more DC links than by chance using binomial probability model (Jiang et al., 2016) 

## 2. Installation and removal

To install this package, start R and enter:
```R
library(devtools)
devtools::install_git("git://github.com/hidelab/diffcoexp.git", branch = "master")
```
The above method does not build and install vignette. To install the package with vignette, enter the following from command line:
```
git clone https://github.com/hidelab/diffcoexp.git
R CMD build diffcoexp
R CMD check diffcoexp_0.0.0.9000.tar.gz
R CMD INSTALL diffcoexp_0.0.0.9000.tar.gz
```
To remove this package, start R and enter:
```R
remove.packages("diffcoexp")
```
## 3. Example

This example illustrates the workflow of downloading gene expression data from GEO and identifying differentially coexpressed links (gene pairs) and enriched genes. 

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
Analysis of all the genes will take about one hour.
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
## 4. References
1. de la Fuente A. From “differential expression” to “differential networking” – identification of dysfunctional regulatory networks in diseases. Trends in Genetics. 2010 Jul;26(7):326–33. 

2. Jiang Z, Dong X, Li Z-G, He F, Zhang Z. Differential Coexpression Analysis Reveals Extensive Rewiring of Arabidopsis Gene Coexpression in Response to Pseudomonas syringae Infection. Scientific Reports [Internet]. 2016 Dec [cited 2017 Sep 20];6(1). Available from: http://www.nature.com/articles/srep35064

3. Yu H, Liu B-H, Ye Z-Q, Li C, Li Y-X, Li Y-Y. Link-based quantitative methods to identify differentially coexpressed genes and gene pairs. BMC bioinformatics. 2011;12(1):315. 
