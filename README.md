diffcoexp
=========
Differential coexpression analysis

## 1. Description

This package is developed for the identification of differentially coexpressed links (gene pairs) and enriched genes.

## 2. Installation and removal

To install this package, start R and enter:
```R
library(devtools)
devtools::install_git("git://github.com/hidelab/diffcoexp.git", branch = "master")
```
The above method does not build and install vignette. To install the package with vignette, enter the following from command line:
```bash
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

```{r message=FALSE, warning = FALSE}
library(GEOquery)
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4158/matrix/GSE4158_series_matrix.txt.gz", "GSE4158_series_matrix.txt.gz")
gse4158 <- getGEO(filename='GSE4158_series_matrix.txt.gz', GSEMatrix=TRUE)
exprs<-exprs(gse4158)
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
Analysis of all the genes will take about one hour. Here only 500 genes were analysed so that the vignette and the package could be built and installed quickly.
```{r message=FALSE, warning = FALSE, results = 'hide'}
exprs<-exprs[c(1:500), ]
exprs.1<-exprs[, c(1:14)]
exprs.2<-exprs[, c(15:26)]
library(diffcoexp)
allowWGCNAThreads()
res=diffcoexp(exprs.1 = exprs.1, exprs.2 = exprs.2, r.method = "spearman" )
```
The results are a list of two data frames, one for differentially co-expressed links (DCLs, gene pairs), one for differentially co-expressed genes (DCGs).
```{r message=FALSE, warning = FALSE}
str(res)
sessionInfo()
```
