#' Compare gene-gene correlation coefficients under two conditions
#'
#' This function calculates correlation coefficients of all gene pairs under two conditions and compare them using Fisher's Z-transformation.
#' @param exprs.1 a data frame or matrix for condition 1, with gene IDs as rownames and sample IDs as column names.
#' @param exprs.2 a data frame or matrix for condition 2, with gene IDs as rownames and sample IDs as column names.
#' @param r.method a character string specifying the method to be used to calculate correlation coefficients.
#' @param q.method a character string specifying the method for adjusting p values.
#' @keywords coexpression
#' @importFrom DiffCorr compcorr
#' @importFrom WGCNA cor
#' @importFrom psych count.pairwise
#' @return a data frame containing the differences between the correlation coefficients under two consitions and their p values. It has the following columns:
#'   \item{\code{Gene.1}}{Gene ID}
#'   \item{\code{Gene.2}}{Gene ID}
#'   \item{\code{cor.1}}{correlation coefficients under condition 1}
#'   \item{\code{cor.2}}{correlation coefficients under condition 2}
#'   \item{\code{cor.diff}}{difference between correlation coefficients under condition 2 and condition 1}
#'   \item{\code{p.1}}{p value under null hypothesis that correlation coefficient under condition 1 equals to zero}
#'   \item{\code{p.2}}{p value under null hypothesis that correlation coefficient under condition 2 equals to zero}
#'   \item{\code{p.diffcor}}{p value under null hypothesis that difference between two correlation coefficients under two conditions equals to zero using Fisher's r-to-Z transformation}
#'   \item{\code{q.1}}{adjusted p value under null hypothesis that correlation coefficient under condition 1 equals to zero}
#'   \item{\code{q.2}}{adjusted p value under null hypothesis that correlation coefficient under condition 2 equals to zero}
#'   \item{\code{q.diffcor}}{adjusted p value under null hypothesis that the difference between two correlation coefficients under two conditions equals to zero using Fisher's r-to-Z transformation}
#' @export
#' @examples
#' data(exprs4158)
#' allowWGCNAThreads()
#' res=comparecor(exprs.1 = exprs.1, exprs.2 = exprs.2, r.method = "spearman")
#' #The result is a data frames.
#' str(res)
"comparecor" <-function(exprs.1, exprs.2, r.method=c('pearson','spearman')[1], q.method=c("BH","holm", "hochberg", "hommel", "bonferroni", "BY","fdr", "none")[1]) {
    exprs.1<-exprs.1[!is.na(rownames(exprs.1)), ]
    exprs.1<-exprs.1[rownames(exprs.1) != "", ]
    exprs.2<-exprs.2[!is.na(rownames(exprs.2)), ]
    exprs.2<-exprs.2[rownames(exprs.2) != "", ]
    if(!all(rownames(exprs.1)==rownames(exprs.2))) {
        stop("rownames of two expression matrices must be the same!")
    }
    genes <- rownames(exprs.1)
    exprs.1 <- as.matrix(exprs.1)
    exprs.2 <- as.matrix(exprs.2)
    if(sum(is.na(exprs.1))==0) {
        cor.1 <- cor(t(exprs.1), method=r.method, use="all.obs")
        n.1 <- ncol(exprs.1)
    } else {
        cor.1 <- cor(t(exprs.1), method=r.method, use="pairwise.complete.obs")
        n.1 <- count.pairwise(t(exprs.1))
        n.1 <- n.1[lower.tri(n.1, diag=F)]
    }

    if(sum(is.na(exprs.2))==0) {
        cor.2 <- cor(t(exprs.2), method=r.method, use="all.obs")
        n.2 <- ncol (exprs.2)
    } else {
        cor.2 <- cor(t(exprs.2), method=r.method, use="pairwise.complete.obs")
        n.2 <- count.pairwise(t(exprs.2))
        n.2 <- n.2[lower.tri(n.2, diag=F)]
    }

    cor.1 <- cor.1[lower.tri(cor.1, diag=F)]
    cor.2 <- cor.2[lower.tri(cor.2, diag=F)]
    rm(exprs.1); rm(exprs.2)

    name.row <- matrix(rep(genes, length(genes)), length(genes), length(genes))
    name.col <- matrix(rep(genes, length(genes)), length(genes), length(genes), byrow=T)
    name.pairs <- matrix(paste(name.row, name.col, sep=','), length(genes), length(genes))
    name.pairs <- name.pairs[lower.tri(name.pairs, diag=F)]
    Gene.1 <- name.row[lower.tri(name.row, diag=F)]
    Gene.2 <- name.col[lower.tri(name.col, diag=F)]
    names(Gene.1)<-names(Gene.2) <- name.pairs
    rm(list=c('name.row', 'name.col'))
    p.1 <- r2p(cor.1, n.1)
    p.2 <- r2p(cor.2, n.2)

    dc<-compcorr(n.1, cor.1, n.2, cor.2)
    res <- data.frame(Gene.1=Gene.1, Gene.2=Gene.2, cor.1 = cor.1, cor.2 = cor.2, cor.diff=cor.2-cor.1, p.1 = p.1, p.2 = p.2, p.diffcor = dc$pval, stringsAsFactors =FALSE)
    res$q.1<-p.adjust(res$p.1, method=q.method)
    res$q.2<-p.adjust(res$p.2, method=q.method)
    res$q.diffcor <- p.adjust(res$p.diffcor, method=q.method)
    return(res)
}

r2p<-function(r, n) {
	t<-r*sqrt((n-2)/(1-r^2))
	p.value <- 2*pt(-abs(t), n-2)
    return(p.value)
}
