#' Compare gene-gene correlation coefficients under two conditions
#'
#' This function calculates correlation coefficients of gene pairs in condition 1 and condition 2 and compare them using Fisher's Z-transformation.
#' @param exprs.1 a data frame or matrix for condition 1, with rows as genes and columns as samples.
#' @param exprs.2 a data frame or matrix for condition 2, with rows as genes and columns as samples.
#' @param r.method a character string specifying the method to be used to calculate correlation coefficients.
#' @keywords coexpression
#' @importFrom DiffCorr compcorr
#' @importFrom WGCNA cor
#' @importFrom psych count.pairwise
#' @return a data frame
#' @export
#' @examples
#' #comparecor()
"comparecor" <-function(exprs.1, exprs.2, r.method=c('pearson','spearman')[1]) {
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
    rm(list=c('name.row', 'name.col'))
    name.pairs <- name.pairs[lower.tri(name.pairs, diag=F)]
    names(cor.1) <- names(cor.2) <- name.pairs

    p.1 <- r2p(cor.1, n.1)
    p.2 <- r2p(cor.2, n.2)

    dc<-compcorr(n.1, cor.1, n.2, cor.2)
    res <- data.frame(cor.1 = cor.1, cor.2 = cor.2, p.1 = p.1, p.2 = p.2, p.diffcor = dc$pval)
    return(res)
}

r2p<-function(r, n) {
	t<-r*sqrt((n-2)/(1-r^2))
	p.value <- 2*pt(-abs(t), n-2)
    return(p.value)
}
