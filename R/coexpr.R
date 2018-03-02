#' Identification of gene pairs coexpressed in at least one of two conditions
#'
#' This function identifies gene pairs coexpressed in at least one of two
#' conditions.
#' @param exprs.1 a SummarizedExperiment, data frame or matrix
#' for condition 1, with gene IDs as rownames and sample IDs as column names.
#' @param exprs.2 a SummarizedExperiment, data frame or matrix
#' for condition 2, with gene IDs as rownames and sample IDs as column names.
#' @param rth the cutoff of r; must be within [0,1].
#' @param qth the cutoff of q-value; must be within [0,1].
#' @param r.method a character string specifying the method to be used to
#' calculate correlation coefficients.
#' @param q.method a character string specifying the method for adjusting p values.
#' @keywords coexpression
#' @importFrom stats p.adjust pbinom pt
#' @export
#' @return a data frame containing gene pairs that are coexpressed in at least
#' one of the conditions with the criteria that the absolute value of
#' correlation coefficient is greater than rth and q value less than qth. It
#' has the following columns:
#'   \item{\code{Gene.1}}{Gene ID}
#'   \item{\code{Gene.2}}{Gene ID}
#'   \item{\code{cor.1}}{correlation coefficients under condition 1}
#'   \item{\code{cor.2}}{correlation coefficients under condition 2}
#'   \item{\code{cor.diff}}{difference between correlation coefficients under
#' condition 2 and condition 1}
#'   \item{\code{p.1}}{p value under null hypothesis that correlation
#' coefficient under condition 1 equals to zero}
#'   \item{\code{p.2}}{p value under null hypothesis that correlation
#' coefficient under condition 2 equals to zero}
#'   \item{\code{p.diffcor}}{p value under null hypothesis that difference
#' between two correlation coefficients under two conditions equals to zero
#' using Fisher's r-to-Z transformation}
#'   \item{\code{q.1}}{adjusted p value under null hypothesis that correlation
#' coefficient under condition 1 equals to zero}
#'   \item{\code{q.2}}{adjusted p value under null hypothesis that correlation
#' coefficient under condition 2 equals to zero}
#'   \item{\code{q.diffcor}}{adjusted p value under null hypothesis that the
#' difference between two correlation coefficients under two conditions equals
#' to zero using Fisher's r-to-Z transformation}
#' @examples
#' data(gse4158part)
#' allowWGCNAThreads()
#' res=coexpr(exprs.1 = exprs.1, exprs.2 = exprs.2, r.method = "spearman")
#' #The result is a data frames.
#' str(res)
"coexpr"<-function(exprs.1, exprs.2, r.method=c('pearson','spearman')[1],
    q.method=c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr",
    "none")[1], rth=0.5, qth=0.1) {
    if (is(exprs.1, "SummarizedExperiment")) {
        exprs.1<- assays(exprs.1)[[1]]
    }
    if (is(exprs.2, "SummarizedExperiment")) {
        exprs.2<- assays(exprs.2)[[1]]
    }
    exprs.1<-exprs.1[!is.na(rownames(exprs.1)), ]
    exprs.1<-exprs.1[rownames(exprs.1) != "", ]
    exprs.2<-exprs.2[!is.na(rownames(exprs.2)), ]
    exprs.2<-exprs.2[rownames(exprs.2) != "", ]
    if(!all(rownames(exprs.1)==rownames(exprs.2))) {
        stop("Rownames of two expression matrices must be the same!")
    }
    x<-comparecor(exprs.1, exprs.2, r.method=r.method)
    if (!is.null(x)) {
        message("Finished running comparecor.")
    }
    x<-subset(x, subset=( (abs(x$cor.1) > rth & x$q.1 < qth) |
        (abs(x$cor.2) > rth & x$q.2 < qth)) )
    return(x)
}
