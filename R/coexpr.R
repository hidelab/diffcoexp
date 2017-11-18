#' Identification of gene pairs coexpressed in at least one of two conditions
#'
#' This function identifies gene pairs coexpressed in at least one of two conditions.
#' @param exprs.1 a data frame or matrix for condition 1, with rows as genes and columns as samples.
#' @param exprs.2 a data frame or matrix for condition 2, with rows as genes and columns as samples.
#' @param rth the cutoff of r; must be within [0,1].
#' @param qth the cutoff of q-value; must be within [0,1].
#' @param r.method a character string specifying the method to be used to calculate correlation coefficients.
#' @param q.method method for adjusting p values.
#' @keywords coexpression
#' @importFrom stats p.adjust pbinom pt
#' @export
#' @return a data frame containing gene pairs that are coexpressed in at least one of the conditions with the criteria that the absolute value of at least one of their correlation coefficients (cor.1 and/or cor.2) is greater that rth with q value (q.1 and/or q.2) less than qth. It has the following columns:
#'   \item{\code{Gene.Pair}}{Gene pair}
#'   \item{\code{cor.1}}{correlation coefficients in condition 1}
#'   \item{\code{cor.2}}{correlation coefficients in condition 2}
#'   \item{\code{p.1}}{p value under null hypothesis that correlation coefficient in condition 1 equals to zero}
#'   \item{\code{p.2}}{p value under null hypothesis that correlation coefficient in condition 2 equals to zero}
#'   \item{\code{p.diffcor}}{p value under null hypothesis that difference between two correlation coefficients under two conditions equals to zero using Fisher's r-to-Z transformation}
#'   \item{\code{q.1}}{adjusted p value under null hypothesis that correlation coefficient in condition 1 equals to zero}
#'   \item{\code{q.2}}{adjusted p value under null hypothesis that correlation coefficient in condition 2 equals to zero}
#'   \item{\code{q.diffcor}}{adjusted p value under null hypothesis that the difference between two correlation coefficients under two conditions equals to zero using Fisher's r-to-Z transformation}
#' @examples
#' data(exprs4158)
#' allowWGCNAThreads()
#' res=coexpr(exprs.1 = exprs.1, exprs.2 = exprs.2, r.method = "spearman")
#' #The result is a data frames.
#' str(res)
"coexpr"<-function(exprs.1, exprs.2, rth=0.5, qth=0.1,
	r.method=c('pearson','spearman')[1],
	q.method=c("BH","holm", "hochberg", "hommel", "bonferroni", "BY","fdr", "none")[1]) {
    if(!all(rownames(exprs.1)==rownames(exprs.2))) {
        stop("rownames of two expression matrices must be the same!")
    }
    x<-comparecor(exprs.1, exprs.2, r.method=r.method)
	if (!is.null(x)) {
		print("Finished running comparecor.")
	}
    x$q.1<-p.adjust(x$p.1, method=q.method)
    x$q.2<-p.adjust(x$p.2, method=q.method)
    x$q.diffcor <- p.adjust(x$p.diffcor, method=q.method)
    x<-subset(x, subset=( (abs(x$cor.1) > rth & x$q.1 < qth) | (abs(x$cor.2) > rth & x$q.2 < qth)) )
    x<-data.frame(Gene.Pair=rownames(x), x,  stringsAsFactors =FALSE)
    return(x)
}
