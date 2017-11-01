#' Identification of co-expressed gene pairs
#'
#' This function is used to identify coexpressed links (gene pairs) in either condition 1 or condition 2.
#' @param exprs.1 a data frame or matrix for condition 1, with rows as genes and columns as samples.
#' @param exprs.2 a data frame or matrix for condition 2, with rows as genes and columns as samples.
#' @param rth the cutoff of r; must be within [0,1].
#' @param qth the cutoff of q-value; must be within [0,1].
#' @param r.method a character string specifying the method to be used to calculate correlation coefficients.
#' @param q.method method for adjusting p values.
#' @keywords coexpression
#' @importFrom stats p.adjust pbinom pt
#' @export
#' @return a data frame
#' @examples
#' #rFilter()
"rFilter"<-function(exprs.1, exprs.2, rth=0.5, qth=0.1,
	r.method=c('pearson','spearman')[1],
	q.method=c("BH","holm", "hochberg", "hommel", "bonferroni", "BY","fdr")[1]) {
    if(!all(rownames(exprs.1)==rownames(exprs.2))) {
        stop("rownames of two expression matrices must be the same!")
    }
    x<-comparecor(exprs.1, exprs.2, r.method=r.method)
	if (!is.null(x)) {
		print("Finished running comparecor.")
	}
    x$q.1<-p.adjust(x$p.1, method=q.method)
    x$q.2<-p.adjust(x$p.2, method=q.method)
    x<-subset(x, subset=( (abs(cor.1) > rth & q.1 < qth) | (abs(cor.2) > rth & q.2 < qth)) )
    return(x)
}


