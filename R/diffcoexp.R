#modified from DCe function of DCGL package
#' Differential co-expression analysis
#'
#' This function is used to identify differentially coexpressed links (gene pairs) and enriched genes.
#' @param exprs.1 a data frame or matrix for condition 1, with rows as genes and columns as samples.
#' @param exprs.2 a data frame or matrix for condition 2, with rows as genes and columns as samples.
#' @param rth the cutoff of r; must be within [0,1].
#' @param qth the cutoff of q-value (adjusted p value); must be within [0,1].
#' @param r.diffth the cutoff of absolute value of the difference between the correlation coefficients of the two conditions; must be within [0,1].
#' @param q.diffth the cutoff of q-value (adjusted p value) of the difference between the correlation coefficients of the two conditions; must be within [0,1].
#' @param q.dcgth the cutoff of q-value (adjusted p value) of the genes enriched in the differentilly correlated gene pairs between the two conditions; must be within [0,1].
#' @param r.method a character string specifying the method to be used to calculate correlation coefficients.
#' @param q.method method for adjusting p values.
#' @keywords coexpression
#' @importFrom  igraph graph.data.frame
#' @export
#' @return a list of two data frames.
#'
#' The DCGs data frame contains genes that contribute to differentially correlated links (gene pairs) with q value less than q.dcgth. It has the following columns:
#'   \item{\code{Gene}}{Gene ID}
#'   \item{\code{All.links}}{Number of links with the absolute correlation coefficients greater than rth and q value less than qth in at least one condition}
#'   \item{\code{DC.links}}{Number of links that passed the criteria for All links and the criteria that the absolute differences between the correlation coefficients in the two condition greater than r.diffth and q value less than q.diffth}
#'   \item{\code{DCL_same}}{Number of subset of DC links with same signed correlation coefficients in both conditions}
#'   \item{\code{DCL_diff}}{Number of subset of DC links with opposite signed correlation coefficients in two conditions but only one of them with the absolute correlation coefficients greater than rth and q value less than qth}
#'   \item{\code{DCL_switch}}{Number of subset of DC links with opposite signed correlation coefficients in two conditions and both of them with the absolute correlation coefficients greater than rth and q value less than qth}
#'   \item{\code{p}}{p value of having >=DC.links given All.links}
#'   \item{\code{q}}{adjusted p value}
#'
#' The DCLs data frame contains the differentially correlated links (gene pairs) that meet the criteria that at least one of their correlation coefficients (cor.1 and/or cor.2) is greater that rth with q value (q.1 and/or q.2) less than qth and the absolute value of the difference between the correlation coefficients under two conditions (cor.diff) is greater than r.diffth with q.diffcor < q.diffth. It has the following columns:
#'   \item{\code{Gene.1}}{Gene ID}
#'   \item{\code{Gene.2}}{Gene ID}
#'   \item{\code{cor.1}}{correlation coefficients in condition 1}
#'   \item{\code{cor.2}}{correlation coefficients in condition 2}
#'   \item{\code{p.1}}{p value of correlation coefficients in condition 1}
#'   \item{\code{p.2}}{p value of correlation coefficients in condition 2}
#'   \item{\code{p.diffcor}}{p value of the test of significance for the difference between two correlation coefficients under two conditions using Fisher's r-to-Z transformation}
#'   \item{\code{q.1}}{adjusted p value of correlation coefficients in condition 1}
#'   \item{\code{q.2}}{adjusted p value of correlation coefficients in condition 2}
#'   \item{\code{q.diffcor}}{adjusted p value of the test of significance for the difference between two correlation coefficients under two conditions using Fisher's r-to-Z transformation}
#'   \item{\code{cor.diff}}{difference between correlation coefficients in condition 2 and condition 1}
#'   \item{\code{type}}{can have value "same signed", "diff signed", or "switched opposites". "same signed" indicates that the gene pair has same signed correlation coefficients under both conditions. "diff signed" indicates that the gene pair has opposite signed correlation coefficients under two conditions and only one of them passed the criteria that the absolute correlation coefficients greater than rth and q value less than qth. "switched opposites" indicates that the gene pair has opposite signed correlation coefficients under two conditions and both of them passed the criteria that the absolute correlation coefficients greater than rth and q value less than qth.}
#' @details diffcoexp function takes two gene expression matrices or data frames under two conditions as input, calculates gene-gene correlations under two conditions and compare them with Fisher's Z transformation, filter the correlation with the rth and qth and the correlation changes with r.diffth and q.diffth.
#' @author Wenbin Wei
#' @examples
#' #diffcoexp()
"diffcoexp" <-
function(exprs.1, exprs.2, rth=0.5, qth=0.1, r.diffth=0.5, q.diffth=0.1, q.dcgth=0.1,
	r.method=c('pearson', 'kendall', 'spearman')[1],
	q.method=c("BH","holm", "hochberg", "hommel", "bonferroni", "BY","fdr")[1]) {
    if(!all(rownames(exprs.1)==rownames(exprs.2))) {
        stop("rownames of two expression matrices must be the same!")
    }
    if (length(rownames(exprs.1))==0 | length(rownames(exprs.2))==0) {
		stop('the expression matrices must have row names specifying the gene names.')
	}
	if ( min(ncol(exprs.1),ncol(exprs.2))<3 ){
		stop('each expression matrix must have at least three or more columns.')
	} else if (min(ncol(exprs.1),ncol(exprs.2))<5 ) {
		warning('the minimum number of columns is less than five and the result may not be reliable.')
	}

 	m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
	genes = rownames(exprs.1)

	cor.filtered = rFilter(exprs.1, exprs.2, r.method=r.method, rth=rth, qth=qth)
	if(!is.null(cor.filtered)) {
		print("Finished running rFilter.")
	}
    cor.filtered$q.diffcor<-p.adjust(cor.filtered$p.diffcor, method=q.method)
    cor.filtered$cor.diff<-cor.filtered$cor.2-cor.filtered$cor.1

    # use strsplit to get two-column edge specification.
    if ( nrow(cor.filtered)==0 ) {
		Result <- emptyresult()
		return(Result)
    } else {
        name.all = strsplit(rownames(cor.filtered), ',')
        name.all = matrix(unlist(name.all), length(name.all),2,byrow=T)
        colnames(name.all) <- c("Gene.1", "Gene.2")
        cor.filtered<-data.frame(name.all, cor.filtered)
    }

  	#############################################################
  	## decide three sets of correlation pairs and organize them into two-columned matrices.
  	#############################################################

  	idx.same = (cor.filtered$cor.1 * cor.filtered$cor.2)>0;
	idx.same[is.na(idx.same)] <- TRUE  ##fixing special cases where cor = NA (caused by at least one constant gene expression vector)
  	idx.diff = (cor.filtered$cor.1 * cor.filtered$cor.2)<0;
	idx.diff[is.na(idx.diff)] <- FALSE
  	idx.switched = (cor.filtered$cor.1 * cor.filtered$cor.2 <0) & ( abs(cor.filtered$cor.1)>=rth & abs(cor.filtered$cor.2)>=rth & cor.filtered$q.1 < qth & cor.filtered$q.2 < qth);
	idx.switched[is.na(idx.switched)] <- FALSE

    cor.same = cor.filtered[idx.same,]
	cor.switched = cor.filtered[idx.switched,]
  	cor.diff = cor.filtered[idx.diff & (!idx.switched), ]

    name.same = NULL
    name.switched = NULL
    name.diff = NULL
  #############################################################
  ## Determine DCLs from same sign correlation pairs
  #############################################################
    n.sameDCL = 0
    if(nrow(cor.same)>1){
		de.s = cor.same$q.diffcor < q.diffth & abs(cor.same$cor.diff) > r.diffth
		DCL.same = cor.same[de.s,]
		name.same = DCL.same[, c("Gene.1","Gene.2")]
		n.sameDCL = nrow(DCL.same)
 	} else {
        DCL.same = NULL
    }

  #############################################################
  ## Determine DCLs from different sign correlation pairs
  #############################################################
    n.diffDCL = 0
    if(nrow(cor.diff)>1){
		de.d = cor.diff$q.diffcor < q.diffth & abs(cor.diff$cor.diff) > r.diffth
		DCL.diff = cor.diff[de.d, ]
		name.diff = DCL.diff[, c("Gene.1","Gene.2")]
  		n.diffDCL = nrow(DCL.diff)
	} else {
        DCL.diff = NULL
    }

################################################################################################
## Determine Switched DCLs if they exist
################################################################################################
	n.switchedDCL = 0
    if(nrow(cor.switched)>1){
        de.switched = cor.switched$q.diffcor < q.diffth & abs(cor.switched$cor.diff) > r.diffth
        DCL.switched = cor.switched[de.switched,]
        name.switched = DCL.switched[, c("Gene.1","Gene.2")]
        n.switchedDCL = nrow(DCL.switched)
    } else {
        DCL.switched = NULL
    }

    n.DCL <- n.sameDCL + n.diffDCL + n.switchedDCL
    print(paste(nrow(cor.filtered), "gene pairs remain after half thresholding."))
    if (n.DCL == 0) {
        print("No DCL meets the thresholds!")
		Result <- emptyresult()
		return(Result)
    } else {
        print(paste(n.DCL, "DCL identified."))
    }
	name.DCL=rbind(name.same, name.diff, name.switched);

####################################
## All links
####################################
	g.all <- graph.data.frame(name.all);
	gene.all <- as.matrix(V(g.all)$name);
	degree.all <- degree(g.all);
#####################################
## DCLs
#####################################
	g.DCL <- graph.data.frame(name.DCL);
	gene.1 <- as.matrix(V(g.DCL)$name);
	degree.DCL <- degree(g.DCL);
######################################
##DCLs of same sign
######################################
    if(n.sameDCL>0) {
        g.same <- graph.data.frame(name.same);
	    g.same.name <- as.matrix(V(g.same)$name);
	    degree.same <- as.matrix(degree(g.same));
    } else {
        degree.same = matrix(0,1,1)
    }

########################################
## DCLs of different sign
########################################
    if(n.diffDCL>0) {
	    g.diff <- graph.data.frame(name.diff);
	    g.diff.name <- as.matrix(V(g.diff)$name);
	    degree.diff <- as.matrix(degree(g.diff));
    } else {
        degree.diff = matrix(0,1,1)
    }

#######################################
## DCLs of switched correlation
#######################################
	if(n.switchedDCL>0) {
		g.switch <- graph.data.frame(name.switched);
		g.switch.name <- as.matrix(V(g.switch)$name);
		degree.switch <- as.matrix(degree(g.switch));
	} else {
        degree.switch = matrix(0,1,1)
	}

#######################################
## Numbers for DCLs of different type.
#######################################
	degree.bind <- matrix(0,m,5)
	row.names(degree.bind) <- genes
	colnames(degree.bind) <- c("All.links", "DC.links", "DCL.same", "DCL.diff", "DCL.switched")

	degree.bind[gene.all,1]=degree.all
	degree.bind[gene.1,2]=degree.DCL
    if(n.sameDCL>0) {
	    degree.bind[g.same.name,3]=degree.same
    }
    if(n.diffDCL>0) {
	    degree.bind[g.diff.name,4]=degree.diff
    }
	if(n.switchedDCL>0) {
		degree.bind[g.switch.name,5]=degree.switch
	}

########################################################
## DCGs Identification
########################################################
 	prob <- nrow(name.DCL)/nrow(name.all)
	p.value <- pbinom(degree.bind[,'DC.links']-1, degree.bind[,'All.links'], prob, lower.tail = F, log.p = FALSE);
 	q.value <- p.adjust(p.value, method=q.method);

 	degree.bind <- cbind(degree.bind, p.value, q.value)
 	colnames(degree.bind) <- c("All.links","DC.links","DCL_same","DCL_diff","DCL_switch","p","q")

 	middle <-sort(as.numeric(degree.bind[,'q']), method = "quick", decreasing=FALSE,index.return=TRUE)$ix
 	DCGs <- degree.bind[middle,]
	DCGs <- as.data.frame(DCGs)
	DCGs <- subset(DCGs, subset= q < q.dcgth)
	DCGs <- cbind(Gene=as.character(rownames(DCGs)), DCGs)
	DCGs$Gene <- as.character(DCGs$Gene)

 #########################################################
    DCLs=data.frame()
    if(n.sameDCL>0) {
        DCLs <- rbind(DCLs, data.frame(DCL.same, type='same signed'))
    }

    if(n.diffDCL>0) {
        DCLs <- rbind(DCLs, data.frame(DCL.diff, type='diff signed'))
    }

    if(n.switchedDCL>0){
        DCLs <- rbind(DCLs, data.frame(DCL.switched, type='switched opposites'))
    }

	DCLs$Gene.1 <-as.character(DCLs$Gene.1)
	DCLs$Gene.2 <-as.character(DCLs$Gene.2)

 	Result <- list(DCGs=DCGs,DCLs=DCLs)
	return(Result)
}

"emptyresult"<-function() {
	DCGs = matrix(0,0, 8)
	colnames(DCGs) = c("Gene", "All.links", "DC.links", "DCL.same", "DCL.diff",
	"DCL.switched", "p", "q")
	DCGs<-as.data.frame (DCGs)
	DCLs = matrix(0,0,12)
	colnames(DCLs) <- c("Gene.1", "Gene.2", "cor.1", "cor.2", "p.1", "p.2",
	"p.diffcor", "q.1", "q.2", "q.diffcor", "cor.diff", "type" )
	DCLs = as.data.frame(DCLs, stringsAsFactors =FALSE)
	Result <- list(DCGs=DCGs,DCLs=DCLs)
	return(Result)
}

#The results of diffcoexp can be further analysed using DRsort fucntion of DCGL package
