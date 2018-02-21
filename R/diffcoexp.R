#modified from DCe function of DCGL package
#' Differential co-expression analysis
#'
#' This function identifies differentially coexpressed links (DCLs) and
#' differentially coexpressed genes (DCGs).
#' @param exprs.1 a SummarizedExperiment, data frame or matrix
#' for condition 1, with gene IDs as rownames and sample IDs as column names.
#' @param exprs.2 a SummarizedExperiment, data frame or matrix
#' for condition 2, with gene IDs as rownames and sample IDs as column names.
#' @param rth the cutoff of r; must be within [0,1].
#' @param qth the cutoff of q-value (adjusted p value); must be within [0,1].
#' @param r.diffth the cutoff of absolute value of the difference between the
#' correlation coefficients of the two conditions; must be within [0,1].
#' @param q.diffth the cutoff of q-value (adjusted p value) of the difference
#' between the correlation coefficients of the two conditions; must be
#' within [0,1].
#' @param q.dcgth the cutoff of q-value (adjusted p value) of the genes
#' enriched in the differentilly correlated gene pairs between the two
#' conditions; must be within [0,1].
#' @param r.method a character string specifying the method to be used to
#' calculate correlation coefficients.
#' @param q.method a character string specifying the method for adjusting p
#' values.
#' @keywords coexpression
#' @importFrom  igraph graph.data.frame
#' @export
#' @return a list of two data frames.
#'
#' The DCGs data frame contains genes that contribute to differentially
#' correlated links (gene pairs) with q value less than q.dcgth.
#' It has the following columns:
#'   \item{\code{Gene}}{Gene ID}
#'   \item{\code{CLs}}{Number of links with the absolute correlation
#' coefficients greater than rth and q value less than qth in at least one
#' condition}
#'   \item{\code{DCLs}}{Number of links that meet the criteria for CLs and the
#' criteria that the absolute differences between the correlation coefficients
#' in the two condition greater than r.diffth and q value less than q.diffth}
#'   \item{\code{DCL.same}}{Number of subset of DCLs with same signed
#' correlation coefficients in both conditions}
#'   \item{\code{DCL.diff}}{Number of subset of DCLs with oppositely signed
#' correlation coefficients under two conditions but only one of them with the
#' absolute correlation coefficients greater than rth and q value less than qth}
#'   \item{\code{DCL.switch}}{Number of subset of DCLs with oppositely signed
#' correlation coefficients under two conditions and both of them with the
#' absolute correlation coefficients greater than rth and q value less than qth}
#'   \item{\code{p}}{p value of having >=DCLs given CLs}
#'   \item{\code{q}}{adjusted p value}
#'
#' The DCLs data frame contains the differentially correlated links (gene pairs)
#' that meet the criteria that at least one of their correlation coefficients
#' (cor.1 and/or cor.2) is greater than rth with q value (q.1 and/or q.2) less
#' than qth and the absolute value of the difference between the correlation
#' coefficients under two conditions (cor.diff) is greater than r.diffth with
#' q.diffcor less than q.diffth. It has the following columns:
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
#'   \item{\code{type}}{can have value "same signed", "diff signed", or
#' "switched opposites". "same signed" indicates that the gene pair has same
#' signed correlation coefficients under both conditions. "diff signed"
#' indicates that the gene pair has oppositely signed correlation coefficients
#' under two conditions and only one of them meets the criteria that the
#' absolute correlation coefficients greater than rth and q value less than qth.
#' "switched opposites" indicates that the gene pair has oppositely signed
#' correlation coefficients under two conditions and both of them meet the
#' criteria that the absolute correlation coefficients greater than rth and q
#' value less than qth.}
#' @details diffcoexp function identifies differentially coexpressed links
#' (DCLs) and differentially coexpressed genes (DCGs). DCLs are gene pairs with
#' significantly different correlation coefficients under two conditions (de la
#' Fuente 2010, Jiang et al., 2016). DCGs are genes with significantly more DCLs
#' than by chance (Yu et al., 2011, Jiang et al., 2016). It takes two gene
#' expression matrices or data frames under two conditions as input, calculates
#' gene-gene correlations under two conditions and compare them with Fisher's Z
#' transformation, filter the correlation with the rth and qth and the
#' correlation changes with r.diffth and q.diffth. It identifies DCGs using
#' binomial probability model (Jiang et al., 2016).
#'
#' The main steps are as follows:
#'
#' a). Correlation coefficients and p values of all gene pairs under two
#' conditions are calculated.
#'
#' b). The difference between the correlation coefficients  under two conditions
#' are calculated and the p value is calculated using Fisher's Z-transformation.
#'
#' c). p values are adjusted.
#'
#' d). Gene pairs (links) coexpressed in at least one condition are identified
#' using the criteria that at least one of the correlation coefficients under
#' two conditions having absolute value greater than the threshold rth and the
#' adjusted p value less than the threshold qth. The links that meet the
#' criteria are included in CLs.
#'
#' e). Differentially coexpressed gene pairs (links) are identified from CLs
#' using the criteria that the absolute value of the difference between the two
#' correlation coefficients is greater than the threshold r.diffth and the
#' adjusted p value is less than the threshold q.diffth. The links that meet
#' the criteria are included in DCLs.
#'
#' f). The DCLs are classified into three categories: "same signed",
#' "diff signed", or "switched opposites". "same signed" indicates that the gene
#' pair has same signed correlation coefficients under both conditions.
#' "diff signed" indicates that the gene pair has oppositely signed correlation
#' coefficients under two conditions and only one of them meets the criteria
#' that the absolute correlation coefficient is greater than the threshold rth
#' and adjusted p value less than the threshold qth. "switched opposites"
#' indicates that the gene pair has oppositely signed correlation coefficients
#' under two conditions and both of them meet the criteria that the absolute
#' correlation coefficients are greater than the threshold rth and adjusted p
#' value less than the threshold qth.
#'
#' g). All the genes in DCLs are tested for their enrichment of DCLs, i.e,
#' whether they have more DCLs than by chance using binomial probability model
#' (Jiang et al., 2016). Those with adjusted p value less than the threshold
#' q.dcgth are included in DCGs.
#' @author Wenbin Wei
#' @references
#' 1. de la Fuente A. From "differential expression" to "differential networking"
#' - identification of dysfunctional regulatory networks in diseases. Trends in
#' Genetics. 2010 Jul;26(7):326-33.
#'
#' 2. Jiang Z, Dong X, Li Z-G, He F, Zhang Z. Differential Coexpression Analysis
#' Reveals Extensive Rewiring of Arabidopsis Gene Coexpression in Response to
#' Pseudomonas syringae Infection. Scientific Reports. 2016 Dec;6(1):35064.
#'
#' 3. Yu H, Liu B-H, Ye Z-Q, Li C, Li Y-X, Li Y-Y. Link-based quantitative
#' methods to identify differentially coexpressed genes and gene pairs. BMC
#' bioinformatics. 2011;12(1):315.
#' @examples
#' data(gse4158part)
#' allowWGCNAThreads()
#' res=diffcoexp(exprs.1 = exprs.1, exprs.2 = exprs.2, r.method = "spearman")
#' #The results are a list of two data frames, one for differentially co-expressed
#' #links (DCLs, gene pairs) and one for differentially co-expressed genes (DCGs).
#' str(res)
"diffcoexp" <-
function(exprs.1, exprs.2, rth=0.5, qth=0.1, r.diffth=0.5, q.diffth=0.1,
    q.dcgth=0.1, r.method=c('pearson', 'kendall', 'spearman')[1],
    q.method=c("BH","holm", "hochberg", "hommel", "bonferroni", "BY", "fdr",
    "none")[1]) {
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
    if (length(rownames(exprs.1))==0 | length(rownames(exprs.2))==0) {
        stop('The expression matrices must have row names specifying the gene names.')
    }
    if ( min(ncol(exprs.1),ncol(exprs.2))<3 ){
        stop('Each expression matrix must have at least three or more columns.')
    } else if (min(ncol(exprs.1),ncol(exprs.2))<5 ) {
        warning('The minimum number of columns is less than five and the result
        may not be reliable.')
    }

    m <- nrow(exprs.1)
    genes = rownames(exprs.1)

    colinks = coexpr(exprs.1, exprs.2, r.method=r.method, rth=rth, qth=qth)
    if(!is.null(colinks)) {
        message("Finished running coexpr.")
    }

    if ( nrow(colinks)==0 ) {
        Result <- emptyresult()
        return(Result)
    }

#colinks$cor.diff<-colinks$cor.2-colinks$cor.1
#############################################################
## decide three sets of correlation pairs and organize them into two-columned matrices.
#############################################################
    idx.same = (colinks$cor.1 * colinks$cor.2)>0;
    idx.same[is.na(idx.same)] <- TRUE
    idx.diff = (colinks$cor.1 * colinks$cor.2)<0;
    idx.diff[is.na(idx.diff)] <- FALSE
    idx.switched = (colinks$cor.1 * colinks$cor.2 <0) &
        ( abs(colinks$cor.1)>=rth & abs(colinks$cor.2)>=rth &
        colinks$q.1 < qth & colinks$q.2 < qth);
    idx.switched[is.na(idx.switched)] <- FALSE

    cor.same = colinks[idx.same,]
    cor.switched = colinks[idx.switched,]
    cor.diff = colinks[idx.diff & (!idx.switched), ]

    name.same = NULL
    name.switched = NULL
    name.diff = NULL

#############################################################
## Determine DCLs from same sign correlation pairs
#############################################################
    n.sameDCL = 0
    if(nrow(cor.same)>1){
        idx.DCL.same = cor.same$q.diffcor < q.diffth &
            abs(cor.same$cor.diff) > r.diffth
        DCL.same = cor.same[idx.DCL.same,]
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
        idx.DCL.diff = cor.diff$q.diffcor < q.diffth &
            abs(cor.diff$cor.diff) > r.diffth
        DCL.diff = cor.diff[idx.DCL.diff, ]
        name.diff = DCL.diff[, c("Gene.1","Gene.2")]
        n.diffDCL = nrow(DCL.diff)
    } else {
        DCL.diff = NULL
    }

#############################################################################
## Determine Switched DCLs if they exist
#############################################################################
    n.switchedDCL = 0
    if(nrow(cor.switched)>1){
        idx.DCL.switched = cor.switched$q.diffcor < q.diffth &
            abs(cor.switched$cor.diff) > r.diffth
        DCL.switched = cor.switched[idx.DCL.switched,]
        name.switched = DCL.switched[, c("Gene.1","Gene.2")]
        n.switchedDCL = nrow(DCL.switched)
    } else {
        DCL.switched = NULL
    }

    n.DCL <- n.sameDCL + n.diffDCL + n.switchedDCL
    message(nrow(colinks), " gene pairs remain after half thresholding.")
    if (n.DCL == 0) {
        message("No DCL meets the thresholds!")
        Result <- emptyresult()
        return(Result)
    } else {
        message(n.DCL, " DCLs identified.")
    }
    name.DCL=rbind(name.same, name.diff, name.switched);

####################################
## colinks
####################################
    name.colinks = colinks[, c("Gene.1", "Gene.2")]
    g.colinks <- graph.data.frame(name.colinks);
    g.colinks.name <- as.matrix(igraph::V(g.colinks)$name);
    degree.colinks <- igraph::degree(g.colinks);

#####################################
## DCLs
#####################################
    g.DCL <- graph.data.frame(name.DCL);
    g.DCL.name <- as.matrix(igraph::V(g.DCL)$name);
    degree.DCL <- igraph::degree(g.DCL);

######################################
##DCLs of same sign
######################################
    if(n.sameDCL>0) {
        g.same <- graph.data.frame(name.same);
        g.same.name <- as.matrix(igraph::V(g.same)$name);
        degree.same <- as.matrix(igraph::degree(g.same));
    } else {
        degree.same = matrix(0,1,1)
    }

########################################
## DCLs of different sign
########################################
    if(n.diffDCL>0) {
        g.diff <- graph.data.frame(name.diff);
        g.diff.name <- as.matrix(igraph::V(g.diff)$name);
        degree.diff <- as.matrix(igraph::degree(g.diff));
    } else {
        degree.diff = matrix(0,1,1)
    }

#######################################
## DCLs of switched correlation
#######################################
    if(n.switchedDCL>0) {
        g.switch <- graph.data.frame(name.switched);
        g.switch.name <- as.matrix(igraph::V(g.switch)$name);
        degree.switch <- as.matrix(igraph::degree(g.switch));
    } else {
        degree.switch = matrix(0,1,1)
    }

#######################################
## Numbers for DCLs of different type.
#######################################
    degree.bind <- matrix(0,m,5)
    row.names(degree.bind) <- genes
    colnames(degree.bind) <- c("CLs", "DCLs", "DCL.same", "DCL.diff", "DCL.switched")

    degree.bind[g.colinks.name,1]=degree.colinks
    degree.bind[g.DCL.name,2]=degree.DCL
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
    prob <- nrow(name.DCL)/nrow(name.colinks)
    p.value <- pbinom(degree.bind[,'DCLs']-1, degree.bind[,'CLs'], prob,
        lower.tail = FALSE, log.p = FALSE);
    q.value <- p.adjust(p.value, method=q.method);
    degree.bind <- cbind(degree.bind, p.value, q.value)
    colnames(degree.bind) <- c("CLs","DCLs","DCL.same","DCL.diff","DCL.switch","p","q")
    DCGs <- degree.bind
    DCGs <- as.data.frame(DCGs)
    DCGs <- subset(DCGs, subset= q < q.dcgth)
    DCGs <- cbind(Gene=as.character(rownames(DCGs)), DCGs)
    DCGs$Gene <- as.character(DCGs$Gene)
    o<-order(DCGs$p)
    DCGs<-DCGs[o,]
    message(length(DCGs$Gene), " DCGs identified.")

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
    colnames(DCGs) = c("Gene", "CLs", "DCLs", "DCL.same", "DCL.diff",
        "DCL.switched", "p", "q")
    DCGs<-as.data.frame (DCGs)
    DCLs = matrix(0,0,12)
    colnames(DCLs) <- c("Gene.1", "Gene.2", "cor.1", "cor.2", "p.1", "p.2",
        "p.diffcor", "q.1", "q.2", "q.diffcor", "cor.diff", "type" )
    DCLs = as.data.frame(DCLs, stringsAsFactors =FALSE)
    Result <- list(DCGs=DCGs, DCLs=DCLs)
    return(Result)
}

#The results of diffcoexp can be further analysed using DRsort fucntion of
#DCGL package
