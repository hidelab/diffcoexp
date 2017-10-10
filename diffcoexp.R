######################################################################################
## exprs.1 a data frame or matrix for condition A, with rows as variables (genes) and columns as samples.
## exprs.2 a data frame or matrix for condition B, with rows as variables (genes) and columns as samples.
## qth the cutoff of q-value; must be within [0,1].
######################################################################################
library(DiffCorr)
library(WGCNA)
library(psych)
r2p<-function(r, n) {
	t<-r*sqrt((n-2)/(1-r^2))
	p.value <- 2*pt(-abs(t), n-2)
    return(p.value)
}

"diffcorcalc" <-function(exprs.1, exprs.2, r.method=c('pearson','spearman')[1]) {
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
        n.1 <-count.pairwise(t(exprs.1))
        n.1 <- n.1 [lower.tri(n.1, diag=F)]
    }

    if(sum(is.na(exprs.2))==0) {
        cor.2 <- cor(t(exprs.2), method=r.method, use="all.obs")
        n.2 <- ncol (exprs.2)
    } else {
        cor.2 <- cor(t(exprs.2), method=r.method, use="pairwise.complete.obs")
        n.2<-count.pairwise(t(exprs.2))
        n.2 <- n.2 [lower.tri(n.2, diag=F)]
    }

    cor.1 <- cor.1[lower.tri(cor.1, diag=F)]
    cor.2 <- cor.2[lower.tri(cor.2, diag=F)]
    rm(exprs.1); rm(exprs.2)

    name.row <- matrix(rep(genes,length(genes)),length(genes),length(genes))
    name.col <- matrix(rep(genes,length(genes)),length(genes),length(genes),byrow=T)
    name.pairs <- matrix(paste(name.row,name.col,sep=','),length(genes),length(genes))
    rm(list=c('name.row','name.col'))
    name.pairs <- name.pairs[lower.tri(name.pairs,diag=F)]
    names(cor.1) <- names(cor.2) <- name.pairs

    p.1 <- r2p(cor.1, n.1)
    p.2 <- r2p(cor.2, n.2)

    dc<-compcorr(n.1, cor.1, n.2, cor.2)
    res <- data.frame(cor.1 = cor.1, cor.2 = cor.2, p.1=p.1, p.2=p.2, p.diffcor=dc$pval)
    return(res)
}

"rFilter"<-function(exprs.1, exprs.2, rth=0.5, qth=0.1,
	r.method=c('pearson','spearman')[1],
	q.method=c("BH","holm", "hochberg", "hommel", "bonferroni", "BY","fdr")[1]) {
    if(!all(rownames(exprs.1)==rownames(exprs.2))) {
        stop("rownames of two expression matrices must be the same!")
    }
    x<-diffcorcalc(exprs.1, exprs.2, r.method=r.method)
	if (!is.null(x)) {
		print("Finished running diffcorcalc.")
	}
    x$q.1<-p.adjust(x$p.1, method=q.method)
    x$q.2<-p.adjust(x$p.2, method=q.method)
    x<-subset(x, subset=(abs(cor.1) > rth & q.1 < qth))
    x<-subset(x, subset=(abs(cor.2) > rth & q.2 < qth))
    return(x)
}

#modified from DCe function of DCGL package
"diffcoexp" <-
function(exprs.1, exprs.2, rth=0.5, qth=0.1, r.diffth=0.5, q.diffth=0.1, q.dcgth=0.1,
	r.method=c('pearson','spearman')[1],
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

	cor.filtered.1 = cor.filtered$cor.1;
	cor.filtered.2 = cor.filtered$cor.2;

  	#############################################################
  	## decide three sets of correlation pairs and organize them into two-columned matrices.
  	#############################################################

  	idx.same = (cor.filtered.1*cor.filtered.2)>0;
	idx.same[is.na(idx.same)] <- TRUE  ##fixing special cases where cor = NA (caused by at least one constant gene expression vector)
  	idx.diff = (cor.filtered.1*cor.filtered.2)<0;
	idx.diff[is.na(idx.diff)] <- FALSE
  	idx.switched = (cor.filtered.1*cor.filtered.2<0) & ( abs(cor.filtered.1)>=rth & abs(cor.filtered.2)>=rth );
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
		DCL.diff = cor.diff[de.d, c("Gene.1","Gene.2")]
		name.diff = DCL.diff[,]
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

 	Result <- list(DCGs=DCGs,DCLs=DCLs)
	return(Result)
}

"emptyresult"<-function() {
	DCGs = matrix(0,0, 7)
	colnames(DCGs) = c("All.links", "DC.links", "DCL.same", "DCL.diff",
	"DCL.switched", "p", "q")
	DCGs<-as.data.frame (DCGs)
	DCLs = matrix(0,0,12)
	colnames(DCLs) <- c("Gene.1", "Gene.2", "cor.1", "cor.2", "p.1", "p.2",
	"p.diffcor", "q.1", "q.2", "q.diffcor", "cor.diff", "type" )
	DCLs = as.data.frame(DCLs)
	Result <- list(DCGs=DCGs,DCLs=DCLs)
	return(Result)
}

#Added DRsort from DCGL package
"DRsort" <- function(DCGs, DCLs, tf2target, expGenes) {
  if (!('DCG' %in% colnames(DCGs))) DCGs <- data.frame(DCG=rownames(DCGs),DCGs)
  if (!('DCG' %in% colnames(DCLs))) {
    DCLs <- DCL.connect.DCG(DCLs,DCGs)
  } ## if input DCLs do not have a "DCG" field, reduce the DCLs to those connecting to DCGs and identify the involved DCGs.
  tf2target <- data.frame(tf2target)
  colnames(tf2target) <- c('TF','target')
  ####### END update 09/12/2013 #####################################################################################################
  ###### TFs with expression data##########
  TFinexpression<-intersect(as.character(tf2target$TF),as.character(expGenes)) ### changed from merge to intersect on 9/16/2013
  TFinexpression<-unique(TFinexpression)
  #######################################     DCG: add one field 'DCGisTF' to indicate whether the DCG is a TF or not        ##################################
  TFs <- unique(tf2target$TF)
  DCGs <- data.frame(DCGisTF=FALSE,DCGs)
  TF.idx <- as.character(DCGs$DCG) %in% TFs
  DCGs$DCGisTF[TF.idx] <- TRUE
  DCGs[,c(1,2)] <- DCGs[,c(2,1)]
  colnames(DCGs)[1:2] <- c('DCG','DCGisTF')
  #######################################   DCG: add one field 'Upstream_TFofDCG' to indicate the TF(s) of the DCG    #############################
  #cat('DCG: Add one field \'Upstream_TFofDCG\' to indicate the TF(s) of the DCG\n')
  DCGistarget<-merge(tf2target,DCGs,by.x="target",by.y="DCG",all.y=T)
  colnames(DCGistarget)[1:2] <- c('DCG','Upstream_TFofDCG')
  #DCGistarget <- DCGistarget[,c('DCG','Upstream_TFofDCG','DCGisTF','dC','DCp.p','All.links.DCe','DC.links','DCL.same','DCL.diff','DCL.switch')]
  DCG2TF <- unique(data.frame(DCG=DCGistarget$DCG,TF=DCGistarget$Upstream_TFofDCG))
  DCG2TF <- DCG2TF[!is.na(DCG2TF[,2]),] ### records with NA value for TF were removed. 10/2/2013 
  if (nrow(DCG2TF)==0) {warning('No DCG2TF is found!\n'); return(NULL)}	    
  ######### this object is outputted for easing TED calculation.
  DCG.factor <- factor(DCG2TF$DCG,levels=unique(DCG2TF$DCG),ordered=T)
  DCGs.TF <- as.character(unlist(by(DCG2TF$TF,DCG.factor,paste,collapse=';')))	
  DCG2TF.byDCG <- data.frame(DCG=levels(DCG.factor),TF=DCGs.TF)
  DCGs <- merge(DCG2TF.byDCG,DCGs,by.x='DCG',by.y='DCG',all.y=T)
  colnames(DCGs)[1:2] <- c('DCG','Upstream_TFofDCG')
  #######################################  DCL: Add one field 'internal.TF' to indicate the internal TF of the DCL pair   ############################
  tf.target <- paste(tf2target[,'TF'],tf2target[,'target'],sep='; ')
  DCL.pair1 <- paste(DCLs$Gene.1,DCLs$Gene.2, sep='; ')
  TF1.idx <- DCL.pair1 %in% tf.target
  DCL.pair2 <- paste(DCLs$Gene.2,DCLs$Gene.1, sep='; ')
  TF2.idx <- DCL.pair2 %in% tf.target
  DCLs <- data.frame(TF=NA,DCLs)
  DCLs[TF1.idx,'TF'] <- as.character(DCLs[TF1.idx,'Gene.1'])
  DCLs[TF2.idx,'TF'] <- as.character(DCLs[TF2.idx,'Gene.2'])
  DCLs[TF1.idx & TF2.idx, 'TF'] <- paste(as.character(DCLs[TF1.idx & TF2.idx,'Gene.1']),as.character(DCLs[TF1.idx & TF2.idx,'Gene.2']),sep='; ')
  colnames(DCLs)[1] <- 'internal.TF'
  #######################################  DCL: Add one field 'common.TF' to indicate the common TF(s) of the DCL pair        ##############
  #cat('DCL: Add one field \'common.TF\' to indicate the common TF(s) of the DCL pair \n')
  tfbridgedDCL<-merge(DCLs,tf2target,by.x="Gene.1",by.y="target") ####nrow of tfbridgedDCL is less than nrow of DCL###
  colnames(tfbridgedDCL)[1:2]<- c('Gene.1','internal.TF')   ## add the 'TF' that regulates 'Gene.1'
  colnames(tfbridgedDCL)[ncol(tfbridgedDCL)] <- 'TF.of.g1' ### added on 9/12/2013
  tfbridgedDCL <- merge(tf2target,tfbridgedDCL,by.x=c('TF','target'),by.y=c('TF.of.g1','Gene.2'))
  colnames(tfbridgedDCL)[1:2] <- c('common.TF','Gene.2')  ## extract the rows in which 'TF' that regulates Gene.1 and Gene.2 both.
  tfbridgedDCL <- unique(tfbridgedDCL)
  if (nrow(tfbridgedDCL)==0) {warning('No TF-bridged-DCL is found!\n'); return(NULL)}
  tfbridgedDCL<-data.frame(common.TFisDCG=FALSE,tfbridgedDCL) 
  common.TFDCG.idx <- as.character(tfbridgedDCL$common.TF) %in% as.character(DCGs[DCGs[,'DCGisTF'],'DCG'])
  tfbridgedDCL$common.TFisDCG[common.TFDCG.idx]<-TRUE
  tfbridgedDCL<-data.frame(common.TFinexprs=FALSE,tfbridgedDCL)
  common.TF.exprsid <- as.character(tfbridgedDCL$common.TF) %in% TFinexpression
  tfbridgedDCL$common.TFinexprs[common.TF.exprsid]<-TRUE
  Pairs <- as.character(paste(tfbridgedDCL$Gene.1,tfbridgedDCL$Gene.2,sep='; '))
  Pairs.factor <- factor(Pairs,levels=unique(Pairs),ordered=T)
  commonTF <- as.character(unlist( by(tfbridgedDCL$common.TF,Pairs.factor,paste,collapse='; ')  ))        
  commonTF <- data.frame(pairID=levels(Pairs.factor),common.TF=commonTF)
  DCL.pairID <- paste(DCLs$Gene.1,DCLs$Gene.2,sep='; ')
  DCLs <- data.frame(pairID=DCL.pairID,DCLs)
  DCLs <- merge(commonTF,DCLs,by.x='pairID',by.y='pairID',all.y=T)
  DRGs<-DCGs[DCGs[,'DCGisTF'] | !is.na(DCGs[,'Upstream_TFofDCG']),]
  DRLs<-DCLs[!is.na(DCLs[,'common.TF']) | !is.na(DCLs[,'internal.TF']),]	
  list(DCGs=DCGs, DCLs=DCLs, DRGs=DRGs, DRLs=DRLs, DCG2TF=DCG2TF, TF_bridged_DCL=tfbridgedDCL)
}