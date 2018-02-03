test_diffcoexp <- function() {
    data(gse4158part)
    allowWGCNAThreads()
    res=diffcoexp(exprs.1 = exprs.1, exprs.2 = exprs.2, r.method = "spearman")
    checkEquals(nrow(res$DCGs), 15)
    checkEquals(ncol(res$DCGs), 8)
    checkEquals(nrow(res$DCLs), 363)
    checkEquals(ncol(res$DCLs), 12)
}

test_coexpr <- function() {
    data(gse4158part)
    allowWGCNAThreads()
    res=coexpr(exprs.1 = exprs.1, exprs.2 = exprs.2, r.method = "spearman")
    checkEquals(nrow(res), 13179)
    checkEquals(ncol(res), 11)
}

test_comparecor <- function() {
    data(gse4158part)
    allowWGCNAThreads()
    res=comparecor(exprs.1 = exprs.1, exprs.2 = exprs.2, r.method = "spearman")
    checkEquals(nrow(res), 79800)
    checkEquals(ncol(res), 11)
}
