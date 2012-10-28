pairtest <- function(case,control,n.perm=0,pheno.perm=NULL) {
  if(!is(case,"SnpMatrix") || !is(control,"SnpMatrix") || !identical(colnames(case),colnames(control)))
    stop("case and control SnpMatrix objects must contain the same SNPs, in the same order\n")
  d <- rbind(case,control) 
  if(!is.null(pheno.perm)) {
    if(nrow(pheno.perm) != nrow(d))
      stop("length of phenotype vector must equal total number of cases and controls.\n")
    n.perm <- ncol(pheno.perm)
  }
  pheno <- rep(c(1,0),times=c(nrow(case),nrow(control)))
  permute <- n.perm > 0
  if(!permute) {
    return(p.value(single.snp.tests(phenotype=pheno,snp.data=d),df=1))
  }
  result <- matrix(NA,ncol(case),n.perm)
  pheno.perm <- genperms(pheno,n.perm)
  for(i in 1:n.perm) {
    cat(".")
    result[,i] <- p.value(single.snp.tests(phenotype=pheno.perm[,i],snp.data=d),df=1)
  }
  return(result)
}
