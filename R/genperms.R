genperms <-
function(pheno,n.perm=0) {
  pheno.perm <- matrix(as.integer(0),length(pheno),n.perm)
  for(j in 1:n.perm) {
    pheno.perm[,j] <- sample(pheno)
  }
  return(pheno.perm)
}
