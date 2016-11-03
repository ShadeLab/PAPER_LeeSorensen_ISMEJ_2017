Sys.getenv("R_LIBS_USER")

library("adephylo")
library("picante")

phylo.ns=read.tree("MASTER_RepSeqs_aligned_ns.tre")
#phylo.d=cophenetic(phylo.ns)
phylo.d=distTips(phylo.ns, method="patristic", useC=TRUE)
phylo.d=as.matrix(phylo.d)
write.table(phylo.d, "Phylod_26Apr16.txt",sep="\t",quote=FALSE)

