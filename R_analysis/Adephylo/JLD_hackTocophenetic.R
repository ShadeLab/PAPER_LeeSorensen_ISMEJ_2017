#From Jack Darcy via Reviewer 2 - a hack to storing large pairwise phylogenetic distances

#######
# use ape pagkage for its "read.tree" function
library(ape)
library(phytools)

# set a seed since I use sample a few times
set.seed(12345)
# read in greengenes tree
# available to download from:
# https://groups.google.com/group/qiime-forum/attach/3b133b0c346110c7/gg_13_5_unannotated.tree.gz?part=4&authuser=0&view=1
# ( just use wget )
#bigtree <- read.tree("gg_13_5_unannotated.tree")

###als use our singletons-removed tree to test code instead of gg tree
bigtree <- read.tree("MASTER_RepSeqs_aligned_ns.tre")

# check how many tips the tree has
ntips <- length(bigtree$tip.label)
ntips
# it has 203452 tips.
# we have 26940 tips

# check tree size in memory
object.size(bigtree) / 1E6
# that's only 18 MB in memory
# 2.78 b memory

# compute pairwise branch-length distance
# this should take some time, and consume a lot of memory
bigmat <- cophenetic.phylo(bigtree)
# DIDN'T WORK!!!! This is because R has a max max vector length of 2147483647 elements. 
# 203452^2 = 41392716304, which is 20 times larger than the max. 
# basically, this means we can't make a square matrix on any tree larger than
# sqrt(2147483647) = 46340 tips. If we found a creative way to make a dist object 
# instead of a matrix, that means our tree can have about twice as many (92K tips).

# Let's downsample to test my function in a reasonable time frame
smalltree_1k <- drop.tip(bigtree, sample(bigtree$tip.label, length(bigtree$tip.label) - 1000))
smallmat_1k <- cophenetic.phylo(smalltree_1k)

cophenetic_jld <- function(tree, ncores=12, percentreport=5){
  require(phytools)
  require(ape)
  require(parallel)
  ntips <- Ntip(tree)
  if(ntips > floor(sqrt(.Machine$integer.max))){stop("tree too big!")}
  tipnames <- tree$tip.label
  
  # pre-allocate output dist object
  distproxy <- rep(0, ntips)
  names(distproxy) <- tipnames
  output <- dist(distproxy)
  
  # make calculations
  writestart <- 1
  lastdone <- 0
  for(i in 1:(ntips - 1)){
    jtips <- tipnames[(i+1):ntips]
    itip <- tipnames[i]
    
    # don't use more cores than ya need
    if(length(jtips) < ncores){
      ncores_i <- length(jtips)
    }else{
      ncores_i <- ncores
    }
    
    idists <- simplify2array(mclapply(jtips, FUN=function(x){fastDist(tree, itip, x)}, mc.cores=ncores_i))
    
    writeend <- writestart + length(idists) - 1
    output[writestart:writeend] <- idists
    writestart <- writeend + 1
    percentdone <- round((writeend / length(output)) * 100)
    if(percentdone >= lastdone + percentreport){
      print(paste(percentdone, "%"))
      lastdone <- percentdone
    }
    
  }
  return(output)
  
}

ptm <- proc.time()
smallmat_1k_JLD <- cophenetic_jld(smalltree_1k, ncores=25)
proc.time() - ptm
# took 760 seconds
(((1000 ^2)-1000)/2) / 760
# 657.23 comparisons per second
((((30000 ^2)-30000)/2) / 657.23)/60/60
# would take 190 hours. 
# VERY CPU-limited though, could probably get a geometric speedup with more cores. 
# throw a few hundred cores at it and you'd probably be OK.


# compare my script with cophenetic.phylo
smallmat_1k <- as.dist(smallmat_1k)
data.frame(labels(smallmat_1k), labels(smallmat_1k_JLD))
# same order
pdf("cophenetic_comparison.pdf")
plot(smallmat_1k_JLD ~ smallmat_1k, xlab="ape", ylab="jld", pch=20, cex=0.5)
abline(a=0, b=1, col="red")
dev.off()





