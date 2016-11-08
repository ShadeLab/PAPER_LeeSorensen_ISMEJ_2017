#MODIFIED to use our dataset (comm.t) instead of "dune" and to only include the abundance-based model.  We also changed the number of patches to by 18 to match with the dataset.
#ORIGINAL scripts available in the appendix of the work below, published in Oikos (Appendix oik.02803)
#Note that beta null models may require ~24 hours to complete
###Read in map file
map=read.table("InputFiles/Centralia_Collapsed_Map_forR.txt", header=TRUE, sep="\t")

###Read in microbial community data
comm=read.table("MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1, check.names=FALSE, sep="\t")

#remove consensus lineage from otu table
rdp=comm[,"ConsensusLineage"]

comm=comm[,-ncol(comm)]

#How many total QCed sequences?
sum(colSums(comm))

#sort community by colnames (to be in the consistent, consecutive order for all analyses)
comm=comm[,order(colnames(comm))]

#who are the singleton OTUs (observed 1 time in an abundance of 1 sequence)?
singletonOTUs=row.names(comm)[rowSums(comm)==1]
length(singletonOTUs)
#total 1374 singleton OTUs
g=grep("_dn", singletonOTUs)
length(g)
#1201 de novo OTUs are singletons

#who are the remaining de novo OTUs?
g=grep("_dn_",row.names(comm))
dn=rdp[g]

rdp.nosigs=rdp[rowSums(comm)>1]

#designate a full dataset
comm.sigs=comm

#remove OTUs with an abundance = 1, across the entire dataset (singleton OTUs)
comm=comm[rowSums(comm)>1,]
sum(colSums(comm))

#transpose matrix
comm.t=t(comm)
###

#######################
### Code for example metacommunity simulation and beta-null deviation calculations
### with "Differentiating between niche and neutral assembly in metacommunities using
### null models of beta-diversity"
### Prepared May 14, 2014 
### Authors Caroline Tucker, Lauren Shoemaker, Brett Melbourne
#######################

## Load required source files and libraries
library(reldist)
library(vegan)
library(bipartite)
source("oik-02803-appendix-to-Tucker2016/MetacommunityDynamicsFctsOikos.R")
source("oik-02803-appendix-to-Tucker2016/PANullDevFctsOikos.R")

##packages for UniFrac Null Model (weighted) #als add
library(GUniFrac)
library(ape)
library(phangorn)
tree <- read.tree("MASTER_RepSeqs_aligned_clean.tre")
is.rooted(tree)

#https://github.com/joey711/phyloseq/issues/235
#FastUniFrac trees are unrooted; calculation is done using mid-point root.
tree <- midpoint(tree)
is.rooted(tree)

#formatting problem with tree tip labels - for some reason tree dn OTUs have extra quotes around them and this needs to be removed
tree$tip.label=gsub("'","", tree$tip.label)


### Prepare and calculate abundance beta-null deviation metric
## Adjusted from Stegen et al 2012 GEB
bbs.sp.site <- comm.t
patches=nrow(bbs.sp.site)
rand <- 999

#note - two randomization runs in < 8 min on my laptop 
null.alphas <- matrix(NA, ncol(comm.t), rand)
null.alpha <- matrix(NA, ncol(comm.t), rand)
expected_beta <- matrix(NA, 1, rand)
null.gamma <- matrix(NA, 1, rand)
null.alpha.comp <- numeric()
bucket_bray_res <- matrix(NA, patches, rand)
bucket_wuf_res <- matrix(NA, patches, rand) #als add

bbs.sp.site = ceiling(bbs.sp.site/max(bbs.sp.site)) 
mean.alpha = sum(bbs.sp.site)/nrow(bbs.sp.site) #mean.alpha
gamma <- ncol(bbs.sp.site) #gamma
obs_beta <- 1-mean.alpha/gamma
obs_beta_all <- 1-rowSums(bbs.sp.site)/gamma

##Generate null patches
for (randomize in 1:rand) {  
  null.dist = comm.t
  for (species in 1:ncol(null.dist)) {
    tot.abund = sum(null.dist[,species])
    null.dist[,species] = 0
    for (individual in 1:tot.abund) {
      sampled.site = sample(c(1:nrow(bbs.sp.site)), 1)
      null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
    }
  }
  
  ##Calculate null deviation for null patches and store
  null.alphas[,randomize] <- apply(null.dist, 2, function(x){sum(ifelse(x > 0, 1, 0))})
  null.gamma[1, randomize] <- sum(ifelse(rowSums(null.dist)>0, 1, 0))
  expected_beta[1, randomize] <- 1 - mean(null.alphas[,randomize]/null.gamma[,randomize])
  null.alpha <- mean(null.alphas[,randomize])
  null.alpha.comp <- c(null.alpha.comp, null.alpha)
  
  bucket_bray <- as.matrix(vegdist(null.dist, "bray"))
  wuf<-(GUniFrac(null.dist, tree, alpha=1)) #als add
  #wuf<-(GUniFrac(comm.t, tree, alpha=1)) #als add test that comparable  values are calculated as with QIIME
  bucket_wuf <- as.matrix(wuf$unifracs[,,"d_1"]) #als add
  diag(bucket_bray) <- NA
  diag(bucket_wuf) <- NA #als add
  bucket_bray_res[,randomize] <- apply(bucket_bray, 2, FUN="mean", na.rm=TRUE)
  bucket_wuf_res[,randomize] <- apply(bucket_wuf, 2, FUN="mean", na.rm=TRUE) #als add
} ## end randomize loop

## Calculate beta-diversity for obs metacommunity
beta_comm_abund <- vegdist(comm.t, "bray")
wuf_comm_abund <- GUniFrac(comm.t, tree, alpha=1) #als add
res_beta_comm_abund <- as.matrix(as.dist(beta_comm_abund))
res_wuf_comm_abund <- as.matrix(as.dist(wuf_comm_abund$unifracs[,,"d_1"])) #als add
diag(res_beta_comm_abund) <- NA
diag(res_wuf_comm_abund) <- NA #als add

# output beta diversity (Bray)
beta_div_abund_stoch <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE)
wuf_div_abund_stoch <- apply(res_wuf_comm_abund, 2, FUN="mean", na.rm=TRUE) #als add

# output abundance beta-null deviation
bray_abund_null_dev <- beta_div_abund_stoch - mean(bucket_bray_res)
wuf_abund_null_dev <- wuf_div_abund_stoch - mean(bucket_wuf_res) #als add


### Outputs:

#beta_div_stoch  - Jaccard beta-diversity for the metacommunity, average value (of all pairwise comparisons) for each patch
#beta_div_abund_stoch - Bray-Curtis beta-diversity for the metacommunity, average value (of all pairwise comparisons) for each patch
#PA_null_dev - presence-absence null deviation values or the metacommunity, average value (of all pairwise comparisons) for each patch
#abund_null_dev - abundance null deviation values or the metacommunity, average value (of all pairwise comparisons) for each patch
###
#END script by Tucker et al.
#######################

#plotting and statistical tests
betanull.out=data.frame(I(beta_div_abund_stoch),I(bray_abund_null_dev),I(wuf_div_abund_stoch),I(wuf_abund_null_dev),I(map[,"SampleID"]),as.character(map[,"Classification"]), as.numeric(map[,"SoilTemperature_to10cm"]), stringsAsFactors=FALSE)
colnames(betanull.out)=c("BRAY_beta_div_abund_stoch", "BRAY_AbundanceNullDeviation", "WUF_div_abund_stoch","WUF_AbundanceNullDeviation","SampleID","Classification", "SoilTemperature_to10cm")
write.table(betanull.out, "bnullout_02nov2016.txt", quote=FALSE, sep="\t")

betanull.out=read.table("bnullout_r1.txt", header=TRUE, sep="\t")
#plottingorder orders samples along a chronosequence and disturbance intensity gradient, by 1) reference samples, 2) fire-affected, sites ranked from hottest to coolest soil temperatures; and 3) recovered sites ranked from hottest to coolest soil temepratures
plottingorder=c(13,15,12,17,14,9,16,1,6,4,11,8,3,7,5,10,2,18)

library("reshape2")
bnull.long=melt(betanull.out, id.vars=c("SampleID", "Classification","SoilTemperature_to10cm"), measure.vars=c("BRAY_AbundanceNullDeviation", "WUF_AbundanceNullDeviation"), col=)

GnYlOrRd=colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)

figX <- ggplot(data=bnull.long, aes(x=Classification, y=as.numeric(value)))+
  geom_boxplot()+ 
  geom_jitter(aes(color=as.numeric(SoilTemperature_to10cm), y=as.numeric(value),cex=1))+
  facet_grid(variable~., scales="free_y")+
  scale_size(guide=FALSE)+
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", guide_legend(title="Temperature (Celsius)"))+
  scale_x_discrete(name="Fire classification", limits=c("Reference", "FireAffected", "Recovered"))+
  scale_y_continuous(name="Abundance Null Deviation")+
  theme_bw(base_size=10)
figX
#ggsave("FigX.eps", width=86, units="mm")

figX2 <- ggplot(data=bnull.long, aes(x=variable, y=as.numeric(value)))+
  geom_boxplot()+ 
  geom_jitter(aes(color=as.numeric(SoilTemperature_to10cm), y=as.numeric(value),cex=1))+
  facet_grid(Classification~., scales="free_y")+
  scale_size(guide=FALSE)+
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", guide_legend(title="Temperature (Celsius)"))+
  scale_x_discrete(name="Beta-null Resemblence", labels=c("Bray-Curtis", "Weighted UniFrac"))+
  scale_y_continuous(name="Abundance Null Deviation")+
  theme_bw(base_size=10)
figX2

bnull.long.bray=bnull.long[bnull.long[,"variable"]=="BRAY_AbundanceNullDeviation",]
figY <- ggplot(data=bnull.long.bray, aes(x=plottingorder, y=as.numeric(value)))+
  geom_point(aes(color=as.numeric(SoilTemperature_to10cm), y=as.numeric(value),cex=1))+
  scale_size(guide=FALSE)+
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", guide_legend(title="Temperature (Celsius)"))+
  scale_x_continuous(name="Disturbance Intensity", breaks=c(1.5,7,15), labels=c("Ref", "FireAffected", "Recovered"))+
  scale_y_continuous(name="Abundance Null Deviation")+
  geom_vline(xintercept=c(2.5,11.5), col="gray", lty="dashed")+
  theme_bw(base_size=10)+
  theme(legend.position="none")
figY


bnull.long.wuf=bnull.long[bnull.long[,"variable"]=="WUF_AbundanceNullDeviation",]
figZ <- ggplot(data=bnull.long.wuf, aes(x=plottingorder, y=as.numeric(value)))+
  geom_point(aes(color=as.numeric(SoilTemperature_to10cm), y=as.numeric(value),cex=1))+
  scale_size(guide=FALSE)+
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", guide_legend(title="Temperature (Celsius)"))+
  scale_x_continuous(name="Disturbance Intensity", breaks=c(1.5,7,15), labels=c("Ref", "FireAffected", "Recovered"))+
  scale_y_continuous(name="Abundance Null Deviation")+
  geom_vline(xintercept=c(2.5,11.5), col="gray", lty="dashed")+
  theme_bw(base_size=10)+
  theme(legend.position="none")
figZ

#Multiplot script written by Winston Chang
load("MiscSourceScripts/multiplot.R")
multiplot(figX, figY, figZ, cols=1)
ggsave("Figures/Fig6AB.eps", width=86, units="mm")

#Pairwise t-tests for Bray Beta Null
t.test(betanull.out[betanull.out[,"Classification"]=="Recovered","BRAY_AbundanceNullDeviation"],betanull.out[betanull.out[,"Classification"]=="FireAffected","BRAY_AbundanceNullDeviation"])
t.test(betanull.out[betanull.out[,"Classification"]=="Recovered","BRAY_AbundanceNullDeviation"],betanull.out[betanull.out[,"Classification"]=="Reference","BRAY_AbundanceNullDeviation"])
t.test(betanull.out[betanull.out[,"Classification"]=="Reference","BRAY_AbundanceNullDeviation"],betanull.out[betanull.out[,"Classification"]=="FireAffected","BRAY_AbundanceNullDeviation"])
#recovered and fire-affected are statistically distinct, p < 0.0006, all other comparisons p > 0.05

#Pairwise t-tests for WUF Beta Null
t.test(betanull.out[betanull.out[,"Classification"]=="Recovered","WUF_AbundanceNullDeviation"],betanull.out[betanull.out[,"Classification"]=="FireAffected","WUF_AbundanceNullDeviation"])
t.test(betanull.out[betanull.out[,"Classification"]=="Recovered","WUF_AbundanceNullDeviation"],betanull.out[betanull.out[,"Classification"]=="Reference","WUF_AbundanceNullDeviation"])
t.test(betanull.out[betanull.out[,"Classification"]=="Reference","WUF_AbundanceNullDeviation"],betanull.out[betanull.out[,"Classification"]=="FireAffected","WUF_AbundanceNullDeviation"])
#recovered and fire-affected are distinct, p < 0.04, all other comparisons p > 0.05

#Are the WUF and Bray beta null correlated?
cor.test(bnull.long.wuf[,"value"], bnull.long.bray[,"value"])
#Pearson's R = 0.71, p = 0.001