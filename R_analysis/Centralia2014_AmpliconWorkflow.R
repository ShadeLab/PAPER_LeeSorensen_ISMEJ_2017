################################
### Code for ecological statistics for
### "Divergent extremes but convergent recovery of bacterial and archaeal soil 
### communities to an ongoing subterranean coal mine fire"
### by SH Lee, JW Sorensen, KL Grady, TC Tobin and A Shade
### Prepared  12 November 2016
### Author: Ashley Shade, Michigan State University; shade.ashley <at> gmail.com
################################
#
# Before you start
# Make sure you are using the latest version of R (and Rstudio)
# The following packages (and their dependencies) are needed to run the whole analysis
# calibrate 1.7.2
# gplots 3.0.1
# ggplot2 2.1.0
# indicspecies 1.7.5
# limma 3.26.9
# mass 7.3-45 (calibrate dependency)
# outliers 0.14
# reshape2 1.4.1
# vegan 2.4-0
# reldist 1.6-6
# bipartite 2.06.1
# GUniFrac 1.0
# ape 3.5
# phangorn 2.0-2
#
################################
### Plotting soil contextual data
################################
#load R libraries for this section
library(ggplot2)
library(reshape2)
library(outliers)

#read in mapping file with soil data
map=read.table("InputFiles/Centralia_Collapsed_Map_forR.txt", header=TRUE, sep="\t")

#plot chemistry v. temperature (Figure 2)
#melt data
map.long=melt(map, id.vars=c("SampleID", "SoilTemperature_to10cm", "Classification"), measure.vars=c("NO3N_ppm","NH4N_ppm","pH","SulfateSulfur_ppm","K_ppm","Ca_ppm","Mg_ppm","OrganicMatter_500","Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per"))

#make a gradient color palette, note bias
GnYlOrRd=colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)

fig2=ggplot(map.long, aes(y=as.numeric(SoilTemperature_to10cm), x=value))+
  
  #add points layer
  geom_point(aes(y=as.numeric(SoilTemperature_to10cm), x=value, shape=Classification, color=as.numeric(SoilTemperature_to10cm)))+
  
  #set facet with 4 columns, make x-axes appropriate for each variable
  facet_wrap(~variable, ncol=4, scales="free_x")+ 
  
  #set gradient for temperature and add gradient colorbar
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", guide_legend(title="Temperature"))+
  
  #omit the legend for the size of the points
  scale_size(guide=FALSE)+
  
  #define the axis labels
  labs(y="Temperature (Celsius)", x="    ")+
  
  #set a simple theme
  theme_bw(base_size=10)

fig2
#ggsave("Figures/Fig2.eps", width=178, units="mm")

##Subset contextual data inclusive of soil quantitative variables
env=map[,c("SoilTemperature_to10cm", "NO3N_ppm", "pH", "K_ppm", "Mg_ppm", "OrganicMatter_500", "NH4N_ppm", "SulfateSulfur_ppm", "Ca_ppm", "Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per","Fire_history")]

##Test for outliers, loop will print all significant outliers and their sampleID - these were not removed from analysis
for (i in 1:ncol(env)){
  x=grubbs.test(env[,i], type=10)
  if(x$p.value < 0.05){
  print(colnames(env)[i])
  print(row.names(env)[env[,i]==max(env[,i])])
    }
  }
#samples 13 (for pH, Ca) and 10 (for NO3N, NH4N,Fe) are common outliers - both have high temps. Sample 3 is also outlier for Mg and OM; this is a recovered site.  Generally this test indicates a lot of variability.

#correlation test between temperature and other soil chemistry
for(i in 1:ncol(env)){
  ct=cor.test(env[,"SoilTemperature_to10cm"],env[,i])
  if (ct$p.value < 0.05){
  print(colnames(env)[i])
  print(ct)
  }
}

#extract means from recovered and reference soils' pH
mean(env[map[,"Classification"]=="Reference","pH"])
mean(env[map[,"Classification"]== "Recovered","pH"])

#plot cell counts and 16S rRNA qPCR data (Supporting Figure 2)
map.long.counts=melt(map, id.vars=c("SampleID", "Classification"), measure.vars=c("rRNA_gene_copies_per_g_dry_soil","CellCounts_per_g_dry_soil"))
labels=c(rRNA_gene_copies_per_g_dry_soil="rRNA gene copies",CellCounts_per_g_dry_soil="Cell counts")

sfig2 <- ggplot(data=map.long.counts, aes(x=Classification, y=value))+
  geom_boxplot() + 
  geom_jitter(aes(shape=Classification))+
  facet_grid(variable~., scales="free_y", labeller=labeller(variable = labels))+
  scale_shape(guide=FALSE)+
  #scale_color_manual(values=colors)+
  scale_x_discrete(name="Fire classification")+
  scale_y_continuous(name="value per g dry soil")+
  theme_bw(base_size=10)
sfig2
#ggsave("Figures/SFig2.eps", width=86, units="mm")

#Pariwise t-tests for cell counts
t.test(map[map[,"Classification"]=="Recovered","CellCounts_per_g_dry_soil"],map[map[,"Classification"]=="FireAffected","CellCounts_per_g_dry_soil"])
t.test(map[map[,"Classification"]=="Recovered","CellCounts_per_g_dry_soil"],map[map[,"Classification"]=="Reference","CellCounts_per_g_dry_soil"])
t.test(map[map[,"Classification"]=="FireAffected","CellCounts_per_g_dry_soil"],map[map[,"Classification"]=="Reference","CellCounts_per_g_dry_soil"])

#Pairwise t-tests for qPCR
t.test(map[map[,"Classification"]=="Recovered","rRNA_gene_copies_per_g_dry_soil"],map[map[,"Classification"]=="FireAffected","rRNA_gene_copies_per_g_dry_soil"])
t.test(map[map[,"Classification"]=="Recovered","rRNA_gene_copies_per_g_dry_soil"],map[map[,"Classification"]=="Reference","rRNA_gene_copies_per_g_dry_soil"])
t.test(map[map[,"Classification"]=="Reference","rRNA_gene_copies_per_g_dry_soil"],map[map[,"Classification"]=="FireAffected","rRNA_gene_copies_per_g_dry_soil"])

################################
### Preparing OTU and distance tables for analysis
################################
#load R libraries for this section
library(ggplot2)
library(reshape2)
library(vegan)

#read in community OTU table, and transpose (rarefied collapsed MASTER table, output from QIIME)
comm=read.table("InputFiles/MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1, check.names=FALSE, sep="\t")

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

### Read in resemblance matrices
#read in weighted unifrac table (output from QIIME)
uf=read.table("InputFiles/weighted_unifrac_MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1)

#sort by rows, columns (so they are in the consecutive order)
uf=uf[order(row.names(uf)),order(colnames(uf))]
uf.d=as.dist(uf)

#read in the unweighted unifrac table (output from QIIME)
uwuf=read.table("InputFiles/unweighted_unifrac_MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1)

#sort by rows, columns (so they are in the consecutive order)
uwuf=uwuf[order(row.names(uwuf)),order(colnames(uwuf))]
uwuf.d=as.dist(uwuf)

#read in the normalized weighted unifrac table (output from QIIME)
nwuf=read.table("InputFiles/weighted_normalized_unifrac_MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1)

#sort by rows, columns (so all tables are in the consecutive order)
nwuf=nwuf[order(row.names(nwuf)),order(colnames(nwuf))]
nwuf.d=as.dist(nwuf)

#assign fire classification
fireclass=map[,"Classification"]

ref.t=comm.t[map$Classification=="Reference",]
rec.t=comm.t[map$Classification=="Recovered",]
fire.t=comm.t[map$Classification=="FireAffected",]

################################
### Calculate and plot within-sample (alpha) diversity
################################
#read in alpha diversity table (output from QIIME)
div=read.table("InputFiles/MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000_alphadiv.txt", header=TRUE)

#sort by sample ID (so that they are in consecutive order)
div=div[order(row.names(div)),]

#calculate pielou's evenness from OTU table
s=specnumber(comm.t)
h=diversity(comm.t,index="shannon")
pielou=h/log(s)

#combine alpha diversity data and fire classification (from map file)
div=cbind(row.names(div),div,pielou, map$Classification)
colnames(div)=c("SampleID", "PD", "Richness", "Pielou", "Classification")
   

#plot (Figure 3)
#reshape the data
div.long=melt(div, id.vars=c("SampleID", "Classification"))

#plot a facet
#comment toggle for color v. bw
colors=c("red", "yellow", "green")

fig3 <- ggplot(data=div.long, aes(x=Classification, y=value))+
  geom_boxplot() + 
  geom_jitter(aes(shape=Classification))+
  #geom_jitter(aes(color=Classification, cex=1.5))+
  facet_grid(variable~., scales="free_y")+
  #scale_shape(guide=FALSE)+
  scale_size(guide=FALSE)+
  scale_color_manual(values=colors)+
  scale_x_discrete(name="Fire classification")+
  scale_y_continuous(name="Diversity value")+
  theme_bw(base_size=10)
fig3
#ggsave("Figures/Fig3.eps", width=86, units="mm")

#ttest
v=c("PD", "Richness", "Pielou")

outdiv=NULL
for(i in 1:length(v)){
  #subset the data to test one phylum at a time
  active=div[div$Classification=="FireAffected",colnames(div)==v[i]]
  recov=div[div$Classification=="Recovered",colnames(div)==v[i]]
  ref=div[div$Classification=="Reference",colnames(div)==v[i]]
  
  #perform the test
  test1=t.test(active, recov, paired=FALSE, var.equal = FALSE)
  test2=t.test(active, ref, paired=FALSE, var.equal = FALSE)
  test3=t.test(ref, recov, paired=FALSE, var.equal = FALSE)
  test1.out=c(v[i],"ActivevRecov",test1$statistic, test1$parameter, test1$p.value)
  test2.out=c(v[i],"ActivevRef",test2$statistic, test2$parameter, test2$p.value)
  test3.out=c(v[i],"RefvRecov",test3$statistic, test3$parameter, test3$p.value)
  outdiv=rbind(outdiv, test1.out, test2.out, test3.out)
  
}
outdiv

################################
### Analysis of technical replicates
################################
#Supporting Table 2 - assessing reproducibility among technical replicates
techdiv=read.table("InputFiles/OTU_hdf5_filteredfailedalignments_rdp_rmCM_even53000_alphadiv.txt") #output from QIIME
techdiv.out=NULL
sampleIDs=c("C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18")
for(i in 1:length(sampleIDs)){
  temp=techdiv[grep(sampleIDs[i], row.names(techdiv)),]
  temp2=c(mapply(mean,temp), mapply(sd,temp))
  techdiv.out=rbind(techdiv.out,temp2)
}
row.names(techdiv.out)=sampleIDs
colnames(techdiv.out)=c("PD_mean", "Richness_mean", "PD_sd", "Richness_sd")
#write.table(techdiv.out, "Results/AlphaDiv_TechnicalReps.txt", quote=FALSE, sep="\t")

#Supporting PCoA (SFig 1)- assessing reproducibility among technical replicates
beta <- read.table("InputFiles/weighted_unifrac_OTU_hdf5_filteredfailedalignments_rdp_rmCM_even53000.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE, row.names=1)

map.f<- read.table("InputFiles/Centralia_Full_Map.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE, row.names=1)
beta <- beta[order(row.names(beta)),order(colnames(beta))]
### Remove Mock
beta <- beta[-55,-55]
library(vegan)
beta.pcoa<- cmdscale(beta, eig=TRUE)
ax1.v.f=beta.pcoa$eig[1]/sum(beta.pcoa$eig)
ax2.v.f=beta.pcoa$eig[2]/sum(beta.pcoa$eig)

coordinates <- as.data.frame(beta.pcoa$points)
Samples <- map$Sample
coordinates$Sample<- map.f$Sample
coordinates_avg_sd <- NULL
for (i in 1:length(Samples)){
  Site <- coordinates[coordinates$Sample==Samples[i],]
  AX1 <- c(mean(Site[,1]),sd(Site[,1]))
  AX2 <- c(mean(Site[,2]),sd(Site[,2]))
  coordinates_avg_sd<- rbind(coordinates_avg_sd,c(AX1,AX2))
  
}
row.names(coordinates_avg_sd)<-Samples

unique(map$Classification)
Class=rep('black',nrow(map))
Class[map$Classification=="FireAffected"]='red'
Class[map$Classification=="Reference"]='green'
Class[map$Classification=="Recovered"]='yellow'

library(calibrate)
#SFig 1
dev.off()
setEPS()
postscript("Figures/SFig1.eps", width = 6, height=6, pointsize=8,paper="special")
plot(coordinates_avg_sd[,1],coordinates_avg_sd[,3] ,cex=1.5,pch=21,bg=Class,main="Averaged Technical Replicates Weighted UniFrac PCoA",xlab= paste("PCoA1: ",100*round(ax1.v.f,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(ax2.v.f,3),"% var. explained",sep=""))
textxy(X=coordinates_avg_sd[,1], Y=coordinates_avg_sd[,3],labs=map$Sample, cex=1)
arrows(coordinates_avg_sd[,1], coordinates_avg_sd[,3]- coordinates_avg_sd[,4], coordinates_avg_sd[,1], coordinates_avg_sd[,3]+ coordinates_avg_sd[,4], length=0.05, angle=90, code=3)
arrows(coordinates_avg_sd[,1]- coordinates_avg_sd[,2], coordinates_avg_sd[,3], coordinates_avg_sd[,1] + coordinates_avg_sd[,2], coordinates_avg_sd[,3], length=0.05, angle=90, code=3)
dev.off()

################################
### Phylum-level responses to fire
################################
#load R libraries for this section
library(ggplot2)

#read in phylum level OTU table (QIIME output)
comm.phylum=read.table("InputFiles/MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000_L2.txt", sep="\t", header=TRUE, row.names=1) #output from QIIME

##sort by sample ID (so that they are in consecutive order)
comm.phylum=comm.phylum[,order(colnames(comm.phylum))]

#combine phyla that contribute less than 0.01 each
below01=comm.phylum[rowSums(comm.phylum)<0.01,]
below01.cs=colSums(below01)

#remove those below01 phyla from the table
comm.phylum=comm.phylum[rowSums(comm.phylum)>0.01,]
#add the summary from the <0.01
comm.phylum=rbind(comm.phylum,below01.cs)
#rename the last row
row.names(comm.phylum)[nrow(comm.phylum)]="Below_0.01"

#for character string trucation in R : http://stackoverflow.com/questions/10883605/truncating-the-end-of-a-string-in-r-after-a-character-that-can-be-present-zero-o
phylumnames=sub(".*p__", "", row.names(comm.phylum))
row.names(comm.phylum)=phylumnames

#assign fire classifications to samples
fireclass=map[,"Classification"]

p.active=comm.phylum[,fireclass=="FireAffected"]
p.recov=comm.phylum[,fireclass=="Recovered"]
p.ref=comm.phylum[,fireclass=="Reference"]

#Calculate a mean phylum rel. abundance across all of the samples that are within each activity group
m.active=apply(p.active,1,mean)
m.recov=apply(p.recov,1,mean)
m.ref=apply(p.ref,1,mean)

m.summary.p=cbind(m.active, m.recov,m.ref)
colnames(m.summary.p)=c("FireAffected", "Recovered", "Reference")

#sort in decreasing total abundance order
m.summary.p=m.summary.p[order(rowSums(m.summary.p),decreasing=TRUE),]

#plot (Figure 5)
m.summary.p.long=melt(m.summary.p, id.vars=row.names(m.summary.p),measure.vars=c("FireAffected", "Recovered", "Reference"))
colors=c("red", "yellow", "green")

fig5=ggplot(m.summary.p.long, aes(x=Var1, y=value, fill=Var2))+
  geom_dotplot(binaxis="y", dotsize = 3)+
  facet_grid(Var2~.)+
  scale_fill_manual(values=colors, guide=FALSE)+
  labs(x="Phylum", y="Mean relative abundance", las=1)+
  theme(axis.text.x = element_text(angle = 90, size = 10, face = "italic"))
fig5
#ggsave("Figures/Fig5.eps", width=178, units="mm")

#Welch's t-test for all phyla
u=row.names(comm.phylum)
out=NULL
for(i in 1:length(u)){
  
  #subset the data to test one phylum at a time
  active=comm.phylum[row.names(comm.phylum)==u[i],fireclass=="FireAffected"]
  recov=comm.phylum[row.names(comm.phylum)==u[i],fireclass=="Recovered"]
  
  #perform the test
  test=t.test(active, recov, paired=FALSE, var.equal = FALSE)
  test.out=c(row.names(comm.phylum)[i],test$statistic, test$parameter, test$p.value)
  out=rbind(out,test.out)
  
}
colnames(out)=c("Phylum", "Tstatistic", "DF", "pvalue")
#all results: Supporting Table 8
out

#extract overrepresented  in fire
out[out[,"pvalue"]<0.05 & out[,"Tstatistic"]>0,]

#extracted overrepresented in recovered
out[out[,"pvalue"]<0.05 & out[,"Tstatistic"]<0,]
#write.table(out, "Results/Phylum_ttest.txt",quote=FALSE, sep="\t")


################################
### Comparative (beta) diversity  
################################
#load R libraries for this section
library(calibrate)
library(ggplot2)
library(vegan)

# use weighted unifrac
uf.pcoa=cmdscale(uf.d, eig=TRUE)
#calculate percent variance explained, then add to plot
ax1.v=uf.pcoa$eig[1]/sum(uf.pcoa$eig)
ax2.v=uf.pcoa$eig[2]/sum(uf.pcoa$eig)

envEF=envfit(uf.pcoa, env)
#Supporting Table 4
envEF

unique(map$Classification)
Class=rep('black',nrow(map))
Class[map$Classification=="FireAffected"]='red'
Class[map$Classification=="Reference"]='green'
Class[map$Classification=="Recovered"]='yellow'

#export figure 4
#textxy is from the calibrate library
dev.off()
setEPS()
postscript("Figures/Fig4.eps", width = 3.385, height=3.385, pointsize=8,paper="special")
plot(uf.pcoa$points[,1],uf.pcoa$points[,2] ,cex=1.5,pch=21,bg=Class,main="Weighted UniFrac PCoA", xlab= paste("PCoA1: ",100*round(ax1.v,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100*round(ax2.v,3),"%var. explained",sep=""))
textxy(X=uf.pcoa$points[,1], Y=uf.pcoa$points[,2],labs=map$SampleID, cex=0.8)
legend('bottomleft',c('Fire Affected','Recovered','Reference'),pch=21,pt.bg=c("red", "yellow", "green"),lty=0)
plot(envEF, p.max=0.10, col="black", cex=1)
dev.off()


#perform hypothesis testing on fire-affected v. recovered+reference sites
#permanova
Class2=sub("green", "yellow", Class)
a=adonis(uf.d~Class2, distance=TRUE, permutations=1000)
a

#multivariate dispersion with Tukey HSD
b=betadisper(uf.d, group=Class2)
TukeyHSD(b, which = "group", ordered = FALSE,conf.level = 0.95)

#mantel w/ spatial distances
space=read.table("InputFiles/spatialdistancematrix.txt", header=TRUE, row.names=1)
space.d=as.dist(space)
mantel(uf.d,space.d)


################################
### Do different resemblances agree in their overarching patterns?
################################
#Supporting Table 3A:  
#the variance explained by each distance (taxonomic/phylogenetic and weighted/unweighted)
bc.d=vegdist(t(comm), method="bray")
sor.d=vegdist(t(comm), method="bray",binary=TRUE)

# PCoA using unweighted unifrac  (QIIME output - unweighted phylogenetic)
uwuf.pcoa=cmdscale(uwuf.d, eig=TRUE)
#calculate percent variance explained, then add to plot
ax1.v.uwuf=uwuf.pcoa$eig[1]/sum(uwuf.pcoa$eig)
ax2.v.uwuf=uwuf.pcoa$eig[2]/sum(uwuf.pcoa$eig)

# PCoA using  bray-curtis  (vegan output - weighted taxonomic)
bc.pcoa=cmdscale(bc.d, eig=TRUE)
#calculate percent variance explained, then add to plot
ax1.v.bc=bc.pcoa$eig[1]/sum(bc.pcoa$eig)
ax2.v.bc=bc.pcoa$eig[2]/sum(bc.pcoa$eig)

#PCoA using sorensen (vegan output - unweighted taxonomic)
sor.pcoa=cmdscale(sor.d, eig=TRUE)
#calculate percent variance explained, then add to plot
ax1.v.sor=sor.pcoa$eig[1]/sum(sor.pcoa$eig)
ax2.v.sor=sor.pcoa$eig[2]/sum(sor.pcoa$eig)

#Mantel and PROTEST tests between all resemblances (Supporting Table 3B)
resem=list(uf.d,uwuf.d,nwuf.d,bc.d,sor.d)
#this loops a bit funny but all pairwise results are available
names=c("weighted_UniFrac", "unweighted_UniFrac", "normalized_weighted_UniFrac", "BrayCurtis", "Sorenson")
m.out=NULL
for (i in 1:length(resem)){
  dist1=resem[[i]]
  print(i)
  
  j=i+1
  
  for(j in 2:length(resem)){
    dist2=resem[[j]]
    print(j)
    
    #Mantel
    m=mantel(dist1,dist2)
    
    #Protest
    pr=protest(dist1,dist2)
    
    #results out
    m.v=c(names[i], names[j],m$statistic, m$signif, pr$t0, pr$ss, pr$signif)  
    m.out=rbind(m.out,m.v)
  }
}

#Supporting Table 3B
colnames(m.out)=c("Dist1", "Dist2", "Mantel_R", "Mantel_p", "PROTEST_R", "PROTEST_m12", "PROTEST_p")
m.out
#write.table(m.out, "Results/MantelDist.txt", quote=FALSE, sep="\t")

################################
### Comparative diversity of fire-affected samples 
################################
#load R libraries for this section
library(vegan)

#reduce uf to fire only
uf.fire=uf[map$Classification=="FireAffected",map$Classification=="FireAffected"]
uf.fire.d=as.dist(uf.fire)
env.fire=env[map$Classification=="FireAffected",]
labels=map[map$Classification=="FireAffected","SampleID"]

#PCoA for fire sites only
uf.fire.pcoa=cmdscale(uf.fire.d, eig=TRUE)
#fit environmental variables
envFIT.fire=envfit(uf.fire.pcoa, env=env.fire)
#print results to screen (Supporting Table 5)
envFIT.fire

#df <- data.frame((envFIT.fire$vectors)$arrows, (envFIT.fire$vectors)$r, (envFIT.fire$vectors)$pvals)
#write.table(df, "Results/ENV_Fire.txt", quote=FALSE, sep="\t")

#calculate %var. explained by each axis
ax1.v.f=uf.fire.pcoa$eig[1]/sum(uf.fire.pcoa$eig)
ax2.v.f=uf.fire.pcoa$eig[2]/sum(uf.fire.pcoa$eig)

#CAP for fire-sites, constrained by temperature
#to determine  explanatory value of abiotic factors for fire-affected sites, after temp is accounted for
#make vector of temperature only
temp=env.fire[,"SoilTemperature_to10cm"]
#CAP
cap1=capscale(uf.fire.d~Condition(temp))
#fit environmental variables
c.ef=envfit(cap1, env.fire)
#print results to screen (Supporting Table 6)
c.ef

#df <- data.frame((c.ef$vectors)$arrows, (c.ef$vectors)$r, (c.ef$vectors)$pvals)
#write.table(df, "Results/CAP.txt", quote=FALSE, sep="\t")

#calculate % var. explained by each axis
ax1.v.f.t=cap1$CA$eig[1]/sum(cap1$CA$eig)
ax2.v.f.t=cap1$CA$eig[2]/sum(cap1$CA$eig)

#Plot:  supporting Figure 4
setEPS()
postscript("Figures/SFig4AB.eps", width = 6.770, height=3.385, pointsize=8,paper="special")
par(mfrow=c(1,2))
plot(uf.fire.pcoa$points[,1],uf.fire.pcoa$points[,2], main= "(A) Fire-affected soils PCoA", type="n",xlab=paste("PCoA1: ",100*round(ax1.v.f,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100*round(ax2.v.f,3),"% var. explained",sep=""))
textxy(X=uf.fire.pcoa$points[,1], Y=uf.fire.pcoa$points[,2],labs=labels, offset=0, cex=0.8)
plot(envFIT.fire, p=0.10)
plot(cap1, cex=0.9,main = "(B) Temperature-constrained \nfire-affected soils PCoA", xlab=paste("CAP Ax1: ",100*round(ax1.v.f.t,3),"%var. explained",sep=""), ylab=paste("CAP Ax2: ",100*round(ax2.v.f.t,3),"%var. explained",sep=""))
plot(c.ef, p= 0.10)
dev.off()

################################
### Sloan neutral model 
################################
#NOTE:  must use full dataset (including singleton OTUs) for this analysis 
#Source for model fits is from Burns et al. ISMEJ 2015, downloaded R code from their supporting materials
#Source code requires:  minpack.lm, Hmisc, stats4 packages - make sure they are installed (and their dependencies)
source("MiscSourceScripts/sncm.fit_function.r")

#assign variables for function
spp=t(comm.sigs)
taxon=as.vector(rdp)

ref.t.sigs=spp[map$Classification=="Reference",]
rec.t.sigs=spp[map$Classification=="Recovered",]
rec.t.sigs.NZ<- rec.t.sigs[,colSums(rec.t.sigs)>0]

fire.t.sigs=spp[map$Classification=="FireAffected",]
fire.t.sigs.NZ<-fire.t.sigs[,colSums(fire.t.sigs)>0]

#Models for the whole community
obs.np=sncm.fit(spp,taxon=rdp, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp,taxon=rdp, stats=TRUE, pool=NULL)

#Models for each classification
#fire affected: total - asks the question:  in itself, do the fire-affected sites follow neutral
obs.fireT=sncm.fit(fire.t.sigs.NZ,taxon=rdp, stats=FALSE, pool=NULL)
sta.fireT=sncm.fit(fire.t.sigs.NZ,taxon=rdp, stats=TRUE, pool=NULL)

#recovered : total - asks the question:  do recovered sites follow neutral expectations?
obs.recT=sncm.fit(rec.t.sigs.NZ,taxon=rdp, stats=FALSE, pool=NULL)
sta.recT=sncm.fit(rec.t.sigs.NZ,taxon=rdp, stats=TRUE, pool=NULL)

results=rbind(sta.np, sta.fireT, sta.recT)
row.names(results)=c("all", "Fire_Affected", "Recovered")


#par(mfrow=c(2,3)) #for plotting in R studio w/out export
l1=list(obs.np, obs.recT, obs.fireT)
l2=list(sta.np, sta.recT, sta.fireT)
names=c("(A) All", "(B) Recovered", "(C) Fire_Affected")
out.sta=NULL

#Plot supporting Fig 5 panels
for(i in 1:length(l1)){
  #define data
  temp=as.data.frame(l1[i])
  sta=as.data.frame(l2[i])
  
  #how many taxa are above their prediction, and below?
  above.pred=sum(temp$freq > (temp$pred.upr), na.rm=TRUE)/sta$Richness
  below.pred=sum(temp$freq < (temp$pred.lwr), na.rm=TRUE)/sta$Richness
  
  out=c(above.pred, below.pred)
  
  ap= temp$freq > (temp$pred.upr)
  bp= temp$freq < (temp$pred.lwr)
  
  #plot figure (SFig5)
  setEPS()
  if(i == 1){
  postscript("Figures/SFig5A.eps", width = 2.33, height=3, pointsize=10,paper="special")
  }
  if (i == 2){
  postscript("Figures/SFig5B.eps", width = 2.33, height=3, pointsize=10,paper="special")
  }
  if (i ==3){
  postscript("Figures/SFig5C.eps", width = 2.33, height=3, pointsize=10,paper="special")
  }
  
  plot(x=log(temp$p), y=temp$freq, main=names[i], xlab="Log Abundance", ylab="Occurrence Frequency")
  points(x=log(temp$p[ap==TRUE]), y=temp$freq[ap==TRUE], col="red", pch=19)
  points(x=log(temp$p[bp==TRUE]), y=temp$freq[bp==TRUE], col="blue", pch=19)
  lines(temp$freq.pred~log(temp$p), col="yellow", lty=1, lwd=6)
  lines(temp$pred.upr~log(temp$p), col="yellow", lty=1, lwd=3)
  lines(temp$pred.lwr~log(temp$p), col="yellow", lty=1, lwd=3)
  
  dev.off()
  
  out.sta=rbind(out.sta, out)
}

colnames(out.sta)=c("%AbovePred", "%BelowPred")
#Supporting Table 7
results=cbind(results, out.sta)
results
#write.table(results, "Results/SloanNeutralModel.txt", quote=FALSE, sep="\t")

################################
### Beta null models
################################
#MODIFIED by als to use our dataset (comm.t) instead of "dune" and to only include the abundance-based model. We also changed the number of patches to by 18 to match with the dataset.
#ORIGINAL scripts available in the appendix of the work below, published in Oikos (Appendix oik.02803, also R_analysis/oik-02803-appendix-to-Tucker2016/)
#Note that beta null models with weighted UniFrac require ~75 hours walltime to complete with 4Gb memory and 1 processing node; beta-null models with Bray-Curtis only require ~30 hours

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
#write.table(betanull.out, "Results/bnullout_r1.txt", quote=FALSE, sep="\t")
#betanull.out=read.table("Results/bnullout_r1.txt", header=TRUE, sep="\t")

##plottingorder orders samples along a chronosequence and disturbance intensity gradient, by 1) reference samples, 2) fire-affected, sites ranked from hottest to coolest soil temperatures; and 3) recovered sites ranked from hottest to coolest soil temepratures
plottingorder=c(13,15,12,17,14,9,16,1,6,4,11,8,3,7,5,10,2,18)

library("reshape2")
bnull.long=melt(betanull.out, id.vars=c("SampleID", "Classification","SoilTemperature_to10cm"), measure.vars=c("BRAY_AbundanceNullDeviation", "WUF_AbundanceNullDeviation"), col=)

GnYlOrRd=colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)

fig6A <- ggplot(data=bnull.long, aes(x=Classification, y=as.numeric(value)))+
  geom_boxplot()+ 
  geom_jitter(aes(color=as.numeric(SoilTemperature_to10cm), y=as.numeric(value)))+
  facet_grid(variable~., scales="free_y")+
  scale_size(guide=FALSE)+
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", guide_legend(title="Temp"))+
  scale_x_discrete(name="Fire classification", limits=c("Reference", "FireAffected", "Recovered"))+
  scale_y_continuous(name="Abundance Null Deviation")+
  theme_bw(base_size=10)
fig6A


bnull.long.bray=bnull.long[bnull.long[,"variable"]=="BRAY_AbundanceNullDeviation",]
fig6B <- ggplot(data=bnull.long.bray, aes(x=plottingorder, y=as.numeric(value)))+
  geom_point(aes(color=as.numeric(SoilTemperature_to10cm), y=as.numeric(value)))+
  scale_size(guide=FALSE)+
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", guide_legend(title="Temperature (Celsius)"))+
  scale_x_continuous(name="Disturbance Intensity", breaks=c(1.5,7,15), labels=c("Ref", "FireAffected", "Recovered"))+
  scale_y_continuous(name="Abundance Null Deviation")+
  geom_vline(xintercept=c(2.5,11.5), col="gray", lty="dashed")+
  theme_bw(base_size=10)+
  theme(legend.position="none")
fig6B


bnull.long.wuf=bnull.long[bnull.long[,"variable"]=="WUF_AbundanceNullDeviation",]
fig6C <- ggplot(data=bnull.long.wuf, aes(x=plottingorder, y=as.numeric(value)))+
  geom_point(aes(color=as.numeric(SoilTemperature_to10cm), y=as.numeric(value)))+
  scale_size(guide=FALSE)+
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", guide_legend(title="Temperature (Celsius)"))+
  scale_x_continuous(name="Disturbance Intensity", breaks=c(1.5,7,15), labels=c("Ref", "FireAffected", "Recovered"))+
  scale_y_continuous(name="Abundance Null Deviation")+
  geom_vline(xintercept=c(2.5,11.5), col="gray", lty="dashed")+
  theme_bw(base_size=10)+
  theme(legend.position="none")
fig6C

#Multiplot script written by Winston Chang
source("MiscSourceScripts/multiplot.R")

dev.off()
setEPS()
postscript("Figures/Fig6ABC.eps", width = 3.385, height=5, pointsize=9,paper="special")
multiplot(fig6A, fig6B, fig6C, cols=1)
dev.off()

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


################################
### Dominant taxa analysis 
################################

#Extract cumulative most abundant OTUs from fire-affected soils - for Table 1
dim(fire)
fire.ordered=fire[order(rowSums(fire),decreasing=TRUE),]
perc=rowSums(fire.ordered)/sum(rowSums(fire.ordered))

#Analysis of the top 10 most prevalent taxa in fire-affected and recovered soils
#libraries needed for this
library(vegan)
library(gplots)

#Do hot soils have consistent dominant membership?
fire=t(fire.t)
fire.new=fire[rowSums(fire)>0,]
rdp.fire=as.vector(rdp.nosigs[rowSums(fire)>0])
dim(fire.new)

rec=t(rec.t)
rec.new=rec[rowSums(rec)>0,]
rdp.rec=as.vector(rdp.nosigs[rowSums(rec)>0])
dim(rec.new)

#Function to provide the OTU numbers and Taxonomic IDs are the top (default=10) in each site.
extractdominant.f<-function(data,rdp,top.no=10){
  out1=NULL
  out2=NULL
  for(i in 1:ncol(data)){
    s=sort(data[,i], decreasing=TRUE, index.return=TRUE)
    otuIDs=names(s$x[1:top.no])
    rdp.out=rdp[s$ix[1:top.no]]
    sampleID=c(rep(colnames(data)[[i]],top.no))
    temp=cbind(sampleID,otuIDs)
    out1=rbind(out1,temp)
    out2=cbind(out2,rdp.out)
  }
  colnames(out2)=colnames(data)
  #write.table(out2, paste("Results/rdp_",top.no,".txt",sep=""), quote=FALSE, sep="\t")
  #who are the top-10 ranked
  u=unique(out1[,2])
  l=length(unique(out1[,2]))
  actual.prop=l/dim(out1)[[1]]
  expected.prop=top.no/dim(out1)[[1]]
  print("Unique OTU IDs within the most abundant")
  print(u)
  print("Number of unique OTUs within the most abundant")
  print(l)
  print("Redundancy index given the number of samples and the top number selected 1.00 means completely nonredundant, every top taxa was observed only 1 time across all samples")
  print(actual.prop)
  print("Expected redundancy index")
  print(expected.prop)
  #print("List of top taxa by sample")
  #print(out2)
  
  return(out2)
}

fire.out=extractdominant.f(fire.new,rdp.fire,10)
rec.out=extractdominant.f(rec.new,rdp.rec,10)

data=NULL
data=fire.new
top.no=10
rdp.in=rdp.fire

subsettop.f=function(data, top.no, rdp.in){
  otuIDs=NULL
  rdpIDs=NULL
  for(i in 1:ncol(data)){
    s=sort(data[,i], decreasing=TRUE, index.return=TRUE)
    otuIDs=c(otuIDs, names(s$x[1:top.no]))
    rdpIDs=c(rdpIDs, rdp.in[s$ix[1:top.no]])
  }
  temp=cbind(otuIDs,rdpIDs)
  #print(temp)
  u.top=unique(otuIDs)
  #temp.u=temp[is.element(temp[,"otuIDs"],u.top),]
  #write.table(temp.u, "Results/OTURDP_Top10.txt", sep="\t", quote=FALSE)
  top10.otu=NULL
  for(j in 1:nrow(data)){
    if(is.element(row.names(data)[j],u.top)){
      top10.otu=rbind(top10.otu,data[j,])
    }
  }
  row.names(top10.otu)=u.top
  colnames(top10.otu)=colnames(data)
  return(top10.otu)
}

topfire=subsettop.f(fire.new,10,rdp.fire)
#how many OTUs are de novo?
length(grep("dn",rownames(topfire)))

#create color pallette; see: http://colorbrewer2.org/ 
hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")
topfire.pa=1*(topfire>0)
sum(rowSums(topfire.pa)==9)

toprec=subsettop.f(rec.new,10, rdp.rec)
#how many OTUs are de novo
length(grep("dn",rownames(toprec)))

#Figure 7
dev.off()
setEPS()
postscript("Figures/Fig7A.eps", width = 3.5, height=7, pointsize=10, paper="special")
heatmap.2(topfire,col=hc(100),scale="column",key=TRUE,symkey=FALSE, trace="none", density.info="none",dendrogram="both", margins=c(5,13), srtCol=90)
dev.off()

setEPS()
postscript("Figures/Fig7B.eps", width = 3.5, height=7, pointsize=10, paper="special")
heatmap.2(toprec,col=hc(100),scale="column",key=TRUE,symkey=FALSE, trace="none", density.info="none",dendrogram="both", margins=c(5,13), srtCol=90)
dev.off()
