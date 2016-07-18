library(vegan)
library(ggplot2)
library(reshape2)

################################
###Plotting soil contextual data - FINISHED
#read in mapping file with soil data
map=read.table("Centralia_Collapsed_Map_forR.txt", header=TRUE, sep="\t")

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
ggsave("Fig2.eps", width=178, units="mm")

##make a dataset just of soil quantitative variables
env=map[,c("SoilTemperature_to10cm", "NO3N_ppm", "pH", "K_ppm", "Mg_ppm", "OrganicMatter_500", "NH4N_ppm", "SulfateSulfur_ppm", "Ca_ppm", "Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per","Fire_history")]

##Test for outliers, loop will print all significant outliers and their sampleID - these were not removed from analysis
library(outliers)
for (i in 1:ncol(env)){
  x=grubbs.test(env[,i], type=10)
  if(x$p.value < 0.05){
  print(colnames(env)[i])
  print(row.names(env)[env[,i]==max(env[,i])])
    }
  }
#samples 13 (for pH, Ca) and 10 (for NO3N, NH4N,Fe) are common outliers - both have high temps. Sample 3 is also outlier for Mg and OM; this is a recovered site.  Generally this test indicates a lot of variability.

#correlation test between temperature ad other soil chemistry
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

#plot cell counts and 16S rRNA qPCR data (S Figure 2)
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
ggsave("SFig2.eps", width=86, units="mm")


#####################################################
###Preparing the OTU and distance tables for analysis - FINISHED
#read in community OTU table, and transpose (rarefied collapsed table, output from QIIME)
comm=read.table("MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1, check.names=FALSE, sep="\t")

#remove consensus lineage from otu table
rdp=comm[,"ConsensusLineage"]

comm=comm[,-ncol(comm)]


#sort community by colnames (to be in the consistent, consecutive order for all analyses)
comm=comm[,order(colnames(comm))]

#who are the singleton OTUs?
singletonOTUs=row.names(comm)[rowSums(comm)==1]
length(singletonOTUs)
#total 1439 singleton OTUs
g=grep("_dn", singletonOTUs)
#1267 de novo OTUs are singletons

#who are the remaining de novo OTUs?
g=grep("_dn_",row.names(comm))
dn=rdp[g]

rdp.nosigs=rdp[rowSums(comm)>1]

#designate a full dataset
comm.sigs=comm

#remove OTUs with an abundance = 1, across the entire dataset (singleton OTUs)
comm=comm[rowSums(comm)>1,]

#transpose matrix
comm.t=t(comm)

#read in weighted unifrac table
uf=read.table("weighted_unifrac_MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1)

#sort by rows, columns (so they are in the consecutive order)
uf=uf[order(row.names(uf)),order(colnames(uf))]
uf.d=as.dist(uf)

#assign fire classification
fireclass=map[,"Classification"]

ref.t=comm.t[map$Classification=="Reference",]
rec.t=comm.t[map$Classification=="Recovered",]
fire.t=comm.t[map$Classification=="FireAffected",]

###
###Calculate and plot alpha diversity - FINISHED
#read in alpha diversity table (output from QIIME)
div=read.table("MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000_alphadiv.txt", header=TRUE)

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
  geom_jitter()+
  #geom_jitter(aes(shape=Classification))+
  #geom_jitter(aes(color=Classification, cex=2))+
  facet_grid(variable~., scales="free_y")+
  #scale_shape(guide=FALSE)+
  scale_size(guide=FALSE)+
  #scale_color_manual(values=colors)+
  scale_x_discrete(name="Fire classification")+
  scale_y_continuous(name="Diversity value")+
  theme_bw(base_size=10)
fig3
ggsave("Fig3.eps", width=86, units="mm")

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


#########################
###Phylum-level responses  - FINISHED
#read in phylum level OTU table (QIIME output)
comm.phylum=read.table("MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000_L2.txt", sep="\t", header=TRUE, row.names=1)

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
row.names(comm.phylum)[30]="Below_0.01"

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

#plot (Figure 5A)
m.summary.p.long=melt(m.summary.p, id.vars=row.names(m.summary.p),measure.vars=c("FireAffected", "Recovered", "Reference"))
colors=c("red", "yellow", "green")

fig5A=ggplot(m.summary.p.long, aes(x=Var1, y=value, fill=Var2))+
  geom_dotplot(binaxis="y", dotsize = 3)+
  facet_grid(Var2~.)+
  scale_fill_manual(values=colors, guide=FALSE)+
  labs(x="Phylum", y="Mean relative abundance", las=1)+
  theme(axis.text.x = element_text(angle = 90, size = 10, face = "italic"))
fig5A
ggsave("Fig5A.eps", width=178, units="mm")

#Welch's t-test for all phyla (Supporting Table X)
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
#all results

out
#extract overrepresented  in fire
out[out[,"pvalue"]<0.05 & out[,"Tstatistic"]>0,]

#xtracted overrepresented in recovered
out[out[,"pvalue"]<0.05 & out[,"Tstatistic"]<0,]


#heatmap of phylum-level responses, standardized (z-score) based on occurrences (oc) - across rows
library(gplots)

#standardize responses
comm.phylum.oc=decostand(comm.phylum, method="standardize", margin=1)
comm.phylum.oc=as.matrix(comm.phylum.oc)


#create color pallette; see: http://colorbrewer2.org/ 
hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")

#more info on heatmaps:  http://www.inside-r.org/packages/cran/gplots/docs/heatmap.2
#dev.off()
fig5B<-heatmap.2(comm.phylum.oc,col=hc(100),scale="row",key=TRUE,symkey=FALSE, trace="none", density.info="none",dendrogram="both", labRow=phylumnames, colCol=c(map$Class_Color), margins=c(5,13), srtCol=90)

#export heatmap 
dev.off()
#178 mm is 7 inches
setEPS()
postscript("Fig5B.eps", width = 7, height=7, pointsize=10, paper="special")
fig5B<-heatmap.2(comm.phylum.oc,
            col=hc(100),
            scale="row",
            key=TRUE,
            symkey=FALSE, 
            trace="none", 
            density.info="none",
            dendrogram="both", 
            labRow=phylumnames, 
            colCol=c(map$Class_Color), 
            margins=c(5,13), 
            srtCol=90)
dev.off()

################################
###Comparative diversity analyses - FINISHED 
# use weighted unifrac table (QIIME output)
library(calibrate)
uf.pcoa=cmdscale(uf.d, eig=TRUE)
#calculate percent variance explained, then add to plot
ax1.v=uf.pcoa$eig[1]/sum(uf.pcoa$eig)
ax2.v=uf.pcoa$eig[2]/sum(uf.pcoa$eig)

envEF=envfit(uf.pcoa, env)
#Supporting Table
envEF

#plot (fig 4)
unique(map$Classification)
Class=rep('black',nrow(map))
Class[map$Classification=="FireAffected"]='red'
Class[map$Classification=="Reference"]='green'
Class[map$Classification=="Recovered"]='yellow'

plot(uf.pcoa$points[,1],uf.pcoa$points[,2] ,cex=1.5,pch=21,bg=Class,main="Weighted UniFrac PCoA",xlab= paste("PCoA1: ",round(ax1.v,3)," var. explained",sep=""), ylab= paste("PCoA2: ",round(ax2.v,3)," var. explained",sep=""))
textxy(X=uf.pcoa$points[,1], Y=uf.pcoa$points[,2],labs=map$SampleID, cex=1)
legend('bottomleft',c('Fire Affected','Recovered','Reference'),pch=21,pt.bg=c("red", "yellow", "green"),lty=0)

#Add env vectors to plot that are significant p < 0.10
plot(envEF, p.max=0.10, col="black")

#export figure 4
dev.off()
setEPS()
postscript("Fig4.eps", width = 3.385, height=3.385, pointsize=8,paper="special")
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
space=read.table("spatialdistancematrix.txt", header=TRUE, row.names=1)
space.d=as.dist(space)
mantel(uf.d,space.d)

#####################
#constrained PCoA (CAP) _ FINISHED
#to ask about explanatory value of abiotic factors for fire-affected sites, after temp is accounted for

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
#print results to screen
envFIT.fire

#df <- data.frame((envFIT.fire$vectors)$arrows, (envFIT.fire$vectors)$r, (envFIT.fire$vectors)$pvals)
#write.table(df, "ENV_Fire.txt", quote=FALSE, sep="\t")

#calculate %var. explained by each axis
ax1.v.f=uf.fire.pcoa$eig[1]/sum(uf.fire.pcoa$eig)
ax2.v.f=uf.fire.pcoa$eig[2]/sum(uf.fire.pcoa$eig)

#CAP for fire-sites, constrained by temperature
#make vector of temperature only
temp=env.fire[,"SoilTemperature_to10cm"]
#CAP
cap1=capscale(uf.fire.d~Condition(temp))
#fit environmental variables
c.ef=envfit(cap1, env.fire)
#print results to screen
c.ef

#df <- data.frame((c.ef$vectors)$arrows, (c.ef$vectors)$r, (c.ef$vectors)$pvals)
#write.table(df, "CAP.txt", quote=FALSE, sep="\t")

#calculate % var. explained by each axis
ax1.v.f.t=cap1$CA$eig[1]/sum(cap1$CA$eig)
ax2.v.f.t=cap1$CA$eig[2]/sum(cap1$CA$eig)

#Plot:  supporting figure
par(mfrow=c(1,2))
plot(uf.fire.pcoa$points[,1],uf.fire.pcoa$points[,2] , main= "Fire-affected soils PCoA", type="n",xlab=paste("PCoA1: ",100*round(ax1.v.f,3)," var. explained",sep=""), ylab= paste("PCoA2: ",100*round(ax2.v.f,3)," var. explained",sep=""))
textxy(X=uf.fire.pcoa$points[,1], Y=uf.fire.pcoa$points[,2],labs=labels, offset=0, cex=0.8)
plot(envFIT.fire, p=0.10)
plot(cap1, cex=0.9,main = "Temperature-constrained fire-affected soils PCoA", xlab=paste("CAP_Ax1: ",100*round(ax1.v.f.t,3),"%var. explained",sep=""), ylab=paste("CAP_Ax2: ",100*round(ax2.v.f.t,3),"%var. explained",sep=""))
plot(c.ef, p= 0.10)


#####################
#Sloan neutral model fitting - FINISHED
#NOTE:  must use full dataset (including singleton OTUs) for this analysis 
#Source for model fits is from Burns et al. ISMEJ 2015, downloaded R code from their supporting materials.
source("sncm.fit_function.R")

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

#Fire Affected (All Sites Pool) - asks the question:  are fire soils a neutral subset of the whole?
obs.fireAP=sncm.fit(fire.t.sigs.NZ,taxon=rdp, stats=FALSE, pool=spp)
sta.fireAP=sncm.fit(fire.t.sigs.NZ,taxon=rdp, stats=TRUE, pool=spp)

# Recovered (All Sites Pool) asks the question:  are rec soils a neutral subset of the whole?
obs.recAP=sncm.fit(rec.t.sigs.NZ,taxon=rdp, stats=FALSE, pool=spp)
sta.recAP=sncm.fit(rec.t.sigs.NZ,taxon=rdp, stats=TRUE, pool=spp)

results=rbind(sta.np, sta.fireT, sta.recT, sta.fireAP, sta.recAP)
row.names(results)=c("all", "fire_total", "recovered_total","Fire_AllPool", "Recovered_AllPool")


#par(mfrow=c(2,3)) #for plotting in R studio w/out export
l1=list(obs.np, obs.recT, obs.fireT, obs.recAP, obs.fireAP)
l2=list(sta.np, sta.recT, sta.fireT, sta.recAP, sta.fireAP)
names=c("All", "Fire_Total", "Recovered_Total", "Recovered_AllPool", "FireAffected_AllPool")
out.sta=NULL

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
  
  #plot figure (SFig3)
  setEPS()
  if(i == 1){
  postscript("SFig3A.eps", width = 2.33, height=3, pointsize=10,paper="special")
  }
  if (i == 2){
  postscript("SFig3B.eps", width = 2.33, height=3, pointsize=10,paper="special")
  }
  if (i ==3){
  postscript("SFig3C.eps", width = 2.33, height=3, pointsize=10,paper="special")
  }
  if (i ==4){
    postscript("SFig3D.eps", width = 2.33, height=3, pointsize=10,paper="special")
  }
  if (i ==5){
    postscript("SFig3E.eps", width = 2.33, height=3, pointsize=10,paper="special")
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
results=cbind(results, out.sta)
write.table(results, "SloanNeutralModel.txt", quote=FALSE, sep="\t")


###########################
###Venn analysis - FINISHED
#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)

#make a binary OTU table for each fire classification
fireclass=map[,"Classification"]

#who are the most prevalent fire-affected OTUs
fire=comm[,fireclass=="FireAffected"]
fire=fire[rowSums(fire)>0,]
order(fire, rowSums(fire), decreasing=TRUE)

active.pa=1*(comm[,fireclass=="FireAffected"]>0)
recov.pa=1*(comm[,fireclass=="Recovered"]>0)
ref.pa=1*(comm[,fireclass=="Reference"]>0)

#summarize and combine
venndata=cbind(1*rowSums(active.pa>0),1*rowSums(recov.pa>0),1*rowSums(ref.pa>0))
colnames(venndata)=c("Fire-affected", "Recovered", "Reference")

#apply venn analysis
v=vennCounts(venndata)
v2=round(v[,"Counts"]/sum(v[,"Counts"]),2)
v[,"Counts"]<-v2
vennDiagram(v)


#Write out the results of venncounts
write.table(v, "VennCounts.txt", quote=FALSE, sep="\t")
par(mfrow=c(1,1))
setEPS()

postscript("Fig6.eps", width = 2.33, height=2.33, pointsize=8,paper="special")
par(mar=c(5,3,2,2)+0.1)
fig6=vennDiagram(v)
dev.off()


#############################
###Indicator species analysis
# resource:  https://cran.r-project.org/web/packages/indicspecies/vignettes/indicspeciesTutorial.pdf
library(indicspecies)

comm.t=as.data.frame(t(comm))
class=as.vector(map$Classification)
ind=multipatt(comm.t, class, control=how(nperm=999), duleg = TRUE, func="IndVal.g")
summary(ind, alpha=0.001, indvalcomp=TRUE, At=0.95, Bt=0.95)


#taxonomic affiliations of the OTUs that are indicators of fire-affected soils (A and B > 0.95, p = 0.001)
rdp.nosigs[row.names(comm)=="704748"]
rdp.nosigs[row.names(comm)=="25116"]
rdp.nosigs[row.names(comm)=="54107"]
#rdp.nosigs[row.names(comm)=="OTU_dn_211"]

#taxonomic affiliations of the top OTUs that are indicators of recovered soils (A and B = 1, p = 0.001)
rdp.nosigs[row.names(comm)=="511701"]
#rdp.nosigs[row.names(comm)=="OTU_dn_3450"]
#rdp.nosigs[row.names(comm)=="OTU_dn_9828"]
#rdp.nosigs[row.names(comm)=="OTU_dn_19940"]
rdp.nosigs[row.names(comm)=="153350"]
#rdp.nosigs[row.names(comm)=="OTU_dn_70"]

which(row.names(comm)=="153350", arr.ind = FALSE)
rdp.nosigs[277]
#which(row.names(comm)=="OTU_dn_70", arr.ind = FALSE)
rdp.nosigs[30]
