################################
### Tutorial for ecological statistics for CMSE 890:Microbial Metagenomics at Michigan State University
### Based on data from
###"Divergent extremes but convergent recovery of bacterial and archaeal soil 
### communities to an ongoing subterranean coal mine fire"
### by SH Lee, JW Sorensen, KL Grady, TC Tobin and A Shade  
### The ISME Journal 2017
### http://www.nature.com/ismej/journal/v11/n6/full/ismej20171a.html 
### Prepared  21 November 2018
### Author: Ashley Shade, Michigan State University; shade.ashley <at> gmail.com
################################

###Packages needed
##
## tidyverse - for formatting data and plotting with ggplot2
## vegan - for ecological statistics
## picante - for calculating phylogenetic diversity
## calibrate - for adding text labels to plots

##Before you start:
# A great intro to R tutorial is linked below. This tutorial is highly recommended for R newbies.
# https://github.com/edamame-course/2018-Tutorials/blob/master/Intro_R_RStudio/R_tutorial.md


################################
### Part 1: Getting to know the mapping file and plotting soil contextual data
################################
library(tidyverse)

#read in mapping file with soil data.  the mapping file contains contextual information that we will use for analysis.
map=read.table("InputFiles/Centralia_Collapsed_Map_forR.txt", header=TRUE, sep="\t")

#check its dimensions.  How many samples are there?  How many taxa are there?
dim(map)

#what does it look like?
head(map)

#what are its column names?  This is all of the information in the table.  Explore!
colnames(map)

#plot the soil chemistry data by temperature.  We will call this sfig3.
#prepare the data for plotting by "melting" to a long format 
map.long=melt(map, id.vars=c("SampleID", "SoilTemperature_to10cm", "Classification"), measure.vars=c("NO3N_ppm","NH4N_ppm","pH","SulfateSulfur_ppm","K_ppm","Ca_ppm","Mg_ppm","OrganicMatter_500","Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per"))

#make a gradient color palette
GnYlOrRd=colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)

#Give the figure a name. Here it is named sfig3 because that is what is was in the paper.
sfig3=ggplot(map.long, aes(y=as.numeric(SoilTemperature_to10cm), x=value))+
  
  #add points layer
  geom_point(aes(x=as.numeric(SoilTemperature_to10cm), y=value, shape=Classification, color=as.numeric(SoilTemperature_to10cm)))+
  
  #set facet with 4 columns, make x-axes appropriate for each variable
  facet_wrap(~variable, ncol=4, scales="free_y")+ 
  
  #set gradient for temperature and add gradient colorbar
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", guide_legend(title="Temperature"))+
  
  #define the axis labels
  labs(x="Temperature (Celsius)", y="    ")+
  
  #set a simple theme
  theme_bw(base_size=10)

sfig3
#If you wanted to save the output, use the code below.
#ggsave("Figures/SFig3.eps", width=178, units="mm")

##Subset contextual data inclusive of soil quantitative variables
env=map[,c("SoilTemperature_to10cm", "NO3N_ppm", "pH", "K_ppm", "Mg_ppm", "OrganicMatter_500", "NH4N_ppm", "SulfateSulfur_ppm", "Ca_ppm", "Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per","Fire_history")]

################################
### Part 2:  Getting to know the OTU table and preparing it for analysis
################################
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

rdp.nosigs=rdp[rowSums(comm)>1]

#designate a full dataset
comm.sigs=comm

#remove OTUs with an abundance = 1, across the entire dataset (singleton OTUs)
comm=comm[rowSums(comm)>1,]
sum(colSums(comm))

#transpose matrix.  
#We do this because some functions anticipate that the data will be samples (rows) x species (columns) instead of species x samples
comm.t=t(comm)


################################
### Part 3:  Calculate and plot within-sample (alpha) diversity
################################
library(vegan)

#calculate richness from OTU table
s=specnumber(comm.t)
#look at the resuelts
s

#calculate pielou's evenness from OTU table
#first, calculate shannon diversity
h=diversity(comm.t, index="shannon")
#look at the results
h
#then, use h (shannon) and s (richness) to calculate Pielou's evennes
pielou=h/log(s)

#calculate faith's phylogenetic diversity
library(picante)
#you need to load in your tree for this (this tree is output from QIIME, but could be from any)
tree <- read.tree("InputFiles/MASTER_RepSeqs_aligned_clean.tre")
pd=pd(comm.t, tree, include.root = FALSE)

library(tidyverse)
#combine alpha diversity data and fire classification (from map file)
div=tibble(map$SampleID,pd$PD,s,pielou,map$Classification)
# bonus!  Here is some info on "tibbles", which are data frames that include different kinds of data (e.g., numeric, factors)
# https://www.statmethods.net/management/reshape.html

#Give your tibble some column names
colnames(div)=c("SampleID", "PD", "Richness", "Pielou", "Classification")

#plot alpha diversity
#first, reshape the data to long format for mapping
div.long=melt(div, id.vars=c("SampleID", "Classification"))

#plot a facet
colors=c("red", "yellow", "green")

fig1 <- ggplot(data=div.long, aes(x=Classification, y=value))+
  geom_boxplot() + 
  geom_jitter(aes(shape=Classification, color=Classification))+
  facet_grid(variable~., scales="free_y") +
  scale_color_manual(values=colors) +
  scale_x_discrete(name="Fire classification") +
  scale_y_continuous(name="Alpha Diversity value") +
  theme_bw(base_size=10)
fig1
#if you want to save the figure output in eps format, use the command below:
#ggsave("Figures/Fig1.eps", width=86, units="mm")

#BONUS!  Check out this tutorial for the basics of data visualization with ggplot2 in R:
# https://github.com/edamame-course/2018-Tutorials/blob/master/data_visualization/data_visualization_tutorial.md 


################################
### Part 4:  Calculate and plot across-sample (beta) diversity
################################
#Calculate Bray-curtis dissimilarity. Remember that Bray-Curtis is a taxonomic, weighted resemblence.
bc.d=vegdist(t(comm), method="bray")

#Calculate Sorensen's dissimilarity.  Remember that Sorensen's is a taxonomic, unweighted resemblence.
sor.d=vegdist(t(comm), method="bray",binary=TRUE)

#Read in: UniFrac distances (output from QIIME).  
#Weighted UniFrac. Remember that weighted UniFrac is a phylogenetic, weighted resemblence.
wuf=read.table("InputFiles/weighted_unifrac_MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1)

#Sort by rows, columns (so they are in the consecutive order, same as the OTU table.  protip:  always check that your samples are in the same order!)
wuf=wuf[order(row.names(wuf)),order(colnames(wuf))]
wuf.d=as.dist(wuf)

#Read in: Unweighted UniFrac (output from QIIME). Remember that unweighted UniFrac is a phylogenetic, unweighted resemblence.
uuf=read.table("InputFiles/unweighted_unifrac_MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1)

#sort by rows, columns (so they are in the consecutive order)
uuf=uuf[order(row.names(uuf)),order(colnames(uuf))]
uuf.d=as.dist(uuf)

##Make PCoAs

#How many fire classifications are there?
unique(map$Classification)
#Map desired colors to the classifications to make our future legend for the PCoA.  
#We start with an empty color vector the same size as the number of samples (18), and then fill them in with a different color for each factor
Class=rep('black',nrow(map))
Class[map$Classification=="FireAffected"]='red'
Class[map$Classification=="Reference"]='green'
Class[map$Classification=="Recovered"]='yellow'

##PCoA A: Bray-Curtis
#We're taking this first ordination plot step by step.  
#First, perform the PCoA
bc.pcoa=cmdscale(bc.d, eig=TRUE)

#Look at the output
bc.pcoa
#Use the command below to interpret each component of the output
?cmdscale

#Second, calculate percent variance explained by each axes 1 and 2, then add to plot later
ax1.v.bc=bc.pcoa$eig[1]/sum(bc.pcoa$eig)
ax2.v.bc=bc.pcoa$eig[2]/sum(bc.pcoa$eig)
#look at the output
ax1.v.bc
ax2.v.bc

#Third, correlate measured environmental variables (env) to the axes 
envEF.bc=envfit(bc.pcoa, env)
#look at the output
envEF.bc

#Finally, bring it all together to plot the ordination.
plot(bc.pcoa$points[,1],bc.pcoa$points[,2],
     cex=1.5,pch=21,
     bg=Class,
     main="Bray Curtis PCoA", 
     xlab= paste("PCoA1: ",100*round(ax1.v.bc,3),"% var. explained",sep=""), 
     ylab= paste("PCoA2: ",100*round(ax2.v.bc,3),"%var. explained",sep=""))

#Add labels to points 
library(calibrate)
textxy(X=bc.pcoa$points[,1], Y=bc.pcoa$points[,2],labs=map$SampleID, cex=0.8)

#Add a legend
legend("topright",c('Fire Affected','Recovered','Reference'), pch=21,pt.bg=c("red", "yellow", "green"),lty=0)
#Protip: use the command below to see options for placement and vertical v. horizontal layout.
?legend 

#Add fitted environmental vectors
plot(envEF.bc, p.max=0.10, col="black", cex=1)

#Exercise:  make the equivalent plot for the other three resemblences.  

##PCoA B: Sorensen
#First, perform the PCoA
sor.pcoa=##FILL IN THE FUNCTIONS AND INPUT
  
#Second, calculate percent variance explained by each axes 1 and 2, then add to plot later
ax1.v.sor=
ax2.v.sor=

#Third, correlate measured environmental variables (env) to the axes 
envEF.sor=

#Finally, bring it all together to plot the ordination. 
  ##Note:  Nothing to fill in here - if you did the top part correctly, this should all work out)
plot(sor.pcoa$points[,1],sor.pcoa$points[,2],
     cex=1.5,pch=21,
     bg=Class,
     main="Sorensen PCoA", 
     xlab= paste("PCoA1: ",100*round(ax1.v.sor,3),"% var. explained",sep=""), 
     ylab= paste("PCoA2: ",100*round(ax2.v.sor,3),"%var. explained",sep=""))
textxy(X=sor.pcoa$points[,1], Y=sor.pcoa$points[,2],labs=map$SampleID, cex=0.8)

legend("topright", c('Fire Affected','Recovered','Reference'),pch=21,pt.bg=c("red", "yellow", "green"),lty=0)
plot(envEF.sor, p.max=0.10, col="black", cex=1)

##PCoA C:  Weighted UniFrac
#First, perform the PCoA
wuf.pcoa=##FILL IN THE FUNCTION AND INPUT
ax1.v.wuf=
ax2.v.wuf=
envEF.wuf=
  
#Plot the ordination:  If everything was filled in propperly above, this should work out.
plot(wuf.pcoa$points[,1],wuf.pcoa$points[,2],
   cex=1.5,pch=21,
   bg=Class,
   main="Sorensen PCoA", 
   xlab= paste("PCoA1: ",100*round(ax1.v.wuf,3),"% var. explained",sep=""), 
   ylab= paste("PCoA2: ",100*round(ax2.v.wuf,3),"%var. explained",sep=""))
textxy(X=wuf.pcoa$points[,1], Y=wuf.pcoa$points[,2],labs=map$SampleID, cex=0.8)
legend("topright", c('Fire Affected','Recovered','Reference'),pch=21,pt.bg=c("red", "yellow", "green"),lty=0)
plot(envEF.wuf, p.max=0.10, col="black", cex=1)

#PCoA D:  Unweighted UniFrac
##Try this one all on your own!  You can do it!!!

##Exercise:  compare the variance explained for each resemblence.  
###Which resemblence has the most explanatory value?  
#Axis 1
ax1.v.bc
ax1.v.sor
ax1.v.wuf
ax1.v.uuf

#Axis 2
ax2.v.bc
ax2.v.sor
ax2.v.wuf
ax2.v.uuf

###Part 5:  Hypothesis testing. 
#perform hypothesis testing to assess differences in community structure between 
#fire-affected, recovered, and reference sites

#permanova - tests for differences in CENTROID and/or DISPERSION
#use weighted UniFrac distance (because it had the most explanatory value of the four tested resemblences)
#The function for permanova in the vegan package is called "adonis" 
a=adonis(wuf.d~Class, distance=TRUE, permutations=1000)
#look at the results.  Is there a difference between these groups?
a

#permdisp - tests only for differences in DISPERSION
#multivariate dispersion
#the function for permdisp in the vegan package is called "betadisper"
b=betadisper(wuf.d, group=Class)
#look at the results.
b

#Are there significant differences between groups? 
#Use a post hoc Tukey's Honestly Significant Difference 
TukeyHSD(b, which = "group", ordered = FALSE,conf.level = 0.95)

#BONUS material!  Collapsing groups
#So, let's focus on differences between the ambient soils v. fire affected comparison.
#The way to do this is to perform the permanova test again, but after collapsing reference and recovered soils together
#perform hypothesis testing on fire-affected v. recovered+reference sites
Class2=sub("green", "yellow", Class)
#permanova
a=adonis(wuf.d~Class2, distance=TRUE, permutations=1000)
#Look at the results.  Notice that the degrees of freedom has decreased because there are now two factors/groups instead of 3
a

#multivariate dispersion with Tukey HSD
b=betadisper(wuf.d, group=Class2)
TukeyHSD(b, which = "group", ordered = FALSE,conf.level = 0.95)

###Congratulations!  You've finished your first ecological analyses!  
###For a taste of the full analysis pipeline and to recreate all of the figures in the paper, check out this public workflow:
### https://github.com/ShadeLab/PAPER_LeeSorensen_ISMEJ_2017/blob/master/R_analysis/Centralia2014_AmpliconWorkflow.R



