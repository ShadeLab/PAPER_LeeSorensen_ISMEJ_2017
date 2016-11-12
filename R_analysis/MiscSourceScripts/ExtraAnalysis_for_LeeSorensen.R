#This panel was omitted from paper
#heatmap of phylum-level responses, standardized (z-score) based on occurrences (oc) - across rows (gplots)
#standardize responses
comm.phylum.oc=decostand(comm.phylum, method="standardize", margin=1)
comm.phylum.oc=as.matrix(comm.phylum.oc)
#create color pallette; see: http://colorbrewer2.org/ 
#hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")

#more info on heatmaps:  http://www.inside-r.org/packages/cran/gplots/docs/heatmap.2
#dev.off()

#fig5B<-heatmap.2(comm.phylum.oc,col=hc(100),scale="row",key=TRUE,symkey=FALSE, trace="none", density.info="none",dendrogram="both", labRow=phylumnames, colCol=c(map$Class_Color), margins=c(5,13), srtCol=90)
#export heatmap 
#dev.off()
#178 mm is 7 inches
#setEPS()
#postscript("Figures/Fig5B.eps", width = 7, height=7, pointsize=10, paper="special")
#fig<-heatmap.2(comm.phylum.oc,
#            col=hc(100),
#            scale="row",
#            key=TRUE,
#            symkey=FALSE, 
#            trace="none", 
#            density.info="none",
#            dendrogram="both", 
#            labRow=phylumnames, 
#            colCol=c(map$Class_Color), 
#            margins=c(5,13), 
#            srtCol=90)
#dev.off()

###Venn analysis - not included in final ms
#Limma must be installed from bioconductor using the commands below
#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
# library(limma)
# 
# #make a binary OTU table for each fire classification
# fireclass=map[,"Classification"]
# 
# active.pa=1*(comm[,fireclass=="FireAffected"]>0)
# recov.pa=1*(comm[,fireclass=="Recovered"]>0)
# ref.pa=1*(comm[,fireclass=="Reference"]>0)
# 
# #summarize and combine
# venndata=cbind(1*rowSums(active.pa>0),1*rowSums(recov.pa>0),1*rowSums(ref.pa>0))
# colnames(venndata)=c("Fire-affected", "Recovered", "Reference")
# 
# #apply venn analysis
# v=vennCounts(venndata)
# v2=round(v[,"Counts"]/sum(v[,"Counts"]),2)
# v[,"Counts"]<-v2
# vennDiagram(v)
# 
# 
# #Write out the results of venncounts
# write.table(v, "VennCounts.txt", quote=FALSE, sep="\t")
# par(mfrow=c(1,1))
# setEPS()
# 
# postscript("Figures/Venn.eps", width = 2.33, height=2.33, pointsize=8,paper="special")
# par(mar=c(5,3,2,2)+0.1)
# fig6=vennDiagram(v)
# dev.off()


#############################
###Indicator species analysis - not included in the final manuscript
# resource:  https://cran.r-project.org/web/packages/indicspecies/vignettes/indicspeciesTutorial.pdf
# library(indicspecies)
# 
# comm.t=as.data.frame(t(comm))
# class=as.vector(map$Classification)
# ind=multipatt(comm.t, class, control=how(nperm=999), duleg = TRUE, func="IndVal.g")
# summary(ind, alpha=0.001, indvalcomp=TRUE, At=0.95, Bt=0.95)
# 
# 
# #taxonomic affiliations of the OTUs that are indicators of fire-affected soils (A and B > 0.95, p = 0.001)
# rdp.nosigs[row.names(comm)=="OTU_dn_34"]
# rdp.nosigs[row.names(comm)=="OTU_dn_2"]
# rdp.nosigs[row.names(comm)=="OTU_dn_125"]
# rdp.nosigs[row.names(comm)=="704748"]
# rdp.nosigs[row.names(comm)=="240440"]
# rdp.nosigs[row.names(comm)=="25116"]
# 
# #taxonomic affiliations of the top OTUs that are indicators of recovered soils (A and B = 1, p = 0.001)
# rdp.nosigs[row.names(comm)=="OTU_dn_5462"]
# rdp.nosigs[row.names(comm)=="511701"]
# rdp.nosigs[row.names(comm)=="OTU_dn_3443"]
# rdp.nosigs[row.names(comm)=="OTU_dn_1439"]
# rdp.nosigs[row.names(comm)=="OTU_dn_2728"]
# rdp.nosigs[row.names(comm)=="OTU_dn_14082"]
# 
# #example to extract RDP taxonomy assignment
# which(row.names(comm)=="704748", arr.ind = FALSE)
# rdp.nosigs[138]
# which(row.names(comm)=="OTU_dn_34", arr.ind = FALSE)
# rdp.nosigs[41]