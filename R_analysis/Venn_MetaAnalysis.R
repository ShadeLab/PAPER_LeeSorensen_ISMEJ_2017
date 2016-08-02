setwd("~/Desktop")
map.f <- read.table("Centralia_Full_Map_Fixed.txt", sep="\t", header= TRUE, row.names=1, stringsAsFactors = FALSE)
map <- read.table("Centralia_Collapsed_Map_forR.txt", sep="\t", header= TRUE, row.names=1, stringsAsFactors = FALSE)
data <- read.table("Carini_RDP_rmCM.txt", sep="\t", header = TRUE, row.names=1, stringsAsFactors = FALSE)

Cen <- data[,grepl("C",colnames(data))]
Car <- data[,grepl("SRR", colnames(data))]

Samples<- unique(map.f$Sample)

collapse <- NULL

for (i in 1:length(Samples)){
  x <- Cen[,grepl(Samples[i], colnames(Cen))]
  collapse <- cbind(collapse,rowSums(x))
}
colnames(collapse) <- Samples

Reference <- collapse[,map$Classification=="Reference"]
FireAffected <- collapse[,map$Classification=="FireAffected"]
Recovered <- collapse[,map$Classification=="Recovered"]

sum(Reference)
sum(FireAffected)
sum(Recovered)
sum(Car)
data.f <-cbind(rowSums(Reference),rowSums(FireAffected), rowSums(Recovered), rowSums(Car))
colnames(data.f) <- c("Reference", "FireAffected", "Recovered", "RelicDNA")
library(vegan)
library(limma)
?rrarefy
data.rare<-rrarefy(t(data.f), 900542)
data.rare <- t(data.rare)
colSums(data.rare)


data.rare.pa <- 1*(data.rare>0)
data.rare.pa.nz <- data.rare.pa[rowSums(data.rare.pa)>0,]
v=vennCounts(data.rare.pa.nz)
v2=round(v[,"Counts"]/sum(v[,"Counts"]),2)
v[,"Counts"]<-v2
vennDiagram(v)
setwd("~/GitHub_Repos/ShadeLab/PAPER_LeeSorensen_inprep/R_analysis/")
data.m<-read.table("MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)

rdp <- data.m[,19]
data.m <- data.m[,-19]
ref.m<- rowSums(data.m[,map$Classification=="Reference"])
rec.m <- rowSums(data.m[,map$Classification=="Recovered"])
fa.m <- rowSums(data.m[,map$Classification=="FireAffected"])

class.m <- cbind(ref.m, rec.m, fa.m)
colnames(class.m) <- c("Reference", "Recovered", "FireAffected")
colSums(class.m)
class.m.rare <- rrarefy(t(class.m), 642000)

rowSums(class.m.rare)
class.m.rare <- t(class.m.rare)
class.m.r.pa <- 1*(class.m.rare>0)
venndata.class <- class.m.r.pa[rowSums(class.m.r.pa)>0,]

v=vennCounts(venndata.class)
v2=round(v[,"Counts"]/sum(v[,"Counts"]),2)
v[,"Counts"]<-v2
vennDiagram(v)