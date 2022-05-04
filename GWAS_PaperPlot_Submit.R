library(ggplot2)
library(cowplot)
library(dplyr)
library(colorRamps)
library(cmocean)
library(viridis)
library(ggnewscale)
library(tidyverse)
library(directlabels)
options(scipen=999)

ThemeSobr=  theme(
  panel.border = element_blank(),  
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  text = element_text(size=12, face="bold"),
  axis.line = element_line(colour = "black", size=0.5),
  legend.spacing.y= unit(0, 'cm'),
  axis.title = element_text(face="bold"),
  plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  axis.text = element_text(size = 10),
  legend.background = element_blank(),
  strip.text = element_text(face="bold", size=14)
)

### Phenotypic PCA ### 
FillCol=c("#3c562be5","#8bb088ff","#593872e5","#648dbee5",
      "#d78938e5","#aa9554e5","#d24437e5","#c586b8e5",
      "#726dd2e5","#9440c9e5","#46b87ee5","#833041e5",
      "#4aaaade5","#b0b43ce5","#a3dc92e5","#d28377e5")
Col=c("#3c562be5","#8bb088ff","#593872e5","#648dbee5",
      "#000000ff","#000000ff","#000000ff","#000000ff",
      "#000000ff","#000000ff","#000000ff","#000000ff",
      "#000000ff","#000000ff","#000000ff","#000000ff")
PhenoOrder=c("silvana","robigus",'illustris',"laura",
             "isabellinus","timaeus","lyrcaeus","seraphion",
             "euphrasius","aurora","euphone","numata",
             "messene","arcuella","bicoloratus","tarapotensis")

#Figure 1C #
# PhenoGeno=read.table("~/Paper/GWAS/V2/CPM_Antpost_sPee_WithGeno_WithPheno", stringsAsFactors = F, header=T)
PhenoGeno=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/CPM_Antpost_sPee_WithGeno_WithPheno", stringsAsFactors = F, header=T)
PhenoGeno$Pheno=factor(PhenoGeno$Pheno, level=PhenoOrder)
base=ggplot(PhenoGeno)
Fig1C=base+geom_point(aes(x=PC_1, y=PC_2, shape=Geno, fill=Pheno, color=Pheno), size=4, alpha=0.9)+
  scale_shape_manual("Genotype", values=c(8,21,22,23,24,25))+
  scale_fill_manual(values=FillCol, guide=F)+
  scale_color_manual(values=Col, guide=F)+
  ThemeSobr+
  xlab("Phenotype PCA 1 (23.48%)")+
  ylab("Phenotype PCA 2 (14.72%)")
  
Fig1C
save_plot("~/Paper/GWAS/V2/Plot/Fig1C.png", Fig1C)
save_plot("~/Paper/GWAS/V2/Plot/Fig1C.svg", Fig1C)

## Figure S8
FigS8=base+geom_point(aes(x=PC_3, y=PC_4, shape=Geno, fill=Pheno, color=Pheno), size=4, alpha=0.9)+
  scale_shape_manual("Genotype", values=c(8,21,22,23,24,25))+
  scale_fill_manual(values=FillCol, guide=F)+
  scale_color_manual(values=Col, guide=F)+
  ThemeSobr+
  xlab("Phenotype PCA 3 (5.16%)")+
  ylab("Phenotype PCA 4 (8.24%)")
save_plot("~/Paper/GWAS/V2/Plot/FigS8.png", FigS8)
save_plot("~/Paper/GWAS/V2/Plot/FigS8.svg", FigS8)


#### Genomic PCA ### Fig S1-S5
data=read.table(paste0("~/Paper/GWAS/V2/PCA/","AllSubSamp_AllRegion.PCA.txt"), header = F, stringsAsFactors = F)
colnames(data)=c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7",
                 "PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16",
                 "PC17","PC18","PC19","PC20","Geno","Pheno","Origin",'Sample','Region')

for (i in unique(data$Region))
{
base=ggplot(data[data$Region==i,])
P1=base+geom_point(aes(x=PC1, y=PC2, color=Origin, shape=Geno), size=4, alpha=0.9)+
  scale_shape_manual("Genotype", values=c(8,21,22,23,24,25))+
  scale_color_manual(values=Col)+
  ThemeSobr+
  xlab("Genotypic PCA 1")+
  ylab("Genotypic PCA 2")+
  facet_wrap(.~paste0("\"",Sample,"\"") , scales = "free")+
  ggtitle("Sample subset:")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
save_plot(paste0("~/Paper/GWAS/V2/Plot/", i,".PCA.png"),ncol=3, P1, base_aspect_ratio = 1.2)
 }


## Fig 3 ##
Col2=scales::viridis_pal(begin=0, end=0.8, option = "A", direction = 1)(2)
ThemeSobr2=  theme(
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  panel.border = element_blank(),  
  panel.grid = element_blank(),
  # panel.background = element_blank(),
  text = element_text(size=14, face="bold"),
  axis.line = element_line(colour = "black", size=0.5),
  legend.spacing.y= unit(0, 'cm'),
  axis.title = element_text(face="bold"),
  plot.title=element_text(size=15, face="bold",hjust=0.5, vjust=2),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  axis.text = element_text(size = 14),
  legend.background = element_blank(),
  strip.text = element_text(face="bold", size=14)
)

## A ## Hn123 multivariate association
data=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/AllSample.SupergeneAndFlank.snps.Hn123.CPM_Antpost_Hn123.PC_1_PC_2_PC_3.Assoc.mqfam.total.reformate.CHRPos.WithEmpP", stringsAsFactors = F, sep="\t", header = T)

data$LogP=-log10(data$P)
data$color="F"
data[data$EMP1==1e-06,]$color="T"
base=ggplot(data)
Assoc_Multi_Hn123=base+
  ylab("-Log10(p)")+
  xlab("")+
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.3) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey30", size=0.3)+
  ThemeSobr2+
  # geom_point( aes(x=CHRPOS, y=LogP, color=color), size=2, alpha=0.8, shape=21, fill="black")+
  # geom_point(data=data[data$color=="T",], aes(x=CHRPOS, y=LogP, color=color), size=2, alpha=0.8, shape=21, fill="black")+
  geom_point( aes(x=CHRPOS, y=LogP, color=color), size=2, alpha=0.8)+
  geom_point(data=data[data$color=="T",], aes(x=CHRPOS, y=LogP, color=color), size=2, alpha=0.8)+
  scale_color_manual(values=Col2, guide=F)+
  xlim(1000000,3500000)+
  ggtitle("Hn123: multivariate association")

## B ## Hn123 multivariate SNP Density
data=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/AllSample.SupergeneAndFlank.snps.Hn123.CPM_Antpost_Hn123.PC_1_PC_2_PC_3.Assoc.mqfam.total.reformate.CHRPos.WithEmpP", stringsAsFactors = F, sep="\t", header = T)
colnames(data)=c("SCAF","SNP","BP","NFAM","NIND","F","P","CHR","CHRPOS","EMP1","Sample", "Partition", "PCs")
SNPfoc=data[data$EMP1==1e-06,]$CHRPOS
UniqueSNPfoc=unique(SNPfoc)
UniqueSNPfoc=UniqueSNPfoc[!is.na(UniqueSNPfoc)]
DataNbSNP=data.frame(Pos=double(), NbSNP=double())
windSize=10000
Increment=100
Start=min(data$CHRPOS)
End=max(data$CHRPOS)
NbWindow=ceiling((End-Start)/Increment)
startWind=Start
for (Window in 1:NbWindow)
{
  EndWind=startWind+windSize
  Mid=startWind+(windSize/2)
  NbSnp=length(UniqueSNPfoc[UniqueSNPfoc>startWind & UniqueSNPfoc<EndWind])
  DataNbSNP[nrow(DataNbSNP)+1,]=c(Mid, NbSnp)
  startWind=startWind+Increment 
}

baseMultiDens=ggplot(DataNbSNP, aes(x=Pos, y=NbSNP))
SNPDens_Step_Multi_Hn123=baseMultiDens+
  ylab("Number of associated SNP")+
  xlab("")+
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.3) +
  geom_line(color=Col2[2], size=1)+
  ThemeSobr2+
  xlim(1000000,3500000)+
  ggtitle("Hn123: multivariate association, associated variant density")

#SNPDens_Step_Multi_Hn123
## C ## Hn123 univariate SNP Density
data=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/AllSample.SupergeneAndFlank.snps.Hn123.CPM_Post_Hn123.AllPartition.Genotypic", stringsAsFactors = F)
colnames(data)=c("SCAF","SNP","BP","A1","TEST","NMISS","BETA","STAT","P","CHR", "CHRPOS","Sample", "Partition", "PC")
data$LogP=-log10(data$P)
data=data[data$TEST!="DOMDEV",]
data=data[data$PC != "PC_3",]
SNPfoc=data[data$LogP>6,]$CHRPOS
UniqueSNPfoc=unique(SNPfoc)
UniqueSNPfoc=UniqueSNPfoc[!is.na(UniqueSNPfoc)]
DataNbSNP=data.frame(Pos=double(), NbSNP=double())
windSize=10000
Increment=100
Start=min(data$CHRPOS)
End=max(data$CHRPOS)
NbWindow=ceiling((End-Start)/Increment)
startWind=Start
for (Window in 1:NbWindow)
{
  EndWind=startWind+windSize
  Mid=startWind+(windSize/2)
  NbSnp=length(UniqueSNPfoc[UniqueSNPfoc>startWind & UniqueSNPfoc<EndWind])
  DataNbSNP[nrow(DataNbSNP)+1,]=c(Mid, NbSnp)
  startWind=startWind+Increment 
}

baseUniDens=ggplot(DataNbSNP, aes(x=Pos, y=NbSNP))
SNPDens_Step_Uni_Hn123=baseUniDens+
  ylab("Number of associated SNP")+
  xlab("Position on Chromosome 15 (bp)")+
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.3) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey30", size=0.3)+
  geom_line(color=Col2[2], size=1)+
  ThemeSobr2+
  xlim(1000000,3500000)+
  ggtitle("Hn123: combined univariate associations, associated variant density")

## D ## Hn0 multivariate association
data=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/OldAssovc/AllHn0_Antpost.PC_1_PC_2_PC_3_PC_4.Assoc.mqfam.total.reformate.CHRPos.WithEmpP", stringsAsFactors = F, sep="\t", header = T)
data$LogP=-log10(data$P)
data$color="F"
data[data$EMP1==1e-06,]$color="T"
base=ggplot(data)
Assoc_Multi_Hn0=base+
  ylab("-Log10(p)")+
  xlab("")+
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.3) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey30", size=0.3)+
  ThemeSobr2+
  # geom_point( aes(x=CHRPOS, y=LogP, color=color), size=2, alpha=0.8, shape=21, fill="black")+
  # geom_point(data=data[data$color=="T",], aes(x=CHRPOS, y=LogP, color=color), size=2, alpha=0.8, shape=21, fill="black")+
  geom_point( aes(x=CHRPOS, y=LogP, color=color), size=2, alpha=0.8)+
  geom_point(data=data[data$color=="T",], aes(x=CHRPOS, y=LogP, color=color), size=2, alpha=0.8)+
  scale_color_manual(values=Col2, guide=F)+
  xlim(1000000,3500000)+
  ggtitle("Hn0: multivariate association")
 Assoc_Multi_Hn0
 
#Assoc_Multi_Hn0
## E ## Hn0 multivariate SNP Dens
data=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/OldAssovc/AllHn0_Antpost.PC_1_PC_2_PC_3_PC_4.Assoc.mqfam.total.reformate.CHRPos.WithEmpP", stringsAsFactors = F, sep="\t", header = T)
colnames(data)=c("SCAF","SNP","BP","NFAM","NIND","F","P","CHR","CHRPOS","EMP1","Sample", "Partition", "PCs")
SNPfoc=data[data$EMP1==1e-06,]$CHRPOS
UniqueSNPfoc=unique(SNPfoc)
UniqueSNPfoc=UniqueSNPfoc[!is.na(UniqueSNPfoc)]
DataNbSNP_MultiHn0=data.frame(Pos=double(), NbSNP=double())
windSize=10000
Increment=100
Start=min(data$CHRPOS)
End=max(data$CHRPOS)
NbWindow=ceiling((End-Start)/Increment)
startWind=Start
for (Window in 1:NbWindow)
{
  EndWind=startWind+windSize
  Mid=startWind+(windSize/2)
  NbSnp=length(UniqueSNPfoc[UniqueSNPfoc>startWind & UniqueSNPfoc<EndWind])
  DataNbSNP_MultiHn0[nrow(DataNbSNP_MultiHn0)+1,]=c(Mid, NbSnp)
  startWind=startWind+Increment 
}

baseMultiDens=ggplot(DataNbSNP_MultiHn0, aes(x=Pos, y=NbSNP))
SNPDens_Step_Multi_Hn0=baseMultiDens+
  ylab("Number of associated SNP")+
  xlab("")+
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.3) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey30", size=0.3)+
  geom_line(color=Col2[2], size=1)+
  ThemeSobr2+
  xlim(1000000,3500000)+
  ggtitle("Hn0: multivariate association, associated variant density")

SNPDens_Step_Multi_Hn0
## F ## Hn0 univariate SNP Dens
data=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/AllSample.SupergeneAndFlank.snps.Hn0.AllPartition.Genotypic", stringsAsFactors = F)
colnames(data)=c("SCAF","SNP","BP","A1","TEST","NMISS","BETA","STAT","P","CHR", "CHRPOS","Sample", "Partition", "PC")
data$LogP=-log10(data$P)
data=data[data$TEST!="DOMDEV",]
data=data[data$PC != "PC_3",]
SNPfoc=data[data$LogP>6,]$CHRPOS
UniqueSNPfoc=unique(SNPfoc)
UniqueSNPfoc=UniqueSNPfoc[!is.na(UniqueSNPfoc)]
DataNbSNP_UniHn0=data.frame(Pos=double(), NbSNP=double())
windSize=10000
Increment=100
Start=min(data$CHRPOS)
End=max(data$CHRPOS)
NbWindow=ceiling((End-Start)/Increment)
startWind=Start
for (Window in 1:NbWindow)
{
  EndWind=startWind+windSize
  Mid=startWind+(windSize/2)
  NbSnp=length(UniqueSNPfoc[UniqueSNPfoc>startWind & UniqueSNPfoc<EndWind])
  DataNbSNP_UniHn0[nrow(DataNbSNP_UniHn0)+1,]=c(Mid, NbSnp)
  startWind=startWind+Increment 
}

baseUniDens=ggplot(DataNbSNP_UniHn0, aes(x=Pos, y=NbSNP))
SNPDens_Step_Uni_Hn0=baseUniDens+
  ylab("Number of associated SNP")+
  xlab("Position on chromosome 15 (bp)")+
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.3) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.3) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey30", size=0.3)+

  geom_line(color=Col2[2], size=1)+
  ThemeSobr2+
  xlim(1000000,3500000)+
  ggtitle("Hn0: combined univariate associations, associated variant density")
SNPDens_Step_Multi_Hn123=SNPDens_Step_Multi_Hn123+
  theme(plot.title = element_text(hjust = 0.5, vjust=10, face="bold.italic"),
                                                         plot.margin = margin(20, 0,5,0))
SNPDens_Step_Uni_Hn123=SNPDens_Step_Uni_Hn123+ 
  theme(plot.title = element_text(hjust = 0.5, vjust=10, face="bold.italic"),
                                                     plot.margin = margin(20, 0,5,0))
SNPDens_Step_Multi_Hn0=SNPDens_Step_Multi_Hn0+ 
  theme(plot.title = element_text(hjust = 0.5, vjust=10, face="bold.italic"),
                                                     plot.margin = margin(20, 0,5,0))
SNPDens_Step_Uni_Hn0=SNPDens_Step_Uni_Hn0+ 
  theme(plot.title = element_text(hjust = 0.5, vjust=10, face="bold.italic"),
                                                 plot.margin = margin(20, 0,5,0))
Assoc_Multi_Hn123=Assoc_Multi_Hn123+theme(plot.title = element_text(face="bold.italic"))
Assoc_Multi_Hn0=Assoc_Multi_Hn0+theme(plot.title = element_text(face="bold.italic"))

plotsHn123 <- align_plots(Assoc_Multi_Hn123,SNPDens_Step_Multi_Hn123, SNPDens_Step_Uni_Hn123, align = 'hv', axis = 'rltp')
MergedPlotHn123=plot_grid(plotsHn123[[1]],plotsHn123[[2]],plotsHn123[[3]],  ncol=1, labels = c('a','b','c'), rel_widths = c(1,1))
plotsHn0 <- align_plots(Assoc_Multi_Hn0,SNPDens_Step_Multi_Hn0, SNPDens_Step_Uni_Hn0, align = 'hv', axis = 'rltp')
MergedPlotHn0=plot_grid(plotsHn0[[1]],plotsHn0[[2]],plotsHn0[[3]],  ncol=1, labels = c('d','e','f'), rel_widths = c(1,1))
MergedPlotAll=MergedPlotAll+ theme(plot.title = element_text(hjust = 0.5, vjust=3))
MergedPlotAll=plot_grid(MergedPlotHn123,MergedPlotHn0,  ncol=2, labels = c('',''), rel_widths = c(1,1))
save_plot("~/Paper/GWAS/V2/Plot/Fig3_V2.png", MergedPlotAll, ncol=2, nrow=3, base_aspect_ratio = 3, bg="transparent")
save_plot("~/Paper/GWAS/V2/Plot/Fig3_V2.svg", MergedPlotAll, ncol=2, nrow=3, base_aspect_ratio = 3)

### Multiplot Hn0 Multi ### Figure S15
Col2=scales::viridis_pal(begin=0, end=0.8, option = "A", direction = 1)(2)
ThemeSobr2=  theme(
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  panel.border = element_blank(),  
  panel.grid = element_blank(),
  # panel.background = element_blank(),
  text = element_text(size=14, face="bold"),
  axis.line = element_line(colour = "black", size=0.5),
  legend.spacing.y= unit(0, 'cm'),
  axis.title = element_text(face="bold"),
  plot.title=element_text(size=15, face="bold",hjust=0.5, vjust=2),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  axis.text = element_text(size = 14),
  legend.background = element_blank(),
  strip.text = element_text(face="bold", size=14)
)
data2=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/AllSample.SupergeneAndFlank.snps.Hn0.CPM_Antpost_Hn0.AllPCs.Assoc.mqfam.total.reformate.CHRPos.WithEmpP", header=T, stringsAsFactors = F, sep="\t")
colnames(data2)=c("SCAF","SNP","BP","NFAM","NIND","F","P", "CHR", "CHRPOS","EMP1", "Sample", "Partition", "PCs")
data2=na.omit(data2)
data2$LogP=-log10(data2$P)
data2$color="F"
data2[data2$EMP1==1e-06,]$color="T"

base=ggplot(data2)
Col2=scales::viridis_pal(begin=0, end=0.8, option = "A", direction = 1)(2)
Plot=base+geom_point(aes(x=CHRPOS, y=LogP, color=color, fill=color), size=3, alpha=0.8, shape=21)+
  xlab("Position on Chromosome 15") + ylab("-log10(P)")+ 
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.8) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey50", size=0.8)+
  theme(plot.background = element_rect(fill = "transparent", colour = NA))+
  facet_wrap(PCs~., ncol=1, scales = "free")+ThemeSobr2+
  xlim(1000000,3500000)+
  scale_color_manual(values=Col2, guide=F)+
  scale_fill_manual(values=Col2, guide=F)
Plot
save_plot("~/Paper/GWAS/V2/Plot/Multiplot_Multivariate_Hn0_Add.png", Plot, nrow=length(unique(data2$PCs)), ncol=1, base_aspect_ratio = 4,limitsize = FALSE)

### Multiplot Hn123 Multi ### Figure S14
data2=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/AllSample.SupergeneAndFlank.snps.Hn123.CPM_Antpost_Hn123.AllPCs.Assoc.mqfam.total.reformate.CHRPos.WithEmpP", header=T, stringsAsFactors = F, sep="\t")
colnames(data2)=c("SCAF","SNP","BP","NFAM","NIND","F","P", "CHR", "CHRPOS","EMP1", "Sample", "Partition", "PCs")
data2=na.omit(data2)
data2$LogP=-log10(data2$P)
data2$color="F"
data2[data2$EMP1==1e-06,]$color="T"

base=ggplot(data2)
Col2=scales::viridis_pal(begin=0, end=0.8, option = "A", direction = 1)(2)
Plot=base+geom_point(data2, aes(x=CHRPOS, y=LogP, color=color, fill=color), size=3, alpha=0.8, shape=21)+
  xlab("Position on Chromosome 15") + ylab("-log10(P)")+ 
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.8) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey50", size=0.8)+
  theme(plot.background = element_rect(fill = "transparent", colour = NA))+
  facet_wrap(PCs~., ncol=1, scales = "free")+ThemeSobr2+
  xlim(1000000,3500000)+
  scale_color_manual(values=Col2, guide=F)+
  scale_fill_manual(values=Col2, guide=F)
Plot
save_plot("~/Paper/GWAS/V2/Plot/Multiplot_Multivariate_Hn123_Add.png", Plot, nrow=length(unique(data2$PCs)), ncol=1, base_aspect_ratio = 4,limitsize = FALSE)

### Multiplot Hn0 Uni ### Figures S18-19
data2=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/AllSample.SupergeneAndFlank.snps.Hn0.AllPartition.Genotypic", stringsAsFactors = F)

colnames(data2)=c("SCAF","SNP","BP","A1","TEST","NMISS","BETA","STAT","P","CHR", "CHRPOS","Sample", "Partition", "PC")
data2=data2[data2$TEST!="DOMDEV",]
data2=data2[(data2$PC != "PC_3" & data2$PC != "PC_2"),]
data2=data2[data2$NMISS > 0.7*max(data2$NMISS),]
data2$Partition=paste0(data2$Partition, ", ", data2$PC)
PartOrder=c("CPM_Antpost_Hn0, PC_1", "CPM_Ant_Hn0, PC_1","CPM_Post_Hn0, PC_1","CPM_Apex_Hn0, PC_1","CPM_Bande_Hn0, PC_1",  
            "CPM_Base_Hn0, PC_1","CPM_Yellow_Hn0, PC_1")
data2$Partition=factor(data2$Partition, level=PartOrder)

data2$LogP=-log10(data2$P)
##Figures S18
baseAdd=ggplot(data2[data2$TEST=="ADD",])
Col2=scales::viridis_pal(begin=0, end=0.8, option = "A", direction = 1)(2)
PlotAdd=baseAdd+geom_point(aes(x=CHRPOS, y=LogP), size=3, alpha=0.8)+
  xlab("Position on Chromosome 15") + ylab("-log10(P)")+ 
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.8) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey50", size=0.8)+
  theme(plot.background = element_rect(fill = "transparent", colour = NA))+
  facet_wrap(Partition~., ncol=1, scale="free")+ThemeSobr2+
  xlim(1000000,3500000)

save_plot("~/Paper/GWAS/V2/Plot/Multiplot_Univariate_Hn0_Add.png", PlotAdd, nrow=length(unique(data2$Partition)), ncol=1, base_aspect_ratio = 4,limitsize = FALSE)

##Figures S19
baseGeno=ggplot(data2[data2$TEST=="GENO_2DF",])
Col2=scales::viridis_pal(begin=0, end=0.8, option = "A", direction = 1)(2)
PlotGeno=baseGeno+geom_point(aes(x=CHRPOS, y=LogP), size=3, alpha=0.8)+
  xlab("Position on Chromosome 15") + ylab("-log10(P)")+ 
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.8) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey50", size=0.8)+
  theme(plot.background = element_rect(fill = "transparent", colour = NA))+
  facet_wrap(Partition~., ncol=1, scale="free")+ThemeSobr2+
  xlim(1000000,3500000)
save_plot("~/Paper/GWAS/V2/Plot/Multiplot_Univariate_Hn0_Geno.png", PlotGeno, nrow=length(unique(data2$Partition)), ncol=1, base_aspect_ratio = 4,limitsize = FALSE)

### Multiplot Hn123 Uni ###
data2=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/AllSample.SupergeneAndFlank.snps.Hn123.CPM_Post_Hn123.AllPartition.Genotypic", stringsAsFactors = F)
colnames(data2)=c("SCAF","SNP","BP","A1","TEST","NMISS","BETA","STAT","P","CHR", "CHRPOS","Sample", "Partition", "PC")
data2$LogP=-log10(data2$P)
Data2Sub=data2[data2$Partition=="CPM_Antpost_Hn123" & data2$PC=="PC_2",]
data2=data2[data2$TEST!="DOMDEV",]
data2=data2[(data2$PC != "PC_3" & data2$PC != "PC_2"),]
data2=rbind(data2, Data2Sub)
data2$Partition=paste0(data2$Partition, ", ", data2$PC)
data2=data2[data2$NMISS > 0.7*max(data2$NMISS),]
PartOrder=c("CPM_Antpost_Hn123, PC_1","CPM_Antpost_Hn123, PC_2", "CPM_Ant_Hn123, PC_1","CPM_Post_Hn123, PC_1","CPM_Apex_Hn123, PC_1","CPM_Bande_Hn123, PC_1",  
              "CPM_Base_Hn123, PC_1","CPM_Yellow_Hn123, PC_1")
  data2$Partition=factor(data2$Partition, level=PartOrder)
  

data2$LogP=-log10(data2$P)
#Figure S16
baseAdd=ggplot(data2[data2$TEST=="ADD",])
Col2=scales::viridis_pal(begin=0, end=0.8, option = "A", direction = 1)(2)
PlotAdd=baseAdd+geom_point(aes(x=CHRPOS, y=LogP), size=3, alpha=0.8)+
  xlab("Position on Chromosome 15") + ylab("-log10(P)")+ 
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.8) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey50", size=0.8)+
  theme(plot.background = element_rect(fill = "transparent", colour = NA))+
  facet_wrap(Partition~., ncol=1, scale="free")+ThemeSobr2+
  xlim(1000000,3500000)

save_plot("~/Paper/GWAS/V2/Plot/Multiplot_Univariate_Hn123_Add.png", PlotAdd, nrow=length(unique(data2$Partition)), ncol=1, base_aspect_ratio = 4,limitsize = FALSE)

#Figure S17
baseGeno=ggplot(data2[data2$TEST=="GENO_2DF",])
Col2=scales::viridis_pal(begin=0, end=0.8, option = "A", direction = 1)(2)
PlotGeno=baseGeno+geom_point(aes(x=CHRPOS, y=LogP), size=3, alpha=0.8)+
  xlab("Position on Chromosome 15") + ylab("-log10(P)")+ 
  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.8) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.8) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey50", size=0.8)+
  theme(plot.background = element_rect(fill = "transparent", colour = NA))+
  facet_wrap(Partition~., ncol=1, scale="free")+ThemeSobr2+
  xlim(1000000,3500000)
save_plot("~/Paper/GWAS/V2/Plot/Multiplot_Univariate_Hn123_Geno.png", PlotGeno, nrow=length(unique(data2$Partition)), ncol=1, base_aspect_ratio = 4,limitsize = FALSE)

## LD plot ## Figure S6
dataLD=read.table("~/Analysis/GWAS/Empirique/LD/AllSample.SupergeneAndFlank.snps.Hn0.Biall.NoMiss.recode.ld.CHRPos.CHRPos", stringsAsFactors = F, header=T)
base=ggplot(dataLD, aes(x=CHRPOS.1, y=CHRPOS, z=R2))
#Hn0
LD_Hn0=base+  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey50", size=0.3) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey50", size=0.3)+
  geom_hline(yintercept = 1416677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_hline(yintercept = 1826677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_hline(yintercept = 3129843, linetype="dashed", color = "grey50", size=0.3) +
  geom_hline(yintercept=2010000, linetype="dashed", color = "grey50", size=0.3)+
  stat_summary_2d(fun=mean, bins=75)+
  ThemeSobr+scale_fill_viridis("LD", option = "B", limits=c(0.2,0.8), breaks=c(0.2,0.5,0.8), labels=c(0.2,0.5,0.8))+
  xlab("Position on chromosome 15")+
  ylab("Position on chromosome 15")+
  xlim(500000,4000000)+
  ylim(500000,4000000)+
  ggtitle("Hn0 samples (only homozygous)")

dataLD=read.table("~/Analysis/GWAS/Empirique/LD/AllSample.SupergeneAndFlank.snps.Hn123.Biall.NoMiss.recode.SNPName.ld.CHRPos.CHRPos", stringsAsFactors = F, header=T)
base=ggplot(dataLD, aes(x=CHRPOS.1, y=CHRPOS, z=R2))
#Hn123
LD_Hn123=base+  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey50", size=0.3) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey50", size=0.3)+
  geom_hline(yintercept = 1416677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_hline(yintercept = 1826677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_hline(yintercept = 3129843, linetype="dashed", color = "grey50", size=0.3) +
  geom_hline(yintercept=2010000, linetype="dashed", color = "grey50", size=0.3)+
  stat_summary_2d(fun=mean, bins=75)+
  ThemeSobr+scale_fill_viridis("LD", option = "B", limits=c(0.2,0.8), breaks=c(0.2,0.5,0.8), labels=c(0.2,0.5,0.8))+
  xlab("Position on chromosome 15")+
  ylab("Position on chromosome 15")+
  xlim(500000,4000000)+
  ylim(500000,4000000)+
  ggtitle("Hn123 samples (homozygous and heterozygous)")

### Homo 123 ###
dataLD=read.table("~/Analysis/GWAS/Empirique/LD/AllSample.SupergeneAndFlank.snps.Hn123_Homo.Biall.NoMiss.recode.SNPName.ld.ld.CHRPos.CHRPos", stringsAsFactors = F, header=T)
base=ggplot(dataLD, aes(x=CHRPOS.1, y=CHRPOS, z=R2))
LD_HomoHn123=base+  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey50", size=0.3) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey50", size=0.3)+
  geom_hline(yintercept = 1416677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_hline(yintercept = 1826677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_hline(yintercept = 3129843, linetype="dashed", color = "grey50", size=0.3) +
  geom_hline(yintercept=2010000, linetype="dashed", color = "grey50", size=0.3)+
  stat_summary_2d(fun=mean, bins=75)+
  ThemeSobr+scale_fill_viridis("LD", option = "B", limits=c(0.2,0.8), breaks=c(0.2,0.5,0.8), labels=c(0.2,0.5,0.8))+
  xlab("Position on chromosome 15")+
  ylab("Position on chromosome 15")+
  xlim(500000,4000000)+
  ylim(500000,4000000)+
  ggtitle("Only Hn123 samples homozygous")

### Hn123 + Hn0 ###
dataLD=read.table("~/Analysis/GWAS/Empirique/LD/AllSample.SupergeneAndFlank.snps.Hn0_Hn123.Biall.NoMiss.recode.SNPName.ld.CHRPos.CHRPos", stringsAsFactors = F, header=T)
base=ggplot(dataLD, aes(x=CHRPOS.1, y=CHRPOS, z=R2))
LD_Hn0Hn123=base+  geom_vline(xintercept = 1416677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_vline(xintercept = 1826677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_vline(xintercept = 3129843, linetype="dashed", color = "grey50", size=0.3) +
  geom_vline(xintercept=2010000, linetype="dashed", color = "grey50", size=0.3)+
  geom_hline(yintercept = 1416677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_hline(yintercept = 1826677, linetype="dashed", color = "grey50", size=0.3) + 
  geom_hline(yintercept = 3129843, linetype="dashed", color = "grey50", size=0.3) +
  geom_hline(yintercept=2010000, linetype="dashed", color = "grey50", size=0.3)+
  stat_summary_2d(fun=mean, bins=75)+
  ThemeSobr+scale_fill_viridis("LD", option = "B", limits=c(0.2,0.8), breaks=c(0.2,0.5,0.8), labels=c(0.2,0.5,0.8))+
  xlab("Position on chromosome 15")+
  ylab("Position on chromosome 15")+
  xlim(500000,4000000)+
  ylim(500000,4000000)+
  ggtitle("Hn0 and Hn123 samples combined")

MergedPlot=plot_grid(LD_Hn123, LD_Hn0, LD_Hn0Hn123, LD_HomoHn123, ncol=2, labels = c('a','b','c','d'))
save_plot("~/Analysis/GWAS/Empirique/LD/CombinedLD.png", MergedPlot,ncol=2, nrow=2, base_aspect_ratio = 1.3)
save_plot("~/Analysis/GWAS/Empirique/LD/CombinedLD.svg", MergedPlot,ncol=2, nrow=2, base_aspect_ratio = 1.3)


### All phenotype PCA ### Figure S9
### Phenotypic PCA ###
FillCol=c("#3c562be5","#8bb088ff","#593872e5","#648dbee5",
          "#d78938e5","#aa9554e5","#d24437e5","#c586b8e5",
          "#726dd2e5","#9440c9e5","#46b87ee5","#833041e5",
          "#4aaaade5","#b0b43ce5","#a3dc92e5","#d28377e5")
Col=c("#3c562be5","#8bb088ff","#593872e5","#648dbee5",
      "#000000ff","#000000ff","#000000ff","#000000ff",
      "#000000ff","#000000ff","#000000ff","#000000ff",
      "#000000ff","#000000ff","#000000ff","#000000ff")
PhenoOrder=c("silvana","robigus",'illustris',"laura",
             "isabellinus","timaeus","lyrcaeus","seraphion",
             "euphrasius","aurora","euphone","numata",
             "messene","arcuella","bicoloratus","tarapotensis")


# PhenoGeno=read.table("~/Paper/GWAS/V2/Phenotypes/UsedPheno/FormatedForPlot/CPM_AllPartition_AnnotedSample.txt", stringsAsFactors = F, header=T)
PhenoGeno=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/CPM_AllPartition_AnnotedSample.txt", stringsAsFactors = F, header=T)
PhenoGeno$Pheno=factor(PhenoGeno$Pheno, level=PhenoOrder)
PhenoGeno=PhenoGeno[!is.na(PhenoGeno$Pheno),]
PartOrder=c("Antpost","Ant","Post","Apex","Bande","Base","Yellow")
PhenoGeno$Partition=factor(PhenoGeno$Partition, level=PartOrder)
base=ggplot(PhenoGeno)
FigAllPheno=base+geom_point(aes(x=PC_1, y=PC_2, shape=Geno, fill=Pheno, color=Pheno), size=5, alpha=0.9)+
  scale_shape_manual("Genotype", values=c(8,21,22,23,24,25))+
  scale_fill_manual(values=FillCol, guide=F)+
  scale_color_manual(values=Col, guide=F)+
  ThemeSobr+
  facet_grid(Partition ~ Sample, scales="free")+
  ggtitle("Sample subset:")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  xlab("Phenotype PCA 1")+
  ylab("Phenotype PCA 2")

 save_plot("~/Paper/GWAS/V2/Plot/FigAllPheno.png", FigAllPheno, ncol = 2, nrow = 7)
 save_plot("~/Paper/GWAS/V2/Plot/FigAllPheno.svg", FigAllPheno, ncol = 2, nrow = 7)

 ### #### Join Plot expression, Association & Annotation ##### Figure S20
 ## A ## Hn123 multivariate association
 data=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/AllSample.SupergeneAndFlank.snps.Hn123.CPM_Antpost_Hn123.PC_1_PC_2_PC_3.Assoc.mqfam.total.reformate.CHRPos.WithEmpP", stringsAsFactors = F, sep="\t", header = T)
 data$LogP=-log10(data$P)
 base=ggplot(data)
 Assoc_Multi_Hn123=base+
   ylab("-Log10(p)")+
   xlab("")+
   geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.3) + 
   geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.3) + 
   geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.3) +
   geom_vline(xintercept=2010000, linetype="dashed", color = "grey30", size=0.3)+
   ThemeSobr2+
   geom_point( aes(x=CHRPOS, y=LogP), size=2, alpha=0.8, shape=21, fill="black")+
   scale_x_continuous(breaks=pretty(1300000:3250000,n=8), limits = c(1300000,3250000))+
      ggtitle("Hn123: multivariate association")
 Assoc_Multi_Hn123
 ## B ## Hn0 multivariate association
 # dataHn0=read.table("~/Analysis/GWAS/Empirique/Multivariate/Hn0/OldAssovc/AllHn0_Antpost.PC_1_PC_2_PC_3_PC_4.Assoc.mqfam.total.reformate.CHRPos.WithEmpP", stringsAsFactors = F, sep="\t", header = T)
 dataHn0=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/AllHn0_Antpost.PC_1_PC_2_PC_3_PC_4.Assoc.mqfam.total.reformate.CHRPos.WithEmpP", stringsAsFactors = F, sep="\t", header = T)
 dataHn0$LogP=-log10(dataHn0$P)
 base=ggplot(dataHn0)
 Assoc_Multi_Hn0=base+
   ylab("-Log10(p)")+
   xlab("")+
   geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.3) + 
   geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.3) + 
   geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.3) +
   geom_vline(xintercept=2010000, linetype="dashed", color = "grey30", size=0.3)+
   ThemeSobr2+
   geom_point( aes(x=CHRPOS, y=LogP), size=2, alpha=0.8, shape=21, fill="black")+
   scale_x_continuous(breaks=pretty(1300000:3250000,n=8), limits = c(1300000,3250000))+
   ggtitle("Hn0: multivariate association")
 
 # C Expression PrePupae ##
  PPex=read.table('~/Paper/GWAS/V2/Submit/CodeAndFiles/ExpressionPPtableWithPos.txt', stringsAsFactors = F, header =T)
  BaseEx=ggplot(PPex)
 PPlot=BaseEx+geom_point(aes(x=PPex$CHRPos, y=PPex$log10FDR, color=as.factor(PPex$test)), size=3)+
   scale_colour_manual(values=c( "blue", "black", "red"))+
   ThemeSobr2+
   geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.3) + 
   geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.3) + 
   geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.3) +
   geom_vline(xintercept=2010000, linetype="dashed", color = "grey30", size=0.3)+
   scale_x_continuous(breaks=pretty(1300000:3250000,n=8), limits = c(1300000,3250000))+
   theme(legend.position = "none")+xlab("")+ylab("Log10(FDR)")
 #D Expression 24h
 ex24=read.table('~/Paper/GWAS/V2/Submit/CodeAndFiles/Expression24HtableWithPos.txt', stringsAsFactors = F, header =T)
 BaseEx24=ggplot(ex24)
 Plot24H=BaseEx24+geom_point(aes(x=ex24$CHRPos, y=ex24$log10FDR, color=as.factor(ex24$test)), size=3)+
   scale_colour_manual(values=c( "blue", "black", "red"))+
   ThemeSobr2+
   geom_vline(xintercept = 1416677, linetype="dashed", color = "grey30", size=0.3) + 
   geom_vline(xintercept = 1826677, linetype="dashed", color = "grey30", size=0.3) + 
   geom_vline(xintercept = 3129843, linetype="dashed", color = "grey30", size=0.3) +
   geom_vline(xintercept=2010000, linetype="dashed", color = "grey30", size=0.3)+
   scale_x_continuous(breaks=pretty(1300000:3250000,n=8), limits = c(1300000,3250000))+
   theme(legend.position = "none")+xlab("Position on chromosome 15")+ylab("Log10(FDR)")
 
 #E Annotation
 gff=read.table("~/Paper/GWAS/V2/Submit/CodeAndFiles/Hmel2_AnnotationGoodName_CHRPos.gff", fill=T, stringsAsFactors = , header=TRUE) # Read the annotation file (cf AnnotationPlot.R)
 CHROM="chr15"# Define the chromosome we are working on
 deb=1300000 # Define the first position of the chromosome we looked at
 fin=3250000  # Define the last position of the chromosome we looked at
 gffCHR=subset(gff, (gff$CHR==CHROM & gff$start >= deb & gff$end <= fin ))
 GenePos=subset(gffCHR, (gffCHR$V3=="gene" & gffCHR$V7=="+"))
 GeneNeg=subset(gffCHR, (gffCHR$V3=="gene" & gffCHR$V7=="-"))
 CDSPos=subset(gffCHR, (gffCHR$V3=="CDS" & gffCHR$V7=="+"))
 CDSNeg=subset(gffCHR, (gffCHR$V3=="CDS" & gffCHR$V7=="-"))
 Gene=subset(gffCHR, gffCHR$V3=="gene")
 Gene$y=-3
 AnnotPlot=ggplot(GenePos)+
   geom_rect(aes(xmin=GenePos$start, xmax=GenePos$end, ymin=1, ymax=2), fill="black")+
   geom_rect(data=GeneNeg, aes(xmin=GeneNeg$start, xmax=GeneNeg$end, ymin=0, ymax=1), fill="black")+
   ylim(0,2)+
   ThemeSobr2+
   theme(legend.position="none",
         axis.text=element_blank(),
         axis.line=element_blank(),
         axis.ticks=element_blank(),
         axis.title=element_blank())
 
 plotChoices <- plot_grid(Assoc_Multi_Hn123,Assoc_Multi_Hn0, PPlot, Plot24H, AnnotPlot, align="v", ncol=1,
                          rel_heights=c(1,1,1,1,0.2))

 save_plot("~/Paper/GWAS/V2/Plot/Assoc_Exp-Annot.svg", plotChoices, nrow= 5, base_aspect_ratio = 4)
 save_plot("~/Paper/GWAS/V2/Plot/Assoc_Exp-Annot.png", plotChoices, nrow= 5, base_aspect_ratio = 4, bg="transparent")
 
 ### Twistt Plot ### Figure S7
 options(scipen=999)
 File="~/Paper/GWAS/V2/Submit/CodeAndFiles/Weight_With_Pos_Supergene_Twistt4.Header.sort.SW25000.txt"
 FileSuff=sub(".txt","", File)
 raw=read.csv(File, sep="\t", header=T, stringsAsFactors = F,check.names=FALSE)
 agp=read.table("~/Paper/GWAS/Hmel2_chromosomes.agp", fill=TRUE, stringsAsFactors = FALSE)
 
 for (i in 1:nrow(raw)) ## Set position along the Chromosome 15 instead that on scaffolds
 { y =  which(agp$V6 == raw[i,1])
 raw$Pos[i]=raw$Pos[i]+agp[y,2]}
 NbTree=sum(raw[1, 3:5])
 Data=raw
 for (i in 3:5){Data[,i]=Data[,i]/NbTree}
 DataWide=Data %>% gather(Topo, Proportion, colnames(Data[,3:5]))
 col=c("#ea0000ff",  "#07008bff",  "#5a5a5aff")
 base=ggplot(DataWide)
 plot=base+geom_step(aes(x=Pos, y=Proportion, colour=Topo), size=1, alpha=0.8)+  
   scale_color_manual(values=col)+ theme_cowplot(14)+
   ylab("Weight")+xlab("Position on chromosome 15")
 
 save_plot(paste(FileSuff,"Full.SW25000.pdf", sep=""), plot, base_aspect_ratio = 3)
 
 