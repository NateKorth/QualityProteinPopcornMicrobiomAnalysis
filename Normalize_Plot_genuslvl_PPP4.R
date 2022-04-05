library(metagenomeSeq)
library(phyloseq)
library(ape)
library(ggplot2)
library(reshape2)
library(vegan)
library(superheat)
library(DESeq2)
library(ggfortify)
library(corrplot)
library(stringr)
library(sommer)
library(data.table)
library(ggpubr)
library(devtools)
library(dplyr)
library(tidyverse)
#another messy set of code with a few valuable pieces mixed in, hope you can find something useful!


setwd("~/PopGhum/Popghum2/FRLysValidation/")

tree=read_tree("exported/PPP4_tree.nwk")
#genuslevel
table.alpha<-read.delim("./exported/PPP4_OTU_genus_Table_silva_taxonomy.tsv")
table.alpha0<-table.alpha[-c(1:6,67)]
#familylevel
#table.alpha<-read.delim("./exported/PPP4_OTU_family_Table_silva_taxonomy.tsv.txt")
#table.alpha0<-table.alpha[-c(1:5,67)]

table.alpha0[1:4,1:4]

#remove ID column and add it as column names
#A.OTU<-table.alpha[-c(1)]
A.OTU<-table.alpha0
OTU_IDs<-as.data.frame(table.alpha$Genus)
#OTU_IDs<-as.data.frame(table.alpha$Family)
rownames(A.OTU)<-OTU_IDs[,1]

A.OTU[1:3,1:37]
A.OTU$Consensus.Lineage<-NULL

#trim low abundance taxa:
B.OTU<-A.OTU
B.OTU$Undetermined<-NULL
B.OTU$AveReads<-rowMeans(B.OTU)
B.OTU<-B.OTU[-which(B.OTU$AveReads<=2),]
B.OTU$AveReads<-NULL

###################################################################
#read in Silva taxonomy | before, in excel seperate taxonomy into columns
taxa<-table.alpha[1:6]

taxa1<-taxa[which(taxa$Genus %in% row.names(B.OTU)),]
#taxa1<-taxa[which(taxa$Family %in% row.names(B.OTU)),]

OTUdata<-AnnotatedDataFrame(taxa1)
rownames(OTUdata)<-taxa1$Genus
#rownames(OTUdata)<-taxa1$Family
#normalize data
#load metadata
Meta<-load_phenoData("PPP4_metadata.txt",sep="\t")
#name the columns in Meta data file
colnames(Meta)<-c("Substrate","Treatment","Time","Rep","Acetate","Propionate","Isobutyrate","Butyrate","Isovalerate","Valerate")
#Match names in meta data to OTU table
ordA<- match(colnames(B.OTU), rownames(Meta))
Ameta <- Meta[ordA,]

#insert subset into phenodataframe:
phenotypeDataA<-AnnotatedDataFrame(Ameta)

#Make MetagenomeSeq object
A.data<-newMRexperiment(B.OTU, phenoData = phenotypeDataA, featureData = OTUdata)

#trim that shit yo: taxa only present in 20% of taxa removed
A.datatrim<-filterData(A.data,present = 10)

#normalize by CSS:
pA<-cumNormStatFast(A.datatrim)
A.datatrim<-cumNorm(A.datatrim,p=pA)

#log=TRUE -log2 transform
PPP_normal_log<-as.data.frame(MRcounts(A.datatrim,norm=TRUE,log=FALSE))
#UFMU_A_normal_log<-log10(UFMU_A_normal_log)
PPP_normal_log[1:5,1:5]
#trimmed relative abundance OTU table:
PPP.relabun.tr<-as.data.frame(MRcounts(A.datatrim,norm=TRUE,log=FALSE))
PPP.relabun.trim<-sweep(PPP.relabun.tr,2,colSums(PPP.relabun.tr),"/")



P.RA<-PPP.relabun.trim
taxa2<-taxa1[which(taxa1$X.OTU.ID %in% row.names(P.RA)),]

P.RA1<-as.data.frame(t(P.RA))
P.RA1$Substrate<-factor(Meta$Substrate,levels=c("Media","Lysine","FL"))
P.RA1$Treatment<-factor(Meta$Treatment,levels=c("Media-0hr","Media-16hr","Media-32hr","Media-48hr","Lysine-16hr","Lysine-32hr","Lysine-48hr","FL-16hr","FL-32hr","FL-48hr"))
P.RA1$Time<-Meta$Time
P.RA1$Rep<-Meta$Rep
P.RA1$Acetate<-as.numeric(Meta$Acetate)
P.RA1$Propionate<-as.numeric(Meta$Propionate)
P.RA1$Isobutyrate<-as.numeric(Meta$Isobutyrate)
P.RA1$Butyrate<-as.numeric(Meta$Butyrate)
P.RA1$Isovaleric<-as.numeric(Meta$Isovalerate)
P.RA1$Valeric<-as.numeric(Meta$Valerate)


#subset by time
P.RA1_32hr<-subset(P.RA1, Time=="32")
#Remove H43:
#P.RA3<-subset(P.RA3, Line!="H43")
#P.RA3<-subset(P.RA3, Line!="H38")
#Remove non-QPP
#P.RA3<-subset(P.RA3, Line2!="H1")
#P.RA3<-subset(P.RA3, Line2!="H2")
#P.RA3<-subset(P.RA3, Line2!="H3")




C<-ggplot(data=P.RA1, aes(x=Treatment,y=Coprococcus3,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Coprococcus 3")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
C+stat_compare_means()

C1<-ggplot(data=P.RA1, aes(x=Treatment,y=Coprococcus1,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Coprococcus1")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
C1+stat_compare_means()

R<-ggplot(data=P.RA1, aes(x=Treatment,y=Roseburia,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Roseburia")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
R+stat_compare_means()

B<-ggplot(data=P.RA1, aes(x=Treatment,y=Butyrate,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Butyrate")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
my_comparisons<-list(c("Media-32hr","Lysine-32hr"),c("Lysine-32hr","FL-32hr"),c("Media-32hr","FL-32hr"))
B+stat_compare_means(comparisons = my_comparisons, label.y=c(4.2,4.2,4.5),method="wilcox.test")

#Propionate/Acetate,NoDifferences
P<-ggplot(data=P.RA1, aes(x=Treatment,y=Propionate,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Propionate")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
P+stat_compare_means()
A<-ggplot(data=P.RA1, aes(x=Treatment,y=Acetate,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Acetate")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
A+stat_compare_means()

#Correlations one line per substrate
ggplot(data=P.RA1,aes(x=Coprococcus1,y=Butyrate,color=Substrate))+scale_color_manual(values = c("grey","firebrick","cornflowerblue"))+geom_point(size=2,aes(color=Substrate,shape=Time))+geom_smooth(method = "lm",se=TRUE,size=1.5)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  theme_classic()+theme(text=element_text(size=20,family = "sans"),plot.margin = margin(10, 20, 14, 5))+ylab("Butyrate (mMol)")+
  xlab("Coprococcus1")+theme(legend.position=c(.9,.12),axis.line = element_line(size=1.5)) 

#single line:
ggplot(data=P.RA1_32hr,aes(x=Roseburia,y=Butyrate))+scale_color_manual(values = c("grey","firebrick","cornflowerblue"))+geom_point(size=2,aes(color=Substrate))+geom_smooth(method = "lm",se=TRUE,size=1.5,colour="black")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  theme_classic()+theme(text=element_text(size=20,family = "sans"),plot.margin = margin(10, 20, 14, 5))+ylab("Butyrate (mMol)")+xlab("Roseburia")+theme(legend.position="none",axis.line = element_line(size=1.5)) 




#OtherBugs
A<-ggplot(data=P.RA1, aes(x=Treatment,y=Bacteroides,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Bacteroides")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
A+stat_compare_means()
A<-ggplot(data=P.RA1, aes(x=Treatment,y=Faecalibacterium,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Faecalibacterium")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
A+stat_compare_means()

A<-ggplot(data=P.RA1, aes(x=Treatment,y=Ruminiclostridium9,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Ruminiclostridium9")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
A+stat_compare_means()

A<-ggplot(data=P.RA1, aes(x=Treatment,y=Phascolarctobacterium,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Phascolarctobacterium")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
A+stat_compare_means()

A<-ggplot(data=P.RA1, aes(x=Treatment,y=Fusobacterium,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Fusobacterium")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
A+stat_compare_means()

A<-ggplot(data=P.RA1, aes(x=Treatment,y=Bilophila,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Bilophila")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
A+stat_compare_means()

A<-ggplot(data=P.RA1, aes(x=Treatment,y=Sutterella,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Sutterella")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
A+stat_compare_means()

A<-ggplot(data=P.RA1, aes(x=Treatment,y=Enterobacteriaceae,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Enterobacteriaceae")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
A+stat_compare_means()



#Family Level
A<-ggplot(data=P.RA1, aes(x=Treatment,y=Lachnospiraceae,fill=Substrate))+  geom_boxplot()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Lachnospiraceae")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
A+stat_compare_means()

#Correlations one line per substrate
ggplot(data=P.RA1,aes(x=Lachnospiraceae,y=Butyrate,color=Substrate))+scale_color_manual(values = c("grey","firebrick","cornflowerblue"))+geom_point(size=2,aes(color=Substrate,shape=Time))+geom_smooth(method = "lm",se=TRUE,size=1.5)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  theme_classic()+theme(text=element_text(size=20,family = "sans"),plot.margin = margin(10, 20, 14, 5))+ylab("Butyrate (mMol)")+xlab("Lachnospiraceae")+theme(legend.position=c(.9,.12),axis.line = element_line(size=1.5)) 

#single line:
ggplot(data=P.RA1_32hr,aes(x=Lachnospiraceae,y=Butyrate))+scale_color_manual(values = c("grey","firebrick","cornflowerblue"))+geom_point(size=2,aes(color=Substrate))+geom_smooth(method = "lm",se=TRUE,size=1.5,colour="black")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  theme_classic()+theme(text=element_text(size=20,family = "sans"),plot.margin = margin(10, 20, 14, 5))+ylab("Butyrate (mMol)")+xlab("Lachnospiraceae")+theme(legend.position="none",axis.line = element_line(size=1.5)) 

P.RA1$Time<-as.character(P.RA1$Time)
P.RA1$Time<-as.numeric(P.RA1$Time)
ggplot(data=P.RA1,aes(x=Time,y=Butyrate,color=Substrate))+scale_color_manual(values = c("grey","firebrick","cornflowerblue"))+geom_point(size=2,aes(color=Substrate))+geom_smooth(se=TRUE,size=1.5)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  theme_classic()+theme(text=element_text(size=20,family = "sans"),plot.margin = margin(10, 20, 14, 5))+ylab("Butyrate (mMol)")+xlab("Time")+theme(legend.position=c(.9,.12),axis.line = element_line(size=1.5)) 


#GenusLevel bar plots
B<-ggplot(data=P.RA1, aes(x=Treatment,y=Butyrate,fill=Substrate))+  geom_bar()+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab("Butyrate")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.9,.12),axis.line = element_line(size=1.5))
my_comparisons<-list(c("Media-32hr","Lysine-32hr"),c("Lysine-32hr","FL-32hr"),c("Media-32hr","FL-32hr"))
B+stat_compare_means(comparisons = my_comparisons, label.y=c(4.2,4.2,4.5),method="wilcox.test")

#remove time 0 and 48:
P.RA2<-subset(P.RA1,Time!="0")
P.RA2<-subset(P.RA2,Time!="48")

my_comparisons<-list(c("Media","Lysine"),c("Lysine","FL"),c("Media","FL"))

B2<-ggplot(data=P.RA2, aes(x=Substrate,y=Butyrate,fill=Substrate))+  geom_bar(stat="summary",fun.y = "mean",position = "dodge",width=1)+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab(NULL)+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_blank())+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position="none",axis.line = element_line(size=1.5))+facet_wrap(~Time,ncol=3)+geom_point()+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),breaks = c(1,2,3,4))+theme(plot.margin=unit(c(0.5,0,0,1.2),"cm"))

B2<-B2+stat_compare_means(comparisons = my_comparisons, label.y=c(4.1,4.1,4.45),method="wilcox.test")
B2
C21<-ggplot(data=P.RA2, aes(x=Substrate,y=Coprococcus1,fill=Substrate))+  geom_bar(stat="summary",fun.y = "mean",position = "dodge",width=1)+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab(NULL)+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_blank())+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position="none",axis.line = element_line(size=1.5))+facet_wrap(~Time,ncol=3)+geom_point()+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),breaks = c(.03,.06,.09))

C21<-C21+stat_compare_means(comparisons = my_comparisons, label.y=c(.102,.102,.111),method="wilcox.test")
C21
C23<-ggplot(data=P.RA2, aes(x=Substrate,y=Roseburia,fill=Substrate))+  geom_bar(stat="summary",fun.y = "mean",position = "dodge",width=1)+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
  ylab(NULL)+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_blank())+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position="none",axis.line = element_line(size=1.5))+facet_wrap(~Time,ncol=3)+geom_point()+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),breaks = c(.03,.06,.09))

C23<-C23+stat_compare_means(comparisons = my_comparisons, label.y=c(.115,.115,.126),method="wilcox.test")
C23




ggarrange(ggarrange(B1,B2,labels = c("B","C")),ggarrange(C1,C21,labels = c("D","E")),ggarrange(C31,C23,labels = c("F","G")),nrow=3)
#########################################################
#C23<-ggplot(data=P.RA2, aes(x=Substrate,y=Coprococcus3,fill=Substrate))+  geom_bar(stat="summary",fun.y = "mean",position = "dodge",width=1)+scale_fill_manual(values = c("grey","firebrick","cornflowerblue"))+
ylab(NULL)+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_blank())+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=NULL,axis.line = element_line(size=1.5))+facet_wrap(~Time,ncol=3)+#geom_point()+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),breaks = c(.01,.02,.03))







#calculate beta diversity:

table.beta<-read.delim("./exported/PPP4_OTU_Table_silva_taxonomy.tsv")
table.beta0<-table.beta[2:61]
rownames(table.beta0)<-table.beta$X.OTU.ID

#betajaccard<-vegdist(t(table.beta0),method="jaccard")
betajaccard<-vegdist(t(B.OTU),method="jaccard")
#betajaccard<-vegdist(t(B.OTU),method="euclidian")
jaccardmds<-metaMDS(betajaccard)
jmdsdata<-as.data.frame(jaccardmds$points)
#jaccardmds$icause

jmdsdata1<-cbind(Meta,jmdsdata)
jmdsdata1$Substrate<-factor(Meta$Substrate,levels=c("Media","Lysine","FL"))
jmdsdata1.1<-subset(jmdsdata1,Time!="48")
jmdsdata2<-subset(jmdsdata1,Time!="0")
jmdsdata3<-subset(jmdsdata2,Time!="48")
jmdsdata4<-subset(jmdsdata3,Time!="16")
jmdsdata4.1<-subset(jmdsdata3,Time!="32")

#alltime
gg <- merge(jmdsdata3,aggregate(cbind(mean.x=MDS1,mean.y=MDS2)~Treatment,jmdsdata3,mean),by="Treatment")


ggplot(gg, aes(x=MDS1,y=MDS2,color=Substrate,shape=Time))+geom_point(size=4)+theme_classic()+scale_color_manual(values = c("grey","firebrick","cornflowerblue"))+
  geom_point(aes(x=mean.x,y=mean.y))+geom_segment(aes(x=mean.x, y=mean.y, xend=MDS1, yend=MDS2),size=1)+theme(text=element_text(size=25,family = "sans"))+ 
  theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=NULL,axis.line = element_line(size=1.5))



#16hr
ggplot(jmdsdata4, aes(x=MDS1,y=MDS2,color=Substrate))+geom_point()

#32hr
ggplot(jmdsdata4.1, aes(x=MDS1,y=MDS2,color=Substrate))+geom_point()


#plot family level barcharts:
library(viridis)
#table.gamma<-read.delim("./exported/PPP4_OTU_family_Table_silva_taxonomy2_.5perTrim.tsv.txt")
table.gamma<-read.delim("./exported/PPP4_OTU_Phylum_Table_silva_taxonomy2_.5perTrim.tsv.txt")
FamTable<-subset(table.gamma,Time!=48)
FamTable$Treatment<-factor(FamTable$Treatment,levels=c("Media-0hr","Media-16hr","Media-32hr","Lysine-16hr","Lysine-32hr","FL-16hr","FL-32hr"))


A<-ggplot(FamTable,aes(fill=Phylum,y=Abundance,x=Treatment))+geom_bar(position="stack", stat="identity")+theme_classic()+
  theme(text=element_text(size=24,family = "sans"),axis.ticks.x=element_blank(),axis.line = element_line(size=1.5),plot.margin = margin(10,0,0,0))+xlab(NULL)+
  ylab(NULL)+scale_y_continuous(expand = c(0,0))

A+scale_fill_viridis(discrete = TRUE,option="H")#+facet_wrap(~Time,ncol=3)
