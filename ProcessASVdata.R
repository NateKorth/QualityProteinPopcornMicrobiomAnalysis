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
library(tidyverse)
##Prepare data for correlation analysis at ASV level, remove negative controls
setwd("~/PopGhum/Popghum2")

tree=read_tree("exported/PPP2_tree.nwk")

table.alpha<-read.delim("./exported/PPP2_OTU_Table_silva_taxonomy.tsv")
table.alpha0<-table.alpha[-c(1:20,382:387)]

table.alpha0[1:4,1:4]

#remove ID column and add it as column names
#A.OTU<-table.alpha[-c(1)]
A.OTU<-table.alpha0
OTU_IDs<-as.data.frame(table.alpha[1])
rownames(A.OTU)<-OTU_IDs[,1]

A.OTU[1:3,1:37]
A.OTU$Consensus.Lineage<-NULL

#trim low abundance taxa:
B.OTU<-A.OTU
B.OTU$Undetermined<-NULL
B.OTU$AveReads<-rowMeans(B.OTU)
B.OTU<-B.OTU[-which(B.OTU$AveReads<=1),]
B.OTU$AveReads<-NULL

###################################################################
#read in Silva taxonomy | before, in excel seperate taxonomy into columns
taxa<-read.delim("./exported/PPP2_silva_taxonomy.tsv",stringsAsFactors = FALSE,sep="\t")
taxa<-taxa[c(1:9)]
ASV.num<-as.data.frame(paste("ASV",seq(1,2425, by=1),sep=""))
names(ASV.num)<-"ASV.num"
ASV.name<-as.data.frame(taxa$Genus)
names(ASV.name)<-"ASV.name"
taxa$ASV.IDs<-paste0(ASV.num$ASV.num,"_",ASV.name$ASV.name,sep="")
taxa$ASV.IDs

taxa1<-taxa[which(taxa$Feature_ID %in% row.names(B.OTU)),]

#ordB<- match(rownames(B.OTU), taxa1$Feature_ID)
#taxa1 <- taxa1[ordB,]

OTUdata<-AnnotatedDataFrame(taxa1)
rownames(OTUdata)<-taxa1$Feature_ID

#normalize data
#load metadata
Meta<-load_phenoData("PPP2_metadata.tsv.txt",sep="\t")
#name the columns in Meta data file
colnames(Meta)<-c("Line","Line2","Group","Plate","Digest","Rep","Subject","Plant","ProteinContent","StarchContent","Prot/Starch","Acetate","Propionate","Isobutyrate","Butyrate","Isovalerate","Valeric","Lysine","Methionine",
                  "Lysine_TotalProtein","InvSimpson","Shannon","Subject2")

#Match names in meta data to OTU table
ordA<- match(colnames(B.OTU), rownames(Meta))
Ameta <- Meta[ordA,]

#insert subset into phenodataframe:
phenotypeDataA<-AnnotatedDataFrame(Ameta)

#Make MetagenomeSeq object
A.data<-newMRexperiment(B.OTU, phenoData = phenotypeDataA, featureData = OTUdata)

#trim that shit yo: taxa only present in 20% of taxa removed--Change this number based on how many samples you have!
A.datatrim<-filterData(A.data,present = 72)

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

ALine<-pData(A.datatrim)$Line
ALine2<-pData(A.datatrim)$Line2
ASubject<-pData(A.datatrim)$Subject
AGroup<-pData(A.datatrim)$Group
APlant<-pData(A.datatrim)$Plant
AProtein<-pData(A.datatrim)$ProteinContent
AStarch<-pData(A.datatrim)$StarchContent
ADigest<-pData(A.datatrim)$Digest
AAcetate<-pData(A.datatrim)$Acetate
APropionate<-pData(A.datatrim)$Propionate
AButyrate<-pData(A.datatrim)$Butyrate
ALysine<-pData(A.datatrim)$Lysine
AMethionine<-pData(A.datatrim)$Methionine
ALysineR<-pData(A.datatrim)$Lysine_TotalProtein
AInvSimpson<-pData(A.datatrim)$InvSimpson
AShannon<-pData(A.datatrim)$Shannon
ASubj2<-pData(A.datatrim)$Subject2
#####################################################################################
#Pick significan OTUs:
#define normalization factor
norm.factor <- normFactors(A.datatrim)
norm.factor <- log2(norm.factor/median(norm.factor)+1)

mod<- model.matrix(~APlant+ASubject+norm.factor)
settings<-zigControl(maxit=10, verbose=TRUE)
fit<-fitZig(obj=A.datatrim,mod=mod,useCSSoffset = TRUE,control=settings)

#table of log fold change coefficients for each OTU-- sort by adjusted p-values
coefs<-MRcoefs(fit,coef=NULL,group=3,number=200)
cc.sig.otus<-coefs[-which(coefs$adjPvalues>=.01),]
sigIDs<-as.data.frame((rownames(cc.sig.otus)))
names(sigIDs)<-"Feature_ID"
sigtaxa<-taxa1[,c("Feature_ID","ASV.IDs")]


sigtaxa<-sigtaxa[which(sigtaxa$Feature_ID %in% sigIDs$Feature_ID),]
ordC<-match(sigtaxa$Feature_ID,row.names(cc.sig.otus))
sigtaxa<-sigtaxa[ordC,]

sigtaxa<-cbind(sigtaxa,cc.sig.otus$adjPvalues)

######################################
#prep data for plotting:
#choose read counts or relative abundance
#P.RA<-PPP_normal_log
P.RA<-PPP.relabun.trim
taxa2<-taxa1[which(taxa1$Feature_ID %in% row.names(P.RA)),]

ASVIDs<-str_replace_all(taxa2$ASV.IDs, "`","")
rownames(P.RA)<-ASVIDs
P.RA1<-as.data.frame(t(P.RA))
ALine1<-as.data.frame(ALine)
P.RA1$Line<-ALine1$ALine
P.RA1$Line<-factor(P.RA1$Line, levels= c( "P1xP2","P2xP1","P1xP3","H20","H25","H28","H38","H43"))

AGroup1<-as.data.frame(AGroup)
P.RA1$Group<-AGroup1$AGroup
P.RA1$Group<-factor(P.RA1$Group, levels= c( "Non-QPP","QPP","QPP-NR"))


ALine21<-as.data.frame(ALine2)
P.RA1$Line2<-ALine21$ALine2
P.RA1$Line2<-factor(P.RA1$Line2, levels= c( "FBB0","FBB16","H1","H2","H3","QPP1A","QPP1B","QPP2","QPP3A","QPP3B"))

ASubject1<-as.data.frame(ASubject)
P.RA1$Subject<-ASubject1$ASubject
APlant1<-as.data.frame(APlant)
P.RA1$Plant<-APlant1$APlant
P.RA1$Plant<-factor(P.RA1$Plant, levels= c("Non-QPP","QPP"))
AProtein1<-as.data.frame(AProtein)
P.RA1$ProteinContent<-AProtein1$AProtein
AStarch1<-as.data.frame(AStarch)
P.RA1$StarchContent<-AStarch1$AStarch
ADigest1<-as.data.frame(ADigest)
P.RA1$Digest<-ADigest1$ADigest
AAcetate1<-as.data.frame(AAcetate)
P.RA1$Acetate<-AAcetate1$AAcetate
AButyrate1<-as.data.frame(AButyrate)
P.RA1$Butyrate<-AButyrate1$AButyrate
APropionate1<-as.data.frame(APropionate)
P.RA1$Propionate<-APropionate1$APropionate
ALysine1<-as.data.frame(ALysine)
P.RA1$Lysine<-ALysine1$ALysine
AMethionine1<-as.data.frame(AMethionine)
P.RA1$Methionine<-AMethionine1$AMethionine
ALysineR1<-as.data.frame(ALysineR)
P.RA1$LysineR<-ALysineR1$ALysineR

AInvSimpson<-as.data.frame(AInvSimpson)
P.RA1$InvSimpson<-AInvSimpson$AInvSimpson
AShannon<-as.data.frame(AShannon)
P.RA1$Shannon<-AShannon$AShannon
ASubj2<-as.data.frame(ASubj2)
P.RA1$Subject2<-ASubj2$ASubj2
P.RA1$norm<-norm.factor
#All data with FBBs
#P.RA2<-subset(P.RA1, Proc=="Boiled" | Proc=="Popped" | Proc=="FBB16" | Proc=="FBB0" | Proc=="FBB24")

#All data without FBBs
P.RA3<-subset(P.RA1, Plant=="QPP"|Plant=="Non-QPP")
#With FBBs:
#P.RA3<-P.RA1
#Remove H43:
#P.RA3<-subset(P.RA3, Line!="H43")
#P.RA3<-subset(P.RA3, Line!="H38")
#Remove non-QPP
#P.RA3<-subset(P.RA3, Line2!="H1")
#P.RA3<-subset(P.RA3, Line2!="H2")
#P.RA3<-subset(P.RA3, Line2!="H3")
#P.RA3$ASV
#Plot taxa:
ggplot(data=P.RA3X, aes(x=Line2,y=ASV13_Roseburia,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("Brown4","cornflowerblue","darkorchid","Red"))+
  ylab("Relative Abundance")+theme_classic()+theme(text=element_text(size=20,family = "sans"))+ggtitle("ASV13_Roseburia")+
  scale_color_manual(values=c("orange", "black","Tan","yellow","red","green4","gold2", "grey"))+facet_wrap(~Subject2, ncol=5)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())


P.RA3$ProteinContent<-as.numeric(P.RA3$ProteinContent)
P.RA3$StarchContent<-as.numeric(P.RA3$StarchContent)
P.RA3$PSRatio<-P.RA3$ProteinContent/P.RA3$StarchContent
P.RA3$PSRatio<-as.character(P.RA3$PSRatio)
P.RA3$Acetate<-as.numeric(P.RA3$Acetate)
P.RA3$Butyrate<-as.numeric(P.RA3$Butyrate)
P.RA3$Propionate<-as.numeric(P.RA3$Propionate)
P.RA3$Lysine<-as.numeric(P.RA3$Lysine)
P.RA3$Methionine<-as.numeric(P.RA3$Methionine)
P.RA3$LysineR<-as.numeric(P.RA3$LysineR)
P.RA3$InvSimpson<-as.numeric(P.RA3$InvSimpson)
P.RA3$Shannon<-as.numeric(P.RA3$Shannon)

P.RA3X<-subset(P.RA3,Subject!="Pool")
#P.RA3X<-subset(P.RA3X,Subject!="S767")
#P.RA3X<-subset(P.RA3X,Subject!="S775")
#P.RA3X<-subset(P.RA3X,Subject!="S776")
#P.RA3X<-subset(P.RA3X,Subject!="S782")

Q<-ggplot(data=P.RA3X, aes(x=Line2,y=Butyrate,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("Butyrate mMol")+theme_classic()+theme(text=element_text(size=20,family = "sans"))+xlab(NULL)+
  scale_color_manual(values=c("orange", "black","Tan","yellow","red","green4","gold2", "grey"))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(0.13,.16))+facet_wrap(~Subject2, ncol=5)#+theme(axis.text.x=element_blank())

my_comparisons<-list(c("H1","QPP1A"),c("H1","QPP1B"),c("H2","QPP2"),c("H3","QPP3A"),c("H3","QPP3B"))
Q+stat_compare_means(comparisons = my_comparisons, label.y=c(11.8,12.3,12.8,13.3,13.8),method="wilcox.test")



#subset by subject:
P.RA3.Pool<-subset(P.RA3, Subject=="Pool")
P.RA3.S767<-subset(P.RA3, Subject=="S767")
P.RA3.S767L<-log(P.RA3.S767[1:149]+0.000001)
P.RA3.S767L<-cbind(P.RA3.S767[150:156],P.RA3.S767L)

P.RA3.S775<-subset(P.RA3, Subject=="S775")
P.RA3.S776<-subset(P.RA3, Subject=="S776")
P.RA3.S782<-subset(P.RA3, Subject=="S782")



#Fit linear mixed model:
fit1<-mmer(ASV21_Bifidobacterium~1, random=~Digest:Line+PSRatio, rcov=~units, data=P.RA3.S767)
summary(fit1)
VC <- summary(fit1)$varcomp
VL <- VC[1,1]
VD <- VC[2,1]
VP <- VC[3,1]
VS <- VC[4,1]
VR <- VC[5,1]
VL1<-VL/(VL+VD+VP+VS+VR)
VD1<-VD/(VL+VD+VP+VS+VR)
VP1<-VP/(VL+VD+VP+VS+VR)
VS1<-VS/(VL+VD+VP+VS+VR)


h2_df <- as.data.frame(cbind("ASV1_Bacteroides","S767",VL1,VD1,VP1,VS1))
names(h2_df)<-c("ASV","Subject","PlantV","DigestBatchV","TotProteinV","TotStarchV")

h2_list<-data.frame(h2_df)
#h2_list<-rbind(h2_list,data.frame(h2_df))

geth2<-function(genus){
  t <- formula(paste0(genus, '~1'))
  #fit2<-mmer(t, random=~Plant+Digest:Line, rcov=~units, data= P.RA3.S782)
  fit2<-mmer(t, random=~Digest:Line+PSRatio, rcov=~units, data= P.RA3.S782)
  VC <- summary(fit2)$varcomp
  #VL <- VC[1,1]
  #VD <- VC[2,1]
  #VR <- VC[3,1]
  VD <- VC[1,1]
  VPSR <- VC[2,1]
  #VS <- VC[3,1]
  VR <- VC[3,1]
  #VL1<-VL/(VL+VD+VR)
  VSP1<-(VPSR)/(VD+VPSR+VR)
  #h2_df <- as.data.frame(cbind(genus,"S782",VL1))
  h2_df <- as.data.frame(cbind(genus,"S782",VSP1))
  #names(h2_df)<-c("ASV","Subject","PlantV")
  names(h2_df)<-c("ASV","Subject","TotalProteinStarchV")
  return(h2_df)
}

h2_listL<-data.frame()
h2_listSP<-data.frame()
#Format Lists of Taxa to get h for:

P.RA3.PoolT<-as.data.frame(t(P.RA3.Pool[c(1:149)]))
P.RA3.PoolT$AveReads<-rowMeans(P.RA3.PoolT)
P.RA3.PoolT<-P.RA3.PoolT[-which(P.RA3.PoolT$AveReads<=0.01),]
P.RA3.PoolT$AveReads<-NULL
PoolASVs<-as.data.frame(row.names(P.RA3.PoolT))
names(PoolASVs)<-c("taxa")

P.RA3.S767T<-as.data.frame(t(P.RA3.S767[c(1:149)]))
P.RA3.S767T$AveReads<-rowMeans(P.RA3.S767T)
P.RA3.S767T<-P.RA3.S767T[-which(P.RA3.S767T$AveReads<=0.01),]
P.RA3.S767T$AveReads<-NULL
S767ASVs<-as.data.frame(row.names(P.RA3.S767T))
names(S767ASVs)<-c("taxa")

P.RA3.S775T<-as.data.frame(t(P.RA3.S775[c(1:149)]))
P.RA3.S775T$AveReads<-rowMeans(P.RA3.S775T)
P.RA3.S775T<-P.RA3.S775T[-which(P.RA3.S775T$AveReads<=0.01),]
P.RA3.S775T$AveReads<-NULL
S775ASVs<-as.data.frame(row.names(P.RA3.S775T))
names(S775ASVs)<-c("taxa")

P.RA3.S776T<-as.data.frame(t(P.RA3.S776[c(1:149)]))
P.RA3.S776T$AveReads<-rowMeans(P.RA3.S776T)
P.RA3.S776T<-P.RA3.S776T[-which(P.RA3.S776T$AveReads<=0.01),]
P.RA3.S776T$AveReads<-NULL
S776ASVs<-as.data.frame(row.names(P.RA3.S776T))
names(S776ASVs)<-c("taxa")

P.RA3.S782T<-as.data.frame(t(P.RA3.S782[c(1:149)]))
P.RA3.S782T$AveReads<-rowMeans(P.RA3.S782T)
P.RA3.S782T<-P.RA3.S782T[-which(P.RA3.S782T$AveReads<=0.01),]
P.RA3.S782T$AveReads<-NULL
S782ASVs<-as.data.frame(row.names(P.RA3.S782T))
names(S782ASVs)<-c("taxa")



#for loop to generate master list of h's
for (i in S782ASVs$taxa){
  h2_listSP<-rbind(h2_listSP,geth2(i))
}

h2_listC<-cbind(h2_listL,h2_listSP$TotalProteinStarchV)
names(h2_listC)<-c("ASV","Subject","QPV","Protein:StarchV")
h2_listC$QPV<-as.numeric(h2_listC$QPV)
h2_listC$`Protein:StarchV`<-as.numeric(h2_listC$`Protein:StarchV`)
#write.csv(h2_listC,"VariationExplained1.csv",quote = FALSE,row.names = FALSE)
h2_listC<-read.csv("VariationExplained1.csv")
h2_listC1<-h2_listC[-which(h2_listC$QPV==0),]


g1<-ggplot(h2_listC1, aes(x=`Protein:StarchV`,y=QPV)) +geom_point()
g1+ xlim(0,0.5)+theme_classic()+theme(text=element_text(size=30,family = "serif"))+ggtitle("Variation Explained")+
  geom_text(label=ifelse(h2_listC1$ASV=="ASV15_Coprococcus3"|h2_listC1$ASV=="ASV44_Butyricicoccus",as.character(c(h2_listC1$ASV,h2_listC1$Subject)),""),hjust=0.1,vjust=0)

#correlation of taxa to scfa
library(corrplot)
#get traits to correlate
P.RA3.S767T1<-cbind(t(P.RA3.S767T),P.RA3.S767[155:162],P.RA3.S767[164:165])
P.RA3.S767T1$Digest<-NULL
P.RA3.S767T1$ASV39_Lachnospiraceae_unknown<-NULL
P.RA3.S767T1$Plant<-P.RA3.S767$Plant
P.RA3.S767TQP<-subset(P.RA3.S767T1,Plant=="Non-QPP")
P.RA3.S767TNQP<-subset(P.RA3.S767T1,Plant!="Non-QPP")  
P.RA3.S767T1$Plant<-NULL
P.RA3.S767TQP$Plant<-NULL
P.RA3.S767TNQP$Plant<-NULL
#S767C<-cor(P.RA3.S767T1)
#S767C<-cor(P.RA3.S767TQP)
S767C<-cor(P.RA3.S767TNQP)
#sigS767C<-cor.mtest(P.RA3.S767T1,conf.level=0.95)
corrplot(S767C,method="color",type="lower",order="hclust",tl.col="black",tl.srt = 45, cex.var=.7 ,col=colorRampPalette(c("mediumorchid4","blueviolet","white","gold","yellow"))(200))#,p.mat=sigS767C$p,sig.level=c(.001,.01,.05),insig="label_sig",pch.cex = .7)

P.RA3.S775T1<-cbind(t(P.RA3.S775T),P.RA3.S775[155:162],P.RA3.S767[164:165])
P.RA3.S775T1$Digest<-NULL

P.RA3.S775T1$Plant<-P.RA3.S775$Plant
P.RA3.S775TQP<-subset(P.RA3.S775T1,Plant=="Non-QPP")
P.RA3.S775TNQP<-subset(P.RA3.S775T1,Plant!="Non-QPP")  
P.RA3.S775T1$Plant<-NULL
P.RA3.S775TQP$Plant<-NULL
P.RA3.S775TNQP$Plant<-NULL
#S775C<-cor(P.RA3.S775T1)
#S775C<-cor(P.RA3.S775TQP)
S775C<-cor(P.RA3.S775TNQP)

corrplot(S775C,method="color",type="lower",order="hclust",tl.col="black",tl.srt = 45,col=colorRampPalette(c("firebrick4", "firebrick3","firebrick1","white","royalblue2","royalblue3","royalblue4"))(200))

P.RA3.S776T1<-cbind(t(P.RA3.S776T),P.RA3.S776[155:162],P.RA3.S767[164:165])
P.RA3.S776T1$Digest<-NULL

P.RA3.S776T1$Plant<-P.RA3.S776$Plant
P.RA3.S776TQP<-subset(P.RA3.S776T1,Plant=="Non-QPP")
P.RA3.S776TNQP<-subset(P.RA3.S776T1,Plant!="Non-QPP")  
P.RA3.S776T1$Plant<-NULL
P.RA3.S776TQP$Plant<-NULL
P.RA3.S776TNQP$Plant<-NULL
#S776C<-cor(P.RA3.S776T1)
#S776C<-cor(P.RA3.S776TQP)
S776C<-cor(P.RA3.S776TNQP)

corrplot(S776C,method="color",type="lower",order="hclust",tl.col="black",tl.srt = 45,col=colorRampPalette(c("darkorange3", "darkorange2","darkorange1","white","green2","green3","darkgreen"))(200))

P.RA3.S782T1<-cbind(t(P.RA3.S782T),P.RA3.S782[155:162],P.RA3.S767[164:165])
P.RA3.S782T1$Digest<-NULL

P.RA3.S782T1$Plant<-P.RA3.S782$Plant
P.RA3.S782TQP<-subset(P.RA3.S782T1,Plant=="Non-QPP")
P.RA3.S782TNQP<-subset(P.RA3.S782T1,Plant!="Non-QPP")  
P.RA3.S782T1$Plant<-NULL
P.RA3.S782TQP$Plant<-NULL
P.RA3.S782TNQP$Plant<-NULL
#S782C<-cor(P.RA3.S782T1)
#S782C<-cor(P.RA3.S782TQP)
S782C<-cor(P.RA3.S782TNQP)

corrplot(S782C,method="color",type="lower",order="hclust",tl.col="black",tl.srt = 45,col=colorRampPalette(c("firebrick4", "firebrick3","firebrick1","white","royalblue2","royalblue3","royalblue4"))(200))

#Correlations of all subjects:
#P.RA3X1<-P.RA3X[-c(150:153)]
#P.RA3X1$PSRatio<-as.numeric(P.RA3X1$PSRatio)
#P.RA3X1$Digest<-as.numeric(P.RA3X1$Digest)
#P.RA3X1<-as.data.frame(t(P.RA3X1))
#P.RA3X1$AveReads<-rowMeans(P.RA3X1)
#P.RA3X1<-P.RA3X1[-which(P.RA3X1$AveReads<=0.01),]
#P.RA3.S782T$AveReads<-NULL
#AllSubC<-cor(as.data.frame(t(P.RA3X1)))
#corrplot(AllSubC,method="color",type="lower",order="hclust",tl.col="black",tl.srt = 45)
#Plot correlated traits:

a<-ggplot(data=P.RA3.S767, aes(x=Line2,y=ASV13_Roseburia,fill=Plant))+  geom_boxplot(lwd=1)+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("ASV13 Roseburia")+theme_classic()+theme(text=element_text(size=20,family = "sans"))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),axis.text.x =element_text(angle = 90,hjust=1),legend.position=c(.2,.8),axis.line = element_line(size=1.5))

b<-ggplot(data=P.RA3.S767, aes(x=Line2,y=Butyrate,fill=Plant))+  geom_boxplot(lwd=1)+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("Butyrate (mMol)")+theme_classic()+theme(text=element_text(size=20,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  theme(axis.ticks.x=element_blank())+xlab(NULL)+ theme(legend.position="none",axis.line = element_line(size=1.5))

c<-ggplot(data=P.RA3.S767, aes(x=Line2,y=ASV15_Coprococcus3,fill=Plant))+  geom_boxplot(lwd=1)+scale_fill_manual(values = c("firebrick2","cornflowerblue"))+
  ylab("ASV15 Coprococcus")+theme_classic()+theme(text=element_text(size=20,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  theme(axis.ticks.x=element_blank())+xlab(NULL)+ theme(legend.position="none",axis.line = element_line(size=1.5)) 

d<-ggplot(data=P.RA3.S767,aes(x=ASV15_Coprococcus3,y=Butyrate,color=Plant))+scale_color_manual(values = c("firebrick2","cornflowerblue"))+geom_point(size=2)+geom_smooth(method = "lm",se=FALSE,size=1.5)+
  theme_classic()+theme(text=element_text(size=20,family = "sans"))+xlab("ASV15 Coprococcus")+ylab(NULL)+
  stat_cor(method="pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+ theme(legend.position="none",axis.line = element_line(size=1.5)) +scale_x_continuous(breaks=c(.01,.02,.03))

e<-ggplot(data=P.RA3.S767,aes(x=ASV13_Roseburia,y=Butyrate,color=Plant))+scale_color_manual(values = c("firebrick2","cornflowerblue"))+geom_point(size=2)+geom_smooth(method = "lm",se=FALSE,size=1.5)+
  theme_classic()+theme(text=element_text(size=20,family = "sans"))+ylab("Butyrate (mMol)")+xlab("ASV13 Roseburia ")+
  stat_cor(method="pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+ theme(legend.position="none",axis.line = element_line(size=1.5),plot.margin = margin(10, 10, 5, 0)) +scale_x_continuous(breaks=c(.02,.05,.08))

g<-ggplot(data=P.RA3.S767,aes(x=LysineR,y=Butyrate))+scale_color_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+geom_point(size=2,aes(color=Plant))+geom_smooth(method = "lm",se=TRUE,size=1.5,colour="black")+#stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  theme_classic()+theme(text=element_text(size=20,family = "sans"),plot.margin = margin(10, 20, 14, 5))+ylab("Butyrate (mMol)")+xlab("pb Lysine/pb amino acids")+theme(legend.position="none",axis.line = element_line(size=1.5)) 

f<-ggplot(data=P.RA3.S767, aes(x=Line2,y=LysineR,fill=Plant))+geom_bar(stat="identity",position=position_dodge(),colour="black")+theme_classic()+theme(text=element_text(size=20,family = "sans"))+ylab("pb Lysine/pb amino acids")+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+ 
  theme(legend.position="none",axis.line = element_line(size=1.5),plot.margin = margin(2, 2, -24, 12),axis.ticks.x=element_blank(),axis.text.x =element_text(angle = 90,hjust=1))+xlab("")+
  scale_y_continuous(expand = c(0,0))


ggarrange(a,c,b,ggarrange(e,d),f,g,labels=c("A","B","C","D","E","F"))



ggarrange(a,b,e,f,g,labels=c("A","B","C","D","E"))


ggarrange(a,b,f,g,labels=c("A","B","C","D"))

#annotate_figure(fig1,top=text_grob("Relationship of Butyrate and taxa in subject 767",face="bold",size=18))

library(phylotools)#assign asvIDs to repseqs:
ASVID_Features<-as.data.frame(cbind(taxa$Feature_ID,taxa$ASV.IDs))
names(ASVID_Features)<-c("seq.name","new.name")
rename.fasta(infile = "./exported/PPP2_rep_seqs.fasta", ASVID_Features, outfile= "./exported/PPP2_rep_seqs2.fasta")
repseqs_trim<-read.fasta("./exported/PPP2_rep_seqs2.fasta")

repseqs_trim2<-repseqs_trim[which(repseqs_trim$seq.name %in% taxa2$ASV.IDs),]
dat2fasta(repseqs_trim2,outfile = "./exported/PPP2_repseqs_trim2")

citation()
#Starch and prot
#single values with sd:
SP_R<-read.csv("StrarchProt_forR2.csv")
#individualvalues for statistics
SP_R2<-read.csv("StrarchProt_forR.csv")

my_comparisons<-list(c("H1","QPP1A"),c("H1","QPP1B"),c("H2","QPP2"),c("H3","QPP3A"),c("H3","QPP3B"))
my_comparisons2<-list(c("H1","H2"),c("H1","H3"),c("H2","H3"),c("QPP1A","QPP1B"),c("QPP1A","QPP2"),c("QPP1A","QPP3A"),c("QPP1A","QPP3B"),c("QPP1B","QPP2"),c("QPP1B","QPP3A"),c("QPP1B","QPP3B"),c("QPP2","QPP3A"),c("QPP2","QPP3B"))

#SP_R<-subset(SP_R, Line!="H43")
#SP_R<-subset(SP_R, Line!="H38")
#my_comparisons<-list(c("H1","QPP1"),c("H2","QPP2"),c("H3","QPP3"))
#my_comparisons2<-list(c("H1","H2"),c("H1","H3"),c("H2","H3"),c("QPP1","QPP2"),c("QPP1","QPP3"),c("QPP2","QPP3"))


S<-ggplot(data=SP_R, aes(x=Line2,y=StarchContent,fill=Plant))+geom_bar(stat="identity",position = "dodge",width=.8)+theme_classic()+
  theme(text=element_text(size=24,family = "sans"))+ylab("Post digestion starch content")+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  geom_errorbar(aes(ymin=StarchContent-StarchSD, ymax=StarchContent+StarchSD),width=.2,position=position_dodge(.9))
S
wilcox.test(SP_R$StarchContent~SP_R$Plant)

Q+stat_compare_means(method="anova")
Q+stat_compare_means(comparisons = my_comparisons)
Q+stat_compare_means(comparisons = my_comparisons2)

pairwise.wilcox.test(SP_R$StarchContent,SP_R$Line2,p.adjust.method = "none")
X<-pairwise.wilcox.test(SP_R$StarchContent,SP_R$Line2,p.adjust.method = "none")
X.1<--log10(X$p.value)
corrplot(X.1,method="color",type="lower",tl.col="black",is.corr=FALSE,col=colorRampPalette(c("white","white","yellow"))(10))


pairwise.wilcox.test(SP_R2$ProteinContent_ALL,SP_R2$Line2,p.adjust.method = "none")
X<-pairwise.wilcox.test(SP_R2$ProteinContent_ALL,SP_R2$Line2,p.adjust.method = "none")
X.0<-X$p.value
X.1<--log10(X$p.value)
C<-corrplot(X.1,method="color",type="lower",tl.col="black",is.corr=FALSE,col=colorRampPalette(c("white","white","blue"))(10))

corrplot(as.matrix(X.1),is.corr=FALSE,type="lower",tl.col="black")


#Q<-ggplot(data=SP_R, aes(x=Line,y=ProteinContent,fill=Plant))+geom_bar(stat="identity",position=position_dodge())+theme_classic()+theme(text=element_text(size=24,family = "sans"))+ylab("Post digestion protein content")+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))
#wilcox.test(SP_R$ProteinContent~SP_R$Plant)

S<-ggplot(data=SP_R, aes(x=Line2,y=StarchContent,fill=Plant))+geom_bar(stat="identity",position=position_dodge())+theme_classic()+theme(text=element_text(size=22,family = "sans"))+ylab("Post digestion starch")+scale_fill_manual(values = c("firebrick2","cornflowerblue"))+ 
  theme(legend.position="none",axis.line = element_line(size=1.5),plot.margin = margin(10, 0, 0, 0),axis.ticks.x=element_blank(),axis.text.x =element_text(angle = 90,hjust=1))+xlab("")+
  scale_y_continuous(expand = c(0,0))+geom_errorbar(aes(ymin=StarchContent-StarchSD, ymax=StarchContent+StarchSD), width=.2,position=position_dodge(.9))#+coord_flip()
S
  
P<-ggplot(data=SP_R2, aes(x=Line2,y=ProteinContent_ALL,fill=Plant))+geom_bar(stat="summary",fun="mean",position=position_dodge())+theme_classic()+theme(text=element_text(size=22,family = "sans"))+ylab("Post digestion protein")+scale_fill_manual(values = c("firebrick2","cornflowerblue"))+ 
  theme(legend.position="none",axis.line = element_line(size=1.5),plot.margin = margin(10, 0, 0, 0),axis.ticks.x=element_blank(),axis.text.x =element_text(angle = 90,hjust=1))+xlab("")+
  scale_y_continuous(expand = c(0,0))+geom_errorbar(aes(ymin=ProteinContent-ProtSD, ymax=ProteinContent+ProtSD), width=.2,position=position_dodge(.9))#+coord_flip()
P

ggplot(data=SP_R, aes(x=Line2,y=ProteinContent,fill=Plant))+theme_classic()+theme(text=element_text(size=20,family = "sans"))+ylab("Post digestion protein content")+scale_fill_manual(values = c("firebrick2","cornflowerblue"))+ 
  theme(legend.position="none",axis.line = element_line(size=1.5),plot.margin = margin(0, 0, 0, 0),axis.ticks.x=element_blank(),axis.text.x =element_text(angle = 90,hjust=1))+xlab("")+
  scale_y_continuous(expand = c(0,0))+stat_compare_means(comparisons = my_comparisons)

ggarrange(P,S,labels = c("A","B"))

P+stat_compare_means(comparisons = my_comparisons)
P+stat_compare_means(comparisons = my_comparisons2)

ggplot(data=P.RA3.S767, aes(x=Line2,y=LysineR,fill=Plant))+geom_bar(stat="identity",position=position_dodge())+theme_classic()+theme(text=element_text(size=24,family = "sans"))+ylab("Post digestion protein content")+scale_fill_manual(values = c("firebrick2","cornflowerblue"))
#plot lysine by plant
ggplot(data=P.RA3.S767, aes(x=Plant,y=LysineR,fill=Plant))+geom_boxplot()+theme_classic()+theme(text=element_text(size=24,family = "sans"),legend.position="None")+ylab("Lysine content")+xlab(NULL)+scale_fill_manual(values = c("firebrick2","cornflowerblue"))

  
 

#Pull Bacteroides/Ruminococcus ASVs:
P.RA3.S767Bact<-P.RA3.S767 %>% select(contains("Bacteroides"))
P.RA3.S767Bact2<-P.RA3.S767 %>% select(contains("Ruminococ"))
P.RA3.S767Bact<-cbind(P.RA3.S767Bact,P.RA3.S767Bact2)
P.RA3.S767Bact$Plant<-P.RA3.S767$Plant
P.RA3.S767Bact.meds<-P.RA3.S767Bact %>% group_by(Plant) %>% summarise_each(funs(mean(.,na.rm = TRUE)))
P.RA3.S767Bact.meds$Subject<-c("S767","S767")

P.RA3.S775Bact<-P.RA3.S775 %>% select(contains("Bacteroides"))
P.RA3.S775Bact2<-P.RA3.S775 %>% select(contains("Ruminococ"))
P.RA3.S775Bact<-cbind(P.RA3.S775Bact,P.RA3.S775Bact2)
P.RA3.S775Bact$Plant<-P.RA3.S775$Plant
P.RA3.S775Bact.meds<-P.RA3.S775Bact %>% group_by(Plant) %>% summarise_each(funs(mean(.,na.rm = TRUE)))
P.RA3.S775Bact.meds$Subject<-c("S775","S775")

P.RA3.S776Bact<-P.RA3.S776 %>% select(contains("Bacteroides"))
P.RA3.S776Bact2<-P.RA3.S776 %>% select(contains("Ruminococ"))
P.RA3.S776Bact<-cbind(P.RA3.S776Bact,P.RA3.S776Bact2)
P.RA3.S776Bact$Plant<-P.RA3.S776$Plant
P.RA3.S776Bact.meds<-P.RA3.S776Bact %>% group_by(Plant) %>% summarise_each(funs(mean(.,na.rm = TRUE)))
P.RA3.S776Bact.meds$Subject<-c("S776","S776")

P.RA3.S782Bact<-P.RA3.S782 %>% select(contains("Bacteroides"))
P.RA3.S782Bact2<-P.RA3.S782 %>% select(contains("Ruminococ"))
P.RA3.S782Bact<-cbind(P.RA3.S782Bact,P.RA3.S782Bact2)
P.RA3.S782Bact$Plant<-P.RA3.S782$Plant
P.RA3.S782Bact.meds<-P.RA3.S782Bact %>% group_by(Plant) %>% summarise_each(funs(mean(.,na.rm = TRUE)))
P.RA3.S782Bact.meds$Subject<-c("S782","S782")


P.RA_Bact_meds<-rbind(P.RA3.S767Bact.meds,P.RA3.S775Bact.meds,P.RA3.S776Bact.meds,P.RA3.S782Bact.meds)
#write.csv(P.RA_Bact_meds,"PPP2_ASVlvl_means1.csv")

#Plot Alpha Diversity:
AlphaMeta<-read.table("PPP2_metadata.tsv.txt",sep="\t",header = TRUE)
AlphaMeta1<-subset(AlphaMeta,Plant=="QPP"|Plant=="Non-QPP")
AlphaMeta2<-subset(AlphaMeta1,Subject!="Pool")
AlphaMeta1.767<-subset(AlphaMeta1,Subject=="S767")
AlphaMeta1.775<-subset(AlphaMeta1,Subject=="S775")
AlphaMeta1.776<-subset(AlphaMeta1,Subject=="S776")
AlphaMeta1.782<-subset(AlphaMeta1,Subject=="S782")

AlphaMeta1<-subset(AlphaMeta,Line2!="QPP"|Plant!="Non-QPP")


a<-ggplot(data=AlphaMeta2, aes(x=Plant,y=InvSimpson,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("Inverted Simpson")+theme_classic()+theme(text=element_text(size=22,family = "sans"))+
xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.2,.8),axis.line = element_line(size=1.5))+facet_wrap(~Subject2,ncol=5)

a<-a+stat_compare_means()

b<-ggplot(data=AlphaMeta1, aes(x=Plant,y=Shannon,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("Shannon")+theme_classic()+theme(text=element_text(size=22,family = "sans"))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.position="None",axis.line = element_line(size=1.5))+facet_wrap(~Subject2,ncol=5)
b<-b+stat_compare_means()
#,axis.text.x =element_text(angle = 90)

ggarrange(a,b)
a
b

#alpha diversity by line:
a<-ggplot(data=AlphaMeta2, aes(x=Line2,y=InvSimpson,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("Inverted Simpson")+theme_classic()+theme(text=element_text(size=22,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.2,.8),axis.line = element_line(size=1.5))+facet_wrap(~Subject2,ncol=5)

#a<-a+stat_compare_means(comparisons = my_comparisons)
a

b<-ggplot(data=AlphaMeta2, aes(x=Line2,y=Shannon,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("Shannon")+theme_classic()+theme(text=element_text(size=22,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.position="None",axis.line = element_line(size=1.5))+facet_wrap(~Subject2,ncol=5)
#b<-b+stat_compare_means()
b

ggarrange(a,b,labels = c("a","b"))


#Plot Butyrate and run ttests:
Q<-ggplot(data=P.RA3X, aes(x=Plant,y=Butyrate,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("Butyrate uMol")+theme_classic()+theme(text=element_text(size=17,family = "sans"))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.87,.62),axis.line = element_line(size=1.5))+facet_wrap(~Subject2,ncol=5)

Q+stat_compare_means()
###
#othertaxa
P.RA3X$ASV8_Bifidobacterium
Q<-ggplot(data=P.RA3X, aes(x=Plant,y=ASV8_Bifidobacterium,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("Butyrate uMol")+theme_classic()+theme(text=element_text(size=17,family = "sans"))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.87,.62),axis.line = element_line(size=1.5))+facet_wrap(~Subject2,ncol=5)

Q+stat_compare_means()
###

wilcox.test(P.RA3.S767$Butyrate~P.RA3.S767$Plant)
wilcox.test(P.RA3.S775$Butyrate~P.RA3.S775$Plant)
wilcox.test(P.RA3.S776$Butyrate~P.RA3.S776$Plant)
wilcox.test(P.RA3.S782$Butyrate~P.RA3.S782$Plant)

fit<-aov(Butyrate~Line2,data=P.RA3.S767)
summary(fit)
TukeyHSD(fit)

#Plot Butyrate and run ttests for line:
Q<-ggplot(data=P.RA3X, aes(x=Line2,y=Butyrate,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("Butyrate")+theme_classic()+theme(text=element_text(size=22,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position=c(.1,.2),axis.line = element_line(size=1.5))+facet_wrap(~Subject2,ncol=5)

my_comparisons<-list(c("H1","QPP1A"),c("H1","QPP1B"),c("H2","QPP2"),c("H3","QPP3A"),c("H3","QPP3B"))

Q+stat_compare_means(comparisons = my_comparisons, label.y = c(12.2,12.6,13,13.4,13.8))#+geom_point(aes(color=Digest))
#Q+stat_compare_means()

#Plot Acetate/Propionate
A<-ggplot(data=P.RA3X, aes(x=Line2,y=Acetate,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("Acetate")+theme_classic()+theme(text=element_text(size=22,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position="none",axis.line = element_line(size=1.5))+facet_wrap(~Subject2,ncol=5)

my_comparisons<-list(c("H1","QPP1A"),c("H1","QPP1B"),c("H2","QPP2"),c("H3","QPP3A"),c("H3","QPP3B"))

A<-A+stat_compare_means(comparisons = my_comparisons, label.y = c(56,57.5,59,60.5,62))

P<-ggplot(data=P.RA3X, aes(x=Line2,y=Propionate,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("Propionate")+theme_classic()+theme(text=element_text(size=22,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position="none",axis.line = element_line(size=1.5))+facet_wrap(~Subject2,ncol=5)

P<-P+stat_compare_means(comparisons = my_comparisons, label.y = c(9,9.4,9.8,10.2,10.6))

ggarrange(A,P,P,A)


#Plot Other taxa and run ttests for line:
Q<-ggplot(data=P.RA3.S767, aes(x=Line2,y=ASV26_Fusobacterium,fill=Plant))+  geom_boxplot()+scale_fill_manual(values = c("firebrick2","cornflowerblue","darkgreen"))+
  ylab("ASV26_Fusobacterium")+theme_classic()+theme(text=element_text(size=25,family = "sans"),axis.text.x =element_text(angle = 90,hjust=1))+
  xlab(NULL)+ theme(axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position='none',axis.line = element_line(size=1.5))+facet_wrap(~Subject2,ncol=5)

my_comparisons<-list(c("H1","QPP1"),c("H1","QPP1.1"),c("H2","QPP2"),c("H3","QPP3"),c("H3","QPP3.1"))
#my_comparisons<-list(c("H1","QPP1"),c("H2","QPP2"),c("H3","QPP3"))

#Q+stat_compare_means(comparisons = my_comparisons, label.y = c(12.2,12.6,13,13.4,13.8))
#Q+stat_compare_means(comparisons = my_comparisons, label.y=c(.033,.034,.035,.036,.037))
Q+stat_compare_means(comparisons = my_comparisons, label.y=c(.078,.081,.084,.087,.09))
P.RA3X$ASV167_ClostridialesFamily8
