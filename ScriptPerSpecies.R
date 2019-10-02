##################################################
## Project: GC Angiosperms
## Date: 04/19
## Author: Maya Schrödl
##################################################
## This script contains the functions that will be
## applied to each species.
##################################################


# Packages ----------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(gridExtra)
library(magrittr)
library(dplyr)
library(minpack.lm)
library(ggrepel)
library(stringr)
library(ggsci)

# Working Directory -------------------------------------------------------
wd="D:/mayasdaten/Documents/IMABEE/Stage_GC/Project_Maya/"

# Species -----------------------------------------------------------------
setwd(paste(wd,"Data",sep=""))
speciesList=read.table("Species_List.txt",h=T,sep="\t")$x # abbreviated species names
speciesInfo=read.table("Species_Info.txt",h=T,sep="\t") # contains taxonomic group, genome size & full species name

# Dataset ---------------------------------------------------------------------
# Loading & correction of CDS without multiple UTR
CDSCorr=function(){
  setwd(paste(wd,"Data/GC_values_NoMultiUTR",sep=""))
  exon=read.table(paste(species,"_exon_NoMultiUTR.txt",sep=""),h=T,sep="\t")
  exon=exon[!exon$nChr=="ChrM" & !exon$nChr=="ChrC",] # no mitochondrial or cytoplasmic DNA
  CDS=exon[exon$type=="CDS",] # only take CDS, not UTR
  return(CDS)
}

# Example Species ---------------------------------------------------------
species="Athaliana" # to test script

#=== SCRIPT ================================================================
# General Summary  ---------------------------------------------------------

# Mean and sd of GC3 (& GCtot) content
MeanSd=function(){ 
  CDS=CDSCorr() # take CDS
  setwd(paste(wd,"Results",sep=""))
  summary=read.table("Summary.txt",h=T)
  index=which(summary$species==species) # index of current species
  summary$GC_total.mean[index]=round(mean(CDS$GC_total),3) # mean GCtot
  summary$GC_total.sd[index]=round(sd(CDS$GC_total),3) # sd Gtot
  summary$GC3.mean[index]=round(mean(CDS$GC3),3) # mean GC3
  summary$GC3.sd[index]=round(sd(CDS$GC3),3) # sd GC3
  write.table(summary, file = "Summary.txt", sep = "\t",quote=F,row.names = F)
}

# DISTRIBUTION ##################################################################
# Distribution of GC content per gene class (gene with a certain exon number)
# weighted by gene length
# Bimodal or unimodal distribtion?

# Calculate weighted mean CDS GC (tot,1,2,3) per gene
GCmean.pergene.weighted=function(){ # GC content weighted by gene length
  CDS=CDSCorr()
  res=CDS %>% # res: weighted mean for each gene
    group_by(keygene) %>% 
    summarize(GC_total= sum(mean(GC_total)*len)/sum(len),GC1 = sum(mean(GC1)*len)/sum(len),GC2 = sum(mean(GC2)*len)/sum(len),GC3=sum(mean(GC3)*len)/sum(len), nbexon=nbexon[1])
   return(res)}

# Distribution of GC (tot,1,2,3) according to mean weighted CDS GCcontent per gene
GCDistribPlot=function(GCtype="GCtot"){
  sequence=GCmean.pergene.weighted()
  if (GCtype=="GCtot"){t=sequence$GC_total}#t: chosen GC content
  if (GCtype=="GC1"){t=sequence$GC1}
  if (GCtype=="GC2"){t=sequence$GC2}
  if (GCtype=="GC3"){t=sequence$GC3}
  
  n=14 #max number of exons per gene
  df=data.frame(cbind(t,sequence$nbexon)) #dataframe with chosen GC content and nbexon for each gene
  df=df[which(df[,2]<=n),] # df[,2]: nbexon - chose only genes with less than n exons
  df$nbexon=as.factor(df[,2])
  
  g=ggplot()+geom_density(data=df,aes(x=df$t,group=df$nbexon,col=df$nbexon))+
    #ggtitle(paste("Distribution of",GCtype,species,sep=" "))+
    xlab("mean GC3 [%]%")+
    # legend:
    theme(legend.title = element_blank())+
    theme(legend.key.size =  unit(0.15, "in"))
  return(g)
}


# GRADIENT ##############################################################

#--- GC gradient as a function of exon number per gene ---------------------

# Gradient from data
# (with median GC content)
Gradient=function(n=14){
  # n: maximum number of exons (per gene)
  data=CDSCorr()
 
  res=data[data$nbexon<=n,] %>% # median GC (tot,1,2,3)+se+ic for each nbexon (i) and num (j)
    group_by(nbexon,num) %>%
    summarize(
      # median GC:
      GC_total.median= median(GC_total, na.rm = TRUE),GC1.median = median(GC1, na.rm = TRUE),GC2.median = median(GC2, na.rm = TRUE),GC3.median = median(GC3, na.rm = TRUE),
      # sd GC:      
      GC_total.sd= sd(GC_total, na.rm = TRUE),GC1.sd = sd(GC1, na.rm = TRUE),GC2.sd = sd(GC2, na.rm = TRUE),GC3.sd = sd(GC3, na.rm = TRUE),lGC=length(GC_total))%>%
      # se GC:
      mutate( GC_total.se=GC_total.sd/sqrt(lGC),GC1.se=GC1.sd/sqrt(lGC),GC2.se=GC2.sd/sqrt(lGC),GC3.se=GC3.sd/sqrt(lGC))  %>%
      # ic GC:
      mutate( GC_total.ic=GC_total.se * qt((1-0.05)/2 + .5, lGC-1),GC1.ic=GC1.se * qt((1-0.05)/2 + .5, lGC-1),GC2.ic=GC2.se * qt((1-0.05)/2 + .5, lGC-1),GC3.ic=GC3.se * qt((1-0.05)/2 + .5, lGC-1))
  
  res$nbexon=as.factor(res$nbexon)
  
  return(res) # res: median GC content for each nbexon (i) and num (j)
}

#--- Gradient Model ----------------------------------------------------------
GradientModel=function(GCtype){ #GCtype: tot,1,2,3
  grad=Gradient() # median GC content for each nbexon (i) and num (j)
  # t: chosen GC content
  # set GC.median to median of chosen GC content
  if (GCtype=="GCtot"){grad$GC.median=grad$GC_total.median}
  if (GCtype=="GC1"){grad$GC.median=grad$GC1.median}
  if (GCtype=="GC2"){grad$GC.median=grad$GC2.median}
  if (GCtype=="GC3"){grad$GC.median=grad$GC3.median}
  
  fn=GC.median~(A-e)*c^(num-1)+(B-e)*d^(as.numeric(nbexon)-num)+e #GC.median: median GC for chosen GC content #Formula (1) in 2.4
  
  A=40; B=40; c=0.5; d=0.5; e=40 # starting values
  mod1=nlsLM(fn,data=grad,start=c(A=A,B=B,c=c,d=d,e=e)) # nls produces non linear model

  return(mod1)
}

# to extract gradient coefficients
GradientModelSummary=function(GCtype="GC3"){
  mod1=GradientModel(GCtype=GCtype)
  
  coef=data.frame(summary(mod1)$coefficients) # parameters A,B,c,d,e
  coef$p=signif(coef$Pr...t..,digits=2) # p values for each parameter
  
  setwd(paste(wd,"Results/Gradient/Model",sep="")) 
  gradSum=read.table(paste("CDSGradientModel_",GCtype,".txt",sep=""),h=T) # to store the parameters for each species
  index=which(gradSum$species==species)
  
  # each parameter:
  gradSum$A[index]=round(coef$Estimate[1],4)
  gradSum$B[index]=round(coef$Estimate[2],4)
  gradSum$c[index]=round(coef$Estimate[3],4)
  gradSum$d[index]=round(coef$Estimate[4],4)
  gradSum$e[index]=round(coef$Estimate[5],4)
  
  # p-values of the model:
  gradSum$p_A[index]=coef$p[1]
  gradSum$p_B[index]=coef$p[2]
  gradSum$p_c[index]=coef$p[3]
  gradSum$p_d[index]=coef$p[4]
  gradSum$p_e[index]=coef$p[5]
  
  # store everything
  write.table(gradSum, file=paste("CDSGradientModel_",GCtype,".txt",sep=""),sep="\t",row.names = F,quote=F)
}

# to apply the model to the nbexon=14
GradientModel_14exons=function(species,GCtype="GC3"){
  setwd(paste(wd,"Results/Gradient/Model",sep=""))
  mod=read.table(paste("CDSGradientModel_",GCtype,".txt",sep=""),h=T)
  index=which(mod$species==species)
  A=mod$A[index];B=mod$B[index];c=mod$c[index];d=mod$d[index];e=mod$e[index]
  fun=function(num){(A-e)*c^(num-1)+(B-e)*d^(14-num)+e}# our gradient model; nbexon=14 #Formula (1) in 2.4
  return(fun)# returns model with parameters of the species
}

# Does the model predict correctly? TESTS ---------------------------------
# Section 2.4.1

# Test summarized genes
# for each num, nbexon:
# mean GC3 weighted by mean exon length
GradientModel_test=function(CDS=CDSCorr()){
  setwd(paste(wd,"Results/Gradient/Model",sep=""))
  mod=read.table(paste("CDSGradientModel_GC3.txt",sep=""),h=T)
  index=which(mod$species==species)
  A=mod$A[index];B=mod$B[index];c=mod$c[index];d=mod$d[index];e=mod$e[index]
  
  fun=function(num,nbexon){(A-e)*c^(num-1)+(B-e)*d^(nbexon-num)+e}
  CDS=CDS[CDS$nbexon<=14,]
  
  MedianLen=CDS%>%
    group_by(nbexon,num)%>%
    summarize(meanlen=mean(len,na.rm=T),GC3.median=median(GC3,na.rm=T)) 
  # mean CDS length for each nbexon and each num
  # median GC3 content for each nbexon and each num (like gradient plot)
  
  GCobs=MedianLen%>% #observed GC per gene class, weighted by (mean) exon length
    group_by(nbexon) %>%
    summarize(GCobs=sum(GC3.median*meanlen)/sum(meanlen)) # weighted mean of GC3 per gene class
  # len: mean CDS length
  
  GCpred=MedianLen%>% #predicted GC per gene class
    group_by(nbexon)%>%
    summarize(GCpred=sum(fun(num,nbexon)*meanlen)/sum(meanlen))# weighted mean of GC3 per gene class
  # len: mean CDS length
  
  GC=GCpred
  GC$GCobs=GCobs$GCobs
  GC$species=rep(species,14)
  
  g=ggplot(data=GC)+geom_point(aes(x=GCobs,y=GCpred,col=nbexon))+ylim(min(GC$GCpred),max(GC$GCpred))+
    geom_abline()+ggtitle(paste(species,"GC3"))
  
  return(list(GC,g))# g: plot
}

# Genes with one exon: overestimated GC for model (for most species)



# bGC? #########################################################################

# Correlation between synonymous (GC3) and non-synonymous (GC1+GC2) positions for each species individually
# along the position on the gene # I did not put this into my report
cor_syn.nonsyn=function(){
CDS=CDSCorr()
CDS=CDS[CDS$nbexon<=14,]
res=CDS%>%
  group_by(nbexon,num)%>%
  summarize(cor=cor.test((GC2+GC1),GC3,method="spearman")$estimate,pvalue=cor.test(GC3,(GC2+GC1),method="spearman")$p.value)

res=res[res$pvalue<0.05,]
res$nbexon=as.factor(res$nbexon)

g=ggplot(data=res)+
  geom_line(aes(x=num,y=cor,col=nbexon,group=nbexon))+
  geom_point(aes(x=num,y=cor,col=nbexon,group=nbexon),shape = 21)+
  # legend:
  theme(legend.title = element_blank())+
  theme(legend.key.size =  unit(0.15, "in"))+  
  theme(legend.text = element_text(size=5))+
  guides(colour = guide_legend(override.aes = list(size=1,linetype=0)))+
  ggtitle(paste(species))+xlab("CDS rank")+ylab("Cor GC3~(GC2+GC1)")
return(g)}


# GENE STRUCTURE ################################################################

# Gene structure distribution
# just to have a quick look on how the genelength is distributed
GeneStructure=function(){
CDS=CDSCorr() # Without MultiUTR
CDS=CDS[CDS$nbexon<15,]
g=ggplot()+geom_histogram(data=CDS,aes(CDS$nbexon))+scale_x_continuous(breaks=seq(1,15,by=1))+
  xlab("number of exons per gene")
return(g)}


# Mean exon number (+var), and number of monoexonic genes
MeanExon_MonoExon=function(){
  CDS=CDSCorr() # Without MultiUTR
  CDS=CDS[CDS$nbexon<14,]
  CDS_keygene=CDS %>%
    group_by(keygene)%>%
    summarize(nbexon=nbexon[1])
  
  nbexon.mean=mean(CDS_keygene$nbexon)
  nbexon.sd=sd(CDS_keygene$nbexon)
  propMonoEx=length(CDS_keygene[CDS_keygene$nbexon==1,])/length(CDS_keygene$nbexon)
  
  setwd(paste(wd,"Results",sep=""))
  SumStr=read.table("Summary_GeneStructure.txt",sep="\t",h=T)
  index=which(SumStr$species==species)
  SumStr$nbexon.mean[index]=nbexon.mean
  SumStr$nbexon.sd[index]=nbexon.sd
  SumStr$prop_MonoEx[index]=propMonoEx
  
  write.table(SumStr,"Summary_GeneStructure.txt",sep="\t",quote=F,row.names=F)}

#MonoExon()

# Correlation Chromosome Lenght & GC ###################################################

#--- Determine Chromosomes ------------------------------------------------
#We want to choose only the biggest sequences (=Chr)
Chromosomes=function(){
  Chr_Scaffolds=exon%>%
    group_by(nChr) %>%
    summarize(lenChr=max(end)-min(start))
  
  Chr_Scaffolds=Chr_Scaffolds[order(Chr_Scaffolds$lenChr,decreasing=T),]
  View(Chr_Scaffolds)
  x=as.numeric(readline("num? "))# Which sequences do we keep as chromosomes?
  print(x)
  
  Chr_Chosen=Chr_Scaffolds[1:x,]
  Chr_Chosen=Chr_Chosen[order(Chr_Chosen$nChr),]
  
  setwd(paste(wd,"Data/Chromosomes",sep=""))
  write.table(Chr_Chosen,paste(species,"_Chromosomes.txt",sep=""),row.names = F,quote=F)}


#--- Calculate GC ------------------------------------------------------------
ChrGC_Summary=function(seq){
  setwd(paste(wd,"Data/Chromosomes",sep=""))
  Chr=read.table(paste(species,"_Chromosomes.txt",sep=""),h=T)
  
  data_Filtered=seq %>% filter(nChr %in% Chr$nChr)
  ChrGC=data_Filtered%>%
    group_by(nChr) %>%
    summarize(GC_total= mean(GC_total, na.rm = TRUE),GC1 = mean(GC1, na.rm = TRUE),GC2 = mean(GC2, na.rm = TRUE),GC3 = mean(GC3, na.rm = TRUE))
  ChrGC$lenChr=Chr$lenChr[match(Chr$nChr,ChrGC$nChr)]
  
  return(ChrGC)
}


#--- Correlations ------------------------------------------------------------

Corr_ChrLength.GC=function(){
  data=ChrGC_Summary(CDSCorr())
  seq2=seq[seq$num==1,] # First CDS (to compare to gradient)
  data2=ChrGC_Summary(seq2)
  
  setwd(paste(wd,"Results/bGC",sep=""))
  ChrSum=read.table(file=paste("Summary_Chromosomes_",sequence,".txt",sep=""),h=T)
  ChrSum$species=as.character(ChrSum$species)
  index=which(ChrSum$species==species)
  
  if (length(data$nChr)>1){
    test=cor.test(data$lenChr,y=data$GC_total,method="spearman")
    ChrSum$cor_GCtotCDS.lenChr[index]=round(as.numeric(test$estimate),3)
    ChrSum$p_GCtotCDS.lenChr[index]=round(as.numeric(test$p.value),5)
    
    test=cor.test(data$lenChr,y=data$GC3,method="spearman")
    ChrSum$cor_GC3CDS.lenChr[index]=round(as.numeric(test$estimate),3)
    ChrSum$p_GC3CDS.lenChr[index]=round(as.numeric(test$p.value),5)
    
    test=cor.test(data2$lenChr,y=data2$GC_total,method="spearman")
    ChrSum$cor_GCtotCDSFirstExon.lenChr[index]=round(as.numeric(test$estimate),3)
    ChrSum$p_GCtotCDSFirstExon.lenChr[index]=round(as.numeric(test$p.value),5)
    
    test=cor.test(data2$lenChr,y=data2$GC3,method="spearman")
    ChrSum$cor_GC3CDSFirstExon.lenChr[index]=round(as.numeric(test$estimate),3)
    ChrSum$p_GC3CDSFirstExon.lenChr[index]=round(as.numeric(test$p.value),5)
  }
  ChrSum$nChr[index]=length(data$nChr)
  write.table(ChrSum,file=paste("Summary_Chromosomes_",sequence,".txt",sep=""),sep="\t",row.names=F,quote=F)
}

#--- Local Structure to determinine centromere (-> Chromosome arms) ----------
CentromTest=function(ChrNmb,div=50){
  Chr=seq[seq$nChr==ChrNmb,]
  numb=max(Chr$start,na.rm=T)/div
  meanGC=data.frame(matrix(rep(NA,div*2),div,2))
  colnames(meanGC)=c("start","GC3")
  
  for (x in 1:div){
    n=x*numb
    meanGC$start[x]=n
    meanGC$GC3[x]=mean(Chr$GC3[Chr$start<n],na.rm=T)
    Chr=Chr[!Chr$start<n,]
  }
  
  ggplot(meanGC, aes(x=meanGC$start,y=meanGC$GC3))+geom_line()+
    ggtitle(paste(species,ChrNmb,"GC3"))+xlab("Position on Chromosome")+ylab("GC3 content")
}
# CentromTest(as.character(Chr$nChr[1]),100)
