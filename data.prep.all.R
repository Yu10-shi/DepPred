library(dplyr)
library(GSA)

setwd("/Users/shiyu/Desktop/UofT/Precticum-projects/data/")

# TCGA data
tcga.exp=read.table(gzfile("TCGA.BRCA.sampleMap_HiSeqV2.gz"), header = T)

# original Broad CCL dependency score
# we use depmap 22Q2 version
score=read.csv("Broad_Achilles_gene_effect.csv") # 17,386*1,086

# broad expression data
exp=read.csv("CCLE_expression.csv") # 1,046*19,222
sample_info=read.csv("sample_info.csv")
# keep CCL with expression data
broad.score=score[score$DepMap_ID%in%exp$X,]
broad.exp=exp[exp$X%in%broad.score$DepMap_ID,]

# rename score snp
snp_name=sapply(2:dim(broad.score)[2], function(i){gsub('[..].*',"", names(broad.score)[i])})
colnames(broad.score)[2:dim(broad.score)[2]]=snp_name

# rename rna seq
exp_name=sapply(2:dim(broad.exp)[2], function(i){gsub('[..].*',"", names(broad.exp)[i])})
colnames(broad.exp)=c("DepMap_ID",exp_name)

# merge with cell line name
sample.info=sample_info[,c(1,3)]

# column is CCL and row is gene feature
broad.final_exp = broad.exp[order(broad.exp$DepMap_ID),]
broad.final_exp=t(broad.final_exp)
broad.final_exp=cbind(rownames(broad.final_exp), broad.final_exp)
broad.final_exp[1,1]="Gene"
broad.final_exp2=data.frame(broad.final_exp)
names(broad.final_exp2)=broad.final_exp2[1,]
broad.final_exp2=broad.final_exp2[-1,]
# Broad expression
broad.final_exp2=broad.final_exp2[order(broad.final_exp2$Gene),]

# calculate sd in original dep score across 1005 CCLs
score_sd=data.frame(apply(broad.score[,c(2:dim(broad.score)[2])], 2, var))

# COMIC cancer gene census
comic=read.csv("COMIC_Census_allSat.csv")

# check sd>0.2
names(score_sd)="sd"
score_sd$Isd=ifelse(score_sd$sd>0.2, 1, 0)
# check gene in comic cancer gene census
score_sd$Icomic=ifelse(rownames(score_sd)%in%comic$Gene.Symbol,1,0)
score_sd$ind=ifelse((score_sd$Isd+score_sd$Icomic)>=1,1,0)

final_broad_gene=rownames(score_sd)[which(score_sd$ind==1)]

# Sanger data
#setwd("/Users/shiyu/Desktop/UofT/Precticum-projects/data/Sanger")
#sanger_score=read.table("01a_qnorm_corrected_logFCs.tsv", header = T)
#sanger_score=read.table("/Users/shiyu/Desktop/UofT/Precticum-projects/data/Sanger/Project_score_combined_Sanger_v1_Broad_20Q2_20210311/Project_score_combined_Sanger_v1_Broad_20Q2_fitness_scores_fold_change_values_20210311.tsv", fill=TRUE)
sanger_score=read.csv("sanger_gene_effect.csv")
#sanger_score=sanger_score[-c(4,5),]

# switch columns and rows of scores
#sanger_score1=as.data.frame(t(sanger_score))

# extract broad and sanger data
#names(sanger_score1)=sanger_score1[1,]
#sanger_score1=sanger_score1[-1,]

#broad.scores=sanger_score1[sanger_score1$source=="Broad",]

#sanger.scores=sanger_score1[sanger_score1$source=="Sanger",]

# sanger expression data

#sanger_exp=read.csv("/Users/shiyu/Desktop/UofT/Precticum-projects/data/Sanger/rnaseq_tpm_20220624.csv", header = T)
#sanger_exp=sanger_exp[,-1]


# take overlap of CCLs between sanger exp and sanger score
#sanger.scores=sanger.scores[sanger.scores$model_id%in%names(sanger_exp),]
#Sanger expression
#sanger_exp=sanger_exp[,names(sanger_exp)%in%c("gene", sanger.scores$model_id)]
# replace NA with 0 Replace NA with mean value no missing value

# keep CCL with expression data
sanger.score=sanger_score[sanger_score$X%in%exp$X,]
sanger.exp=exp[exp$X%in%sanger.score$X,]

# rename score snp
sanger_snp_name=sapply(2:dim(sanger.score)[2], function(i){gsub('[..].*',"", names(sanger.score)[i])})
colnames(sanger.score)[2:dim(sanger.score)[2]]=sanger_snp_name
names(sanger.score)[1]="DepMap_ID"

# rename rna seq
sanger_exp_name=sapply(2:dim(sanger.exp)[2], function(i){gsub('[..].*',"", names(sanger.exp)[i])})
colnames(sanger.exp)=c("DepMap_ID",sanger_exp_name)

# column is CCL and row is gene feature
sanger.final_exp = sanger.exp[order(sanger.exp$DepMap_ID),]
sanger.final_exp=t(sanger.final_exp)
sanger.final_exp=cbind(rownames(sanger.final_exp), sanger.final_exp)
sanger.final_exp[1,1]="Gene"
sanger.final_exp2=data.frame(sanger.final_exp)
names(sanger.final_exp2)=sanger.final_exp2[1,]
sanger.final_exp2=sanger.final_exp2[-1,]
# Sanger expression
sanger.final_exp2=sanger.final_exp2[order(sanger.final_exp2$Gene),]

# calculate sd in dep score across sanger data
sanger_score_sd=data.frame(apply(sanger.score[,c(2:dim(sanger.score)[2])], 2, var))

# check sd>0.2
names(sanger_score_sd)="sd"
sanger_score_sd$Isd=ifelse(sanger_score_sd$sd>0.2, 1, 0)
# check gene in comic cancer gene census
sanger_score_sd$Icomic=ifelse(rownames(sanger_score_sd)%in%comic$Gene.Symbol,1,0)
sanger_score_sd$ind=ifelse((sanger_score_sd$Isd+sanger_score_sd$Icomic)>=1,1,0)

final_sanger_gene=rownames(sanger_score_sd)[which(sanger_score_sd$ind==1)]

# common sanger and broad data
common.sanger.broad = sanger.score[sanger.score$DepMap_ID%in%broad.score$DepMap_ID,]
common.broad.sanger = broad.score[broad.score$DepMap_ID%in%sanger.score$DepMap_ID,]

# unique broad data
unique.broad.score.to.sanger = broad.score[!broad.score$DepMap_ID%in%common.broad.sanger$DepMap_ID,]
unique.broad.exp.to.sanger = broad.final_exp2[,names(broad.final_exp2)%in%c("Gene", unique.broad.score.to.sanger$DepMap_ID)]

# unique sanger data
unique.sanger.score.to.broad = sanger.score[!sanger.score$DepMap_ID%in%common.broad.sanger$DepMap_ID,]
unique.sanger.exp.to.broad = sanger.final_exp2[,names(sanger.final_exp2)%in%c("Gene", unique.sanger.score.to.broad$DepMap_ID)]

# RNAi data
rnai.lfc = read.csv("Marcotte_LFC_matrix.csv")
rnai_exp = read.csv("RNAseq_lRPKM_data.csv")
rnai_exp$X = sub(" .*", "", rnai_exp$X)
#names(rnai.exp)[1]='Barcode.Sequence'
#gene.mapping = read.csv("shRNAmapping.csv")

#rnai.score=read.csv("RNAi_(Achilles+DRIVE+Marcotte,_DEMETER2).csv")
rnai_score = read.csv("D2_Achilles_gene_dep_scores.csv")
rnai_score$X = sub(" .*", "", rnai_score$X)

# keep CCL with expression data
rnai_score=rnai_score[,names(rnai_score)%in%names(rnai_exp)]
rnai_exp=rnai_exp[,names(rnai_exp)%in%names(rnai_score)]

# assume different gene expresses similarly in the same cell line
# replace NA with column mean
for (i in 2:dim(rnai_score)[2]) {
  rnai_score[,i][is.na(rnai_score[,i])]=mean(rnai_score[,i], na.rm=TRUE)
}


# extract RNAi data with expression information
#rnai.score=rnai.score[rnai.score$DepMap_ID%in%rnai_exp$X,]
#rnai.exp=rnai_exp[rnai_exp$X%in%rnai.score$DepMap_ID,]

# rename RNAi exp data
#rnai.exp_name=sapply(2:dim(rnai.exp)[2], function(i){gsub('[..].*',"", names(rnai.exp)[i])})
#colnames(rnai.exp)=c("DepMap_ID",rnai.exp_name)

# rename score snp
#rnai.snp_name=sapply(2:dim(rnai.score)[2], function(i){gsub('[..].*',"", names(rnai.score)[i])})
#colnames(rnai.score)[2:dim(rnai.score)[2]]=rnai.snp_name

# extract unique RNAi CCL
#unique.rnai.score=rnai.score[!rnai.score$DepMap_ID%in%score$DepMap_ID,]
#unique.rnai.score=unique.rnai.score[order(unique.rnai.score$DepMap_ID),]
#unique.rnai.exp=rnai.exp[!rnai.exp$DepMap_ID%in%exp$DepMap_ID,]

# order rnai exp by dep_id
rnai_exp=rnai_exp[order(rnai_exp$X),]

# delete NA rows
rnai.exp=rnai_exp[!rnai_exp$X=="",]

# calculate sd in original dep score across 487 CCLs
rnai_score_sd=data.frame(apply(rnai_score[,c(2:dim(rnai_score)[2])], 1, var))

# check sd>0.2
names(rnai_score_sd)="sd"
rownames(rnai_score_sd) = rnai_score$X
rnai_score_sd$Isd=ifelse(rnai_score_sd$sd>0.2, 1, 0)
# check gene in comic cancer gene census
rnai_score_sd$Icomic=ifelse(rownames(rnai_score_sd)%in%comic$Gene.Symbol,1,0)
rnai_score_sd$ind=ifelse((rnai_score_sd$Isd+rnai_score_sd$Icomic)>=1,1,0)

# final gene for RNAi dataset
final_rnai_gene=rownames(rnai_score_sd)[which(rnai_score_sd$ind==1)]

# column is CCL and row is gene feature
#rnai.final_exp=t(unique.rnai.exp)
#rnai.final_exp=cbind(rownames(rnai.final_exp), rnai.final_exp)
#rnai.final_exp[1,1]="Gene"
#rnai.final_exp2=data.frame(rnai.final_exp)
#names(rnai.final_exp2)=rnai.final_exp2[1,]
#rnai.final_exp2=rnai.final_exp2[-1,]

# RNAi expression
#rnai.final_exp2=rnai.final_exp2[order(rnai.final_exp2$Gene),]
names(rnai.exp)[1] = "Gene"

# take overlap of final gene
final_gene=Reduce(intersect, list(final_broad_gene, final_rnai_gene, final_sanger_gene))

# breast cancer patients
brca.code = read.csv("Breast cancer code.csv")
names(tcga.exp)=gsub(".01", "", names(tcga.exp))
names(tcga.exp)=gsub(".11", "", names(tcga.exp))
names(tcga.exp)=gsub(".06", "", names(tcga.exp))
names(tcga.exp)=gsub("\\.", "-", names(tcga.exp))
tcga.exp=tcga.exp[, names(tcga.exp)%in%c("sample", brca.code$bcr_patient_barcode)]
names(tcga.exp)[1]="Gene"
tcga.exp1=tcga.exp[order(tcga.exp$Gene),]

# take overlap of gene features among broad, sanger, rnai and tcga
ovl.fea = Reduce(intersect, list(unique.broad.exp.to.sanger$Gene, unique.sanger.exp.to.broad$Gene, rnai.exp$Gene, tcga.exp1$Gene))

# Check if final gene in expression features and define final depoi
final_depoi = final_gene[final_gene%in%ovl.fea]
write.table(final_depoi,"/Users/shiyu/Desktop/UofT/DepMAP/SL Data/gene_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# extract exp data of broad, sanger, rnai and tcga
final.broad.exp = unique.broad.exp.to.sanger[unique.broad.exp.to.sanger$Gene%in%ovl.fea,]
final.broad.exp = final.broad.exp[!duplicated(final.broad.exp$Gene),]
write.table(final.broad.exp, "broad_exp0424.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t") # 18,926*1,006

final.sanger.exp = unique.sanger.exp.to.broad[unique.sanger.exp.to.broad$Gene%in%ovl.fea,]
final.sanger.exp = final.sanger.exp[!duplicated(final.sanger.exp$Gene),]
#final.sanger.exp[is.na(final.sanger.exp)]=0
write.table(final.sanger.exp, "sanger_exp0424.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

final.rnai.exp = rnai.exp[rnai.exp$Gene%in%ovl.fea,]
final.rnai.exp = final.rnai.exp[!duplicated(final.rnai.exp$Gene),]
write.table(final.rnai.exp, "rnai_exp0424.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

tcga.exp.final = tcga.exp1[tcga.exp1$Gene%in%ovl.fea,]
write.table(tcga.exp.final, "tcga_exp0424.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# extract depscore data of broad, sanger, rnai and tcga
# final broad score file
final.broad.score=unique.broad.score.to.sanger[,(names(unique.broad.score.to.sanger)%in%c("DepMap_ID", final_depoi))] #1005*136
# keep matrix format for score file
rownames(final.broad.score)=final.broad.score[,1]
final.broad.score=final.broad.score[,-1]
final.broad.score=as.data.frame(t(final.broad.score))
# order score by gene
final.broad.score=final.broad.score[order(rownames(final.broad.score)),]
write.table(final.broad.score,"broad_depscore0424.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# final Sanger score file
#sanger.final_score = sanger.scores %>% select(-c(model_name, source))
#names(sanger.final_score)[1] = "DepMap_ID"
#sanger.final_score=sanger.final_score[,(names(sanger.final_score)%in%c("DepMap_ID", final_depoi))] #316*136
# keep matrix format for score file
#rownames(sanger.final_score)=sanger.final_score[,1]
#sanger.final_score=sanger.final_score[,-1]
#sanger.final_score=sanger.final_score[order(rownames(sanger.final_score)),]
#sanger.final_score=as.data.frame(t(sanger.final_score))
# order rnai score by gene
#sanger.final_score=sanger.final_score[order(rownames(sanger.final_score)),]
#write.table(sanger.final_score,"sanger_depscore0410.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

final.sanger.score=unique.sanger.score.to.broad[,(names(unique.sanger.score.to.broad)%in%c("DepMap_ID", final_depoi))] #1005*136
# keep matrix format for score file
rownames(final.sanger.score)=final.sanger.score[,1]
final.sanger.score=final.sanger.score[,-1]
final.sanger.score=as.data.frame(t(final.sanger.score))
# order score by gene
final.sanger.score=final.sanger.score[order(rownames(final.sanger.score)),]
write.table(final.sanger.score,"sanger_depscore0424.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


# final RNAi score file
#rnai.final_score=unique.rnai.score[,(names(unique.rnai.score)%in%c("DepMap_ID", final_depoi))] #128*136
# keep matrix format for score file
#rownames(rnai.final_score)=rnai.final_score[,1]
#rnai.final_score=rnai.final_score[,-1]
#rnai.final_score=as.data.frame(t(rnai.final_score))
# order rnai score by gene
#rnai.final_score=rnai.final_score[order(rownames(rnai.final_score)),]

names(rnai_score)[1]='DepMap_ID'
final.rnai.score=rnai_score[rnai_score$DepMap_ID%in%final_depoi,] #1005*136
final.rnai.score=final.rnai.score[order(final.rnai.score$DepMap_ID),]
write.table(final.rnai.score,"rnai_depscore0424.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

gsea=GSA.read.gmt("c2.cgp.v6.2.symbols.gmt")
#tumor_fprint=read.delim("crispr_gene_fingerprint_cgp.txt", header = T)
gsea_name=gsea$geneset.names
#gsea_name=gsea$geneset.names[gsea$geneset.names%in%tumor_fprint$GeneSet]
#gsea_num=which(gsea$geneset.names%in%tumor_fprint$GeneSet)
fprint_data=matrix(0, 3433, length(final_depoi))
for(i in 1:length(final_depoi)){
  for(j in 1:3433){
    if(final_depoi[i]%in%gsea$genesets[[j]]){
      fprint_data[j,i]=1
    }else{
      fprint_data[j,i]=0
    }
  }
}
fprint_data=data.frame(gsea_name,fprint_data)
colnames(fprint_data)=c("GeneSet",final_depoi)
fprint_data=fprint_data[rowSums(fprint_data[,-1])>0,] #2930*602
write.table(fprint_data,"all_ccl_fprint0424.txt", row.names = FALSE, quote = FALSE, sep = "\t")
#final_gene[!final_gene%in%colnames(tumor_fprint)]
#tumor_fprint=data.frame(tumor_fprint$GeneSet,tumor_fprint[,colnames(tumor_fprint)%in%final_gene])

# C6 oncogenic 
gsea=GSA.read.gmt("c6.all.v2023.1.Hs.symbols.gmt")
gsea_name=gsea$geneset.names
# fingerprint data
fprint=matrix(0, 189, length(final_depoi))
for(i in 1:length(final_depoi)){
  for(j in 1:189){
    #num=gsea_num[j]
    if(final_depoi[i]%in%gsea$genesets[[j]]){
      fprint[j,i]=1
    }else{
      fprint[j,i]=0
    }
  }
}
fprint=data.frame(gsea_name,fprint)
colnames(fprint)=c("GeneSet",final_depoi)
fprint=fprint[rowSums(fprint[,-1])>0,] #185*602
write.table(fprint,"all_ccl_c6_fprint0424.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# C7 Immune
gsea=GSA.read.gmt("c7.immunesigdb.v2023.1.Hs.symbols.gmt")
gsea_name=gsea$geneset.names
# fingerprint data
fprint=matrix(0, 4872, length(final_depoi))
for(i in 1:length(final_depoi)){
  for(j in 1:4872){
    #num=gsea_num[j]
    if(final_depoi[i]%in%gsea$genesets[[j]]){
      fprint[j,i]=1
    }else{
      fprint[j,i]=0
    }
  }
}
fprint=data.frame(gsea_name,fprint) #4872*136
colnames(fprint)=c("GeneSet",final_depoi)
fprint=fprint[rowSums(fprint[,-1])>0,] # 4,866*602
write.table(fprint,"all_ccl_c7_fprint0424.txt", row.names = FALSE, quote = FALSE, sep = "\t")


# distribution plot
library(ggplot2)

broad.exp.plot = as.data.frame(as.numeric(as.vector(as.matrix(final.broad.exp[,-1]))))
broad.exp.plot$type = "Broad"
names(broad.exp.plot)[1] = "value"

sanger.exp.plot = as.data.frame(as.numeric(as.vector(as.matrix(final.sanger.exp[,-1]))))
sanger.exp.plot$type = "Sanger"
names(sanger.exp.plot)[1] = "value"

rnai.exp.plot = as.data.frame(as.numeric(as.vector(as.matrix(final.rnai.exp[,-1]))))
rnai.exp.plot$type = "RNAi"
names(rnai.exp.plot)[1] = "value"

tcga.exp.plot = as.data.frame(as.numeric(as.vector(as.matrix(tcga.exp.final[,-1]))))
tcga.exp.plot$type = "TCGA"
names(tcga.exp.plot)[1] = "value"

exp.plot = rbind(broad.exp.plot, sanger.exp.plot, rnai.exp.plot, tcga.exp.plot)

p = ggplot(exp.plot, aes(x=value))+
  geom_density(aes(colour=factor(type))) +
  labs(colour='Data Type', x = "Gene expression values", y="Density")+
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 12), legend.text = element_text(size = 10))
p

# Our algorithm result
#broad.true = read.csv("/Users/shiyu/Desktop/UofT/Precticum-projects/data/broad_true_20240423.csv", header = F)
broad.true = read.csv("/Users/shiyu/Desktop/UofT/DepMap/broad_true_20240608.csv", header = F)
broad.true.y = c(t(broad.true))
#broad.predict = read.csv("/Users/shiyu/Desktop/UofT/Precticum-projects/data/broad_predict_20240423.csv", header = F)
broad.predict = read.csv("/Users/shiyu/Desktop/UofT/DepMap/broad_predict_20240608.csv", header = F)
broad.predict.y = c(t(broad.predict))
broad = data.frame(broad.true.y, broad.predict.y)
names(broad) = c("Y", "X")

plot(broad$X, broad$Y, main = "Broad CCL data testing", xlab = "DepPred-Predicted dependency score", ylab = "Original dependency score", pch = 15, col="blue",
     cex.lab=1.5, cex.axis=1.5, ylim = c(-3, 1), xlim = c(-3, 1))
abline(lm(Y ~ X, data = broad), col = "black", lwd=4)
text(-1.75, 0.5, expression(rho == 0.74), col = "black", font = 4, cex = 1.5)
text(0.5, -0.5, expression(y==1.11*x - 0.04), col = "red", font = 4, cex = 1.5)

# dependency score density plot
broad.true.y.density = as.data.frame(broad.true.y)
broad.true.y.density$type = "CCL original"
names(broad.true.y.density)[1] = "Dependency scores"

broad.predict.y.density = as.data.frame(broad.predict.y)
broad.predict.y.density$type = "CCL predicted"
names(broad.predict.y.density)[1] = "Dependency scores"

#predicted TCGA
tumor_matrix = read.csv("/Users/shiyu/Desktop/UofT/DepMAP/tcga_broad_predict_20240608.csv",header=FALSE)
tcga.predict.y = c(t(tumor_matrix))
tcga.predict.y.density = as.data.frame(tcga.predict.y)
tcga.predict.y.density$type = "Tumor predicted"
names(tcga.predict.y.density)[1] = "Dependency scores"

all.y.density = rbind(broad.true.y.density, broad.predict.y.density, tcga.predict.y.density)

p = ggplot(all.y.density, aes(x=`Dependency scores`))+
  geom_density(aes(colour=factor(type))) +
  labs(colour=' ', x = "Dependency scores", y="Density")+
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 12), legend.text = element_text(size = 10))
p

# DeepDEP result
deepdep.broad.true.y = read.csv("/Users/shiyu/Desktop/UofT/Precticum-projects/data/deepdep_broad_true_20240502.csv", header = F)
#deepdep.broad.true.y = c(t(broad.true))
deepdep.broad.predict.y = read.csv("/Users/shiyu/Desktop/UofT/Precticum-projects/data/deepdep_broad_step3_20240502.csv", header = F)
#deepdep.broad.predict.y = c(t(broad.predict))
deepdep.broad = data.frame(deepdep.broad.true.y, deepdep.broad.predict.y)
names(deepdep.broad) = c("Y", "X")

plot(deepdep.broad$X, deepdep.broad$Y, main = "Broad CCL data testing", xlab = "Predicted dependency score", ylab = "Original dependency score", pch = 15, col="pink",
     cex.lab=1.5, cex.axis=1.5, ylim = c(-3, 2), xlim = c(-3, 1))
abline(lm(Y ~ X, data = deepdep.broad), col = "lightblue", lwd=4)
text(-1.75, 1, expression(rho == 0.7246), col = "black", font = 4, cex = 1.5)


