library(DESeq2)
library(tidyverse)

##Lets load the data
infile_feature_counts<-"Rdata/vehicle_drug_feature_counts.txt"
df_counts<-read.delim(infile_feature_counts,comment.char = "#")

##Lets remove the columns we don't need
df_counts<-df_counts %>% select(-Chr, -Start,-End,-Strand,-Length)

##To makes samples easier to read remove file extension from sample name
colnames(df_counts)<-gsub(pattern=".bam",replacement = "",
                                  x=colnames(df_counts),fixed=T)

##Lets make a meta-data table (Required for DESeq2)
df_meta<-data.frame(Sample=colnames(df_counts)[2:5],
                    Info=colnames(df_counts)[2:5])

df_meta<-df_meta %>% separate(col=Info,into=c("Condition","Replicate"),sep="_")

##Put the count data into a matrix format so it can be used
##by DESeq2

m_counts<-as.matrix(df_counts[,2:ncol(df_counts)])
rownames(m_counts)<-df_counts$Geneid

##Make sure the colnames of the count data matches the sample
##names in the meta-data. Important because matching between
##count data and meta-data is position based
all(colnames(m_counts)==df_meta$Sample)

##Create the DESeq2 data object

dds<-DESeqDataSetFromMatrix(countData = m_counts,
                            colData = df_meta,
                            design=~Condition)

##Perform sample ordination to assess sample quality
##PCA's generally don't work will with RNA-seq data
##The rlog function helps normalize the data 
##in a PCA friendly manner
rld<-rlog(dds)

DESeq2::plotPCA(rld,intgroup="Condition")

##Perform DESeq2 analysis to detect differences based on Condition
dds<-DESeq(dds)

##Extract the results
##Use contrast to specify the comparison to show and what the 
##reference is. drug is in the numerator and vehicle the
##denominator
##No filtering performed. Specify alpha for 
##optimizing independent filtering
res<-results(dds,contrast = c("Condition","drug","vehicle"),
             alpha=0.05,pAdjustMethod = "fdr")

df_res<-as.data.frame(res)
df_res$Transcript_ID<-rownames(df_res)
col_order<-c("Transcript_ID",colnames(df_res)[1:6])
df_res<-df_res[,col_order]
View(df_res)

##Split into Significantly up or down-regulated
df_res_sig_up<-subset(df_res,log2FoldChange>0 & padj<=0.1)
df_res_sig_down<-subset(df_res,log2FoldChange<0 & padj<=0.1)

##Lets also get the normalized data for graphing 
df_counts<-as.data.frame(counts(dds,normalized=T))
df_counts$Transcript_ID<-rownames(df_counts)
col_order<-c("Transcript_ID",colnames(df_counts[1:4]))
df_counts<-df_counts[,col_order]
##Write out all of the results
write_tsv(df_res,"Routput/all_deseq2_results.tsv")
write_tsv(df_res_sig_down,"Routput/sig_down_deseq2_results.tsv")
write_tsv(df_res_sig_up,"Routput/sig_up_deseq2_results.tsv")
write_tsv(df_counts,"Routput/deseq_normalized_counts.tsv")

