##Workshop 1
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("recount3")
install.packages("tidyverse")
library("recount3")
library("edgeR")
library("tidyverse")
available_projects()
indata<-read.delim("Session1_GSE81518_NormalizedCounts.txt", sep="\t")
head(indata)
dim(indata)
sample_names<-colnames(indata[-1])
split_sample_names<-str_split_fixed(sample_names, "_",2)
tissue<-split_sample_names[,2]
patient<-str_sub(sample_names,9,9)
sample.sums<-colSums(indata[ ,-1])
counts_per_gene_per_category<-data.frame("Gene"=indata[,1], "BM"=rowSums(indata[,c(2,4,6)]), "CNS"=rowSums(indata[,c(3,5,7)]))
filter(counts_per_gene_per_category, BM==0 | CNS==0) |> nrow()
boxplot(indata[,-1])
highly_expressed_genes<- which(indata[,-1] > 3000000, arr.ind = TRUE)
indata[highly_expressed_genes[,1],1]
long_format_indata <- pivot_longer(indata,-Gene, names_to = "Sample", values_to="counts")|> 
  mutate("Patient"= str_sub(Sample,9,9),
         "Tissue"= str_split_fixed(Sample, "_",2)[,2])
install.packages("ggplot2")
library(ggplot2)
ggplot(data = long_format_indata) + geom_boxplot(aes(x=Sample, y=counts))
long_format_indata <- long_format_indata |> 
  mutate("Log10_Counts"=log10(counts))
ggplot(data = long_format_indata) + geom_boxplot(aes(x=Sample, y=Log10_Counts))
VEGFA_long_format_indata <- filter(long_format_indata, Gene =="VEGFA")
ggplot(data = VEGFA_long_format_indata) + geom_point(aes(x=Tissue, y=counts, col=Patient))

##RNAseq pre-processing