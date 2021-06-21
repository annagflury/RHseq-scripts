library(matrixTests)
library(dplyr)

#import N2 combined effect size file (83)
x=read.table("/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/ed_combined_reps_n2_overlap.txt", header=FALSE, sep=",", col.names=paste0("V",seq_len(132)), fill = TRUE, stringsAsFactors = FALSE)
rownames(x) <-x[,1]
x[,1] <- NULL
is.data.frame(x) && all(sapply(x, is.numeric))


#import ED combined effect size file (24)
y=read.table("/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/n2_combined_reps_ed_overlap.txt", header=FALSE, sep=",", col.names=paste0("V",seq_len(345)), fill = TRUE, stringsAsFactors = FALSE)
rownames(y) <-y[,1]
y[,1] <- NULL
is.data.frame(y) && all(sapply(y, is.numeric))

#find genes with biggest difference between N2 and Ed
#means=data.frame(n2_ave=rowMeans(x, na.rm=TRUE, dims=1), ed_ave=rowMeans(y, na.rm=TRUE, dims=1))
#means$difference <- c(means$n2_ave - means$ed_ave)
#means=means[order(means$difference, decreasing = TRUE),]
#means

#perform rowwise wilcoxon
#my_data <- row_wilcoxon_twosample(x, y, alternative="two.sided",mu=0,exact=NA,correct=TRUE)
#my_data <- tibble::rownames_to_column(my_data, "gene")
#my_data <- data.frame(my_data$gene,my_data$pvalue)
#names(my_data)[1] <- "gene"
#names(my_data)[2] <- "pvalue"

my_data <- row_t_welch(x, y)
my_data <- tibble::rownames_to_column(my_data, "gene")
my_data <- data.frame(my_data$gene,my_data$pvalue)
names(my_data)[1] <- "gene"
names(my_data)[2] <- "pvalue"

#multiple testing correction
adjusted = c(p.adjust(my_data$pvalue,method="BH"))
my_data <- cbind(my_data, adjusted)
my_data <- my_data[order(my_data$pvalue),]
my_data

#outfile
write.csv(my_data, file="/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/pvalues_ttest_cutoff_5.txt", row.names=FALSE, quote=FALSE)



