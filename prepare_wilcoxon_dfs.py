import pandas as pd
import numpy as np
import csv
from functools import reduce


genes1 = []
genes2 = []
with open('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/ed_combined_reps.txt') as file:
 	for line in file.readlines():
 		genes1.append(line.split(","))

with open('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/n2_combined_reps.txt') as file:
 	for line in file.readlines():
 		genes2.append(line.split(","))
		

genes_in_1 = set()
genes_in_2 = set()
to_remove_from_both = set()

#Create a set of the genes that are in input_file #1
for gene_log_ratio in genes1:
	gene = gene_log_ratio[0]
	#save genes with just 1 insert for removal from both output sets
	if (len(gene_log_ratio) <= 2):
		to_remove_from_both.add(gene)
	else:
		genes_in_1.add(gene)

#Create a set of the genes that are in input_file #2
for gene_log_ratio in genes2:
	gene = gene_log_ratio[0]
	#save genes with just 1 insert for removal from both output sets
	if (len(gene_log_ratio) <= 2):
		to_remove_from_both.add(gene)
	else:	
		genes_in_2.add(gene)

#remove the genes that only have 1 insertion in either input-file
genes_in_1 -= to_remove_from_both
genes_in_2 -= to_remove_from_both

#Get total Count of insertions for each file
total_insertions_in_1 = 0
total_insertions_in_2 = 0
for gene_log_ratios in genes1:
	total_insertions_in_1 += len(gene_log_ratios) - 1
for gene_log_ratios in genes2:
	total_insertions_in_2 += len(gene_log_ratios) - 1

#Remove the genes that don't exist in the other input-set
genes1 = list(map(lambda line: ",".join(line), list(filter(lambda row : row[0] in genes_in_2, genes1))))
genes2 = list(map(lambda line: ",".join(line), list(filter(lambda row : row[0] in genes_in_1, genes2))))


output_file1 = "ed_combined_reps_n2_overlap.txt"
output_file2 = "n2_combined_reps_ed_overlap.txt"

with open('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/' + output_file1, 'w') as out_file:
	out_file.writelines(genes1)

with open('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/' + output_file2, 'w') as out_file:
	out_file.writelines(genes2)

print("Total insertions in " + output_file1 + ":" + str(total_insertions_in_1))
print("Total insertions in " + output_file2 + ":" + str(total_insertions_in_2))



