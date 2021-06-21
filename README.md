# RHseq-scripts

RH-seq Pipeline

1. Cut adapter out of end of read: 3’ adapter sequence: AGGCGCGCCACAGTACGAGT
	⁃	fastx_clipper -a AGGCGCGCCACAGTACGAGT -i input.fq -o output_no_adapter.fastq
2. Trim first 6 bp of file, so every read begins with barcode
	⁃	fastx_trimmer -f 7 -i input_no_adapter.fastq -o no_adapter_trimmed.fastq
3. Separate barcodes (indexes)
	⁃	cat /Users/anna/no_adapter_trimmed.fastq | fastx_barcode_splitter.pl --bcfile /Users/anna/Documents/buck/sequences/barcodes.txt --bol --mismatches 0 --prefix /Users/anna/ --suffix ".fq"
4. Map reads with map_and_pool_elegans.py
	⁃	conda activate py2
	⁃	python map_and_pool_elegans.py /directory/BC.fq /output_directory/
5. Get normalized log ratios with log_ratios.py
6. Combine replicates with 5_combine_replicates.py
7. Remove quotes with:
	⁃	sed ’s/“//g’ input.txt > output.txt
8. Remove n2/ed3077 before WBGene with:
	- sed’s/n2//g’ input.txt > output.txt
9. Filter ED/N2 genes in common with prepare_wilcoxon_dfs.py
10. Get longest insertion count to input into R script for each file (+1):
	awk -F '[,]' '{print $1, NF-1}' /Users/anna/Documents/buck/rhseq/ed_combined_reps_n2_overlap.txt | sort -k2 -n
11. In wilcox_matrix.R:
	replace col.names=paste0("V",seq_len(71)) with longest insertion count +1 from step 10 for ed and n2 respectively
	run wilcox matrix
