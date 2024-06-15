#/data/rep2_bismark_yao/methylation_extractor/CHG_context_3days107-LFK10433_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt
#/data/rep2_bismark_yao/methylation_extractor/CHG_context_3daysWT-LFK10433_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt
#/data/rep2_bismark_yao/methylation_extractor/CHG_context_7days107-LFK10433_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt
#/data/rep2_bismark_yao/methylation_extractor/CHG_context_7daysWT-LFK10433_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt
#/data/rep2_bismark_yao/methylation_extractor/CHG_context_S0day107-LFK10433_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt
#/data/rep2_bismark_yao/methylation_extractor/CHG_context_S0dayWT-LFK10433_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt
# looks like:
#Bismark methylation extractor version v0.24.1
#A00599:494:HHVH2DSX7:1:1101:10411:1016_1:N:0:ACTAGTGG+TGGCTAAT  -       1       9179512 x
#A00599:494:HHVH2DSX7:1:1101:10411:1016_1:N:0:ACTAGTGG+TGGCTAAT  -       1       9179503 x
#A00599:494:HHVH2DSX7:1:1101:10411:1016_1:N:0:ACTAGTGG+TGGCTAAT  -       1       9179438 x
#A00599:494:HHVH2DSX7:1:1101:10411:1016_1:N:0:ACTAGTGG+TGGCTAAT  -       1       9179431 x
#A00599:494:HHVH2DSX7:1:1101:10411:1016_1:N:0:ACTAGTGG+TGGCTAAT  -       1       9179422 x

# Pair seq result to genome name based on coordinates
## Path to Araport11 protein_coding_position
#/home/vcm/araport_reference/protein_coding.Araport11.position
###
#1# STEP 1: ### pairing genome
###

zcat Araport11.gtf.gz |perl -alne '{next unless $F[2] eq "gene" ;/gene_id \"(.*?)\";/; print "$F[0]\t$F[3]\t$F[4]\t$1" }' >protein_coding.Araport11.position

### up500
## prepare promoter position by minus 500bp of the START
# output path: /data/rep2_bismark_yao/methylation_extractor/promoter/up500
# ref_genome path:  /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_up500.Araport11.position
zcat /home/vcm/araport_reference/Araport11.gtf.gz |perl -alne '{next unless $F[2] eq "gene" ;/gene_id \"(.*?)\";/; print "$F[0]\t$F[3]\t$F[3]\t$1" }' | awk '{ $2 = $3 - 500; print }' | awk '$2 >= 0' | sed 's/ /\t/g' >promoter_up500.Araport11.position
# Notice that there're 3 lines with negative value:
# less promoter_up500.Araport11.position | grep '-'
# Chr5 -498 2 AT5G00730
# ChrC -496 4 ATCG00010
# ChrC -117 383 ATCG00020
# Remove these three genes by adding: awk '$2 >= 0' 

### down500
# output path: /data/rep2_bismark_yao/methylation_extractor/promoter/up500
# ref_genome path:  /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_up500.Araport11.position
zcat /home/vcm/araport_reference/Araport11.gtf.gz |perl -alne '{next unless $F[2] eq "gene" ;/gene_id \"(.*?)\";/; print "$F[0]\t$F[3]\t$F[3]\t$1" }' | awk '{ $2 = $3; $3 = $3 + 500; print }' | awk '$2 >= 0' | sed 's/ /\t/g' >promoter_down500.Araport11.position

### up1k
## prepare promoter position by minus 1kbp of the START
# output path: /data/rep2_bismark_yao/methylation_extractor/promoter/up1k
# ref_genome path:  /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_up1k.Araport11.position
# Remove these three genes by adding: awk '$2 >= 0'
zcat /home/vcm/araport_reference/Araport11.gtf.gz |perl -alne '{next unless $F[2] eq "gene" ;/gene_id \"(.*?)\";/; print "$F[0]\t$F[3]\t$F[3]\t$1" }' | awk '{ $2 = $3 - 1000; print }' | awk '$2 >= 0' | sed 's/ /\t/g' >promoter_up1k.Araport11.position


### delimite the csv file (wait for matching for gene_id) by tab in tsv format
for file in CHG_context_*; do sed 's/,/\t/g' $file > ${file%_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt}.tsv; echo ${file%_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt}.tsv; done
for file in CHH_context_*; do sed 's/,/\t/g' $file > ${file%_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt}.tsv; echo ${file%_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt}.tsv; done
for file in CpG_context_*; do sed 's/,/\t/g' $file > ${file%_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt}.tsv; echo ${file%_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt}.tsv; done
#${file%_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.txt}.tsv

#change tsv to bed format
for file in CHG*.tsv; do tail -n +2 $file | awk -v OFS='\t' '{print $3, $4, $4, $2, $5}' > ./bed/${file%.tsv}.bed; echo ${file%.tsv}.bed; done
for file in CHH*.tsv; do tail -n +2 $file | awk -v OFS='\t' '{print $3, $4, $4, $2, $5}' > ./bed/${file%.tsv}.bed; echo ${file%.tsv}.bed; done
for file in CpG*.tsv; do tail -n +2 $file | awk -v OFS='\t' '{print $3, $4, $4, $2, $5}' > ./bed/${file%.tsv}.bed; echo ${file%.tsv}.bed; done

# make first column to be 'Chr' + num
for file in CHG_context_*.bed; do sed -i 's/^/Chr/' $file; echo $file; done
for file in CHH_context_*.bed; do sed -i 's/^/Chr/' $file; echo $file; done
for file in CpG_context_*.bed; do sed -i 's/^/Chr/' $file; echo $file; done

###
#3# Step 3: bedtools intersect
###

### use bedtools to intersect overlapped column
#bedtools intersect -a /data/rep2_bismark_yao/methylation_extractor/bed/CHG_context_3days107-LFK10433.bed -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb
for file in CHG_context_*.bed; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > ${file%.bed}_matched.bed; echo ${file%.bed}_matched.bed; done
for file in CHH_context_*.bed; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > ${file%.bed}_matched.bed; echo ${file%.bed}_matched.bed; done
for file in CpG_context_*.bed; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > ${file%.bed}_matched.bed; echo ${file%.bed}_matched.bed; done

#up500 promoter
for file in CHG_context_*.bed; do bedtools intersect -a $file -b /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_up500.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > /data/rep2_bismark_yao/methylation_extractor/promoter/up500/${file%.bed}_up500_matched.bed; echo ${file%.bed}_up500_matched.bed; done &
for file in CHH_context_*.bed; do bedtools intersect -a $file -b /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_up500.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > /data/rep2_bismark_yao/methylation_extractor/promoter/up500/${file%.bed}_up500_matched.bed; echo ${file%.bed}_up500_matched.bed; done &
for file in CpG_context_*.bed; do bedtools intersect -a $file -b /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_up500.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > /data/rep2_bismark_yao/methylation_extractor/promoter/up500/${file%.bed}_up500_matched.bed; echo ${file%.bed}_up500_matched.bed; done &

#down500 promoter
for file in ./CHG_context_*.bed; do bedtools intersect -a $file -b /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_down500.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > /data/rep2_bismark_yao/methylation_extractor/promoter/down500/${file%.bed}_down500_matched.bed; echo ${file%.bed}_down500_matched.bed; done &
for file in ./CHH_context_*.bed; do bedtools intersect -a $file -b /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_down500.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > /data/rep2_bismark_yao/methylation_extractor/promoter/down500/${file%.bed}_down500_matched.bed; echo ${file%.bed}_down500_matched.bed; done &
for file in ./CpG_context_*.bed; do bedtools intersect -a $file -b /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_down500.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > /data/rep2_bismark_yao/methylation_extractor/promoter/down500/${file%.bed}_down500_matched.bed; echo ${file%.bed}_down500_matched.bed; done &

#up1k promoter
for file in CHG_context_*.bed; do bedtools intersect -a $file -b /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_up1k.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > /data/rep2_bismark_yao/methylation_extractor/promoter/up1k/${file%.bed}_up1k_matched.bed; echo ${file%.bed}_up1k_matched.bed; done &
for file in CHH_context_*.bed; do bedtools intersect -a $file -b /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_up1k.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > /data/rep2_bismark_yao/methylation_extractor/promoter/up1k/${file%.bed}_up1k_matched.bed; echo ${file%.bed}_up1k_matched.bed; done &
for file in CpG_context_*.bed; do bedtools intersect -a $file -b /data/rep2_bismark_yao/methylation_extractor/promoter/promoter_up1k.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > /data/rep2_bismark_yao/methylation_extractor/promoter/up1k/${file%.bed}_up1k_matched.bed; echo ${file%.bed}_up1k_matched.bed; done &

# Remove the last duplicated column
# output path: /data/rep2_bismark_yao/methylation_extractor/annotated_output
for file in *_matched.bed; do cut -f 1-7 $file > /data/rep2_bismark_yao/methylation_extractor/annotated_output/$file; echo $file; done
# /data/rep2_bismark_yao/methylation_extractor/promoter/up500
for file in *_matched.bed; do cut -f 1-7 $file > /data/rep2_bismark_yao/methylation_extractor/promoter/up500/annotated_output/$file; echo $file; done
# /data/rep2_bismark_yao/methylation_extractor/promoter/down500
for file in *_matched.bed; do cut -f 1-7 $file > /data/rep2_bismark_yao/methylation_extractor/promoter/down500/annotated_output/$file; echo $file; done
# /data/rep2_bismark_yao/methylation_extractor/promoter/up1k
for file in *_matched.bed; do cut -f 1-7 $file > /data/rep2_bismark_yao/methylation_extractor/promoter/up1k/annotated_output/$file; echo $file; done

####### ignore
# ### For IGV visualization, make a bedGraph file
# # Need a value for each methylated site!
# # in annotated output
# # for file in *.deduplicated.bismark.cov.bed; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,3,4,10 -c 10 -o collapse | awk -v OFS='\t' '{print $1, $2, $3, $4, $5}' > ${file%.bed}_matched.bed; echo ${file%.bed}_matched.bed; done
# # for file in *_matched.bed; do  awk -v OFS='\t' '{print $1, $2, $2, $2, $5}' $file
# bedtools intersect -a 3days107-LFK10433_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.bismark.cov.bed -b CHG_context_3days107-LFK10433.bed -wa -wb
# for file1 in testuniq_batmeth_*D107X*; do echo $file1; for file2 in testuniq_batmeth_*DWTX*; do echo $file2; join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.5,2.6 <(sort -k2 $file1) <(sort -k2 $file2) > ${file1%107X.tsv}X_commonlines.tsv; echo ${file1%107X.tsv}X_commonlines.tsv; done; done &
# join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,2.4 <(sort -k2 CHG_context_3days107-LFK10433.bed) <(sort -k2 3days107-LFK10433_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.bismark.cov.bed) > CHG_3days107.count.bed &
# join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,2.4 <(sort -k2 CHH_context_3days107-LFK10433.bed) <(sort -k2 3days107-LFK10433_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.bismark.cov.bed) > CHH_3days107.count.bed &
# join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,2.4 <(sort -k2 CpG_context_3days107-LFK10433.bed) <(sort -k2 3days107-LFK10433_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.bismark.cov.bed) > CpG_3days107.count.bed &

# join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,2.4 <(sort -k2 CHG_context_S0day107-LFK10431.bed) <(sort -k2 S0day107-LFK10431_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.bismark.cov.bed) > CHG_S0day107.count.bed &
# join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,2.4 <(sort -k2 CHH_context_S0day107-LFK10431.bed) <(sort -k2 S0day107-LFK10431_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.bismark.cov.bed) > CHH_S0day107.count.bed &
# join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,2.4 <(sort -k2 CpG_context_S0day107-LFK10431.bed) <(sort -k2 S0day107-LFK10431_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.bismark.cov.bed) > CpG_S0day107.count.bed &

# # for file in *_context_3days107*.bed; do new_name=$(basename "$file" | cut -d '-' -f 1); join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,2.4 <(sort -k2 $file) <(sort -k2 3days107-LFK10433_L1_1_val_1_val_1_bismark_hisat2_pe.deduplicated.bismark.cov.bed) > ${new_name}.count.bed; echo ${new_name}.count.bed; done &
#####################
#####################
# Separate upmethylated or downmethylated
# use the -v option to pass shell variables to awk
for file in CHG_context_*_matched.bed; do awk -v outfile="${file%.bed}_X.bed" '$4 == "X" {print > outfile}' "$file"; awk -v outfile="${file%.bed}_x.bed" '$4 == "x" {print > outfile}' "$file"; echo ${file%.bed}_X.bed; echo ${file%.bed}_x.bed; done
for file in CHH_context_*_matched.bed; do awk -v outfile="${file%.bed}_H.bed" '$4 == "H" {print > outfile}' "$file"; awk -v outfile="${file%.bed}_h.bed" '$4 == "h" {print > outfile}' "$file"; echo ${file%.bed}_H.bed; echo ${file%.bed}_h.bed; done
for file in CpG_context_*_matched.bed; do awk -v outfile="${file%.bed}_Z.bed" '$4 == "Z" {print > outfile}' "$file"; awk -v outfile="${file%.bed}_z.bed" '$4 == "z" {print > outfile}' "$file"; echo ${file%.bed}_Z.bed; echo ${file%.bed}_z.bed; done

# first discard column 2, 3, 4. # a separate step for record, reorganize and output to bed5 file
### This step is essential for batmeth batDMR!!! ###
for file in *matched_*.bed; do awk -v OFS='\t' '{print $1, $2, $5, $6, $7}' $file > ${file%.bed}.bed5.tsv; echo ${file%.bed}.bed5.tsv; done

# second, based on column 7, count how many times does this gene occurs and add the count number as column 8.
# Please note that column 3 of the output file is meaningless as many sites are combined into 1
# sample: awk -F'\t' -v OFS='\t' '{count[$7]++; lines[$7]=$0} END {for (gene in count) print lines[gene], count[gene]}' CHG_context_3days107-LFK10433_matched_x.bed > test.bed
for file in *_context_*_matched_*.bed; do awk -F'\t' -v OFS='\t' '{count[$7]++; lines[$7]=$0} END {for (gene in count) print lines[gene], count[gene]}' $file > ${file%.bed}_counted8.tsv; echo ${file%.bed}_counted8.tsv; done 
# check file size: tree -h | grep counted8

# Discard unnecessary col 2 3 4
for file in *counted8*; do awk -v OFS='\t' '{print $1, $5, $6, $7, $8}' $file > ${file%_counted8.tsv}_counted5.tsv; echo ${file%_counted8.tsv}_counted5.tsv; done
# check file size: tree -h | grep counted5
######################

###
#4# Step 4: Prepare batmeth input by R
###

### BATMETH_PROCESS !!!

# ### Download to local for R analysis
# #path: /data/rep2_bismark_yao/methylation_extractor/annotated_output/*counted5*
# #scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/annotated_output/*counted5* "C:\Users\Jerry\OneDrive - Duke University\Prof. Lee Lab\2023Fall_meeting\met_dist_rep2"

# # gene body bed5.tsv
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/annotated_output/*matched_x.bed5.tsv "D:\data\BS_CXX_coverage\bed5\unmet"
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/annotated_output/*matched_X.bed5.tsv "D:\data\BS_CXX_coverage\bed5\met"
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/annotated_output/*matched_h.bed5.tsv "D:\data\BS_CXX_coverage\bed5\unmet"
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/annotated_output/*matched_H.bed5.tsv "D:\data\BS_CXX_coverage\bed5\met"
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/annotated_output/*matched_z.bed5.tsv "D:\data\BS_CXX_coverage\bed5\unmet"
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/annotated_output/*matched_Z.bed5.tsv "D:\data\BS_CXX_coverage\bed5\met"

# # up500 promoter bed5.tsv
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/promoter/up500/bed/*matched_x.bed5.tsv "D:\data\BS_CXX_coverage\up500\unmet"
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/promoter/up500/bed/*matched_X.bed5.tsv "D:\data\BS_CXX_coverage\up500\met"
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/promoter/up500/bed/*matched_h.bed5.tsv "D:\data\BS_CXX_coverage\up500\unmet"
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/promoter/up500/bed/*matched_H.bed5.tsv "D:\data\BS_CXX_coverage\up500\met"
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/promoter/up500/bed/*matched_z.bed5.tsv "D:\data\BS_CXX_coverage\up500\unmet"
# scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/methylation_extractor/promoter/up500/bed/*matched_Z.bed5.tsv "D:\data\BS_CXX_coverage\up500\met"

# # Bismark coverage file
# # add 'Chr' before first col, change name to .csv.bed
# for file in *.cov; do cp $file $file.bed; echo $file.bed; done 
# for file in *.cov.bed; do sed -i 's/^/Chr/' $file; echo $file; done

# #for row in *.cov.bed; awk 'name=$0','index=$1' row for nrow in  do awk -F , '$4 == "something4"

### Run R processing function: batmeth_processing
# After R preparing for Batmeth2 input
# ### upload file to vcm@10.200.13.52 & cm@dku-vcm-3031.vm.duke.edu
# # CHG
# scp "D:/data/BS_CXX_coverage/up500/batmeth_*x.tsv" vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/Batmeth2/input
# # CHH
# scp "D:/data/BS_CXX_coverage/up500/batmeth_*h.tsv" vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/Batmeth2/input
# # CpG
# scp "D:/data/BS_CXX_coverage/up500/batmeth_*z.tsv" vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/Batmeth2/input

###
#5# Step 5: Run batmeth batDMR
###

# batmeth2 path: /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2
# batmeth prepare ref genome in path ref_genome, TAIR10_chr_all.fas
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 index -g ./TAIR10_chr_all.fas
# Genome path: /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas
# Diff meth analysis
for file in batmeth_*D*.tsv; do less $file | sort | uniq > testuniq_$file; echo testuniq_$file; done

# CHG
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CHG_batmeth_3Dx.dmraw.txt -1 batmeth_3D107x.tsv -2 batmeth_3DWTx.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CHG_batmeth_7Dx.dmraw.txt -1 batmeth_7D107x.tsv -2 batmeth_7DWTx.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CHG_batmeth_0Dx.dmraw.txt -1 batmeth_0D107x.tsv -2 batmeth_0DWTx.tsv &


for file1 in testuniq_batmeth_*D107X*; do echo $file1; for file2 in testuniq_batmeth_*DWTX*; do echo $file2; join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.5,2.6 <(sort -k2 $file1) <(sort -k2 $file2) > ${file1%107X.tsv}X_commonlines.tsv; echo ${file1%107X.tsv}X_commonlines.tsv; done; done &
for file in testuniq_batmeth_*DX_commonlines.tsv; do awk '{print $1, $2, $3, $4, $5, $6}' $file | sort > ${file%X_commonlines.tsv}107X.tsv; echo ${file%X_commonlines.tsv}107X.tsv; done &
for file in testuniq_batmeth_*DX_commonlines.tsv; do awk '{print $1, $2, $3, $4, $7, $8}' $file | sort > ${file%X_commonlines.tsv}WTX.tsv; echo ${file%X_commonlines.tsv}WTX.tsv; done &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm testuniq_batmeth_0DX.dmraw.txt -1 testuniq_batmeth_0D107X.tsv -2 testuniq_batmeth_0DWTX.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm testuniq_batmeth_3DX.dmraw.txt -1 testuniq_batmeth_3D107X.tsv -2 testuniq_batmeth_3DWTX.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm testuniq_batmeth_7DX.dmraw.txt -1 testuniq_batmeth_7D107X.tsv -2 testuniq_batmeth_7DWTX.tsv &
# rename output fine
for file in testuniq_batmeth_*DX.dmraw.txt; do new_name="CHG_${file#*_}"; mv "$file" "$new_name"; done

# CHH
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CHH_batmeth_3Dh.dmraw.txt -1 batmeth_3D107h.tsv -2 batmeth_3DWTh.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CHH_batmeth_7Dh.dmraw.txt -1 batmeth_7D107h.tsv -2 batmeth_7DWTh.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CHH_batmeth_0Dh.dmraw.txt -1 batmeth_0D107h.tsv -2 batmeth_0DWTh.tsv &

/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CHH_batmeth_3DH.dmraw.txt -1 batmeth_3D107H.tsv -2 batmeth_3DWTH.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CHH_batmeth_7DH.dmraw.txt -1 batmeth_7D107H.tsv -2 batmeth_7DWTH.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CHH_batmeth_0DH.dmraw.txt -1 batmeth_0D107H.tsv -2 batmeth_0DWTH.tsv &
# CpG
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CpG_batmeth_3Dz.dmraw.txt -1 batmeth_3D107z.tsv -2 batmeth_3DWTz.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CpG_batmeth_7Dz.dmraw.txt -1 batmeth_7D107z.tsv -2 batmeth_7DWTz.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CpG_batmeth_0Dz.dmraw.txt -1 batmeth_0D107z.tsv -2 batmeth_0DWTz.tsv &

for file1 in testuniq_batmeth_*D107Z*; do echo $file1; for file2 in testuniq_batmeth_*DWTZ*; do echo $file2; join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.5,2.6 <(sort -k2 $file1) <(sort -k2 $file2) > ${file1%107Z.tsv}Z_commonlines.tsv; echo ${file1%107Z.tsv}Z_commonlines.tsv; done; done &
for file in testuniq_batmeth_*DZ_commonlines.tsv; do awk '{print $1, $2, $3, $4, $5, $6}' $file | sort > ${file%Z_commonlines.tsv}107Z.tsv; echo ${file%Z_commonlines.tsv}107Z.tsv; done &
for file in testuniq_batmeth_*DZ_commonlines.tsv; do awk '{print $1, $2, $3, $4, $7, $8}' $file | sort > ${file%Z_commonlines.tsv}WTZ.tsv; echo ${file%Z_commonlines.tsv}WTZ.tsv; done &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm testuniq_batmeth_0D107Z.dmraw.txt -1 testuniq_batmeth_0D107Z.tsv -2 testuniq_batmeth_0DWTZ.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm testuniq_batmeth_3D107Z.dmraw.txt -1 testuniq_batmeth_3D107Z.tsv -2 testuniq_batmeth_3DWTZ.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm testuniq_batmeth_7D107Z.dmraw.txt -1 testuniq_batmeth_7D107Z.tsv -2 testuniq_batmeth_7DWTZ.tsv &
# rename output fine
for file in testuniq_batmeth_*DZ.dmraw.txt; do new_name="CpG_${file#*_}"; mv "$file" "$new_name"; done


### up1k
# CHG
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/up1k_batmeth2/input/dmr_output/batmeth_0D_CHG_1k.dmraw.txt -1 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_0D107_CHG_1k.tsv -2 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_0DWT_CHG_1k.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/up1k_batmeth2/input/dmr_output/batmeth_3D_CHG_1k.dmraw.txt -1 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_3D107_CHG_1k.tsv -2 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_3DWT_CHG_1k.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/up1k_batmeth2/input/dmr_output/batmeth_7D_CHG_1k.dmraw.txt -1 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_7D107_CHG_1k.tsv -2 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_7DWT_CHG_1k.tsv &
# CHH
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/up1k_batmeth2/input/dmr_output/batmeth_0D_CHH_1k.dmraw.txt -1 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_0D107_CHH_1k.tsv -2 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_0DWT_CHH_1k.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/up1k_batmeth2/input/dmr_output/batmeth_3D_CHH_1k.dmraw.txt -1 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_3D107_CHH_1k.tsv -2 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_3DWT_CHH_1k.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/up1k_batmeth2/input/dmr_output/batmeth_7D_CHH_1k.dmraw.txt -1 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_7D107_CHH_1k.tsv -2 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_7DWT_CHH_1k.tsv &
# CpG
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/up1k_batmeth2/input/dmr_output/batmeth_0D_CpG_1k.dmraw.txt -1 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_0D107_CpG_1k.tsv -2 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_0DWT_CpG_1k.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/up1k_batmeth2/input/dmr_output/batmeth_3D_CpG_1k.dmraw.txt -1 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_3D107_CpG_1k.tsv -2 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_3DWT_CpG_1k.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/up1k_batmeth2/input/dmr_output/batmeth_7D_CpG_1k.dmraw.txt -1 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_7D107_CpG_1k.tsv -2 /data/rep2_bismark_yao/up1k_batmeth2/input/batmeth_7DWT_CpG_1k.tsv &
# terminate called after throwing an instance of 'std::logic_error'
# what():  Couldn't parse a line "Chr3  14199215        +       CpG     1225678621      NA".
#!# If contains 'NA', run these to remove those lines:
#!# awk '!/NA/' input_file.txt > output_file.txt {linux}
#!# {/r} 
# # Print the row numbers with NA values
# out_df <- in_df[complete.cases(in_df), ]

### down500
# CHG
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/down500_batmeth2/input/dmr_output/batmeth_0D_CHG_down500.dmraw.txt -1 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_0D107_CHG_down500.tsv -2 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_0DWT_CHG_down500.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/down500_batmeth2/input/dmr_output/batmeth_3D_CHG_down500.dmraw.txt -1 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_3D107_CHG_down500.tsv -2 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_3DWT_CHG_down500.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/down500_batmeth2/input/dmr_output/batmeth_7D_CHG_down500.dmraw.txt -1 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_7D107_CHG_down500.tsv -2 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_7DWT_CHG_down500.tsv &
# CHH
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/down500_batmeth2/input/dmr_output/batmeth_0D_CHH_down500.dmraw.txt -1 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_0D107_CHH_down500.tsv -2 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_0DWT_CHH_down500.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/down500_batmeth2/input/dmr_output/batmeth_3D_CHH_down500.dmraw.txt -1 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_3D107_CHH_down500.tsv -2 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_3DWT_CHH_down500.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/down500_batmeth2/input/dmr_output/batmeth_7D_CHH_down500.dmraw.txt -1 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_7D107_CHH_down500.tsv -2 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_7DWT_CHH_down500.tsv &
# CpG
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/down500_batmeth2/input/dmr_output/batmeth_0D_CpG_down500.dmraw.txt -1 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_0D107_CpG_down500.tsv -2 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_0DWT_CpG_down500.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/down500_batmeth2/input/dmr_output/batmeth_3D_CpG_down500.dmraw.txt -1 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_3D107_CpG_down500.tsv -2 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_3DWT_CpG_down500.tsv &
/home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/down500_batmeth2/input/dmr_output/batmeth_7D_CpG_down500.dmraw.txt -1 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_7D107_CpG_down500.tsv -2 /data/rep2_bismark_yao/down500_batmeth2/input/batmeth_7DWT_CpG_down500.tsv &


# test
# /data/rep2_bismark_yao/up500_batmeth2/input/test/batmeth_3D107X.tsv /data/rep2_bismark_yao/up500_batmeth2/input/test/batmeth_3DWTX.tsv
# /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/up500_batmeth2/input/test/testCHG_batmeth_3DX.dmraw.txt -1 /data/rep2_bismark_yao/up500_batmeth2/input/test/testbatmeth_3D107X.tsv -2 /data/rep2_bismark_yao/up500_batmeth2/input/test/testbatmeth_3DWTX.tsv &
# join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.5,2.6 <(sort -k2 testdedup_batmeth_3D107X.tsv) <(sort -k2 testdedup_batmeth_3DWTX.tsv) > common_lines_3D.tsv
# awk '{print $1, $2, $3, $4, $5, $6}' common_lines_3D.tsv | sort > "testuniq_batmeth_3D107X.tsv" 
# awk '{print $1, $2, $3, $4, $7, $8}' common_lines_3D.tsv | sort > "testuniq_batmeth_3DWTX.tsv" 
# /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm /data/rep2_bismark_yao/up500_batmeth2/input/test/testuniq_batmeth_3DX.dmraw.tsv -1 testuniq_batmeth_3D107X.tsv -2 testuniq_batmeth_3DWTX.tsv &
# It works!!!!!!!!!!
# for file in batmeth_*D*.tsv; do less $file | sort | uniq > testuniq_$file; echo testuniq_$file; done # before all the steps!
# less testbatmeth_3D107X.tsv | sort | uniq > testdedup_batmeth_3D107X.tsv
# testuniq_batmeth_0D107X.tsv  testuniq_batmeth_3D107X.tsv  testuniq_batmeth_7D107X.tsv
# join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.5,2.6 <(sort -k2 testuniq_batmeth_0D107X.tsv) <(sort -k2 testuniq_batmeth_0DWTX.tsv) > common_lines_0DX.tsv &
# awk '{print $1, $2, $3, $4, $5, $6}' common_lines_0DX.tsv | sort > "testuniq_batmeth_0D107X.tsv" &
# awk '{print $1, $2, $3, $4, $7, $8}' common_lines_0DX.tsv | sort > "testuniq_batmeth_0DWTX.tsv" &
# join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.5,2.6 <(sort -k2 testuniq_batmeth_3D107X.tsv) <(sort -k2 testuniq_batmeth_3DWTX.tsv) > common_lines_3DX.tsv &
# awk '{print $1, $2, $3, $4, $5, $6}' common_lines_3DX.tsv | sort > "testuniq_batmeth_3D107X.tsv" &
# awk '{print $1, $2, $3, $4, $7, $8}' common_lines_3DX.tsv | sort > "testuniq_batmeth_3DWTX.tsv" &
# join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.5,2.6 <(sort -k2 testuniq_batmeth_7D107X.tsv) <(sort -k2 testuniq_batmeth_7DWTX.tsv) > common_lines_7DX.tsv &
# awk '{print $1, $2, $3, $4, $5, $6}' common_lines_7DX.tsv | sort > "testuniq_batmeth_7D107X.tsv" &
# awk '{print $1, $2, $3, $4, $7, $8}' common_lines_7DX.tsv | sort > "testuniq_batmeth_7DWTX.tsv" &

# # for file1 in *batmeth_*D107X*; do echo $file1; for file2 in *batmeth_*DWTX*; do echo $file2; done; done
# for file1 in *batmeth_*D107X*; do echo $file1; for file2 in *batmeth_*DWTX*; do echo $file2; join -1 2 -2 2 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.5,2.6 <(sort -k2 $file1) <(sort -k2 $file2) > ${file1%107X.tsv}X_commonlines.tsv; echo ${file1%107X.tsv}X_commonlines.tsv; done; done &
# for file in batmeth_*DX_commonlines.tsv; do awk '{print $1, $2, $3, $4, $5, $6}' $file | sort > testuniq_${file%X_commonlines.tsv}107X.tsv; echo testuniq_${file%X_commonlines.tsv}107X.tsv; done &
# for file in batmeth_*DX_commonlines.tsv; do awk '{print $1, $2, $3, $4, $7, $8}' $file | sort > testuniq_${file%X_commonlines.tsv}WTX.tsv; echo testuniq_${file%X_commonlines.tsv}WTX.tsv; done &
# for file1 in testuniq_batmeth_*D107X.tsv; do echo $file1; for file2 in testuniq_batmeth_*DWTX.tsv; do echo $file2; /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm CHG_batmeth_3DX.dmraw.txt -1 $file1 -2 $file2; echo CHG_batmeth_3DX.dmraw.txt; done; done
# # change all filtration to the lowest
# # /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2 batDMR -pvalue 1 -FDR 1 -methdiff 0.01 -g /data/rep2_bismark_yao/Batmeth2/ref_genome/TAIR10_chr_all.fas -o_dm batmeth_3D.dmraw.txt -1 batmeth_3D107x.tsv -2 batmeth_3DWTx.tsv
# # in case methylated count > total count, filter (CHGX, CpGZ)
# for file in batmeth_*D*.tsv; do 

###
#6# Step 6: Rearrange batDMR output [dmraw] to [dm5] [unmatched_dm5] for annotation
###

#####
### up500
# the output file batmeth_*D.dm.txt:
# Chr1    32727   +       CHG     3.6033e-11      7.57843e-11     74      199     3       130     0.348782
# Select the needed field (1,2,5,6,11): Chr1    32727     3.6033e-11      7.57843e-11     0.348782
# CHG
for file in CHG*x.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &
for file in CHG*X.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &

# CHH
for file in CHH*h.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &
for file in CHH*H.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &

# CpG
for file in CpG*z.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &
for file in CpG*Z.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &

# Output: Chr1    32727 (Start-500)   3.6033e-11      7.57843e-11     0.348782
# pair again with Araport11.position to get gene name back
# get start position back by col2 + 500
# CHG
for file in CHG*.dm5.tsv; do awk -v OFS='\t' '{ $2 = $2 + 500; print $1, $2, $2, $3, $4, $5}' $file > ${file%.dm5.tsv}.unmatched.dm5.tsv; echo ${file%.dm5.tsv}.unmatched.dm5.tsv; done &
# CHH
for file in CHH*.dm5.tsv; do awk -v OFS='\t' '{ $2 = $2 + 500; print $1, $2, $2, $3, $4, $5}' $file > ${file%.dm5.tsv}.unmatched.dm5.tsv; echo ${file%.dm5.tsv}.unmatched.dm5.tsv; done &
# CpG
for file in CpG*.dm5.tsv; do awk -v OFS='\t' '{ $2 = $2 + 500; print $1, $2, $2, $3, $4, $5}' $file > ${file%.dm5.tsv}.unmatched.dm5.tsv; echo ${file%.dm5.tsv}.unmatched.dm5.tsv; done &
#######
###  up1k
# CHG
for file in batmeth_*D_CHG_1k.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &
# CHH
for file in batmeth_*D_CHH_1k.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &
# CpG
for file in batmeth_*D_CpG_1k.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &

# get start position back by col2 + 1000
# CHG
for file in batmeth_*D_CHG_1k.dm5.tsv; do awk -v OFS='\t' '{ $2 = $2 + 1000; print $1, $2, $2, $3, $4, $5}' $file > ${file%.dm5.tsv}.unmatched.dm5.tsv; echo ${file%.dm5.tsv}.unmatched.dm5.tsv; done &
# CHH
for file in batmeth_*D_CHH_1k.dm5.tsv; do awk -v OFS='\t' '{ $2 = $2 + 1000; print $1, $2, $2, $3, $4, $5}' $file > ${file%.dm5.tsv}.unmatched.dm5.tsv; echo ${file%.dm5.tsv}.unmatched.dm5.tsv; done &
# CpG
for file in batmeth_*D_CpG_1k.dm5.tsv; do awk -v OFS='\t' '{ $2 = $2 + 1000; print $1, $2, $2, $3, $4, $5}' $file > ${file%.dm5.tsv}.unmatched.dm5.tsv; echo ${file%.dm5.tsv}.unmatched.dm5.tsv; done &
######
###  down500
# CHG
for file in batmeth_*D_CHG_down500.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &
# CHH
for file in batmeth_*D_CHH_down500.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &
# CpG
for file in batmeth_*D_CpG_down500.dmraw.txt; do awk -v OFS='\t' '{print $1, $2, $5, $6, $11}' $file > ${file%raw.txt}5.tsv; echo ${file%raw.txt}5.tsv; done &

# get start position back by col2 - 500
# CHG
for file in batmeth_*D_CHG_down500.dm5.tsv; do awk -v OFS='\t' '{ $2 = $2 - 500; print $1, $2, $2, $3, $4, $5}' $file > ${file%.dm5.tsv}.unmatched.dm5.tsv; echo ${file%.dm5.tsv}.unmatched.dm5.tsv; done &
# CHH
for file in batmeth_*D_CHH_down500.dm5.tsv; do awk -v OFS='\t' '{ $2 = $2 - 500; print $1, $2, $2, $3, $4, $5}' $file > ${file%.dm5.tsv}.unmatched.dm5.tsv; echo ${file%.dm5.tsv}.unmatched.dm5.tsv; done &
# CpG
for file in batmeth_*D_CpG_down500.dm5.tsv; do awk -v OFS='\t' '{ $2 = $2 - 500; print $1, $2, $2, $3, $4, $5}' $file > ${file%.dm5.tsv}.unmatched.dm5.tsv; echo ${file%.dm5.tsv}.unmatched.dm5.tsv; done &
######

###
#7# Step 7: Annotate [dm] output with geneID into final output
###
#####
### up500
# for file in *batmeth_*D.unmatched.dm5.tsv; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb; done
# Output: Chr5    26831185        26831185        1.11264e-61     2.73127e-60     -0.393073       Chr5    26831185        26833903        AT5G67250
# Final Output: Chr1    1728658 1729789 AT1G05783       3.40384e-180    3.42944e-178    -0.269964
# CHG
for file in CHG_batmeth_*D*.unmatched.dm5.tsv; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,8,9,10,4,5,6 -c 7 -o collapse | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' > ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; echo ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; done &
# CHH
for file in CHH_batmeth_*D*.unmatched.dm5.tsv; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,8,9,10,4,5,6 -c 7 -o collapse | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' > ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; echo ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; done &
# CpG
for file in CpG_batmeth_*D*.unmatched.dm5.tsv; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,8,9,10,4,5,6 -c 7 -o collapse | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' > ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; echo ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; done &
#####

### up1k
# CHG
for file in batmeth_*D_CHG_1k.unmatched.dm5.tsv; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,8,9,10,4,5,6 -c 7 -o collapse | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' > ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; echo ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; done &
# CHH
for file in batmeth_*D_CHH_1k.unmatched.dm5.tsv; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,8,9,10,4,5,6 -c 7 -o collapse | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' > ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; echo ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; done &
# CpG
for file in batmeth_*D_CpG_1k.unmatched.dm5.tsv; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,8,9,10,4,5,6 -c 7 -o collapse | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' > ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; echo ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; done &
#####

### down500
# CHG
for file in batmeth_*D_CHG_down500.unmatched.dm5.tsv; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,8,9,10,4,5,6 -c 7 -o collapse | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' > ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; echo ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; done &
# CHH
for file in batmeth_*D_CHH_down500.unmatched.dm5.tsv; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,8,9,10,4,5,6 -c 7 -o collapse | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' > ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; echo ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; done &
# CpG
for file in batmeth_*D_CpG_down500.unmatched.dm5.tsv; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,8,9,10,4,5,6 -c 7 -o collapse | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' > ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; echo ${file%.unmatched.dm5.tsv}.matched.dm5.tsv; done &
#####

# Now download to local computer for R analysis Venn plot
# CHG
scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/up500_batmeth2/input/CHG_batmeth_*Dx.matched.dm5.tsv "D:/data/BS_CXX_coverage/up500/output/down"
scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/up500_batmeth2/input/CHG_batmeth_*DX.matched.dm5.tsv "D:/data/BS_CXX_coverage/up500/output/up"
# CHH
scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/up500_batmeth2/input/CHH_batmeth_*Dh.matched.dm5.tsv "D:/data/BS_CXX_coverage/up500/output/down"
scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/up500_batmeth2/input/CHH_batmeth_*DH.matched.dm5.tsv "D:/data/BS_CXX_coverage/up500/output/up"
# CpG
scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/up500_batmeth2/input/CpG_batmeth_*Dz.matched.dm5.tsv "D:/data/BS_CXX_coverage/up500/output/down"
scp vcm@dku-vcm-3031.vm.duke.edu:/data/rep2_bismark_yao/up500_batmeth2/input/CpG_batmeth_*DZ.matched.dm5.tsv "D:/data/BS_CXX_coverage/up500/output/up"


# less RIPChIP | grep ',' | wc -l
# less /home/vcm/araport_reference/protein_coding.Araport11.position | grep 'AT1G01010'
#########################################

# For Zhujun proportion heatmap
tail -n +2 CpG_batmeth_3DZ.dmraw.txt_combined.file.temp.txt | tr ':' '\t' | awk -v OFS='\t' '{ $2 = $2 + 500; print $1, $2, $2, $4, $5, $6, $7, $8}' > ./test/CpG_3DWT_proportion.test.tsv
for file in *.dmraw.txt_combined.file.temp.txt; do tail -n +2 $file | tr ':' '\t' | awk -v OFS='\t' '{ $2 = $2 + 500; print $1, $2, $2, $4, $5, $6, $7, $8}' > ./proportion_heatmap/${file%.dmraw.txt_combined.file.temp.txt}_proportion.test.tsv; echo ${file%.dmraw.txt_combined.file.temp.txt}_proportion.test.tsv; done

for file in *.dmraw.txt; do awk -v OFS='\t' '{ $2 = $2 + 500; print $1, $2, $2, $4, $7, $8, $9, $10}' $file > ./proportion_heatmap/${file%.dmraw.txt}_proportion.test.tsv; echo ${file%.dmraw.txt}_proportion.test.tsv; done
# # pair with genome
# head CpG_3DWT_proportion.test.tsv
# Chr1    6788    CpG     698     1885    27      96
# Chr1    23121   CpG     3       102     1       18
# Chr1    28500   CpG     13346   14165   3552    3794
# Chr1    31170   CpG     8025    9063    1628    1938

# $1, $2, $4, $5, $6, $7, $8, $12
# Chr5    26969516        26969516        CpG     5       126                     Chr5    26969516        26970668        AT5G67640
for file in *proportion*; do bedtools intersect -a $file -b /home/vcm/araport_reference/protein_coding.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,6,7,8,12 -c 8 -o collapse | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8}' > ./proportion_heatmap/${file%.tsv}.matched.tsv; echo ${file%.tsv}.matched.tsv; done &
