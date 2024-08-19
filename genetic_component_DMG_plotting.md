# Step1: construct GTF files with only 1 type of genetic component

> genetic component in GTF:
gene
exon
five_prime_UTR
three_prime_UTR
transposabal_element
CDS
mRNA
protein
miRNA_primary_transcript
miRNA

> path:

```
(base) vcm@test:~/araport_reference$ pwd
/home/vcm/araport_reference
├── [2.9M]  Araport11_gene.gtf
├── [ 84M]  Araport11.gtf
├── [6.3M]  Araport11.gtf.gz
├── [139M]  Araport11_seq_20220914.fa
├── [ 28M]  Araport11_seq_20220914.fa.gz
├── [3.7M]  Araport11_transposon.gtf
├── [6.3M]  araport_reference
├── [ 11M]  exon.Araport11.position
├── [2.5M]  five_prime_UTR.Araport11.position
├── [1.2M]  gene.Araport11.position
├── [   0]  intron.Araport11.position
├── [1.0M]  protein_coding.Araport11.position
└── [2.2M]  three_prime_UTR.Araport11.position
```


> input file - Araport11.gtf
```
Chr1    Araport11       transposable_element_gene       433031  433819  .       -       .       transcript_id "AT1G02228"; gene_id "AT1G02228";
Chr1    Araport11       transposable_element_gene       846664  847739  .       +       .       transcript_id "AT1G03420"; gene_id "AT1G03420";
Chr1    Araport11       transposable_element_gene       2415041 2415970 .       +       .       transcript_id "AT1G07800"; gene_id "AT1G07800";
Chr1    Araport11       transposable_element_gene       2531695 2534786 .       -       .       transcript_id "AT1G08105"; gene_id "AT1G08105";
Chr1    Araport11       transposable_element_gene       2790290 2793641 .       +       .       transcript_id "AT1G08735"; gene_id "AT1G08735";
...
```

> output file - protein_coding.Araport11.position
```
Chr1    3631    5899    AT1G01010
Chr1    6788    9130    AT1G01020
Chr1    11649   13714   AT1G01030
Chr1    23121   31227   AT1G01040
Chr1    28500   28706   AT1G01046
Chr1    31170   33171   AT1G01050
Chr1    33365   37871   AT1G01060
...
```

> output file - exon.Araport11.position
```
Chr1    3631    3913    AT1G01010       exon
Chr1    3996    4276    AT1G01010       exon
Chr1    4486    4605    AT1G01010       exon
Chr1    4706    5095    AT1G01010       exon
Chr1    5174    5326    AT1G01010       exon
...
```

> output file - five_prime_UTR.Araport11.position
```
Chr1    3631    3759    AT1G01010       five_prime_UTR
Chr1    8667    9130    AT1G01020       five_prime_UTR
Chr1    8667    8737    AT1G01020       five_prime_UTR
Chr1    8443    8464    AT1G01020       five_prime_UTR
Chr1    8571    9130    AT1G01020       five_prime_UTR
Chr1    8443    8464    AT1G01020       five_prime_UTR
...
```

> output file - three_prime_UTR.Araport11.position
```
Chr1    5631    5899    AT1G01010       three_prime_UTR
Chr1    6788    6914    AT1G01020       three_prime_UTR
Chr1    6788    7069    AT1G01020       three_prime_UTR
Chr1    7157    7314    AT1G01020       three_prime_UTR
Chr1    6788    6914    AT1G01020       three_prime_UTR
...
```

> output file - transposable_element_gene.Araport11.position
```
working!!!
```

> Code
```
(protein_coding)
zcat Araport11.gtf.gz |perl -alne '{next unless $F[2] eq "gene" ;/gene_id \"(.*?)\";/; print "$F[0]\t$F[3]\t$F[4]\t$1" }' > protein_coding.Araport11.position
(exon)
zcat Araport11.gtf.gz |perl -alne '{next unless $F[2] eq "exon" ;/gene_id \"(.*?)\";/; print "$F[0]\t$F[3]\t$F[4]\t$1\t$F[2]" }' > exon.Araport11.position
(five_prime_UTR)
zcat Araport11.gtf.gz |perl -alne '{next unless $F[2] eq "five_prime_UTR" ;/gene_id \"(.*?)\";/; print "$F[0]\t$F[3]\t$F[4]\t$1\t$F[2]" }' > five_prime_UTR.Araport11.position
(three_prime_UTR)
zcat Araport11.gtf.gz |perl -alne '{next unless $F[2] eq "three_prime_UTR" ;/gene_id \"(.*?)\";/; print "$F[0]\t$F[3]\t$F[4]\t$1\t$F[2]" }' > three_prime_UTR.Araport11.position
zcat Araport11.gtf.gz |perl -alne '{next unless $F[2] eq "transposable_element_gene" ;/gene_id \"(.*?)\";/; print "$F[0]\t$F[3]\t$F[4]\t$1\t$F[2]" }' > transposable_element_gene.Araport11.position
```

# Step2: bed intersect

> input file - methylation.bedgraph from bismark, modified for IGV
> path:
```
/home/vcm/IGV/Methylation/*.bedGraph
/home/vcm/araport_reference
/home/vcm/BSseq_rep2_batmeth/genetic_component
```
> sample
```
track type=bedGraph
1       0       1       0
1       1       2       10
1       2       3       10
1       7       8       11.1111111111111
1       8       9       11.1111111111111
1       9       10      20
1       14      15      0
...
```

### 2.1 change the bedgraph file to bed file for bedtools

> Code
```
# make first column to be 'Chr' + num
```
```
for file in /home/vcm/IGV/Methylation/*.bedGraph; do tail -n +2 $file > ${file%.bedGraph}_unannotated.bed; echo ${file%.bedGraph}_unannotated.bed; done

mv ~/IGV/Methylation/*_unannotated.bed ~/BSseq_rep2_batmeth/genetic_component

for file in *_unannotated.bed; do sed -i 's/^/Chr/' $file; echo $file; done
```

> output file - *_unannotated.bed
```
Chr1    0       1       0
Chr1    1       2       10
Chr1    2       3       10
Chr1    7       8       11.1111111111111
Chr1    8       9       11.1111111111111
Chr1    9       10      20
...
```

### 2.2 bedtools intersect
*change to an env installed with bedtools*

```
bedtools intersect -a A.bed -b B.Araport11.position -wa -wb

for file in CHG_context_*.bed; do bedtools intersect -a A.bed -b B.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,4,5,7,8,9 -c 7 -o collapse > C.bed; echo $file; done &
```

> Trial sample ... completed!
bedtools intersect -a 3days107_unannotated.bed -b /home/vcm/araport_reference/exon.Araport11.position -wa -wb
> output
```
Chr1    191183  191184  0       Chr1    191138  191475  AT1G01520       exon
Chr1    191184  191185  0       Chr1    191138  191258  AT1G01520       exon
Chr1    191184  191185  0       Chr1    191138  191258  AT1G01520       exon
...
```
> So, we only want to keep column 1,2,3,4,8,9
> (Chr1    191183  191184  0       AT1G01520       exon)
bedtools intersect -a 3days107_unannotated.bed -b /home/vcm/araport_reference/exon.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,3,4,8 -c 9 -o collapse > test.bed

> Code
#### Exon
```
for file in *_unannotated.bed; do bedtools intersect -a $file -b /home/vcm/araport_reference/exon.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,3,4,8 -c 9 -o collapse > ${file%_unannotated.bed}_exon_annotated.bed; echo ${file%_unannotated.bed}_exon_annotated.bed; done &
```
#### five_prime_UTR
```
for file in *_unannotated.bed; do bedtools intersect -a $file -b /home/vcm/araport_reference/five_prime_UTR.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,3,4,8 -c 9 -o collapse > ${file%_unannotated.bed}_5UTR_annotated.bed; echo ${file%_unannotated.bed}_5UTR_annotated.bed; done &
```
#### three_prime_UTR
```
for file in *_unannotated.bed; do bedtools intersect -a $file -b /home/vcm/araport_reference/three_prime_UTR.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,3,4,8 -c 9 -o collapse > ${file%_unannotated.bed}_3UTR_annotated.bed; echo ${file%_unannotated.bed}_3UTR_annotated.bed; done &
```
#### transposable_element_gene
```
for file in *_unannotated.bed; do bedtools intersect -a $file -b /home/vcm/araport_reference/transposable_element_gene.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,3,4,8 -c 9 -o collapse > ${file%_unannotated.bed}_transposable_element_gene_annotated.bed; echo ${file%_unannotated.bed}_transposable_element_gene_annotated.bed; done &
```

# Step 3: run batmeth2 plotting

### 3.1 index reference genome
`batMeth2 index -g genome_name.fa`
#### Araport11 sequence:
```
batmeth2 index -g Araport11.fas
```
> genome path: /home/vcm/araport_reference/araport11_ref/Araport11.fas
> However, Araport11 sequence file does not work. It returns with error:
```
[MM] /home/vcm/Batmeth2_download/BatMeth2/bin/batmeth2zation of DNA methylation data
        calculate differential DNA methylation cytosine or regiocross gene/TE/peak/bwame index Araport11.fas.batmeth2.fa
sh: 1: /home/vcm/Batmeth2_download/BatMeth2/bin/batmeth2zation: not found
sh: 2: calculate: not found
```
#### TAIR10 sequence
```
wget https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
mv download\?filePath\=Genes%2FTAIR10_genome_release%2FTAIR10_chromosome_files%2FTAIR10_chr_all.fas.gz TAIR10_chr_all.fas.gz
gunzip --keep TAIR10_chr_all.fas.gz

batmeth2 index TAIR10_chr_all.fas
```
> We can add prefix to the genome by parameter (--gp genome_path_prefix), e.g., TAIR10
```
[MM] /home/vcm/Batmeth2_download/BatMeth2/bin/genome2cg -g TAIR10_chr_all.fas
Converting C->T
Converting G->A
[MM] /home/vcm/Batmeth2_download/BatMeth2/bin/genomebinLen TAIR10_chr_all.fas
[MM] /home/vcm/Batmeth2_download/BatMeth2/bin/bwame index TAIR10_chr_all.fas.batmeth2.fa
[MM index] Pack FASTA... %.0/2f sec
[MM index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=478674536, availableWord=45680880
[BWTIncConstructFromPacked] 10 iterations done. 74387544 characters processed.
[BWTIncConstructFromPacked] 20 iterations done. 138353096 characters processed.
[BWTIncConstructFromPacked] 30 iterations done. 195201800 characters processed.
[BWTIncConstructFromPacked] 40 iterations done. 245725048 characters processed.
[BWTIncConstructFromPacked] 50 iterations done. 290626248 characters processed.
[BWTIncConstructFromPacked] 60 iterations done. 330530504 characters processed.
[BWTIncConstructFromPacked] 70 iterations done. 365993464 characters processed.
[BWTIncConstructFromPacked] 80 iterations done. 397509016 characters processed.
[BWTIncConstructFromPacked] 90 iterations done. 425516104 characters processed.
[BWTIncConstructFromPacked] 100 iterations done. 450404856 characters processed.
[BWTIncConstructFromPacked] 110 iterations done. 472522008 characters processed.
[bwt_gen] Finished constructing BWT in 114 iterations.
[MM index] 81.19 seconds elapse.
[MM index] Update BWT... 0.76 sec
[MM index] Pack forward-only FASTA... 0.58 sec
[MM index] Construct SA from BWT and Occ... 42.46 sec
[main] Version: 0.7.17-r1198-dirty
[main] CMD: /home/vcm/Batmeth2_download/BatMeth2/bin/bwame index TAIR10_chr_all.fas.batmeth2.fa
[main] Real time: 126.007 sec; CPU: 125.982 sec
```
> Genome path: /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas


### 3.2 Calulate mC across predefined regions
> Format:
```
 with bed file:
    methyGff -B -o gene.meth -G genome.fa -b gene.bed -m output.methrario.txt

with multiple bed file:
    methyGff -B -o expressed.gene.meth unexpressed.gene.meth \
        -G genome.fa -b expressed.gene.bed unexpressed.gene.bed -m output.methrario.txt
```
> output.methratio.txt
```
==> calmeth_7D107.methratio.txt <==
#chromsome      loci    strand  context C_count CT_count        methRatio       eff_CT_count    rev_G_count     rev_GA_count    MethContext       5context
Chr1    1       +       CHH     0       17      0.000000        17.0    50      50      U       NNCCC
Chr1    2       +       CHH     0       17      0.000000        17.0    50      50      U       NCCCT
Chr1    3       +       CHH     0       17      0.000000        17.0    50      50      U       CCCTA
Chr1    8       +       CHH     1       19      0.052632        19.0    53      53      U       AACCC
Chr1    9       +       CHH     0       18      0.000000        18.0    54      54      U       ACCCT
Chr1    10      +       CHH     0       23      0.000000        23.0    56      56      U       CCCTA
Chr1    15      +       CHH     0       23      0.000000        23.0    61      61      U       AACCC
Chr1    16      +       CHH     0       25      0.000000        25.0    61      61      U       ACCCT
Chr1    17      +       CHH     0       25      0.000000        25.0    61      61      U       CCCTA
...
```

> Trial on 1 file -- finished. Now working on draft!

#### Exon
> 7D
```
methyGff -B -o exon_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff -P -o exon_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff -B -o exon_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &

methyGff -P -o exon_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &
```

> 3D
```
methyGff -B -o exon_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff -P -o exon_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff -B -o exon_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &

methyGff -P -o exon_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &
```

> 0D
```
methyGff -B -o exon_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff -P -o exon_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff -B -o exon_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &

methyGff -P -o exon_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &
```

#### five_prime_UTR (5UTR)
> 7D
```
methyGff -B -o 5UTR_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff -P -o 5UTR_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff --TSS --TTS --GENE -o 5UTR_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff -B -o 5UTR_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &

methyGff -P -o 5UTR_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o 5UTR_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &
```

> 3D
```
methyGff -B -o 5UTR_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff -P -o 5UTR_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff --TSS --TTS --GENE -o 5UTR_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff -B -o 5UTR_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &

methyGff -P -o 5UTR_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o 5UTR_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &
```

> 0D
```
methyGff -B -o 5UTR_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff -P -o 5UTR_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff --TSS --TTS --GENE -o 5UTR_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff -B -o 5UTR_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &

methyGff -P -o 5UTR_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o 5UTR_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_5UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &
```

#### three_prime_UTR (3UTR) ...working
> 7D
```
methyGff -B -o 3UTR_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff -P -o 3UTR_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff --TSS --TTS --GENE -o 3UTR_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff -B -o 3UTR_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &

methyGff -P -o 3UTR_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o 3UTR_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &
```

> 3D
```
methyGff -B -o 3UTR_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff -P -o 3UTR_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff --TSS --TTS --GENE -o 3UTR_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff -B -o 3UTR_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &

methyGff -P -o 3UTR_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o 3UTR_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &
```

> 0D
```
methyGff -B -o 3UTR_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff -P -o 3UTR_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff --TSS --TTS --GENE -o 3UTR_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff -B -o 3UTR_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &

methyGff -P -o 3UTR_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o 3UTR_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_3UTR_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &
```

#### 
> 7D
```
methyGff -B -o exon_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff -P -o exon_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt &

methyGff -B -o exon_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &

methyGff -P -o exon_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt &
```

> 3D
```
methyGff -B -o exon_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff -P -o exon_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_3D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3D107.methratio.txt &

methyGff -B -o exon_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &

methyGff -P -o exon_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_3DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/3daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_3DWT.methratio.txt &
```

> 0D
```
methyGff -B -o exon_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff -P -o exon_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_0D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0day107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0D107.methratio.txt &

methyGff -B -o exon_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &

methyGff -P -o exon_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &

methyGff --TSS --TTS --GENE -o exon_0DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/S0dayWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_0DWT.methratio.txt &
```

> methGff:
> We can modify this parameter '-distance' to remove the up/down stream region:
> -d/--distance
> DNA methylation level distributions in body and <INT>-bp flanking sequences. The distance of upstream and downstream. default:2000
> can not be set to 0! Try -d 1 10 100 1000
> strange... all -d of TSSprofile does not work? It's that because I only used BED file specifying exon region???? But -d larger than 100 works???
> The TSS profile is flooded with nan, but the centerprofile seems workable. Try -center!


> don't include -d and -bl at the same time!
> Include -bl will return errors??? why???

> we can also change the gene length as well, since exon usually doesn't each 2K:
> -bl/--bodyLen
> Body length to which all regions will be fit. (default: same as -d, which is 2k bp)
> *by this mean, we can plot gene at different length category - but first we have to filter the GTF file / BED file based on length feature in prior.
> Command Format :   methyGff [options] -o <OUT_PREFIX> -G GENOME -gff <GFF file>/-gtf <GTF file>/-b <bed file>/-b4 <bed4 file> -m <from Split methratio outfile> [-B] [-P]


Usage:
	-o|--out         Output file prefix
	-G|--genome      Genome
	-m|--methratio   Methratio output file.
	-c|--coverage    >= <INT> coverage. default:4
	-C               <= <INT> coverage. default 600.
	-nC              >= <INT> Cs per bins or genes. default:1
	-gtf|-gff        Gtf/gff file
	-b|--BED         Bed file, chrom start end
	-b4              Bed file, chrom start end strand
	-b5              Bed file, chrom start end geneid strand
	-d|--distance    DNA methylation level distributions in body and <INT>-bp flanking sequences. The distance of upstream and downstream. default:2000
	-B|--body        For different analysis input format, gene/TEs body methylation level. [Different Methylation Gene(DMG/DMT...)]
	-P|--promoter    For different analysis input format.[Different Methylation Promoter(DMP)]
	--TSS            Caculate matrix for TSS. [Outfile: outPrefix.TSS.cg.n.txt]
	--TTS            Caculate matrix for TTS. [Outfile: outPrefix.TTS.cg.n.txt] 
	--GENE           Caculate matrix for GENE and flank 2k. [Outfile: outPrefix.GENE.cg.n.txt] 
	-s|--step        Gene body and their flanking sequences using an overlapping sliding window of 2% of the sequence length at a step of 1% of the sequence length. So default step: 0.01 (1%), used for profile.
	-hs|--heatmapstep        Gene body and their flanking sequences using an overlapping sliding window of 20% of the sequence length at a step of 10% of the sequence length. So default step: 0.1 (10%), used for heatmap.
	-bl|--bodyLen    Body length to which all regions will be fit. (default: same as -d)
	-S|--chromStep   Caculate the density of genes/TEs in chromsome using an overlapping sliding window of 100000bp at a step of 50000bp, must equal "-s" in Split.. default step: 50000(bp)
	-h|--help 

> Error message:
```
==> test1exon1_7DWT.meth.centerprofile.txt <==
CG	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan	-nan
CHG     ...
```

> Try step length - default as 0.01, for finer resolution to 0.001
*...Not working, it will shrink off the middle part of the gene... I guess it's smaller the size to 1/10 compared than before?*
> -s 0.1 doesn't work. Any changes on -s won't work!

methyGff -B -o test1exon_s0.1_7D107.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7days107_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7D107.methratio.txt -s 0.1 &

methyGff -B -o test1exon_s0.1_7DWT.meth -G /home/vcm/araport_reference/TAIR10_ref/TAIR10_chr_all.fas -b ~/BSseq_rep2_batmeth/genetic_component/7daysWT_exon_annotated.bed -m ~/BSseq_rep2_batmeth/genetic_component/methratio/calmeth_7DWT.methratio.txt -s 0.1 &


### 3.3 plot meth profile (landscape)
> not working, all goes to nan... maybe the issue of predefined region by BED???
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f test1exon1_7DWT.meth.TSSprofile.txt test1exon1_7D107.meth.TSSprofile.txt -l exon_7DWT exon_7D107 --outFileName test1exon1_plot_profile_mCG_TSSTES_7D.pdf -s 1 1 1 -xl up TSS TES down --context CG

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f test1exon10_7DWT.meth.TSSprofile.txt test1exon10_7D107.meth.TSSprofile.txt -l exon_7DWT exon_7D107 --outFileName test1exon10_plot_profile_mCG_TSSTES_7D.pdf -s 1 1 1 -xl up TSS TES down --context CG

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f test1exon100_7DWT.meth.TSSprofile.txt test1exon100_7D107.meth.TSSprofile.txt -l exon_7DWT exon_7D107 --outFileName test1exon100_plot_profile_mCG_TSSTES_7D.pdf -s 1 1 1 -xl up TSS TES down --context CG

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f test1exon1000_7DWT.meth.TSSprofile.txt test1exon1000_7D107.meth.TSSprofile.txt -l exon_7DWT exon_7D107 --outFileName test1exon1000_plot_profile_mCG_TSSTES_7D.pdf -s 1 1 1 -xl up TSS TES down --context CG


python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f testexon7DWT.meth.TSSprofile.txt testexon7D107.meth.TSSprofile.txt -l gene_7DWT gene_7D107 --outFileName testexon_plot_profile_CHG_TSSTES_7D.pdf -s 1 1 1 -xl up2k TSS TES down2k --context CHG &

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f testexon7DWT.meth.TSSprofile.txt testexon7D107.meth.TSSprofile.txt -l gene_7DWT gene_7D107 --outFileName testexon_plot_profile_CHH_TSSTES_7D.pdf -s 1 1 1 -xl up2k TSS TES down2k --context CHH &
```



> usage: bt2profile.py [-f MRFILE [MRFILE ...]] [-l LABEL [LABEL ...]] --outFileName
                     FILENAME [--sample SAMPLE [SAMPLE ...]] [-s SCALE [SCALE ...]]
                     [-xl XLABEL [XLABEL ...]] [-yl YLABEL] [-t TITLE [TITLE ...]]
                     [--yMin YMIN [YMIN ...]] [--yMax YMAX [YMAX ...]]
                     [--color COLOR [COLOR ...]]
                     [--legend {0,1,2,3,4,5,6,7,8,9,10,11,12}] [--lastlegend LASTLEGEND]
                     [--legendsize LEGENDSIZE] [--context CONTEXT] [--pergroup PERGROUP]
                     [-ft IMAGE_FORMAT] [--dpi DPI] [--help]

*-d 100 + centerprofile : works*
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f test1exon100_7DWT.meth.centerprofile.txt test1exon100_7D107.meth.centerprofile.txt \
-l exon_7DWT exon_7D107 --outFileName test1exon_d100_plot_profile_mCG_center_7D.pdf \
-s 1 1 -xl up center down --context CG &

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f test1exon100_7DWT.meth.centerprofile.txt test1exon100_7D107.meth.centerprofile.txt \
-l exon_7DWT exon_7D107 --outFileName test1exon_d100_plot_profile_CHG_center_7D.pdf \
-s 1 1 -xl up center down --context CHG &

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f test1exon100_7DWT.meth.centerprofile.txt test1exon100_7D107.meth.centerprofile.txt \
-l exon_7DWT exon_7D107 --outFileName test1exon_d100_plot_profile_CHH_center_7D.pdf \
-s 1 1 -xl up center down --context CHH &
```

*-s 0.001 -- does not work*
*-s 0.1 -- does not work*
*Just don't edit -s!!!!!*

#### plot center profile for genetic components!!!
> 7D
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f exon_7DWT.meth.centerprofile.txt exon_7D107.meth.centerprofile.txt \
-l exon_7DWT exon_7D107 --outFileName exon_plot_centerprofile_mCG_7D.pdf \
-s 1 1 -xl up center down --context CG &

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f exon_7DWT.meth.centerprofile.txt exon_7D107.meth.centerprofile.txt \
-l exon_7DWT exon_7D107 --outFileName exon_plot_centerprofile_CHG_7D.pdf \
-s 1 1 -xl up center down --context CHG &

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f exon_7DWT.meth.centerprofile.txt exon_7D107.meth.centerprofile.txt \
-l exon_7DWT exon_7D107 --outFileName exon_plot_centerprofile_CHH_7D.pdf \
-s 1 1 -xl up center down --context CHH &
```
> 3D
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f exon_3DWT.meth.centerprofile.txt exon_3D107.meth.centerprofile.txt \
-l exon_3DWT exon_7D107 --outFileName exon_plot_centerprofile_mCG_3D.pdf \
-s 1 1 -xl up center down --context CG &

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f exon_3DWT.meth.centerprofile.txt exon_3D107.meth.centerprofile.txt \
-l exon_3DWT exon_3D107 --outFileName exon_plot_centerprofile_CHG_3D.pdf \
-s 1 1 -xl up center down --context CHG &

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f exon_3DWT.meth.centerprofile.txt exon_3D107.meth.centerprofile.txt \
-l exon_3DWT exon_3D107 --outFileName exon_plot_centerprofile_CHH_3D.pdf \
-s 1 1 -xl up center down --context CHH &
```
> 0D
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f exon_0DWT.meth.centerprofile.txt exon_0D107.meth.centerprofile.txt \
-l exon_0DWT exon_0D107 --outFileName exon_plot_centerprofile_mCG_0D.pdf \
-s 1 1 -xl up center down --context CG &

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f exon_0DWT.meth.centerprofile.txt exon_0D107.meth.centerprofile.txt \
-l exon_0DWT exon_0D107 --outFileName exon_plot_centerprofile_CHG_0D.pdf \
-s 1 1 -xl up center down --context CHG &

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f exon_0DWT.meth.centerprofile.txt exon_0D107.meth.centerprofile.txt \
-l exon_0DWT exon_0D107 --outFileName exon_plot_centerprofile_CHH_0D.pdf \
-s 1 1 -xl up center down --context CHH &
```

#### Use this!!! AverMeth!!! for whole genome!!!
*try AverMethylevel...*
*Missing middle part; don't know what happen! It's the issue of exon file.*
*If we try the whole file like below, it works:) and it's BEAUTIFUL!!!!!!! Try this later!!!*
*Working*
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f methyGff_for-heatmap_7DWT.genomewide.gene.meth.AverMethylevel.txt methyGff_for-heatmap_7D107.genomewide.gene.meth.AverMethylevel.txt \
-l gene_7DWT gene_7D107 --outFileName testAverM_plot_profile_mCG_7D.pdf \
-s 1 1 -xl up center down --context CG &

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f methyGff_for-heatmap_7DWT.genomewide.gene.meth.AverMethylevel.txt methyGff_for-heatmap_7D107.genomewide.gene.meth.AverMethylevel.txt \
-l gene_7DWT gene_7D107 --outFileName testAverM_plot_profile_CHG_7D.pdf \
-s 1 1 -xl up center down --context CHG &

python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py \
-f methyGff_for-heatmap_7DWT.genomewide.gene.meth.AverMethylevel.txt methyGff_for-heatmap_7D107.genomewide.gene.meth.AverMethylevel.txt \
-l gene_7DWT gene_7D107 --outFileName testAverM_plot_profile_CHH_7D.pdf \
-s 1 1 -xl up center down --context CHH &
```

### 3.4 plot meth heatmap
change needed
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.cg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.cg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.chg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.chg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.chh.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.chh.txt -o bt2heatmap_0D107.genomewide.gene.meth.TSSTTS.cgchgchh.pdf --plotmatrix 3x2 --colorMap vlag --centerlabel center -z mCG CHG CHH --zMax 0.3 0.05 0.03 -t "0D107 TSS TTS" &
```


