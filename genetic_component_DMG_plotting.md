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
├── [3.7M]  Araport11_transposon.gtf
├── [6.3M]  araport_reference
├── [ 11M]  exon.Araport11.position
├── [1.2M]  gene.Araport11.position
├── [   0]  intron.Araport11.position
└── [1.0M]  protein_coding.Araport11.position
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

> Trial sample
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
```
for file in *_unannotated.bed; do bedtools intersect -a $file -b /home/vcm/araport_reference/exon.Araport11.position -wa -wb | bedtools groupby -i - -g 1,2,3,4,8 -c 9 -o collapse > ${file%_unannotated.bed}_annotated.bed; echo ${file%_unannotated.bed}_annotated.bed; done &
```

# Step 3: run batmeth2 plotting

### 3.1 prepare batmeth2 input file
?

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

> Trial on 1 file
```
methyGff -B -o test7D107.meth -G 7days107-LFK10435_L1_1.fq.clean.gz 7days107-LFK10435_L1_2.fq.clean.gz -b 7days107_unannotated.bed -m calmeth_7D107.methratio.txt
```

