# Step1: construct GTF files with only 1 type of genetic component

"path:"
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

"input file - Araport11.gtf"
```
Chr1    Araport11       transposable_element_gene       433031  433819  .       -       .       transcript_id "AT1G02228"; gene_id "AT1G02228";
Chr1    Araport11       transposable_element_gene       846664  847739  .       +       .       transcript_id "AT1G03420"; gene_id "AT1G03420";
Chr1    Araport11       transposable_element_gene       2415041 2415970 .       +       .       transcript_id "AT1G07800"; gene_id "AT1G07800";
Chr1    Araport11       transposable_element_gene       2531695 2534786 .       -       .       transcript_id "AT1G08105"; gene_id "AT1G08105";
Chr1    Araport11       transposable_element_gene       2790290 2793641 .       +       .       transcript_id "AT1G08735"; gene_id "AT1G08735";
...
```

"output file - protein_coding.Araport11.position"
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

"output file - exon.Araport11.position"
```
Chr1    3631    3913    AT1G01010       exon
Chr1    3996    4276    AT1G01010       exon
Chr1    4486    4605    AT1G01010       exon
Chr1    4706    5095    AT1G01010       exon
Chr1    5174    5326    AT1G01010       exon
...
```

"output file - five_prime_UTR.Araport11.position"
```
Chr1    3631    3759    AT1G01010       five_prime_UTR
Chr1    8667    9130    AT1G01020       five_prime_UTR
Chr1    8667    8737    AT1G01020       five_prime_UTR
Chr1    8443    8464    AT1G01020       five_prime_UTR
Chr1    8571    9130    AT1G01020       five_prime_UTR
Chr1    8443    8464    AT1G01020       five_prime_UTR
...
```

"output file - three_prime_UTR.Araport11.position"
```
Chr1    5631    5899    AT1G01010       three_prime_UTR
Chr1    6788    6914    AT1G01020       three_prime_UTR
Chr1    6788    7069    AT1G01020       three_prime_UTR
Chr1    7157    7314    AT1G01020       three_prime_UTR
Chr1    6788    6914    AT1G01020       three_prime_UTR
...
```

"Code"
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





