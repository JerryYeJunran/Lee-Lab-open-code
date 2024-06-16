## 0. preparation
### 0.1. download dependencies according to the webpage

```
https://batmeth2-docs.readthedocs.io/en/latest/index.html
```

#### 0.2. fastp
```
conda install fastp
```
#### 0.3. samtools
```
conda install samtools
```
#### 0.4. libpng

```
https://askubuntu.com/questions/895897/error-while-loading-shared-libraries-libpng12-so-0
```
**install autogen**
```
sudo apt-get update

sudo apt-get upgrade 

sudo apt-get install libtool autoconf build-essential pkg-config automake tcsh

wget http://archive.ubuntu.com/ubuntu/pool/main/libp/libpng/libpng_1.2.54.orig.tar.xz
```

**download libpng**
```
tar xvf  libpng_1.2.54.orig.tar.xz

cd libpng-1.2.54

./autogen.sh

./configure

make -j8 

sudo make install

sudo ldconfig
```

## 1. build index referemce genome
### 1.1 create path: ./batmeth_index and download ref_genome.fa

> download fasta file, all chromosome sequences, from TAIR10
```
cd ./batmeth

wget <"link to TAIR10 fasta file">
```

### 1.2 index reference genome

```
batMeth2 index -g genome_name.fa
```
> batmeth2 index -g TAIR10_chr_all.fas
* now genome path is: ./batmeth_index/TAIR10_chr_all.fas *
back to main path: cd ..

## 2. align gzipped raw reads

```
Paired-end:
batmeth2 align\
-g ./batmeth_index/TAIR10_chr_all.fas\
-1 Read1.fq.gz -2 Read2.fq.gz\
-o <prefix_name_added before output file name>\
-p <threads number using>
--fastp <path to fastp>
*fastp: /home/vcm/miniconda3/envs/batmeth/bin/fastp*
```
> batmeth2 align --fastp /home/vcm/miniconda3/envs/batmeth/bin/fastp -g ./batmeth2_index/TAIR10_chr_all.fas -1 7days107-LFK10435_L1_1.fq.gz -2 7days107-LFK10435_L1_2.fq.gz -o align_7D107 -p 8 &
>
> batmeth2 align --fastp /home/vcm/miniconda3/envs/batmeth/bin/fastp -g ./batmeth2_index/TAIR10_chr_all.fas -1 7daysWT-LFK10434_L1_1.fq.gz -2 7daysWT-LFK10434_L1_2.fq.gz -o align_7DWT -p 8 &
>
> batmeth2 align --fastp /home/vcm/miniconda3/envs/batmeth/bin/fastp -g ./batmeth2_index/TAIR10_chr_all.fas -1 3days107-LFK10433_L1_1.fq.gz -2 3days107-LFK10433_L1_2.fq.gz -o align_3D107 -p 8 &
>
> batmeth2 align --fastp /home/vcm/miniconda3/envs/batmeth/bin/fastp -g ./batmeth2_index/TAIR10_chr_all.fas -1 3daysWT-LFK10432_L1_1.fq.gz -2 3daysWT-LFK10432_L1_2.fq.gz -o align_3DWT -p 8 &
>
> batmeth2 align --fastp /home/vcm/miniconda3/envs/batmeth/bin/fastp -g ./batmeth2_index/TAIR10_chr_all.fas -1 S0day107-LFK10431_L1_1.fq.gz -2 S0day107-LFK10431_L1_2.fq.gz -o align_0D107 -p 8 &
>
> batmeth2 align --fastp /home/vcm/miniconda3/envs/batmeth/bin/fastp -g ./batmeth2_index/TAIR10_chr_all.fas -1 S0dayWT-LFK10430_L1_1.fq.gz -2 S0dayWT-LFK10430_L1_2.fq.gz -o align_0DWT -p 8 &
>
```
Or you can pause by pressing Ctrl+z, then type 'bg' in the shell.
```
## 3. Calculate DNA methylation level

> with bam file:
```
batmeth calmeth [options] -g genome.fa  -b alignment.sort.bam -m output.methrario.txt
```
> with sam file:
```
batmeth calmeth [options] -g genome.fa  -i alignment.sort.sam -m output.methrario.txt
```

> batmeth2 calmeth -g ./batmeth2_index/TAIR10_chr_all.fas -b ./align_7D107.sort.bam -m calmeth_7D107 &
>
> batmeth2 calmeth -g ./batmeth2_index/TAIR10_chr_all.fas -b ./align_7DWT.sort.bam -m calmeth_7DWT &
>
> batmeth2 calmeth -g ./batmeth2_index/TAIR10_chr_all.fas -b ./align_3D107.sort.bam -m calmeth_3D107 &
>
> batmeth2 calmeth -g ./batmeth2_index/TAIR10_chr_all.fas -b ./align_3DWT.sort.bam -m calmeth_3DWT &
>
> batmeth2 calmeth -g ./batmeth2_index/TAIR10_chr_all.fas -b ./align_0D107.sort.bam -m calmeth_0D107 &
>
> batmeth2 calmeth -g ./batmeth2_index/TAIR10_chr_all.fas -b ./align_0DWT.sort.bam -m calmeth_0DWT &

```
1. methratio
    Chromosome Loci Strand Context C_count CT_count methlevel eff_CT_count rev_G_count rev_GA_count MethContext 5context
    # ex. Chr1    61      +       CHH     3       11      0.286364        10.5    20      21      hU      ATCTT
    # C_count      The number of C in this base pair.
    # CT_count     The number of coverage in this base pair.
    # eff_CT_count Adjust read coverage based on opposite strand.
    # rev_G_count  The number of G in the reverse strand.
    # rev_GA_count The number of coverage in the reverse strand.
    # MethContext  M/Mh/H/hU/U, M means the methylation level ≥ 80%, etc
2. methBins
    Chrom BinIndex methlevel context
    # ex. Chr1    1       0.113674        CG
    # The BinIndex is defined by -s paramater in calmeth.
    # This file can be used for visualization the DNA methylation level acorss the chromosome.
3. Region
    chrom regionStart strand context c_count ct_count
    # ex. Chr1    1001    +       CG      1       227
    # The bins methylation level output file (BS.mr_Region.C*.txt) can be used to do DMR detection.
4. mCdensity
    CG/CHG/CHH C count in [0, 1%) [1%, 2%) ... [49%, 50%) ... [99%, 100%]
    # According to the DNA methylation level, the number of cytosine sites at different methylation levels was counted from 0 to 100.
5. mCcatero
    Average DNA methylation level including mC, mCG and other states.
```

## 4. Convert methratio file to bigwig for IGV visualization
#### Devided by 3 methylation types mCG, CHG, CHH, and total C
#### 4.1 index reference genome fasta file
```
samtools faidx ./batmeth2_index/TAIR10_chr_all.fas
```
#### 4.2 transform methratio.txt to .bw

```
python your_path_to_batmeth/batmeth2_to_bigwig.py -sort genome.fa.fai prefix.methratio.txt
```
*you can view your batmeth2 path by 'whereis batmeth2'*
*remember to delete the last path*
```
# my example
batmeth2: /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2
```

> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2_to_bigwig.py -sort ./batmeth2_index/TAIR10_chr_all.fas.fai calmeth_7D107.methratio.txt
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2_to_bigwig.py -sort ./batmeth2_index/TAIR10_chr_all.fas.fai calmeth_7DWT.methratio.txt
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2_to_bigwig.py -sort ./batmeth2_index/TAIR10_chr_all.fas.fai calmeth_3D107.methratio.txt
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2_to_bigwig.py -sort ./batmeth2_index/TAIR10_chr_all.fas.fai calmeth_3DWT.methratio.txt
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2_to_bigwig.py -sort ./batmeth2_index/TAIR10_chr_all.fas.fai calmeth_0D107.methratio.txt
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/batmeth2_to_bigwig.py -sort ./batmeth2_index/TAIR10_chr_all.fas.fai calmeth_0DWT.methratio.txt

#### 4.3 download bigwig file to own computer and load into IGV

> scp jerry@10.200.13.122:your_path_to_file path_on_your_own_computer

## 5. Diff methylation calling


[ Main paramaters ]
	-o_dm
	output file
	-o_dmr
	dmr output file when use auto detect by dmc
	-g|--genome
	Genome files
	-1
	sample1 methy files, sperate by space.
	-2
	sample2 methy files, sperate by space.
	-mindmc	
	min dmc sites in dmr region. [default : 4]
	-minstep
	min step in bp [default : 100]
	-maxdis
	max length of dmr [default : 0]
	-pvalue
	pvalue cutoff, default: 0.01
	-FDR
	adjust pvalue cutoff default : 1.0
	-methdiff
	the cutoff of methylation differention. default: 0.25 [CpG]
	-element
	caculate predefinded region, input file with id.
	-context
	Context for DM. [CG/CHG/CHH/ALL]
	-L
	predefinded regions or loci.
	-gz
	gzip input file.
	-h|--help


[Output file]
	1. DMC
	# format
	Chrom position starnd context pvalue adjust_pvalue combine_pvalue corrected_pvalue \
	cover_sample1 meth_sample1 cover_sample2 cover_sample2 meth.diff
	2. DMR
	# format
	Chrom start end methlevelInSample1 methlevelInSample2 NdmcInRegion hypermdc,hypodmc


#### 5.1 get dmc and dmr results

```
batDMR -g genome.fa -o_dm mutant.output.dmc -o_dmr mutant.output.dmr \
-1 mutant.methratio.txt -2 WT.methratio.txt \
-methdiff 0.2 -minstep 100 -mindmc 5 -pval 0.01
```

**1. Pre-definded regions (Gene/TE/UTR/CDS or other regions)**
```
BatMeth2 batDMR -g genome -L -o_dm dm.output.txt -1 [sample1.methC.txt replicates ..] \
-2 [sample2.methC.txt replicates ..]
```
**2. Auto define DMR region according the dmc**

***Recommended***
```
BatMeth2 batDMR -g genome -o_dm dm.output.txt -o_dmr dmr.output.txt -1 [sample1.methC.txt replicates ..] \
-2 [sample2.methC.txt replicates ..]
```

> batmeth2 batDMR -g ./batmeth2_index/TAIR10_chr_all.fas -o_dm batdmr_dm_7D.txt -o_dmr batdmr_dmr_7D.txt -1 calmeth_7DWT.methratio.txt -2 calmeth_7D107.methratio.txt &
>
> batmeth2 batDMR -g ./batmeth2_index/TAIR10_chr_all.fas -o_dm batdmr_dm_3D.txt -o_dmr batdmr_dmr_3D.txt -1 calmeth_0DWT.methratio.txt -2 calmeth_3D107.methratio.txt &
>
> batmeth2 batDMR -g ./batmeth2_index/TAIR10_chr_all.fas -o_dm batdmr_dm_0D.txt -o_dmr batdmr_dmr_0D.txt -1 calmeth_0DWT.methratio.txt -2 calmeth_0D107.methratio.txt &

#### 5.2 obtained hyper、hypo dmc/dmr from dmc/dmr results

```
awk -v OFS="\t" 'gsub(/\,/,"\t",$NF)' mutant.output.dmr | awk '$(NF-2)>4 && $NF<=1'  > mutant.output.hyper.dmr

awk -v OFS="\t" 'gsub(/\,/,"\t",$NF)' mutant.output.dmr | awk '!($(NF-2)>4 && $NF<=1)'  > mutant.output.hypo.dmr

awk '$NF>0' mutant.output.dmc | awk '{print $1"\t"$2"\t"$2}' > mutant.output.hyper.dmc

awk '$NF<0' mutant.output.dmc | awk '{print $1"\t"$2"\t"$2}' > mutant.output.hypo.dmc
```

## 6. Calulate mC across predefined regions


Example:
```
with gtf file:
    	methyGff -B -o gene.meth -G genome.fa -gtf gene.gtf -m output.methrario.txt

with multiple gtf file:
	methyGff -B -o expressed.gene.meth unexpressed.gene.meth \
        -G genome.fa  -gtf expressed.gene.gtf unexpressed.gene.gtf -m output.methrario.txt

with bed file:
	methyGff -B -o gene.meth -G genome.fa -b gene.bed -m output.methrario.txt

 with multiple bed file:
	methyGff -B -o expressed.gene.meth unexpressed.gene.meth \
        -G genome.fa -b expressed.gene.bed unexpressed.gene.bed -m output.methrario.txt
```
```
> methyGff [options] -o <OUT_PREFIX> -G GENOME -gff <GFF file>/-gtf <GTF file>/-b <bed file> -m <from Split methratio outfile> [-B][-P]
```

## 7. Plot diff meth with batmeth

#### 7.0 required installation
```
pip install numpy
pip install pandas
pip install matplotlib
pip install seaborn
```
#### 7.1 Plot methylation profile

> [option] -B -P --TSS --TTS --GENE

`-B`
```
methyGff -B -o methyGff_7DWT.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_7DWT.methratio.txt &

methyGff -B -o methyGff_7D107.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_7D107.methratio.txt &

methyGff -B -o methyGff_3DWT.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_3DWT.methratio.txt &

methyGff -B -o methyGff_3D107.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_3D107.methratio.txt &

methyGff -B -o methyGff_0DWT.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_0DWT.methratio.txt &

methyGff -B -o methyGff_0D107.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_0D107.methratio.txt &
```
`-P`
```
methyGff -P -o methyGff_promoter_7DWT.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_7DWT.methratio.txt &

methyGff -P -o methyGff_promoter_7D107.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_7D107.methratio.txt &

methyGff -P -o methyGff_promoter_3DWT.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_3DWT.methratio.txt &

methyGff -P -o methyGff_promoter_3D107.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_3D107.methratio.txt &

methyGff -P -o methyGff_promoter_0DWT.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_0DWT.methratio.txt &

methyGff -P -o methyGff_promoter_0D107.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gff batmeth2_index/TAIR10_GFF3_genes_transposons.gff -m calmeth_0D107.methratio.txt &
```

#### 7.2 Visualization

#### 7.2.1 plot profile - methylation line chart

```
#!# In file bt2profile.py, line 111:
Remove the quotation mark at line 45 around the number 45. 
```
**A. Gene Body (-B)**

**1. TSS, TES**

> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_7DWT.meth.TSSprofile.txt methyGff_7D107.meth.TSSprofile.txt -l gene_7DWT gene_7D107 --outFileName plot_profile_TSSTES_7D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_3DWT.meth.TSSprofile.txt methyGff_3D107.meth.TSSprofile.txt -l gene_3DWT gene_3D107 --outFileName plot_profile_TSSTES_3D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_0DWT.meth.TSSprofile.txt methyGff_0D107.meth.TSSprofile.txt -l gene_0DWT gene_0D107 --outFileName plot_profile_TSSTES_0D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
>

**2. Center**

> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_7DWT.meth.TSSprofile.txt methyGff_7D107.meth.TSSprofile.txt -l gene_7DWT gene_7D107 --outFileName plot_profile_center_7D.pdf -s 1 1 -xl up2k center down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_3DWT.meth.TSSprofile.txt methyGff_3D107.meth.TSSprofile.txt -l gene_3DWT gene_3D107 --outFileName plot_profile_center_3D.pdf -s 1 1 -xl up2k center down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_0DWT.meth.TSSprofile.txt methyGff_0D107.meth.TSSprofile.txt -l gene_0DWT gene_0D107 --outFileName plot_profile_center_0D.pdf -s 1 1 -xl up2k center down2k &
>

**B. Promoter**

**1. TSS, TES**

> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_promoter_7DWT.meth.TSSprofile.txt methyGff_7D107.meth.TSSprofile.txt -l promoter_7DWT promoter_7D107 --outFileName plot_promoter_profile_TSSTES_7D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_promoter_3DWT.meth.TSSprofile.txt methyGff_3D107.meth.TSSprofile.txt -l promoter_3DWT promoter_3D107 --outFileName plot_promoter_profile_TSSTES_3D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_promoter_0DWT.meth.TSSprofile.txt methyGff_0D107.meth.TSSprofile.txt -l promoter_0DWT promoter_0D107 --outFileName plot_promoter_profile_TSSTES_0D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
>

**2. Center**

> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_promoter_7DWT.meth.TSSprofile.txt methyGff_7D107.meth.TSSprofile.txt -l promoter_7DWT promoter_7D107 --outFileName plot_promoter_profile_center_7D.pdf -s 1 1 -xl up2k center down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_promoter_3DWT.meth.TSSprofile.txt methyGff_3D107.meth.TSSprofile.txt -l promoter_3DWT promoter_3D107 --outFileName plot_promoter_profile_center_3D.pdf -s 1 1 -xl up2k center down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_promoter_0DWT.meth.TSSprofile.txt methyGff_0D107.meth.TSSprofile.txt -l promoter_0DWT promoter_0D107 --outFileName plot_promoter_profile_center_0D.pdf -s 1 1 -xl up2k center down2k &
>

#### 7.2.2 Plot methylation profile compared to random gene list 

## (pending test)

$ BatMeth2 methyGff -o active random \
    -G genome.fa -m methratio.txt \
    -b active.bed random.bed -B

$ bt2profile.py -f active.centerprofile.txt \
    random.centerprofile.txt \
    -l active random \
    --outFileName active_random.output.meth.pdf \
    -s 1 1 -xl up2k center down2k

#### 7.2.3 Plot methylation profile compared to average gene list 

## (pending test)

$ bt2profile.py -f H3K27me3.bdgene.AverMethylevel.txt \
    H3K27me3.unbdgene.AverMethylevel.txt \
    -l H3K27me3.bdgene H3K27me3.unbdgene \
    --outFileName H3K27me3.output.meth.pdf \
    -s 1 1 1 -xl up2k TSS TES down2k

## 7.3 Plot bt2heatmap

#### 7.3.1 whole gene + up/down stream

$ python bt2heatmap.py -m H3K4me3.bdgene.GENE.cg.txt -l bg \
-o test0.pdf -z k43 -sl TSS -el TTS

#### 7.3.2 TSS & TES

$ python bt2heatmap.py -m H3K4me3.bdgene.TSS.cg.txt H3K4me3.bdgene.TTS.cg.txt \
    -l tss tts -o test.pdf --zMax 0.1 --colorMap vlag --centerlabel center -z bd

#### 7.3.3 TSS & TES, with distinguishing mCG, CHG, CHH methylation

$ python bt2heatmap.py -m H3K4me3.bdgene.TSS.cg.txt H3K4me3.bdgene.TTS.cg.txt \
    H3K4me3.bdgene.TSS.chg.txt H3K4me3.bdgene.TTS.chg.txt \
    H3K4me3.bdgene.TSS.chh.txt H3K4me3.bdgene.TTS.chh.txt \
    -l H3K4me3.bdgene-tss H3K4me3.bdgene-tts \
    -o H3K4me3.bdgene.TSS_TTS.heatmap.pdf --plotmatrix 3x2 \
    --centerlabel center -z cg chg chh --zMax 0.3 1 0.01

# ADDING (1)
## 7. Plot diff meth with batmeth
### Using DEG gene list gtf file to see if DE genes have distinct methylation pattern from all gene:

> DE gene list GTF path:
>
> /data/rep2_bismark_yao/Batmeth2/DEG
>
> 0D_log2genetable_GTF.tsv  3D_log2genetable_GTF.tsv  7D_log2genetable_GTF.tsv

> log2genetable_GTF path:
>
> ~/annotate/0D_log2genetable_GTF.tsv  3D_log2genetable_GTF.tsv  7D_log2genetable_GTF.tsv

> with gtf file:
>
> target path: /data/rep2_bismark_yao/Batmeth2
```
methyGff -B -o methyGff_7DWT_DEG.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/annotate/7D_log2genetable_GTF.tsv -m calmeth_7DWT.methratio.txt &
methyGff -B -o methyGff_7D107_DEG.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/annotate/7D_log2genetable_GTF.tsv -m calmeth_7D107.methratio.txt &
methyGff -B -o methyGff_3DWT_DEG.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/annotate/7D_log2genetable_GTF.tsv -m calmeth_3DWT.methratio.txt &
methyGff -B -o methyGff_3D107_DEG.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/annotate/7D_log2genetable_GTF.tsv -m calmeth_3D107.methratio.txt &
methyGff -B -o methyGff_0DWT_DEG.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/annotate/7D_log2genetable_GTF.tsv -m calmeth_0DWT.methratio.txt &
methyGff -B -o methyGff_0D107_DEG.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/annotate/7D_log2genetable_GTF.tsv -m calmeth_0D107.methratio.txt &
```
#### 7.2.1 plot profile - methylation line chart
**A. Gene Body (-B)**

**1. TSS, TES**

> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_7DWT_DEG.gene.meth.TSSprofile.txt methyGff_7D107_DEG.gene.meth.TSSprofile.txt -l DEgene_7DWT DEgene_7D107 --outFileName plot_profile_DEG_TSSTES_7D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_3DWT_DEG.gene.meth.TSSprofile.txt methyGff_3D107_DEG.gene.meth.TSSprofile.txt -l DEgene_3DWT DEgene_3D107 --outFileName plot_profile_DEG_TSSTES_3D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_0DWT_DEG.gene.meth.TSSprofile.txt methyGff_0D107_DEG.gene.meth.TSSprofile.txt -l DEgene_0DWT DEgene_0D107 --outFileName plot_profile_DEG_TSSTES_0D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
>

**2. Center**

> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_7DWT_DEG.gene.meth.TSSprofile.txt methyGff_7D107_DEG.gene.meth.TSSprofile.txt -l DEgene_7DWT DEgene_7D107 --outFileName plot_profile_DEG_center_7D.pdf -s 1 1 -xl up2k center down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_3DWT_DEG.gene.meth.TSSprofile.txt methyGff_3D107_DEG.gene.meth.TSSprofile.txt -l DEgene_3DWT DEgene_3D107 --outFileName plot_profile_DEG_center_3D.pdf -s 1 1 -xl up2k center down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_0DWT_DEG.gene.meth.TSSprofile.txt methyGff_0D107_DEG.gene.meth.TSSprofile.txt -l DEgene_0DWT DEgene_0D107 --outFileName plot_profile_DEG_center_0D.pdf -s 1 1 -xl up2k center down2k &
>

**Let's compare DEG with all_gene in one plot together!**

**A. Gene Body (-B)**

**1. TSS, TES**

> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_7DWT_DEG.gene.meth.TSSprofile.txt methyGff_7D107_DEG.gene.meth.TSSprofile.txt methyGff_7DWT.meth.TSSprofile.txt methyGff_7D107.meth.TSSprofile.txt -l DEgene_7DWT DEgene_7D107 gene_7DWT gene_7D107 --outFileName plot_profile_DEG_ALL_TSSTES_7D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_3DWT_DEG.gene.meth.TSSprofile.txt methyGff_3D107_DEG.gene.meth.TSSprofile.txt methyGff_3DWT.meth.TSSprofile.txt methyGff_3D107.meth.TSSprofile.txt -l DEgene_3DWT DEgene_3D107 gene_3DWT gene_3D107 --outFileName plot_profile_DEG_ALL_TSSTES_3D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
>
> python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_0DWT_DEG.gene.meth.TSSprofile.txt methyGff_0D107_DEG.gene.meth.TSSprofile.txt methyGff_0DWT.meth.TSSprofile.txt methyGff_0D107.meth.TSSprofile.txt -l DEgene_0DWT DEgene_0D107 gene_0DWT gene_0D107 --outFileName plot_profile_DEG_ALL_TSSTES_0D.pdf -s 1 1 1 -xl up2k TSS TES down2k &
> 


