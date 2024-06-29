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

#### 0.5 gsl-2.4
download_path = /home/vcm/Batmeth2_download/gsl/gsl-2.4

```
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/home/vcm/Batmeth2_download/BatMeth2/gsl-2.4/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/home/vcm/Batmeth2_download/BatMeth2/gsl-2.4/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH::/home/vcm/Batmeth2_download/BatMeth2/gsl-2.4/lib
export LIBRARY_PATH=$LIBRARY_PATH::/home/vcm/Batmeth2_download/BatMeth2/gsl-2.4/lib

source ~/.bashrc
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
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f methyGff_7DWT.meth.TSSprofile.txt methyGff_7D107.meth.TSSprofile.txt -l gene_7DWT gene_7D107 --outFileName plot_profile_mCG_TSSTES_7D.pdf -s 1 1 1 -xl up2k TSS TES down2k --context CG &
```
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f methyGff_7DWT.meth.TSSprofile.txt methyGff_7D107.meth.TSSprofile.txt -l gene_7DWT gene_7D107 --outFileName plot_profile_CHG_TSSTES_7D.pdf -s 1 1 1 -xl up2k TSS TES down2k --context CHG &
```
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f methyGff_7DWT.meth.TSSprofile.txt methyGff_7D107.meth.TSSprofile.txt -l gene_7DWT gene_7D107 --outFileName plot_profile_CHH_TSSTES_7D.pdf -s 1 1 1 -xl up2k TSS TES down2k --context CHH &
```

```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f methyGff_3DWT.meth.TSSprofile.txt methyGff_3D107.meth.TSSprofile.txt -l gene_3DWT gene_3D107 --outFileName plot_profile_mCG_TSSTES_3D.pdf -s 1 1 1 -xl up2k TSS TES down2k --context CG &
```
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f methyGff_3DWT.meth.TSSprofile.txt methyGff_3D107.meth.TSSprofile.txt -l gene_3DWT gene_3D107 --outFileName plot_profile_CHG_TSSTES_3D.pdf -s 1 1 1 -xl up2k TSS TES down2k --context CHG &
```
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f methyGff_3DWT.meth.TSSprofile.txt methyGff_3D107.meth.TSSprofile.txt -l gene_3DWT gene_3D107 --outFileName plot_profile_CHH_TSSTES_3D.pdf -s 1 1 1 -xl up2k TSS TES down2k --context CHH &
```

```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f methyGff_0DWT.meth.TSSprofile.txt methyGff_0D107.meth.TSSprofile.txt -l gene_0DWT gene_0D107 --outFileName plot_profile_mCG_TSSTES_0D.pdf -s 1 1 1 -xl up2k TSS TES down2k --context CG &
```
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f methyGff_0DWT.meth.TSSprofile.txt methyGff_0D107.meth.TSSprofile.txt -l gene_0DWT gene_0D107 --outFileName plot_profile_CHG_TSSTES_0D.pdf -s 1 1 1 -xl up2k TSS TES down2k --context CHG &
```
```
python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2profile.py -f methyGff_0DWT.meth.TSSprofile.txt methyGff_0D107.meth.TSSprofile.txt -l gene_0DWT gene_0D107 --outFileName plot_profile_CHH_TSSTES_0D.pdf -s 1 1 1 -xl up2k TSS TES down2k --context CHH &
```

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

#### 7.2.2 Ploting methylation profile with DEG list
usage: bt2profile.py [-f MRFILE [MRFILE ...]] [-l LABEL [LABEL ...]] --outFileName FILENAME [--sample SAMPLE [SAMPLE ...]] [-s SCALE [SCALE ...]] [-xl XLABEL [XLABEL ...]] [-yl YLABEL]
                     [-t TITLE [TITLE ...]] [--yMin YMIN [YMIN ...]] [--yMax YMAX [YMAX ...]] [--color COLOR [COLOR ...]] [--legend {0,1,2,3,4,5,6,7,8,9,10,11,12}] [--lastlegend LASTLEGEND]
                     [--legendsize LEGENDSIZE] [--context CONTEXT] [--pergroup PERGROUP] [-ft IMAGE_FORMAT] [--dpi DPI] [--help]

options:
  -f MRFILE [MRFILE ...], --mrfile MRFILE [MRFILE ...]
                        DNA AverMethylevel files, seperate by space. eg. wildtype.AverMethy.txt
  -l LABEL [LABEL ...], --label LABEL [LABEL ...]
                        Labels of samples, sperate by space. eg. -l widetype
  --outFileName FILENAME, -o FILENAME
                        Output file name.
  --sample SAMPLE [SAMPLE ...]
                        The interval of N data is a group, and the average value is taken as the representative. Please note that the number of labels should correspondto the number of
                        samples after averaging.
  -s SCALE [SCALE ...], --scale SCALE [SCALE ...]
                        Visual X-axis spacing, default upsteam:body:downstream is 1:1:1 (-s 1 1 1), which should be consistent with -b and -bl parameters in BatMeth2:methyGff,and separated
                        by spaces
  -xl XLABEL [XLABEL ...], --xlabel XLABEL [XLABEL ...]
                        Consistent with the -s parameter, if the -s parameter is 1 1 1, i.e. 1:1:1,then the corresponding X-axis label is UP TSS TES Down
  -yl YLABEL, --ylabel YLABEL
                        y-axis label
  -t TITLE [TITLE ...], --title TITLE [TITLE ...]
                        Title of the plot, to be printed on top of the generated image. Leave blank for no title.
  --yMin YMIN [YMIN ...]
                        Minimum value for the Y-axis. Multiple values, separated by spaces can be set for each profile. If the number of yMin values is smaller thanthe number of plots, the
                        values are recycled.
  --yMax YMAX [YMAX ...]
                        Maximum value for the Y-axis. Multiple values, separated by spaces can be set for each profile. If the number of yMin values is smaller thanthe number of plots, the
                        values are recycled.
  --color COLOR [COLOR ...]
                        List of colors to use, should same as the number of samples,Color names and html hex strings (e.g., #eeff22) are accepted. The color names should be space separated.
                        For example, --color red blue green
  --legend {0,1,2,3,4,5,6,7,8,9,10,11,12}
                        The location of the legend. best : 0, upper right : 1, upper left : 2, lower left : 3, lower right : 4, right : 5, center left : 6, center right: 7, lower center: 8,
                        upper center: 9, center : 10, out : 11, none : 12
  --lastlegend LASTLEGEND
                        Only show the last figure's legend.
  --legendsize LEGENDSIZE
                        the text size of the legend.
  --context CONTEXT     
  			choices=["ALL", "C", "CG", "CHG","CHH"],
  --color COLOR
  			List of colors to use, should same as the number of samples,Color names and html hex strings (e.g., #eeff22) are accepted. The color names should be space separated.
                        For example, --color red blue green
  --pergroup PERGROUP   plot cg/ch of the same sample in one fig,only useful when have more than 1 sample input file
  -ft IMAGE_FORMAT, --image_format IMAGE_FORMAT
                        The file format, e.g. 'png', 'pdf', 'svg', ... The behavior when this is unset is documented under fname.
  --dpi DPI             Set the DPI to save the figure. default: 200

```
methyGff -B -o methyGff_7DWT_DEG.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/annotate/7D_log2genetable_GTF.tsv -m calmeth_7DWT.methratio.txt
```

```
python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2profile.py -f methyGff_3DWT_DEG.gene.meth.TSSprofile.txt methyGff_3D107_DEG.gene.meth.TSSprofile.txt -l DEgene_3DWT DEgene_3D107 --outFileName plot_profile_DEG_center_3D.pdf -s 1 1 -xl up2k center down2k &
```
#### 7.2.3 Plot methylation profile compared to random gene list 

> ## (pending test)
> 
> $ BatMeth2 methyGff -o active random \
>     -G genome.fa -m methratio.txt \
>     -b active.bed random.bed -B
> 
> *where is the random bed file? How can I create such file?*
> BatMeth2 methyGff -o active random -G ./batmeth2_index/TAIR10_chr_all.fas -m calmeth_7DWT.methratio.txt -b methyGff_7DWT_DEG.bed Pb_Random.bed -B &
> 
> Processing 1 out of 2. InFile: methyGff_7DWT_DEG.bed
> File methyGff_7DWT_DEG.bed Cannot be opened ....
> 
> $ bt2profile.py -f active.centerprofile.txt \
>     random.centerprofile.txt \
>     -l active random \
>     --outFileName active_random.output.meth.pdf \
>     -s 1 1 -xl up2k center down2k

> #### 7.2.3 Plot methylation profile compared to average gene list 
> 
> ## (pending test)
> 
> $ bt2profile.py -f H3K27me3.bdgene.AverMethylevel.txt \
>     H3K27me3.unbdgene.AverMethylevel.txt \
>     -l H3K27me3.bdgene H3K27me3.unbdgene \
>     --outFileName H3K27me3.output.meth.pdf \
>     -s 1 1 1 -xl up2k TSS TES down2k

## 7.3 Plot bt2heatmap

**options:**
  -f MRFILE [MRFILE ...], --mrfile MRFILE [MRFILE ...]
                        input methylevel files, wildtype.body.c*.txt
  -m MATRIXFILE [MATRIXFILE ...], --matrixfile MATRIXFILE [MATRIXFILE ...]
                        input methylevel matrix files, wildtype.GENE.cg.txt
  -l SAMPLESLABEL [SAMPLESLABEL ...], --samplesLabel SAMPLESLABEL [SAMPLESLABEL ...]
                        the label of the samples
  -z GROUPLABELS [GROUPLABELS ...], --groupLabels GROUPLABELS [GROUPLABELS ...]
                        Labels for the regions plotted in the heatmap. If more than one region is being plotted, a list of labels separated by spaces is required. If a label itself contains a
                        space, then quotes are needed. For example, --groupLabels label_1, "label 2".
  -sl STARTLABEL, --startlabel STARTLABEL
                        the start label of the samples
  -el ENDLABEL, --endlabel ENDLABEL
                        the end label of the samples
  -pl CENTERLABEL, --centerlabel CENTERLABEL
                        the center label of the samples
  --plotmatrix PLOTMATRIX
                        1x1, default, row x col, order by columun, for exsample, 2x3 :file1 file2 file3file4 file5 file6
  --outFileName FILENAME, -o FILENAME
                        Output file name.
  -c COLORMAP [COLORMAP ...], --colorMap COLORMAP [COLORMAP ...]
                        Color map to use for the heatmap. If more than one heatmap is being plotted the color of each heatmap can be enter individually (e.g. `--colorMap Reds Blues`).The
                        available options are: 'magma', 'inferno', 'plasma', 'viridis', 'cividis', 'twilight', 'twilight_shifted', 'turbo', 'Blues', 'BrBG', 'BuGn', 'BuPu', 'CMRmap', 'GnBu',
                        'Greens', 'Greys', 'OrRd', 'Oranges', 'PRGn', 'PiYG', 'PuBu', 'PuBuGn', 'PuOr', 'PuRd', 'Purples', 'RdBu', 'RdGy', 'RdPu', 'RdYlBu', 'RdYlGn', 'Reds', 'Spectral',
                        'Wistia', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'afmhot', 'autumn', 'binary', 'bone', 'brg', 'bwr', 'cool', 'coolwarm', 'copper', 'cubehelix', 'flag', 'gist_earth',
                        'gist_gray', 'gist_heat', 'gist_ncar', 'gist_rainbow', 'gist_stern', 'gist_yarg', 'gnuplot', 'gnuplot2', 'gray', 'hot', 'hsv', 'jet', 'nipy_spectral', 'ocean', 'pink',
                        'prism', 'rainbow', 'seismic', 'spring', 'summer', 'terrain', 'winter', 'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3', 'tab10', 'tab20',
                        'tab20b', 'tab20c', 'grey', 'gist_grey', 'gist_yerg', 'Grays', 'rocket', 'mako', 'icefire', 'vlag', 'flare', 'crest'
  --alpha ALPHA         The alpha channel (transparency) to use for the heatmaps. The default is 1.0 and values must be between 0 and 1.
  --colorList COLORLIST [COLORLIST ...]
                        List of colors to use to create a colormap. For example, if `--colorList black,yellow,blue` is set (colors separated by comas) then a color map that starts with black,
                        continues to yellow and finishes in blue is created. If this option is selected, it overrides the --colorMap chosen. The list of valid color names can be seen here:
                        http://matplotlib.org/examples/color/named_colors.html The number of transitions is defined by the --colorNumber option.
  --colorNumber COLORNUMBER
                        --colorList is required for an effect. This controls the number of transitions from one color to the other. If --colorNumber is the number of colors in --colorList
                        then there will be no transitions between the colors.
  --missingDataColor MISSINGDATACOLOR
                        If --missingDataAsZero was not set, such cases will be colored in white by default.
  --sortRegions {descend,ascend,no}
                        Whether the heatmap should present the regions sorted. The default is to sort in descending order based on the mean value per region.
  --sortUsing {mean,median,max,min,sum}
                        Indicate which method should be used for sorting. For each row the method is computed.
  --sortUsingSamples SORTUSINGSAMPLES [SORTUSINGSAMPLES ...]
                        List of sample numbers (order as in matrix), which are used by --sortUsing for sorting. If no value is set, it uses all samples. Example: --sortUsingSamples 1 3
  --linesAtTickMarks    Draw dashed lines from all tick marks through the heatmap. This is then similar to the dashed line draw at region bounds when using a reference point and --sortUsing
                        region_length
  --clusterUsingSamples CLUSTERUSINGSAMPLES [CLUSTERUSINGSAMPLES ...]
                        List of sample numbers (order as in matrix), that are used for clustering by --kmeans or --hclust if not given, all samples are taken into account for clustering.
                        Example: --ClusterUsingSamples 1 3
  --kmeans KMEANS       Number of clusters to compute. When this option is set, the matrix is split into clusters using the k-means algorithm. Only works for data that is not grouped,
                        otherwise only the first group will be clustered.
  --hclust HCLUST       Number of clusters to compute. When this option is set, then the matrix is split into clusters using the hierarchical clustering algorithm, using "ward linkage". Only
                        works for data that is not grouped, otherwise only the first group will be clustered. --hclust could be very slow if you have >1000 regions. In those cases, you might
                        prefer --kmeans or if more clustering methods are required you can save the underlying matrix and run the clustering using other software. The plotting of the
                        clustering may fail with an error if a cluster has very few members compared to the total number of regions.
  -s SCALE [SCALE ...], --scale SCALE [SCALE ...]
                        Maximum value for the Y-axis. Multiple values, separated by spaces can be set for each profile. If the number of yMin values is smaller thanthe number of plots, the
                        values are recycled.
  -t TITLE [TITLE ...], --title TITLE [TITLE ...]
                        Title of the plot, to be printed on top of the generated image. Leave blank for no title.
  --zMin ZMIN [ZMIN ...]
                        Values to anchor the colormap
  --zMax ZMAX [ZMAX ...]
                        Values to anchor the colormap, Maximum value for the heatmap.
  -ft IMAGE_FORMAT, --image_format IMAGE_FORMAT
                        The file format, e.g. 'png', 'pdf', 'svg', ... The behavior when this is unset is documented under fname.
  --perGroup            The default is to plot all groups of regions by sample. Using this option instead plots all samples by group of regions. Note that this is only useful if you have
                        multiple groups of regions. by sample rather than group.
  --dpi DPI             Set the DPI to save the figure. default: 100
  --figsize FIGSIZE     Set the figure size to save the figure. [with]x[height], default: 1.5x11
  --boxAroundHeatmaps BOXAROUNDHEATMAPS
                        By default black boxes are plot around heatmaps. This can be turned off by setting --boxAroundHeatmaps no
  --help, -h            show this help message and exit


#### 7.3.1 whole gene + up/down stream

Inout Files:
methyGff_*D*.meth.body.cg.txt
methyGff_*D*.meth.body.chg.txt
methyGff_*D*.meth.body.chh.txt
Refering to:
bt2heatmap_methyGff_0D107_DEG.gene_body_CHH.pdf
bt2heatmap_methyGff_0D107_body_CHH.pdf
bt2heatmap_methyGff_0DWT_DEG.gene_body_CHH.pdf
bt2heatmap_methyGff_0DWT_body_CHH.pdf
bt2heatmap_methyGff_3D107_DEG.gene_body_CHH.pdf
bt2heatmap_methyGff_3D107_body_CHH.pdf
bt2heatmap_methyGff_3DWT_DEG.gene_body_CHH.pdf
bt2heatmap_methyGff_3DWT_body_CHH.pdf
bt2heatmap_methyGff_7D107_DEG.gene_body_CHH.pdf
bt2heatmap_methyGff_7D107_body_CHH.pdf
bt2heatmap_methyGff_7DWT_DEG.gene_body_CHH.pdf
bt2heatmap_methyGff_7DWT_body_CHH.pdf

####
# WORKING!!!!
####

#### 7.3.0

ref_gtf path = /home/vcm/annotate/Araport11.gtf
`for-heatmap: used to run bt2heatmap`
`genome wide: using the whole genome gtf ref_genome`
`gene: gene region, other than promoter region`

methyGff --TSS --TTS --GENE -o methyGff_for-heatmap_0DWT.genomewide.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/araport_reference/Araport11.gtf -m calmeth_0DWT.methratio.txt &
methyGff --TSS --TTS --GENE -o methyGff_for-heatmap_3DWT.genomewide.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/araport_reference/Araport11.gtf -m calmeth_3DWT.methratio.txt &
methyGff --TSS --TTS --GENE -o methyGff_for-heatmap_7DWT.genomewide.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/araport_reference/Araport11.gtf -m calmeth_7DWT.methratio.txt &
methyGff --TSS --TTS --GENE -o methyGff_for-heatmap_0D107.genomewide.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/araport_reference/Araport11.gtf -m calmeth_0D107.methratio.txt &
methyGff --TSS --TTS --GENE -o methyGff_for-heatmap_3D107.genomewide.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/araport_reference/Araport11.gtf -m calmeth_3D107.methratio.txt &
methyGff --TSS --TTS --GENE -o methyGff_for-heatmap_7D107.genomewide.gene.meth -G ./batmeth2_index/TAIR10_chr_all.fas -gtf /home/vcm/araport_reference/Araport11.gtf -m calmeth_7D107.methratio.txt &


#### 7.3.1 bt2heatmap

python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m test_methyGff_7DWT.gene.meth.TSS.cg.txt test_methyGff_7DWT.gene.meth.TTS.cg.txt -l TSS TTS -o test_methyGff_7DWT.gene.meth.TSSTTS.cg.pdf --colorMap vlag --centerlabel center -z mCG

`-TSS + -TTS` `mCG`
python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.cg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.cg.txt -l TSS TTS -o methyGff_0D107.genomewide.gene.meth.TSSTTS.cg.pdf --colorMap vlag --centerlabel center -z mCG &

python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.chg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.chg.txt -l TSS TTS -o methyGff_0D107.genomewide.gene.meth.TSSTTS.chg.pdf --colorMap vlag --centerlabel center -z CHG &

python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.chh.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.chh.txt -l TSS TTS -o methyGff_0D107.genomewide.gene.meth.TSSTTS.chh.pdf --colorMap vlag --centerlabel center -z CHH &

python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_3D107.genomewide.gene.meth.TSS.cg.txt methyGff_for-heatmap_3D107.genomewide.gene.meth.TTS.cg.txt -l TSS TTS -o methyGff_3D107.genomewide.gene.meth.TSSTTS.cg.pdf --colorMap vlag --centerlabel center -z mCG -l TSS TTS -t 3D107 --zMax 0.4 &

python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_3D107.genomewide.gene.meth.TSS.chg.txt methyGff_for-heatmap_3D107.genomewide.gene.meth.TTS.chg.txt -l TSS TTS -o methyGff_3D107.genomewide.gene.meth.TSSTTS.chg.pdf --colorMap vlag --centerlabel center -z CHG -l TSS TTS -t 3D107 --zMax 0.3 &

python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_3D107.genomewide.gene.meth.TSS.chh.txt methyGff_for-heatmap_3D107.genomewide.gene.meth.TTS.chh.txt -l TSS TTS -o methyGff_3D107.genomewide.gene.meth.TSSTTS.chh.pdf --colorMap vlag --centerlabel center -z CHH -l TSS TTS -t 3D107 --zMax 0.1 &

?????

#### ????? Need to work in real server for more memory usage!
 python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.cg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.cg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.chg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.chg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.chh.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.chh.txt -o bt2heatmap_0D107.genomewide.gene.meth.TSSTTS.cgchgchh.pdf --plotmatrix 3x2 --colorMap vlag --centerlabel center -z mCG CHG CHH --zMax 0.3 0.05 0.03 -t "0D107 TSS TTS" &

  python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_3D107.genomewide.gene.meth.TSS.cg.txt methyGff_for-heatmap_3D107.genomewide.gene.meth.TTS.cg.txt methyGff_for-heatmap_3D107.genomewide.gene.meth.TSS.chg.txt methyGff_for-heatmap_3D107.genomewide.gene.meth.TTS.chg.txt methyGff_for-heatmap_3D107.genomewide.gene.meth.TSS.chh.txt methyGff_for-heatmap_3D107.genomewide.gene.meth.TTS.chh.txt -o bt2heatmap_3D107.genomewide.gene.meth.TSSTTS.cgchgchh.pdf --plotmatrix 3x2 --colorMap vlag --centerlabel center -z mCG CHG CHH --zMax 0.3 0.05 0.03 -t "0D107 TSS TTS" &

  python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_7D107.genomewide.gene.meth.TSS.cg.txt methyGff_for-heatmap_7D107.genomewide.gene.meth.TTS.cg.txt methyGff_for-heatmap_7D107.genomewide.gene.meth.TSS.chg.txt methyGff_for-heatmap_7D107.genomewide.gene.meth.TTS.chg.txt methyGff_for-heatmap_7D107.genomewide.gene.meth.TSS.chh.txt methyGff_for-heatmap_7D107.genomewide.gene.meth.TTS.chh.txt -o bt2heatmap_7D107.genomewide.gene.meth.TSSTTS.cgchgchh.pdf --plotmatrix 3x2 --colorMap vlag --centerlabel center -z mCG CHG CHH --zMax 0.3 0.05 0.03 -t "0D107 TSS TTS" &

 python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_0DWT.genomewide.gene.meth.TSS.cg.txt methyGff_for-heatmap_0DWT.genomewide.gene.meth.TTS.cg.txt methyGff_for-heatmap_0DWT.genomewide.gene.meth.TSS.chg.txt methyGff_for-heatmap_0DWT.genomewide.gene.meth.TTS.chg.txt methyGff_for-heatmap_0DWT.genomewide.gene.meth.TSS.chh.txt methyGff_for-heatmap_0DWT.genomewide.gene.meth.TTS.chh.txt -o bt2heatmap_0DWT.genomewide.gene.meth.TSSTTS.cgchgchh.pdf --plotmatrix 3x2 --colorMap vlag --centerlabel center -z mCG CHG CHH --zMax 0.3 0.05 0.03 -t "0D107 TSS TTS" &

  python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_3DWT.genomewide.gene.meth.TSS.cg.txt methyGff_for-heatmap_3DWT.genomewide.gene.meth.TTS.cg.txt methyGff_for-heatmap_3DWT.genomewide.gene.meth.TSS.chg.txt methyGff_for-heatmap_3DWT.genomewide.gene.meth.TTS.chg.txt methyGff_for-heatmap_3DWT.genomewide.gene.meth.TSS.chh.txt methyGff_for-heatmap_3DWT.genomewide.gene.meth.TTS.chh.txt -o bt2heatmap_3DWT.genomewide.gene.meth.TSSTTS.cgchgchh.pdf --plotmatrix 3x2 --colorMap vlag --centerlabel center -z mCG CHG CHH --zMax 0.3 0.05 0.03 -t "0D107 TSS TTS" &

  python /home/vcm/Batmeth2_download/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_7DWT.genomewide.gene.meth.TSS.cg.txt methyGff_for-heatmap_7DWT.genomewide.gene.meth.TTS.cg.txt methyGff_for-heatmap_7DWT.genomewide.gene.meth.TSS.chg.txt methyGff_for-heatmap_7DWT.genomewide.gene.meth.TTS.chg.txt methyGff_for-heatmap_7DWT.genomewide.gene.meth.TSS.chh.txt methyGff_for-heatmap_7DWT.genomewide.gene.meth.TTS.chh.txt -o bt2heatmap_7DWT.genomewide.gene.meth.TSSTTS.cgchgchh.pdf --plotmatrix 3x2 --colorMap vlag --centerlabel center -z cg chg chh -z mCG CHG CHH --zMax 0.3 0.05 0.03 -t "0D107 TSS TTS" &

?????
python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.cg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.cg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.chg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.chg.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TSS.chh.txt methyGff_for-heatmap_0D107.genomewide.gene.meth.TTS.chh.txt -l TSS TTS TSS TTS TSS TTS -o methyGff_0D107.genomewide.gene.meth.TSSTTS.cgchgchh.pdf --colorMap vlag --centerlabel center &

python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m \
methyGff_for-heatmap_3D107.genomewide.gene.meth.TSS.cg.txt \
methyGff_for-heatmap_3D107.genomewide.gene.meth.TTS.cg.txt \
methyGff_for-heatmap_3D107.genomewide.gene.meth.TSS.chg.txt \
methyGff_for-heatmap_3D107.genomewide.gene.meth.TTS.chg.txt \
methyGff_for-heatmap_3D107.genomewide.gene.meth.TSS.chh.txt \
methyGff_for-heatmap_3D107.genomewide.gene.meth.TTS.chh.txt \
-l mCG_TSS mCG_TTS CHG_TSS CHG_TTS CHH_TSS CHH_TTS \
-o test_methyGff_3D107.genomewide.gene.meth.TSSTTS.cgchgchh.pdf \
--plotmatrix 3x2 --colorMap vlag --centerlabel center \
-t whole_genome_3D107 \
-z mCG CHG CHH &

python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m \
methyGff_for-heatmap_7D107.genomewide.gene.meth.TSS.cg.txt \
methyGff_for-heatmap_7D107.genomewide.gene.meth.TTS.cg.txt \
methyGff_for-heatmap_7D107.genomewide.gene.meth.TSS.chg.txt \
methyGff_for-heatmap_7D107.genomewide.gene.meth.TTS.chg.txt \
methyGff_for-heatmap_7D107.genomewide.gene.meth.TSS.chh.txt \
methyGff_for-heatmap_7D107.genomewide.gene.meth.TTS.chh.txt \
-l mCG_TSS mCG_TTS CHG_TSS CHG_TTS CHH_TSS CHH_TTS \
-o test_methyGff_7D107.genomewide.gene.meth.TSSTTS.cgchgchh.pdf \
--plotmatrix 3x2 --colorMap vlag --centerlabel center \
-t whole_genome_7D107 \
-z mCG CHG CHH &

python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m \
methyGff_for-heatmap_0DWT.genomewide.gene.meth.TSS.cg.txt \
methyGff_for-heatmap_0DWT.genomewide.gene.meth.TTS.cg.txt \
methyGff_for-heatmap_0DWT.genomewide.gene.meth.TSS.chg.txt \
methyGff_for-heatmap_0DWT.genomewide.gene.meth.TTS.chg.txt \
methyGff_for-heatmap_0DWT.genomewide.gene.meth.TSS.chh.txt \
methyGff_for-heatmap_0DWT.genomewide.gene.meth.TTS.chh.txt \
-l mCG_TSS mCG_TTS CHG_TSS CHG_TTS CHH_TSS CHH_TTS \
-o test_methyGff_0DWT.genomewide.gene.meth.TSSTTS.cgchgchh.pdf \
--plotmatrix 3x2 --colorMap vlag --centerlabel center \
-t whole_genome_0DWT \
-z mCG CHG CHH &

python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m \
methyGff_for-heatmap_3DWT.genomewide.gene.meth.TSS.cg.txt \
methyGff_for-heatmap_3DWT.genomewide.gene.meth.TTS.cg.txt \
methyGff_for-heatmap_3DWT.genomewide.gene.meth.TSS.chg.txt \
methyGff_for-heatmap_3DWT.genomewide.gene.meth.TTS.chg.txt \
methyGff_for-heatmap_3DWT.genomewide.gene.meth.TSS.chh.txt \
methyGff_for-heatmap_3DWT.genomewide.gene.meth.TTS.chh.txt \
-l mCG_TSS mCG_TTS CHG_TSS CHG_TTS CHH_TSS CHH_TTS \
-o test_methyGff_3DWT.genomewide.gene.meth.TSSTTS.cgchgchh.pdf \
--plotmatrix 3x2 --colorMap vlag --centerlabel center \
-t whole_genome_3DWT \
-z mCG CHG CHH &

python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m \
methyGff_for-heatmap_7DWT.genomewide.gene.meth.TSS.cg.txt \
methyGff_for-heatmap_7DWT.genomewide.gene.meth.TTS.cg.txt \
methyGff_for-heatmap_7DWT.genomewide.gene.meth.TSS.chg.txt \
methyGff_for-heatmap_7DWT.genomewide.gene.meth.TTS.chg.txt \
methyGff_for-heatmap_7DWT.genomewide.gene.meth.TSS.chh.txt \
methyGff_for-heatmap_7DWT.genomewide.gene.meth.TTS.chh.txt \
-l mCG_TSS mCG_TTS CHG_TSS CHG_TTS CHH_TSS CHH_TTS \
-o test_methyGff_7DWT.genomewide.gene.meth.TSSTTS.cgchgchh.pdf \
--plotmatrix 3x2 --colorMap vlag --centerlabel center \
-t whole_genome_7DWT \
-z mCG CHG CHH &


`-GENE` `mCG`
python /home/vcm/miniconda3/envs/batmeth/BatMeth2/bin/bt2heatmap.py -m test_methyGff_7DWT.gene.meth.GENE.cg.txt -l Geme -o test_methyGff_7DWT.gene.meth.GENE.cg.pdf --colorMap vlag --centerlabel center -z mCG

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


