## 0 Files and preparation
#### 0.1. ref_gtf file
**Path to Araport11 protein_coding_position**
```
/home/vcm/araport_reference/protein_coding.Araport11.position
```
> $ head /home/vcm/araport_reference/protein_coding.Araport11.position
> 
> Chr1	3631	5899	AT1G01010
> 
> Chr1	6788	9130	AT1G01020
>
> Chr1	11649	13714	AT1G01030
>
> Chr1	23121	31227	AT1G01040
>
> Chr1	28500	28706	AT1G01046
>
> Chr1	31170	33171	AT1G01050
>
> Chr1	33365	37871	AT1G01060
>
> Chr1	38444	41017	AT1G01070
>
> Chr1	44970	47059	AT1G01080
>
>  Chr1	47234	49304	AT1G01090
> 

#### 0.2. target file path
*target file must be tab-delimited (\t)*
*target file must be in bed format: Chr_number  Start  End  ...*
```
~/annotate/0D_log2genetable.txt
~/annotate/3D_log2genetable.txt
~/annotate/7D_log2genetable.txt
```
> $ head ~/annotate/0D_log2genetable.txt
>
> ID	PPDE	PostFC	C1Mean	C2Mean	log2FC	log2_FDR
>
> AT1G07887	1	0.007387157	2.019619049	374.1787752	-7.080765044940396	0
>
> AT1G53480	1	0.030704193	28.91644564	965.4529907	-5.025420504627368	0
>
> AT2G01422	1	0.094681511	13.87824472	153.7498566	-3.400773459887111	0
>
> AT2G33175	1	69.56867144	332.6631522	4.042537859	6.120365864561457	0
>
> AT2G36490	1	0.082607247	285.7158436	3467.055651	-3.5975878372784478	0
>
> AT3G05945	1	9.051509061	348.6079391	37.8466168	3.1781583373149056	0
>
> AT3G19350	1	10.95219257	120.4369134	10.31504715	3.4531478134093248	0
>
> AT3G30720	1	0.130960161	46.46318128	359.7658596	-2.932800094338966	0
>
> AT1G51640	1	8.212735522	44.27969005	4.732875962	3.0378628391395925	0
> 
#### 0.3. Environment with bedtools
`conda activate bed`

## 1. use operation code to pair up two files based on gene names
```
awk 'NR==FNR {a[$1]; next} ($2 in a) {print}' fileA.tsv fileB.tsv > fileC.tsv
```
> Explanation:
>
> NR==FNR {a[$1]; next}: This part processes fileA.tsv (the first file).
>
> NR is the total number of records read so far, and FNR is the number of records read from the current file.
>
> When NR==FNR, it means awk is processing the first file (fileA.tsv).
>
> Here, a[$1] creates an associative array a where the key is the first column of fileA.tsv.
*Where fileA is the list of interested genes, fileB is the ref_genome gtf file, fileC is the new output file*

**Sample**
> for file in ~/annotate/*D_log2genetable.txt; do awk 'NR==FNR {a[$1]; next} ($4 in a) {print}' $file /home/vcm/araport_reference/protein_coding.Araport11.position > ${file%.txt}_GTF.tsv; echo ${file%.txt}_GTF.tsv; done



