## 0 Files and preparation
#### 0.1. ref_gtf file
**Path to Araport11 protein_coding_position**
```
/home/vcm/araport_reference/protein_coding.Araport11.position
```
#### 0.2. target file path
*target file must be tab-delimited (\t)*
*target file must be in bed format: Chr_number  Start  End  ...*
```
~/annotate/0D_log2genetable.txt
~/annotate/3D_log2genetable.txt
~/annotate/7D_log2genetable.txt
```
#### 0.3. Environment with bedtools
`conda activate bed`

## 1. use operation code to pair up two files based on gene names

