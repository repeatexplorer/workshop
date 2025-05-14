# Running RepeatExplorer2 in command line using Singularity container

## Prerequisites 
 
- A working **mamba/conda**  setup (e.g. Mambaforge).
- Singularity ≥ 3.6 installed in a Conda env.
- The RepeatExplorer2 image already present at

```txt
/mnt/data/repex_tarean_0.3.12-7a7dc9e.sif
```

## Install Singularity (do not run, singularity is already installed) 

If you haven’t yet, create a Singularity env:
```bash
mamba create -n singularity -c conda-forge singularity>3.6.3
```

Activate it when you need Singularity:
```bash
conda activate singularity
```

## Set Up Your Analysis Directory 

Everything happens under `~/repeatexplorer`:
```bash
cd
mkdir -p repeatexplorer
cd repeatexplorer
```

## Singularity image 

Singulary image can be obtained using command:
```bash
# DO NOT RUN THIS
singularity pull library://repeatexplorer/default/repex_tarean:0.3.12-7a7dc9e
```
> **Note** : You do *not* need to `pull` or the container—just use already downloaded image at `/mnt/data`.



Getting help for `seqclust` (the main repeat clustering tool) inside the container:

```bash
singularity exec -e -B $HOME /mnt/data/repex_tarean_0.3.12-7a7dc9e.sif seqclust --help
```



## Running Clustering on Test Data 

 
**Download a small paired-end test set** :

```bash
cd ~/repeatexplorer
wget https://bitbucket.org/petrnovak/repex_tarean/raw/devel/test_data/LAS_paired_10k.fas
```
 
**Inspect it** :

```bash
head LAS_paired_10k.fas
seqkit stat LAS_paired_10k.fas
```
 
**Run RepeatExplorer’s clustering** :


```bash
cd ~/repeatexplorer
singularity exec -e \
  -B $HOME \
  /mnt/data/repex_tarean_0.3.12-7a7dc9e.sif \
  seqclust -p -v re_test LAS_paired_10k.fas
```
## Example of RepeatExplorer2 clustering including data pre-processing

### Download data:

Get data and check the quality of reads using FastQC program:
```bash
# Clustering example 1
cd ~/repeatexplorer
mkdir single_species
cd single_species
# Sample of T.cacao pair-end Illumina reads
wget https://github.com/kavonrtep/example_data/raw/master/SRR089356_1.fastq.gz
wget https://github.com/kavonrtep/example_data/raw/master/SRR089356_2.fastq.gz
fastqc *.fastq.gz
````

### Data Pre-processing with fastp 

 
**Run fastp**  in pair-end mode with HTML and JSON report

```bash

fastp \
  -i SRR089356_1.fastq.gz \
  -I SRR089356_2.fastq.gz \
  -o SRR089356_1_clean.fastq.gz \
  -O SRR089356_2_clean.fastq.gz \
  --trim_front1 10 --trim_front2 10 \
  --length_required 90 \
  --html fastp_report.html \
  --thread 4
```
 
**Check output stats** :


```bash
seqkit stats *clean*.fastq.gz
fastqc *clean*.fastq.gz
```

### Sampling & Interleaving 

 
**Sample to target coverage**  (e.g. 5 000 reads):


```bash
seqtk sample -s 10 SRR089356_1_clean.fastq.gz 5000 > sample_1.fastq
seqtk sample -s 10 SRR089356_2_clean.fastq.gz 5000 > sample_2.fastq
```
 
**Interleave and convert to FASTA** :


```bash
seqtk mergepe sample_1.fastq sample_2.fastq > merged.fastq
seqtk seq -A merged.fastq > merged.fasta
```

## Clustering with Default Settings 



```bash
cd ~/repeatexplorer
singularity exec -e \
  -B $HOME \
  /mnt/data/repex_tarean:0.3.12-7a7dc9e.sif \
  seqclust -p -v re_output_run1 \
  single_species/merged.fasta
```

## Comparative Analysis Workflow 

**Prepare two datasets**  in `~/repeatexplorer/comparative/`:


```bash
cd ~/repeatexplorer
# Get data from comparative analysis
cd ~/repeatexplorer
mkdir comparative
cd comparative
wget  https://github.com/kavonrtep/example_data/raw/master/SRR9938304_1.fastq.gz
wget  https://github.com/kavonrtep/example_data/raw/master/SRR9938304_2.fastq.gz
wget  https://github.com/kavonrtep/example_data/raw/master/SRR089356_1.fastq.gz
wget  https://github.com/kavonrtep/example_data/raw/master/SRR089356_2.fastq.gz

```
 
**QC & filter with fastp**  (repeat section 6 per dataset).

```bash
fastp \
  -i SRR9938304_1.fastq.gz \
  -I SRR9938304_2.fastq.gz \
  -o SRR9938304_1_clean.fastq.gz \
  -O SRR9938304_2_clean.fastq.gz \
  --trim_front1 10 --trim_front2 10 \
  --length_required 90 \
  --max_len1 90 --max_len2 90 \
  --html fastp_report.html \
  --thread 4
  
fastp \
    -i SRR089356_1.fastq.gz \
    -I SRR089356_2.fastq.gz \
    -o SRR089356_1_clean.fastq.gz \
    -O SRR089356_2_clean.fastq.gz \
    --trim_front1 10 --trim_front2 10 \
    --max_len1 90 --max_len2 90 \
    --length_required 90 \
    --html fastp_report.html \
    --thread 4
```
 
**Sample & interleave**  

```bash
seqtk sample -s 10 SRR9938304_1_clean.fastq.gz 5000 > sample_CA_1.fastq
seqtk sample -s 10 SRR9938304_2_clean.fastq.gz 5000 > sample_CA_2.fastq

seqtk sample -s 10 SRR089356_1_clean.fastq.gz 5000 > sample_CB_1.fastq
seqtk sample -s 10 SRR089356_2_clean.fastq.gz 5000 > sample_CB_2.fastq
 ```
Interlave and convert to FASTA:

```bash
seqtk mergepe sample_CA_1.fastq sample_CA_2.fastq > merged_CA.fastq
seqtk mergepe sample_CB_1.fastq sample_CB_2.fastq > merged_CB.fastq
# convert to FASTA
seqtk seq -A merged_CA.fastq > merged_CA.fasta
seqtk seq -A merged_CB.fastq > merged_CB.fasta
```

**Add prefixes & concatenate** :

```bash
# Rename to add CA/CB prefixes
seqtk rename merged_CA.fasta CA > merged_CA_prefix.fastq
seqtk rename merged_CB.fasta CB > merged_CB_prefix.fastq
cat merged_CA_prefix.fastq merged_CB_prefix.fastq > CA_CB_final.fasta
```
 
**Run comparative clustering** :


```bash
singularity exec -e \
  -B $HOME \
  /mnt/data/repex_tarean_0.3.12-7a7dc9e.sif \
  seqclust -p \
           --prefix_length 2 \
           -v re_output_comparative \
           CA_CB_final.fasta
```


## Custom TEMP Directory (Do not run)


If you need to redirect Singularity’s temp files:

```bash
singularity exec -e \
  --env TEMP=/mnt/tmp \
  -B /mnt/tmp:/mnt/tmp \
  -B $HOME \
  /mnt/data/repex_tarean_0.3.12-7a7dc9e.sif \
  seqclust -p -v re_test LAS_paired_10k.fas
```

Ensure `/mnt/tmp` (or your chosen path) has enough free space and is writable.

