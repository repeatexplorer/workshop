# TEtrimmer


## Introduction of TEtrimmer
Many tools have been developed for the discovery and annotation of transposable elements (TEs). 
However, the high-quality TE consensus library construction still requires manual curation of TEs, 
which is time-consuming and needs experts with an in-depth understanding of TE biology.

TEtrimmer (https://github.com/qjiangzhao/TEtrimmer) is a software designed to automate the manual curation of TEs. The input can be a TE library from 
*de novo* TE discovery tools, such as EDTA and RepeatModeler2, or a TE library from closely related species. 
For each input consensus sequence, TEtrimmer automatically performs BLASTN search,  sequence extraction, extension, 
multiple sequence alignment (MSA), MSA clustering, MSA cleaning, TE boundary definition, and TE classification. 
TEtrimmer also provides a graphical user interface (GUI) to inspect and improve predicted TEs, which can assist 
achieving manual curation-level TE consensus libraries easily. 

## Run TEtrimmer
### Copy TEtrimmer source code to your local folder
```bash
cp -r /mnt/data/TEtrimmer /home/helix/Documents
```

### Check TEtrimmer help messages

```bash
# Activate the conda "singularity" environment

conda activate singularity

# Test TEtrimmer with Singularity

singularity exec --writable-tmpfs \
--bind /home/helix/Documents/TEtrimmer:/bind_path \
--bind /mnt/data/pfam_database:/pfam \
/mnt/data/tetrimmer_1.4.0--hdfd78af_0.sif \
python /bind_path/tetrimmer/TEtrimmer.py --help
```


### Run TEtrimmer with Singularity
```bash
singularity exec --writable-tmpfs \
--bind /home/helix/Documents/TEtrimmer:/bind_path \
--bind /mnt/data/pfam_database:/pfam \
/mnt/data/tetrimmer_1.4.0--hdfd78af_0.sif \
python /bind_path/tetrimmer/TEtrimmer.py \
-i /bind_path/tests/test_input.fa \
-g /bind_path/tests/test_genome.fasta \
-o /bind_path/TEtrimmer_output \
--pfam_dir /pfam \
-t 5 --classify_all
```
### Output structure

- 📁**Classification** - *This folder is used for TE classifications.*  
- 📁**Multiple_sequence_alignment** - *All raw files will be stored in this folder if < --debug > is enabled.*
  - 📄**error_file.txt** - *Error file to store all error messages, only visible if errors were found.*
- 📁**Single_fasta_files** - *All sequences in the input file will be separated into single FASTA files and stored here.*
- 📁**TEtrimmer_for_proof_curation** - *This folder contains files used for manual inspection of TEtrimmer annotations.* 
  - 📁**Annotation_perfect** - *Four files are associated with each sequence as shown below.*
    - 📄**TE_name.raw.fa** - *Multiple sequence alignment file before TE boundary definition.*
    - 📄**TE_name.fa** - *Multiple sequence alignment file after TE boundary definition, which is used to generate the consensus sequence.*
    - 📄**TE_name.pdf** - *Plot file used to evaluate output.*
    - 📄**TE_name.cluster.fa** - *Multiple sequence alignment file before clustering.*
  - 📁**Annotation_good** 
  - 📁**Annotation_check_recommended**
  - 📁**Annotation_check_required**
  - 📁**Clustered_proof_curation** - *This folder contains all the output files from folder "Annotation_perfect", "Annotation_good", "Annotation_check_recommended", and "Annotation_check_required". The difference is TEtrimmer group similar output TEs into one "Cluster", which can make it easier to compare similar outputs.*
  - 📁**TE_low_copy** - *This folder contains low copy TEs.*
  - 📁**TE_skipped** - *Contains TE_Aid plots for all skipped TEs.*
- 📄**Sequence_name_mapping.txt** - *This file connects the input sequence names with the modified names from TEtrimmer.*
- 📄**summary.txt** - *Summary file.* 
- 📄**TEtrimmer_consensus.fasta** - *TE consensus library file before de-duplication.*
- 📄**TEtrimmer_consensus_merged.fasta** - *TE consensus library file after de-duplication.*


## Inspect outputs with TEtrimmerGUI
### Start TEtrimmerGUI
```bash
# Use the conda "base" environment for TEtrimmerGUI

conda deactivate

# Provide execute permission 

chmod +x /home/helix/Documents/TEtrimmer/tetrimmerGUI/blast/*

# Check TEtrimmerGUI help messages

python /home/helix/Documents/TEtrimmer/tetrimmerGUI/TEtrimmerGUI.py --help

# Run TEtrimmerGUI
python /home/helix/Documents/TEtrimmer/tetrimmerGUI/TEtrimmerGUI.py \
-i /home/helix/Documents/TEtrimmer/TEtrimmer_output/TEtrimmer_for_proof_curation \
-g /home/helix/Documents/TEtrimmer/tests/test_genome.fasta \
-o /home/helix/Documents
```


# Repeat Annotations with REPET

https://forgemia.inra.fr/urgi-anagen/wiki-repet/-/wikis/REPET_Snakemake-:-REPET-V4

# Repeat Annotation Pipeline (Viridiplantae)

The Assembly annotation pipeline is a comprehensive genome annotation workflow that
integrates multiple specialized tools to identify and classify various types of repetitive
elements in genomic sequences.

- The pipeline uses DANTE/DANTE_LTR to identify intact LTR retrotransposons.
- TideCluster is employed to identify tandem repeats, with separate processes for
  default-length and short monomer repeats.
- The pipeline then creates custom libraries of repeat sequences, including those from LTR
  retrotransposons and tandem repeats. This library can be supplemented with user-provided
  custom repeat databases.
- After building these repeat libraries, the pipeline uses RepeatMasker to annotate the
  genome comprehensively.
- The workflow produces detailed GFF3 annotation files for different repeat classes 
  (mobile elements, simple repeats, low complexity regions, rDNA), density visualizations
  as bigWig files, and summary statistics and plots.

## Input Data

| Description                            | Type  | File Name / Location                                      |
|----------------------------------------|-------|-----------------------------------------------------------|
| Genome assembly for repeat annotation  | FASTA | `/mnt/data/tiny_pea.fasta`                                |
| Tandem repeat library                  | FASTA | `/mnt/data/Tandem_repeat_library.fasta`                   |
| RepeatMasker custom library            | FASTA | `/mnt/data/RM_custom_library.fasta`                       |
| Configuration file for annotation pipeline | YAML  | `/mnt/data/config.yaml`                                   |
| Singularity container for annotation pipeline | SIF   | `/mnt/data/assembly_repeat_annotation_pipeline_0.6.7.sif` |


The pipeline supports supplying an additional custom repeat library via the
optional `custom_library`
parameter. When provided, RepeatMasker will use this library for similarity-based
annotation. Sequences must be in FASTA format, with headers following the convention:
`>repeatname#class/subclass`

The pipeline recognizes the following classification categories:

```txt
Class_I/LINE
Class_I/LTR/Ty1_copia
Class_I/LTR/Ty1_copia/Ale
Class_I/LTR/Ty1_copia/Angela
Class_I/LTR/Ty1_copia/Bianca
Class_I/LTR/Ty1_copia/Ikeros
Class_I/LTR/Ty1_copia/Ivana
Class_I/LTR/Ty1_copia/SIRE
Class_I/LTR/Ty1_copia/TAR
Class_I/LTR/Ty1_copia/Tork
Class_I/LTR/Ty1_copia/Alexandra
Class_I/LTR/Ty1_copia/Bryana
Class_I/LTR/Ty1_copia/Ferco
Class_I/LTR/Ty3_gypsy
Class_I/LTR/Ty3_gypsy/chromovirus
Class_I/LTR/Ty3_gypsy/chromovirus/Ferney
Class_I/LTR/Ty3_gypsy/chromovirus/CRM
Class_I/LTR/Ty3_gypsy/chromovirus/Reina
Class_I/LTR/Ty3_gypsy/chromovirus/Tekay
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tatius
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Athila
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/Ogre
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/Retand
Class_II/Subclass_1/TIR/EnSpm_CACTA
Class_II/Subclass_1/TIR/hAT
Class_II/Subclass_1/TIR/MITE
Class_II/Subclass_1/TIR/MITE/Stowaway
Class_II/Subclass_1/TIR/MuDR_Mutator
Class_II/Subclass_1/TIR/PIF_Harbinger
Class_II/Subclass_1/TIR/Tc1_Mariner
Class_II/Subclass_2/Helitron
rDNA_45S/18S
rDNA_45S/25S
rDNA_45S/5.8S
rDNA_5S/5S
```

Currently, the pipeline bundles DANTE, DANTE_LTR, TideCluster and RepeatMasker—but not
DANTE_TIR—because DANTE_TIR is still experimental. If you also need to annotate DNA
transposons, you must supply your own custom library. We recommend generating that library
with the standalone DANTE_TIR pipeline and then manually curating its entries to minimize
false-positive annotations. Alternatively, you can use Class II transposon clusters
obtained via RepeatExplorer2. Library of rDNA sequences is built-in in the pipeline. If Class II
custom library is provided, the pipeline uses this library to identify LTR-RT elements
with potential Class II insertions. Such elements are not used in similarity-based
annotation.

## Pipeline steps

- DANTE annotation, DANTE output filtering
- DANTE_LTR annotation
- Annotation of tandem repeats using TideCluster
- Annotation of tandem repeats using custom library.
- Building custom libraries for RepeatMasker using DANTE_LTR and TideCluster outputs and custom
  libraries
- Library reduction (optional)
- RepeatMasker annotation of dispersed repeats
- RepeatMasker annotation of tandem repeats
- Reconciliation of RepeatMasker outputs and tandem repeat annotations
- Coverage track generation


## Running the Repeat Annotation Pipeline

### Prerequisites
- Singularity installed (available in `singularity` conda environment)
- Pipeline singularity container (DOI: https://doi.org/10.5281/zenodo.15234515) can be downloaded
from: https://zenodo.org/records/15234516/files/assembly_repeat_annotation_pipeline_0.6.7.sif?download=1

### Setup
To run the pipeline, you need a configuration file named `config.yaml`. This file should contain the following parameters:

```yaml
genome_fasta: /mnt/data/test_data/tiny_pea.fasta
output_dir: output
custom_library: /mnt/data/test_data/RM_custom_library.fasta
tandem_repeat_library: /mnt/data/test_data/Tandem_repeat_library.fasta
```

Additional optional parameters are:
```yaml
repeatmasker_sensitivity: default # possible values are : sensitive, default, quick,
reduce_library: True # possible values are: True, False, if missing, True is used
```


> **Note**: Use absolute paths for all files in the configuration file. The pipeline will
> create a directory named `output` in the current working directory to store
> all the output files.

### Running the pipeline

```bash
conda activate singularity
mkdir ~/Repeat_Annotations
cd ~/Repeat_Annotations
singularity run -B /mnt/data -B $PWD /mnt/data/assembly_repeat_annotation_pipeline_0.6.7.sif -c /mnt/data/config.yaml -t 10
# the running time is about 1 hour 
```
> **Details**: 
> - The `-B` flag binds host directories /mnt/data and the current working directory
> - $PWD is shell variable set to the current working directory
> - -t 10 sets the number of threads to use during execution
> - -c specifies the path to your configuration file
> - All data must be either in the current working directory or in the /mnt/data directory to be accessible by the container


### Expected running times

| Species                      | Genome Size | Runtime | CPU Cores | RAM    |
|------------------------------|-------------|---------|-----------|--------|
| Filipendula ulmaria          | 275 Mb      | 56 min  | 20        | 256 GB | 
| Ailanthus altissima assembly | 0.94 Gb     | 5 hrs   | 20        | 256 GB |
| Ballota nigra assembly       | 1.2 Gb      | 9 hrs   | 20        | 256 GB |
| Pisum sativum assembly       | 4.2 Gb      | 4 days  | 24        | 256 GB |
| Vicia faba                   | 12.1 Gb     | 22 days | 36        | 384 GB | 



## Output data structure

```text
<output_dir>/
├── DANTE/
│   ├── DANTE.gff3                 # Raw DANTE annotation of retrotransposon domains
│   └── DANTE_filtered.gff3        # Filtered DANTE output (high-confidence domains)
├── DANTE_LTR/
│   ├── DANTE_LTR.gff3             # LTR retrotransposon coordinates from DANTE_LTR
│   ├── DANTE_LTR_summary.html     # Graphical summary report of LTR-RT detection
│   └── LTR_RTs_library.fasta      # Extracted LTR-RT sequences library
├── TideCluster/
│   ├── default/                                  # TideCluster with default setting monopmer size  
│   │   ├── TideCluster_clustering.gff3           # GFF3 file with tandem repeats detected by TideHunter
│   │   ├── TideCluster_tidehunter_short.gff3     # GFF3 file with tandem repeat clusters
│   │   ├── TideCluster_annotation.gff3           # GFF3 file with tandem repeat clusters annotated by RepeatMasker.
│   │   └── *_10k.bw, *_100k.bw                   # Tandem repeat density tracks (10 Kb and 100 Kb windows)
│   └── short_monomer/
│       ├── TideCluster_clustering.gff3           # Short-monomer repeat clusters
│       ├── TideCluster_annotation.gff3           # Short-monomer repeat clusters
│       └── TideCluster_tidehunter_short.gff3     # Short-monomer predictions
├── Libraries/
│   ├── class_ii_library.fasta         # DNA transposon library (Class II) if provided or empty
│   ├── LTR_RTs_library_clean.fasta    # Filtered LTR-RT library (no subclass 2 contamination)
│   ├── combined_library.fasta         # Concatenated library for RepeatMasker (LTR + custom + rDNA)
│   └── combined_library_reduced.fasta # Reduced-size library (optional clustering)
├── RepeatMasker/
│   ├── RM_on_combined_library.out             # Raw RepeatMasker output
│   ├── RM_on_combined_library.gff3            # RepeatMasker annotations converted to GFF3 format
│   ├── RM_on_combined_library_plus_DANTE.gff3 # Merged RepeatMasker + DANTE filtered domains
│   └── Repeat_Annotation_NoSat.gff3           # Final annotation with satellites subtracted
├── Repeat_Annotation_NoSat_split_by_class_gff3/
│         ├── 
│         ├──                                       # GFF3 files for each repeat class         
│
├── Repeat_Annotation_NoSat_split_by_class_bigwig/
│         ├──                                       # BigWig files for each repeat class
│         ├──          
│
├── All_Ty1_Copia_RepeatMasker.gff3              |
├── All_Ty3_Gypsy_RepeatMasker.gff3              |
├── Low_complexity_RepeatMasker.gff3             | Similarity basesd annotations
├── Mobile_elements_RepeatMasker.gff3            |  (with RepeatMasker)
├── Simple_repeats_RepeatMasker.gff3             |
├── Tandem_repeats_RepeatMasker.gff3             |
├── rDNA_RepeatMasker.gff3                       |
│
├── DANTE_filtered.gff3                            | 
├── DANTE_LTR.gff3                                 | Structure based annotations
├── Tandem_repeats_TideCluster_annotated.gff3      |
├── Tandem_repeats_TideCluster.gff3                |
│
├── all_repeats_for_masking.bed       # Merged BED of all repeats for masking steps
├── gaps_10plus.bed                   # Coordinates of gaps ≥10 Ns in the assembly
├── summary_statistics.csv            # Per-repeat class genome proportions
├── summary_plots.pdf                 # PDF report with plots of repeat distribution along chromosomes
├── TideCluster_report.html           # Graphical HTML summary of TideCluster results
└── DANTE_LTR_report.html             # Graphical summary report of LTR-RT detection
```

> **Note**
> - BigWig files (*.bw) are generated for both RepeatMasker and TideCluster annotations at 10Kb and 100Kb resolution.
> - Per-class GFF3 splits (e.g., All_Ty1_Copia.gff3, Simple_repeats.gff3) reside under Repeat_Annotation_NoSat_split_by_class_gff3/
> - Top level directory contains symbolic links to the most important annotation files
> - Tandem repeat annotations are available as structure-based (TideCluster) and similarity-based (RepeatMasker) annotations.


# Make IGV visualization with Example Analysis on Complete Genome Assembly of *P. sativum* Cameor v2

Output from repeat annotation pipeline for the complete genome assembly is stored in
`/mnt/data/example_analyses/RepeatAnnotations`

- run IGV
- Open Genome Assembly: **Genomes** -> ***Load Genome from File...***
  - Select the genome assembly file from `/mnt/data/Pisum_assembly_cameor_ver_2.fasta`
- Import following annotation tracks located in `/mnt/data/example_analyses/RepeatAnnotations` directory with **File** -> ***Load from File...***:

  
  - `Repeat_Annotation_NoSat_split_by_class_bigwig/10k/All_Ty1_Copia_10k.bw`
  - `Repeat_Annotation_NoSat_split_by_class_bigwig/10k/All_Ty3_Gypsy_10k.bw`
  - `Repeat_Annotation_NoSat_split_by_class_bigwig/10k/Mobile_elements_10k.bw` 
  - `Repeat_Annotation_NoSat_split_by_class_bigwig/10k/rDNA_10k.bw`
  - `Repeat_Annotation_NoSat_split_by_class_bigwig/10k/Simple_repeats_10k.bw`
  - `TideCluster/default/TideCluster_clustering_10k.bw`
  - `TideCluster/default/TideCluster_clustering_split_files_bigwig/10k/TRC_1_10k.bw`
  - `TideCluster/default/TideCluster_clustering_split_files_bigwig/10k/TRC_2_10k.bw`
  - `TideCluster/default/TideCluster_clustering_split_files_bigwig/10k/TRC_3_10k.bw`
  - `TideCluster/default/TideCluster_clustering_split_files_bigwig/10k/TRC_4_10k.bw`
  - `TideCluster/default/TideCluster_clustering_split_files_bigwig/10k/TRC_5_10k.bw`
  - `Tandem_repeats_TideCluster.gff3`
  - `Tandem_repeats_TideCluster_annotated.gff3`
  - `Tandem_repeats_RepeatMasker.gff3`
  - `TideCluster/TideCluster_clustering_default_and_short_merged.gff3`
  - `DANTE_LTR.gff3`
  - `all_repeats_for_masking.bed`
  - `gaps_10plus.bed`
  - `Low_complexity_RepeatMasker.gff3`
  - `Simple_repeats_RepeatMasker.gff3`
  - `Mobile_elements_RepeatMasker.gff3`


Regions of interest: 

- `chr5:205,846,880-208,075,181` : Example of structure based TideCluster annotation and similarity based RepeatMasker annotation based TideCluster library
- `chr6:165,888,582-165,958,151` : Example of DANTE_LTR annotation and RepeatMasker annotation