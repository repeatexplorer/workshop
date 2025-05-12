# Input Data


| Description                                   | Type  | File Name / Location                                      |
| --------------------------------------------- | ----- | --------------------------------------------------------- |
| Genome assembly for tandem repeat annotation  | FASTA | `/mnt/data/tiny_pea.fasta`                                |
| Tandem repeat library                         | FASTA | `/mnt/data/Tandem_repeat_library.fasta`                   |
| RepeatMasker custom library                   | FASTA | `/mnt/data/RM_custom_library.fasta`                       |
| Configuration file for annotation pipeline    | YAML  | `/mnt/data/config.yaml`                                   |
| Singularity container for annotation pipeline | SIF   | `/mnt/data/assembly_repeat_annotation_pipeline_0.6.7.sif` |

# Example Analysis on Full Data (*P. sativum* Cameor v2 Assembly)


| Description                     | Directory                                      |
| ------------------------------- | ---------------------------------------------- |
| TideCluster analysis            | `/mnt/data/example_analyses/TideCluster`       |
| DANTE analysis                  | `/mnt/data/example_analyses/DANTE`             |
| DANTE_LTR analysis              | `/mnt/data/example_analyses/DANTE_LTR`         |
| DANTE_TIR analysis              | `/mnt/data/example_analyses/DANTE_TIR`         |
| Full repeat annotation pipeline | `/mnt/data/example_analyses/RepeatAnnotations` |

# Conda Environments

## Available Environments

```bash
conda env list
```

```txt
# conda environments:
#
base                 * /home/helix/miniforge3
cd-hit                 /home/helix/miniforge3/envs/cd-hit
dante_ltr              /home/helix/miniforge3/envs/dante_ltr
dante_tir              /home/helix/miniforge3/envs/dante_tir
singularity            /home/helix/miniforge3/envs/singularity
tidecluster            /home/helix/miniforge3/envs/tidecluster
```

## Example Installations (do not run)

All programs can be installed using Mamba or Conda; examples use Mamba:

```bash
mamba create -n dante_ltr   -c conda-forge -c bioconda -c petrnovak dante dante_ltr
mamba create -n dante_tir   -c conda-forge -c r        -c bioconda -c petrnovak dante_tir
mamba create -n singularity -c conda-forge -c bioconda singularity
mamba create -n tidecluster -c conda-forge -c bioconda -c petrnovak tidecluster
```

# Singularity Container for Repeat Annotation

The Singularity container (DOI:
[https://doi.org/10.5281/zenodo.15234515](https://doi.org/10.5281/zenodo.15234515)) can be downloaded from:

[https://zenodo.org/records/15234516/files/assembly_repeat_annotation_pipeline_0.6.7.sif?download=1](https://zenodo.org/records/15234516/files/assembly_repeat_annotation_pipeline_0.6.7.sif?download=1)

# TideCluster

> **Note:** The tool is under development; see [TideCluster
> GitHub](https://github.com/kavonrtep/TideCluster) for updates.
>
> **Credits**: TideCluster uses TideHunter for tandem repeat detection
> (Gao et al., *Bioinformatics* 35(14):i200--i207;
> [https://doi.org/10.1093/bioinformatics/btz376](https://doi.org/10.1093/bioinformatics/btz376))

## Installation (do not run; already installed)

```bash
mamba create -n tidecluster -c conda-forge -c bioconda -c petrnovak tidecluster
```

### Prerequisites

- Activate the `tidecluster` Conda environment.
- Working directory: `~/tidecluster/`.
- Input FASTA: `/mnt/data/tiny_pea.fasta`.
- Reference library (RepeatMasker format):
  `/mnt/data/Tandem_repeat_library.fasta`.
- CPU threads: adjust `-c` as needed (example uses 10).

## Running Individual Analysis Steps with Default Parameters

1. Open a terminal, switch to the TideCluster directory, and activate
   the environment:

   ```bash
   cd ~/tidecluster/
   conda activate tidecluster
   ```
2. Execute TideHunter:

   ```bash
   TideCluster.py tidehunter \
     -f /mnt/data/tiny_pea.fasta \
     -pr tiny_pea_default \
     -c 10
   ```

## Visualization in IGV

1. Launch IGV.
2. **Genomes → Load Genome from File** → select
   `/mnt/data/tiny_pea.fasta`.
3. **File → Load from File** → select:
   - `~/tidecluster/tiny_pea_default_tidehunter.gff3`
   - `~/tidecluster/tiny_pea_default_chunks.bed`
4. Right-click tracks to set:
   - `*_chunks.bed` → **Squished**
   - `*_tidehunter.gff3` → **Expanded**
5. **File → Save Session...**

## Clustering of Tandem Repeats

```bash
TideCluster.py clustering \
  -f /mnt/data/tiny_pea.fasta \
  -pr tiny_pea_default \
  -c 10
```

- Inspect results in IGV:
  - **File → Load from File** → select
    `~/tidecluster/tiny_pea_default_clustering.gff3`.

## Annotation Using Custom Reference Library

```bash
TideCluster.py annotation \
  -pr tiny_pea_default \
  -l /mnt/data/Tandem_repeat_library.fasta \
  -c 10
```

- Inspect results in IGV:
  - **File → Load from File** → select
    `~/tidecluster/tiny_pea_default_annotation.gff3`.

## Update GFF3 to Show Annotations from the Custom Library

1. Open `tiny_pea_default_annotation.tsv` in LibreOffice Calc (or
   Excel).
2. Sort by **Column C** (descending).
3. Remove rows with score \< 0.5.
4. Delete **Column C**.
5. Save as `tiny_pea_default_annotation_refDB.csv`.

```bash
tc_update_gff3.py \
  -g tiny_pea_default_annotation.gff3 \
  -t tiny_pea_default_annotation_refDB.csv \
  -o tiny_pea_default_annotation_refDB.gff3
```

- Inspect results in IGV:
  - **File → Load from File** → select
    `~/tidecluster/tiny_pea_default_annotation_refDB.gff3`.
  - **File → Save Session...**

## Consensus Sequence Generation with TAREAN

```bash
TideCluster.py tarean \
  -f /mnt/data/tiny_pea.fasta \
  -pr tiny_pea_default \
  -c 10
```

- Inspect `tiny_pea_default_tarean_index.html` in your web browser.

## Automatic Pipeline Execution (Short Monomers)

*In this example, we demonstrate automatic execution of all analysis
steps with custom settings to detect short monomer repeats.*

```bash
TideCluster.py run_all \
  -f /mnt/data/tiny_pea.fasta \
  -pr tiny_pea_short_monomers \
  -l /mnt/data/Tandem_repeat_library.fasta \
  -c 10 \
  -T "-p 10 -P 39 -c 5 -e 0.25"
```

- Inspect results in IGV:
  - **File → Load from File** → select:
    - `~/tidecluster/tiny_pea_short_monomers_tidehunter.gff3`
    - `~/tidecluster/tiny_pea_short_monomers_annotation.gff3`
  - **File → Save Session...**

> **Note:** No TAREAN output is generated if no repeat family passes the
> 50 kb total length cutoff.


## TideCluster output files

#### Tidehunter Step

the `prefix` is the prefix used in the command line so here it is `tiny_pea_default` and `tiny_pea_short_monomers`

- `prefix_tidehunter.gff3` - GFF3 file with tandem repeats detected by TideHunter.
- `prefix_chunks.bed` - BED file showing how the reference sequence was split into chunks for parallel processing.

#### Clustering Step

- `prefix__tidehunter_short.gff3` GFF3 file with tandem repeats shorter than the minimum length threshold used in the clustering step.
- `prefix_clustering.gff3` - GFF3 file with tandem repeats identified by `mmseqs2` and `BLASTN`.
  Tandem repeat regions in the GFF3 file are labeled by **T**andem **R**epeat **C**luster ID (TRC1, TRC2, etc.). Each TRC is described by the `repeat_type` attribute. `repeat_type` can be either TR (Tandem Repeat) or SSR (Simple Sequence Repeat).
- `prefix_clustering.gff3_1.gff3` - Intermediate file with tandem repeats clustered by `mmseqs2`.
- `prefix_consensus` - Directory with consensus sequences for each cluster as identified by TideHunter. There is one FASTA file per cluster. Each FASTA file contains all consensus sequences identified by TideHunter for a given cluster.
- `prefix_consensus_1` - Intermediate directory with consensus sequences for each cluster as identified by `mmseqs2`.
- `prefix_clustering_split_files` - Directory with GFF3 files, one for each TRC cluster. Each GFF3 file contains tandem repeat regions for a single TRC cluster.

#### Annotation Step

- `prefix_annotation.gff3` - GFF3 file with tandem repeats annotated by RepeatMasker.
  Annotations are shown as additional attributes in the GFF3 file.
- `prefix_annotation.tsv` - Summarized annotation for each TRC cluster in a tab-delimited format.
- `prefix_annotation_split_files` - Directory with GFF3 files, one for each TRC cluster.
  Each GFF3 file contains tandem repeat annotations for a single TRC cluster.

#### TAREAN Step

- `prefix_index.html` - Main HTML report, other reports are linked from this file.
- `prefix_tarean_report.html` - HTML report with tandem repeat annotations.
- `prefix_tarean_report.tsv` - File with tandem repeat annotations in a tab-delimited format.
- `prefix_kite_report.html` - HTML report with KITE analysis.
- `prefix_trc_superfamilies.html`  HTML report with TRC superfamilies.
- `prefix_trc_superfamilies.tsv`  File with TRC superfamilies in a tab-delimited format.
- `prefix_tarean` - Directory containing subdirectories with detailed TAREAN output for each TRC cluster.
- `prefix_consensus_dimer_library.fasta` - FASTA file with consensus sequences for
  each TRC cluster. This sequences can be used as a library for similarity based
  annotation using RepeatMasker. This file is created only for TRC clusters that
  pass the minimum combined length threshold.


# Annotation of Protein Domains with DANTE

## Prerequisites

- Activate the `dante_ltr` Conda environment.
- Working directory: `~/dante_ltr/`.
- Input FASTA: `/mnt/data/tiny_pea.fasta`.

## Running DANTE

```bash
conda activate dante_ltr
# Check versions
dante --version
dante_ltr --version
```

> **Note:** DANTE_LTR versions up to 0.3.5.3 are compatible with REXdb
> Viridiplantae v3.0. Versions ≥4.0.1 support REXdb Viridiplantae v3.0
> and v4.0. DANTE versions up to 0.1.9 include REXdb Viridiplantae 3.0;
> versions ≥0.2.0 include REXdb Viridiplantae 4.0.

```bash
mkdir ~/te_annotation
cd ~/te_annotation/
GENOME=/mnt/data/tiny_pea.fasta
dante -q $GENOME -o DANTE.gff3 -c 15
# This step takes about 6 minutes
```

> **Note:** DANTE is optimized for genome assemblies. For many short
> sequences, use `-S (--short_reads)`.

DANTE.gff3 contains detected protein domain annotations and can be used
as input for DANTE_LTR and DANTE_TIR. For phylogenetic analyses, filter
domains using `dante_gff_output_filtering.py`.

### Example of Filtering DANTE Output

```bash
# Default filtering
dante_gff_output_filtering.py --dom_gff DANTE.gff3 \
  -ouf DANTE_filtered_default.gff3 \
  -dps DANTE_filtered_default.fasta
```

Default thresholds:

- Max interruptions: 3 (per 100 AA)
- Min alignment length proportion: 0.8
- Max alignment length proportion: 1.2
- Min alignment identity: 0.35
- Min alignment similarity: 0.45

```bash
# Default thresholds with additional filters for Ty3/gypsy and RT domains
dante_gff_output_filtering.py --dom_gff DANTE.gff3 \
  -el gypsy -sd RT \
  -ouf DANTE_filtered_Gypsy_RT.gff3 \
  -dps DANTE_filtered_Gypsy_RT.fasta
```

## Extracting DNA Sequences from Filtered DANTE Output

```bash
dante_gff_to_dna.py \
  -i $GENOME \
  -d DANTE_filtered_Gypsy_RT.gff3 \
  -out DANTE_filtered_Gypsy_RT_dna \
  -ex
```

DNA sequences are output to the specified directory, one FASTA file per
lineage.

# Identification of Full-Length LTR-Retrotransposons with DANTE_LTR

## Prerequisites

- Activate the `dante_ltr` Conda environment.
- Input FASTA: `/mnt/data/tiny_pea.fasta`.
- Unfiltered DANTE output: `DANTE.gff3`.

## Running DANTE_LTR

```bash
cd ~/te_annotation/
dante_ltr -g DANTE.gff3 -s $GENOME -o DANTE_LTR -c 5 -M 1
```

> **Details:** - Input must be unfiltered DANTE GFF3. - `-M 1` sets the
> maximum missing protein domains. - `--te_constrains` specifies TE
> search constraints. - `--no_ambiguous_domains` removes ambiguous
> domains from analysis.

**Output files:**

- `DANTE_LTR.gff3`: Annotated elements.
- Extracted DNA sequences:
  - `DANTE_LTR_D.fasta` (partial elements without LTRs)
  - `DANTE_LTR_DL.fasta` (elements with LTRs and protein domains)
  - `DANTE_LTR_DLP.fasta` (elements with LTRs, domains, PBS)
  - `DANTE_LTR_DLT.fasta` (elements with LTRs, domains, TSD)
  - `DANTE_LTR_DLTP.fasta` (complete elements with LTRs, domains,
    PBS, TSD)
- Summaries:
  - `DANTE_LTR_statistics.csv`
  - `DANTE_LTR_summary.csv`
  - `DANTE_LTR_summary.html`

## Creating an LTR-RT Library for RepeatMasker

```bash
dante_ltr_to_library -g DANTE_LTR.gff3 -s $GENOME -o DANTE_LTR_library -c 5 -m 5
```

> **Options:** - `-m 5`: minimum cluster coverage (default: 3). -
> `--proportion_min PROPORTION_MIN` (default: 0.95).

**Output (`DANTE_LTR_library/`):**

```
DANTE_LTR_library/
├── TE_DL.fasta
├── TE_DLP.fasta
├── TE_DLT.fasta
├── TE_DLTP.fasta
├── TE_DLplus.fasta
├── TE_all.fasta
└── mmseqs2/
    ├── mmseqs_all_seqs.fasta
    ├── mmseqs_cluster.tsv
    ├── mmseqs_rep_seq.fasta
    ├── mmseqs_representative_seq_clean.fasta
    ├── mmseqs_representative_seq_clean_rm_compatible.fasta  <- Custom library for RepeatMasker
    └── partitioned_s900_w1000.fasta
```

# Annotation of DNA Transposons with DANTE_TIR

> **Note:** The tool is under development; see [DANTE_TIR
> GitHub](https://github.com/kavonrtep/dante_tir) for updates.

## Prerequisites

- Activate the `dante_tir` Conda environment.
- Input FASTA: `/mnt/data/tiny_pea.fasta`.
- Unfiltered DANTE output: `DANTE.gff3`.

```bash
cd ~/te_annotation/
conda activate dante_tir
dante_tir.py -g DANTE.gff3 -f $GENOME -o DANTE_TIR -c 10
```

**Output (`DANTE_TIR/`):**

- `DANTE_TIR_final.gff3` (final annotations)
- Extracted DNA sequences:
  - `DANTE_TIR_EnSpm_CACTA.fasta`
  - `DANTE_TIR_MuDR_Mutator.fasta`
  - `DANTE_TIR_hAT.fasta`
  - `DANTE_TIR_final.fasta` (all elements)
- `TIR_classification_summary.txt` (tabular summary)

# Exploration of Results in IGV

*(TODO: use example data)*
