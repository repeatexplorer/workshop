# Input Data

| Description                                   | Type  | File name / Location                                      |
|-----------------------------------------------|-------|-----------------------------------------------------------|
| Genome assembly for tandem repeat annotation  | FASTA | `/mnt/data/tiny_pea.fasta`                                |
| Tandem repeat library                         | FASTA | `/mnt/data/Tandem_repeat_library.fasta`                   |
| RepeatMasker custom library                   | FASTA | `/mnt/data/RM_custom_library.fasta`                       |
| Configuration file for annotation pipeline    | YAML  | `/mnt/data/config.yaml`                                   |
| Singularity container for annotation pipeline | SIF   | `/mnt/data/assembly_repeat_annotation_pipeline_0.6.7.sif` |

# Example Analysis on Full Data (*P. sativum* Cameor v2 Assembly)

| Description                     | Directory                                      |
|---------------------------------|------------------------------------------------|
| DANTE analysis                  | `/mnt/data/example_analyses/DANTE`             |
| DANTE_LTR analysis              | `/mnt/data/example_analyses/DANTE_LTR`         |
| DANTE_TIR analysis              | `/mnt/data/example_analyses/DANTE_TIR`         |
| Full repeat annotation pipeline | `/mnt/data/example_analyses/RepeatAnnotations` |

# Conda Environments

## Available environments

```bash
conda env list
```

```txt
# conda environments:
#
base                 * /home/helix/miniforge3
dante_ltr              /home/helix/miniforge3/envs/dante_ltr
dante_tir              /home/helix/miniforge3/envs/dante_tir
singularity            /home/helix/miniforge3/envs/singularity
tidecluster            /home/helix/miniforge3/envs/tidecluster
```

## Example Installations (do not run)

All programs can be installed using `mamba`:

```bash
mamba create -n dante_ltr -c conda-forge -c bioconda -c petrnovak dante dante_ltr
mamba create -n dante_tir -c conda-forge -c r -c bioconda -c petrnovak dante_tir
mamba create -n singularity -c conda-forge -c bioconda singularity
mamba create -n tidecluster -c conda-forge -c bioconda -c petrnovak tidecluster
```

# Singularity Container for Repeat Annotation

Singularity container ([DOI: 10.5281/zenodo.15234515](https://doi.org/10.5281/zenodo.15234515)) can be downloaded from:
[Zenodo Link](https://zenodo.org/records/15234516/files/assembly_repeat_annotation_pipeline_0.6.7.sif?download=1)

# TideCluster


## Installation (already installed; do not run)

```bash
mamba create -n tidecluster -c conda-forge -c bioconda -c petrnovak tidecluster
```

### Prerequisites

- Conda environment `tidecluster` activated
- Working directory: `~/tidecluster/`
- Input FASTA: `/mnt/data/tiny_pea.fasta`
- Reference library (RepeatMasker format): `/mnt/data/Tandem_repeat_library.fasta`
- CPU threads: adjust `-c` as needed (example uses 4)

## Running Individual Analysis Steps

Activate environment and switch directory:

```bash
cd ~/tidecluster/
conda activate tidecluster
```

Run TideHunter:

```bash
TideCluster.py tidehunter -f /mnt/data/tiny_pea.fasta -pr tiny_pea_default -c 4
```

## Visualization in IGV

- Launch IGV
- **Genomes → Load Genome from File** → `/mnt/data/tiny_pea.fasta`
- **File → Load from File** → select:
  - `~/tidecluster/tiny_pea_default_tidehunter.gff3`
  - `~/tidecluster/tiny_pea_default_chunks.bed`
- Set tracks:
  - `*_chunks.bed` → Squished
  - `*_tidehunter.gff3` → Expanded
- **File → Save Session**

## Clustering Tandem Repeats

```bash
TideCluster.py clustering -f /mnt/data/tiny_pea.fasta -pr tiny_pea_default -c 4
```
- Inspect in IGV:
  - Load: `~/tidecluster/tiny_pea_default_clustering.gff3`

## Annotation with Custom Reference Library

```bash
TideCluster.py annotation -pr tiny_pea_default -l /mnt/data/Tandem_repeat_library.fasta -c 4
```
- Inspect in IGV: Load `~/tidecluster/tiny_pea_default_annotation.gff3`

## Updating GFF3 with Selected Clusters

- Open `tiny_pea_default_annotation.tsv` in LibreOffice/Excel
- Sort by Column C (descending), remove rows with scores < 0.5
- Delete Column C, save as `tiny_pea_default_annotation_refDB.csv`

Update GFF3:

```bash
update_gff3.py -g tiny_pea_default_annotation.gff3 -t tiny_pea_default_annotation_refDB.csv -o tiny_pea_default_annotation_refDB.gff3
```

- Inspect in IGV: Load `~/tidecluster/tiny_pea_default_annotation_refDB.gff3`

## Consensus Sequence Generation (TAREAN)

```bash
TideCluster.py tarean -f /mnt/data/tiny_pea.fasta -pr tiny_pea_default -c 4
```
- Inspect HTML report in browser.

## Automatic Pipeline Execution (Short Monomers)

```bash
TideCluster.py run_all -f /mnt/data/tiny_pea.fasta -pr tiny_pea_short_monomers -l /mnt/data/Tandem_repeat_library.fasta -c 4 -T "-p 10 -P 39 -c 5 -e 0.25"
```
- Inspect in IGV:
  - Load:
    - `~/tidecluster/tiny_pea_short_monomers_tidehunter.gff3`
    - `~/tidecluster/tiny_pea_short_monomers_annotation.gff3`

> **NOTE:** No TAREAN output if no repeat family exceeds 50 kb length.

# Annotation of Protein Domains with DANTE

### Running DANTE

```bash
mkdir ~/te_annotation
cd ~/te_annotation/
conda activate dante_ltr
GENOME=/mnt/data/tiny_pea.fasta
dante -q $GENOME -o DANTE.gff3 -c 15
```

> **NOTE:** Use `-S` option for short sequences.

### Filtering DANTE Output

Default filtering:

```bash
dante_gff_output_filtering.py --dom_gff DANTE.gff3 -ouf DANTE_filtered_default.gff3 -dps DANTE_filtered_default.fasta
```

Ty3/Gypsy RT filtering:

```bash
dante_gff_output_filtering.py --dom_gff DANTE.gff3 -el gypsy -sd RT -ouf DANTE_filtered_Gypsy_RT.gff3 -dps DANTE_filtered_Gypsy_RT.fasta
```

### Extract DNA Sequences

```bash
dante_gff_to_dna.py -i $GENOME -d DANTE_filtered_Gypsy_RT.gff3 -out DANTE_filtered_Gypsy_RT_dna -ex
```

# Identification of Full-length LTR Retrotransposons (DANTE_LTR)

```bash
dante_ltr -g DANTE.gff3 -s $GENOME -o DANTE_LTR -c 5 -M 1
```

### Create RepeatMasker Library

```bash
dante_ltr_to_library -g DANTE_LTR.gff3 -s $GENOME -o DANTE_LTR_library -c 5 -m 5
```

# Annotation of DNA Transposons (DANTE_TIR)

```bash
cd ~/te_annotation/
conda activate dante_tir
dante_tir.py -g DANTE.gff3 -f $GENOME -o DANTE_TIR -c 10
```

# Exploration of Results in IGV

TODO - use example data

