# Repeat Annotaion Pipeline (Viridiplantae)

The Assembly annotation pipeline is a comprehensive genome annotation workflow that
integrates multiple specialized tools to identify and classify various types of repetitive
elements in genomic sequences.

- The pipeline uses DANTE/DANTE_LTR to identify intact LTR retrotransposons.
- TideCluster is employed to identify tandem repeats, with separate processes for
  default-length and short monomer repeats.
- The pipeline then creates custom libraries of repeat sequences, including those from LTR
  retrotransposons and Tandem repeats. This library can be supplements with user-provided
  custom repeat databases.
- After building these repeat libraries, the pipeline uses RepeatMasker to annotate the
  genome comprehensively.
- The workflow produces detailed GFF3 annotation files for different repeat classes (
  mobile elements, simple repeats, low complexity regions, rDNA), density visualizations
  as bigWig files, and summary statistics and plots.

## Input Data

| Description                                   | Type  | File Name / Location                                      |
|-----------------------------------------------|-------|-----------------------------------------------------------|
| Genome assembly for tandem repeat annotation  | FASTA | `/mnt/data/tiny_pea.fasta`                                |
| Tandem repeat library                         | FASTA | `/mnt/data/Tandem_repeat_library.fasta`                   |
| RepeatMasker custom library                   | FASTA | `/mnt/data/RM_custom_library.fasta`                       |
| Configuration file for annotation pipeline    | YAML  | `/mnt/data/config.yaml`                                   |
| Singularity container for annotation pipeline | SIF   | `/mnt/data/assembly_repeat_annotation_pipeline_0.6.7.sif` |

## Configuration file structure

```yaml
genome_fasta: /mnt/data/test_data/tiny_pea.fasta
output_dir: Repeat_Annotations
custom_library: /mnt/data/test_data/RM_custom_library.fasta
tandem_repeat_library: /mnt/data/test_data/Tandem_repeat_library.fasta
```

> **Note**: Use absolute paths for all files in the configuration file. The pipeline will
> create a directory named `Repeat_Annotations` in the current working directory to store
> the output files.

The pipeline supports supplying a additional custom repeat library via the
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
obtained via RepeatExplorer2. Library of rDNA sequences is built-in pipeline. 
