# wgd-nextflow

A Nextflow DSL2 pipeline for whole genome duplication (WGD), synteny-aware orthogroup construction, gene family evolution, tree cleaning, reconciliation, and rediploidisation analyses.

The pipeline integrates:

- Optional ab initio annotation with ANNEVO
- Genome preprocessing, including Tiberius annotation normalisation
- Primary transcript selection and repository ID standardisation
- Annotation parsing and validation for GENESPACE
- Optional external OrthoFinder, or internal OrthoFinder through GENESPACE
- GENESPACE synteny-aware orthogroup construction
- Orthogroup filtering and optional tandem duplicate collapse
- Alternative alignment and gene-tree branches:
  - MACSE coding-sequence alignment followed by nucleotide IQ-TREE inference
  - MAFFT amino-acid alignment followed by amino-acid IQ-TREE inference
- Optional gene-tree cleaning
- Optional AleRax gene tree/species tree reconciliation
- Optional rediploidisation classification, Circos input generation, Circos plotting, and species-tree plotting

---

## Table of contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Optional: pre-pull containers](#optional-pre-pull-containers)
- [Input files](#input-files)
  - [1. `genomes.tsv`](#1-genomestsv)
  - [2. Genome annotation and sequence files](#2-genome-annotation-and-sequence-files)
  - [3. Genome FASTA files for ANNEVO](#3-genome-fasta-files-for-annevo)
  - [4. Chromosome dictionary / chromosome-length files](#4-chromosome-dictionary--chromosome-length-files)
  - [5. Species tree](#5-species-tree)
- [Running the pipeline](#running-the-pipeline)
- [Start modes](#start-modes)
  - [`full`](#full)
  - [`parsed`](#parsed)
  - [`genespace`](#genespace)
  - [`alerax`](#alerax)
  - [`redip`](#redip)
  - [`redip_rooted`](#redip_rooted)
  - [`treeclean`](#treeclean)
- [Alignment and gene-tree branches](#alignment-and-gene-tree-branches)
  - [MACSE nucleotide branch](#macse-nucleotide-branch)
  - [MAFFT amino-acid branch](#mafft-amino-acid-branch)
- [Tree cleaning](#tree-cleaning)
- [Example parameter file](#example-parameter-file)
- [Optional ANNEVO annotation](#optional-annevo-annotation)
- [AleRax](#alerax-1)
- [Rediploidisation analysis](#rediploidisation-analysis)
  - [Position sources](#position-sources)
  - [Gene-tree rooting and classification](#gene-tree-rooting-and-classification)
  - [Already-rooted gene trees](#already-rooted-gene-trees)
  - [Rediploidisation species-tree plots](#rediploidisation-species-tree-plots)
- [Modular execution](#modular-execution)
- [Outputs](#outputs)
  - [Genome repository output](#genome-repository-output)
  - [GENESPACE output](#genespace-output)
  - [Post-processing output](#post-processing-output)
  - [Tree-cleaning output](#tree-cleaning-output)
  - [AleRax output](#alerax-output)
  - [Rediploidisation output](#rediploidisation-output)
- [Common issues](#common-issues)
- [Typical workflows](#typical-workflows)
- [Citation](#citation)

---

## Overview

This pipeline can perform:

1. Optional genome annotation with ANNEVO
2. Tiberius annotation normalisation, when `genome_source` is `tiberius`
3. Primary transcript selection
4. Chromosome naming standardisation using chromosome dictionaries or chromosome-length BED files
5. GENESPACE genome repository staging
6. GENESPACE annotation parsing and parse-output validation
7. Orthogroup inference with external OrthoFinder or GENESPACE-internal OrthoFinder
8. Synteny-aware orthogroup construction with GENESPACE
9. GENESPACE result validation
10. Orthogroup filtering, including optional outgroup requirements
11. Optional tandem duplicate collapse
12. Orthogroup FASTA writing for CDS and peptide sequences
13. One or both downstream alignment and tree-building branches:
    - CDS alignment with MACSE followed by nucleotide IQ-TREE analyses
    - peptide alignment with MAFFT followed by amino-acid IQ-TREE analyses
14. Optional gene-tree cleaning, either during a standard downstream run or as a standalone start mode
15. Optional gene tree/species tree reconciliation with AleRax, including multi-model runs
16. Optional rediploidisation analysis from unrooted IQ-TREE outputs or already-rooted gene trees
17. Optional Circos link, Circos input, Circos plot, and rediploidisation species-tree plot generation

The pipeline supports multiple entry points, so users can run the full workflow from genome inputs or restart from existing GENESPACE, IQ-TREE, tree-cleaning, AleRax, or rediploidisation inputs.

---

## Requirements

- Nextflow >= 23.10
- Apptainer, Singularity, or Docker
- HPC recommended for large datasets
- GPU access required only when running ANNEVO, depending on the selected profile and ANNEVO configuration

---

## Optional: pre-pull containers

To avoid repeated downloads and filesystem latency issues on HPC systems, containers can be pre-pulled:

```bash
bash utils/pull_containers.sh
```

If `.sif` files are present in the `apptainer/` directory, the pipeline will use them automatically.

Container selection order is:

1. Explicit path in `params.sif_paths`
2. Local SIF file in `params.apptainer_dir`
3. GHCR container image from `params.containers`

For example, the pipeline may look for local SIFs matching patterns such as:

```text
apptainer/genespace_*.sif
apptainer/macse_*.sif
apptainer/mafft_*.sif
apptainer/iqtree_*.sif
apptainer/alerax_*.sif
apptainer/annevo_*.sif
apptainer/gffutils_*.sif
apptainer/redip_*.sif
```

---

## Input files

### 1. `genomes.tsv`

A tab-separated file with required columns:

```text
genome_id    genome_source    ploidy
```

Required columns:

- `genome_id`: genome identifier; must match input file names
- `genome_source`: annotation/source label used by preprocessing and GENESPACE parsing
- `ploidy`: integer ploidy value

Common `genome_source` values include:

```text
ncbi
ensembl
phytozome
helixer
tiberius
```

When `genome_source` is `tiberius`, the pipeline runs the Tiberius normalisation step before primary transcript selection.

Optional columns used by downstream analyses:

- `outgroup`: outgroup status or outgroup tier used by rediploidisation rooting and, when enabled, outgroup-aware orthogroup filtering
- `redip`: marks species included as focal rediploidisation/WGD ingroup species

Example:

```text
genome_id	genome_source	ploidy	outgroup	redip
Eric_Calluna_vul	ensembl	2	0	1
Sarr_Sarracenia_pur	ncbi	2	0	1
Eben_Diospyros_lot	ncbi	2	0	1
Cornales_Cornus_flo	ncbi	2	1	0
Gunnerales_Myrothamnus_fla	phytozome	1	2	0
```

Providing an `outgroup` column allows the post-GENESPACE filter to require an outgroup gene copy when `require_outgroup_og: true`. Outgroup information is also used for rooting gene trees prior to rediploidisation analysis.

When integer outgroup tiers are used, the lowest tier present in a gene tree is tried first during rooting. This allows a preferred outgroup tier to be used before secondary outgroups.

Providing a `redip` column is required when running the rediploidisation subworkflow. Each species marked as rediploidisation-positive is treated as part of the WGD ingroup and as a focal species for classification.

---

### 2. Genome annotation and sequence files

For each genome, the default expected files are:

```text
gff/<genome_id>.gff3
protein/<genome_id>.pep
cds/<genome_id>.cds
```

The default extensions can be changed with:

```yaml
ext:
  gff: "gff3"
  pep: "pep"
  cds: "cds"
  fasta: "fa"
```

---

### 3. Genome FASTA files for ANNEVO

Genome FASTA files are required only when running the optional ANNEVO annotation workflow:

```yaml
annotation:
  run: true
  tool: "annevo"
```

When ANNEVO is enabled, the pipeline searches for genome FASTA files in `fasta_dir`.

Supported FASTA extensions include:

```text
<genome_id>.fa
<genome_id>.fa.gz
<genome_id>.fasta
<genome_id>.fasta.gz
<genome_id>.fna
<genome_id>.fna.gz
```

Example:

```text
fasta/Eric_Calluna_vul.fa.gz
fasta/Sarr_Sarracenia_pur.fa.gz
```

---

### 4. Chromosome dictionary / chromosome-length files

For each genome, provide either:

```text
chr_dict/<genome_id>.tsv
```

or:

```text
chr_dict/<genome_id>_chr_lengths.bed
```

These files are used to standardise or simplify chromosome names for downstream analyses and plotting.

If running rediploidisation analysis, chromosome-length BED files are required for focal redip species because they are used to prepare Circos plots.

#### Style A: two-column mapping

```text
old_seqid	new_seqid
NC_001234.1	chr1
NC_001235.1	chr2
```

#### Style B: BED-like chromosome-length file

```text
#chrom	chromStart	chromEnd	name
chr1	0	88077321	NC_001234.1
chr2	0	76011234	NC_001235.1
```

For this style, the mapping is interpreted as:

```text
old_seqid = column 4
new_seqid = column 1
```

When building chromosome-length BED files from NCBI assembly reports, the sequence identifier in column 4 must match the identifiers used in the GFF. For RefSeq-derived annotations this will usually be the RefSeq accession, for example `NC_...`, rather than the GenBank chromosome accession, for example `CM...`.

---

### 5. Species tree

A species tree in Newick format is required when using any of the following:

```yaml
use_species_tree_for_orthofinder: true
use_species_tree_for_alerax: true
run_rediploidisation: true
start_mode: "redip"
start_mode: "redip_rooted"
```

Example:

```yaml
species_tree: "path/to/species_tree.nwk"
```

---

## Running the pipeline

Basic run:

```bash
nextflow run main.nf \
  -profile <profile_name> \
  -params-file params.yaml \
  -with-apptainer \
  -output-dir results
```

Resume a previous run:

```bash
nextflow run main.nf \
  -profile <profile_name> \
  -params-file params.yaml \
  -with-apptainer \
  -output-dir results \
  -resume
```

Available profiles may include:

```text
docker
eddie
cropdiversity
```

---

## Start modes

The pipeline supports seven start modes:

| Mode | Description |
|---|---|
| `full` | Run the full pipeline from genome inputs |
| `parsed` | Start from a parsed GENESPACE working directory and run GENESPACE plus downstream steps |
| `genespace` | Start from a completed GENESPACE working directory and run downstream steps |
| `alerax` | Run only AleRax from existing IQ-TREE gene trees and matching nucleotide alignments |
| `redip` | Run only rediploidisation from existing unrooted IQ-TREE outputs |
| `redip_rooted` | Run only rediploidisation from already-rooted gene trees |
| `treeclean` | Run only tree cleaning from existing OG FASTAs, IQ-TREE outputs, and NT alignments, with optional AleRax/rediploidisation afterwards |

---

### `full`

Runs the complete workflow from genome annotations and sequence files.

This mode can optionally run ANNEVO before preprocessing if:

```yaml
annotation:
  run: true
  tool: "annevo"
```

If ANNEVO is disabled, the pipeline uses the supplied GFF, CDS, peptide, and chromosome dictionary files directly. Tiberius annotations are normalised automatically when `genome_source` is `tiberius`.

---

### `parsed`

Starts from an existing parsed GENESPACE working directory.

Requires:

```yaml
start_mode: "parsed"
existing_genespace_wd: "path/to/workingDir"
```

The pipeline creates a parse-completion marker, validates that the expected parsed inputs exist, then runs OrthoFinder or GENESPACE as configured.

---

### `genespace`

Starts from a completed GENESPACE working directory.

Requires:

```yaml
start_mode: "genespace"
existing_genespace_wd: "path/to/workingDir"
```

The pipeline validates that the directory contains completed GENESPACE output and then reruns downstream analyses only. This is useful for rerunning alignment, tree-building, tree cleaning, AleRax, or rediploidisation analyses without rerunning GENESPACE.

---

### `alerax`

Runs only the AleRax workflow from existing IQ-TREE gene trees and matching nucleotide alignments.

Requires:

```yaml
start_mode: "alerax"
run_alerax: true
```

If `use_species_tree_for_alerax: true`, also provide:

```yaml
species_tree: "path/to/species_tree.nwk"
```

Existing inputs are specified under `alerax`:

```yaml
alerax:
  gene_trees_dir: "path/to/iqtree_nt"
  nt_alignments_dir: "path/to/macse"
```

The `gene_trees_dir` should contain IQ-TREE directories named like:

```text
og_<orthogroup>/og_<orthogroup>_iqtree.treefile
```

The nucleotide alignment directory should contain files named like:

```text
og_<orthogroup>_NT.fasta
```

---

### `redip`

Runs only the rediploidisation workflow from existing unrooted IQ-TREE outputs.

Requires:

```yaml
start_mode: "redip"
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"
```

Also requires existing IQ-TREE outputs:

```yaml
rediploidisation:
  gene_trees_dir: "path/to/iqtree_aa_or_iqtree_nt"
```

The `gene_trees_dir` should contain directories named like:

```text
og_<orthogroup>/og_<orthogroup>_iqtree.treefile
og_<orthogroup>/og_<orthogroup>.iqtree.status
```

Only IQ-TREE directories with status `OK` are passed into the rediploidisation workflow.

If positions are taken from GENESPACE-derived BED files, provide:

```yaml
rediploidisation:
  genespace_wd: "path/to/workingDir"
```

If using a standalone positions file instead, set:

```yaml
rediploidisation:
  positions_source: "positions"
  positions: "path/to/positions.tsv"
```

---

### `redip_rooted`

Runs only the rediploidisation workflow from already-rooted gene trees. This mode skips the rooting step inside the rediploidisation subworkflow.

Requires:

```yaml
start_mode: "redip_rooted"
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"
```

Also provide a directory of rooted trees:

```yaml
rediploidisation:
  rooted_gene_trees_dir: "path/to/rooted_trees"
```

Rooted tree filenames must end in:

```text
.rooted.treefile
```

Examples:

```text
og_8819.rooted.treefile
og_8819_subtree_001.rooted.treefile
```

The complete filename prefix before `.rooted.treefile` is preserved as the tree ID.

As in `redip` mode, provide either a GENESPACE working directory for BED-derived positions:

```yaml
rediploidisation:
  genespace_wd: "path/to/workingDir"
```

or a standalone positions file:

```yaml
rediploidisation:
  positions_source: "positions"
  positions: "path/to/positions.tsv"
```

---

### `treeclean`

Runs the tree-cleaning workflow from existing orthogroup FASTAs, IQ-TREE directories, and nucleotide alignments.

Requires:

```yaml
start_mode: "treeclean"
run_tree_cleaning: true
```

Also provide:

```yaml
treeclean:
  og_fasta_dir: "path/to/og_fasta_cds"
  gene_trees_dir: "path/to/iqtree_nt"
  nt_alignments_dir: "path/to/macse"
```

Expected input patterns are:

```text
og_fasta_dir/og_<orthogroup>.fasta
gene_trees_dir/og_<orthogroup>/og_<orthogroup>_iqtree.treefile
gene_trees_dir/og_<orthogroup>/og_<orthogroup>.iqtree.status
nt_alignments_dir/og_<orthogroup>_NT.fasta
```

Only IQ-TREE directories with status `OK` are used.

`treeclean` mode can optionally run AleRax and/or rediploidisation after cleaning:

```yaml
run_alerax: true
run_rediploidisation: true
```

If rediploidisation is enabled after standalone tree cleaning, provide either:

```yaml
rediploidisation:
  genespace_wd: "path/to/workingDir"
```

or:

```yaml
rediploidisation:
  positions_source: "positions"
  positions: "path/to/positions.tsv"
```

---

## Alignment and gene-tree branches

The pipeline can run one or both downstream alignment and tree-building branches. Select the branch with:

```yaml
alignment_method: "macse_nt"
```

Allowed values are:

| Value | Alignment branch run | Gene trees inferred from |
|---|---|---|
| `macse_nt` | MACSE only | nucleotide alignments |
| `mafft_aa` | MAFFT only | amino-acid alignments |
| `both` | MACSE and MAFFT | both nucleotide and amino-acid alignments |

No unnecessary alignment branch is run unless it is explicitly requested.

Relevant parameters include:

```yaml
alignment_method: "macse_nt"
mafft_bin: "mafft"
mafft_opts: "--auto --amino"
iqtree_nt_model: "MFP"
iqtree_aa_model: "MFP"
```

### MACSE nucleotide branch

The MACSE branch:

1. Takes per-orthogroup CDS FASTAs
2. Writes amino-acid and nucleotide MACSE alignments
3. Infers nucleotide gene trees with IQ-TREE
4. Produces MACSE and IQ-TREE reports

Outputs are written under:

```text
post_genespace/macse/
post_genespace/iqtree_nt/
```

### MAFFT amino-acid branch

The MAFFT branch:

1. Takes per-orthogroup protein FASTAs
2. Writes amino-acid MAFFT alignments
3. Infers amino-acid gene trees with IQ-TREE
4. Produces MAFFT and IQ-TREE reports

Outputs are written under:

```text
post_genespace/mafft_aa/
post_genespace/iqtree_aa/
```

When both branches are run, AleRax and rediploidisation prefer the amino-acid IQ-TREE branch unless cleaned trees or an explicit tree directory are supplied.

---

## Tree cleaning

Tree cleaning can run in two ways.

### Inline tree cleaning during `full`, `parsed`, or `genespace` runs

Enable with:

```yaml
run_tree_cleaning: true
```

Inline tree cleaning currently requires the MACSE nucleotide branch:

```yaml
alignment_method: "macse_nt"
```

or:

```yaml
alignment_method: "both"
```

This is because the tree-cleaning workflow uses CDS FASTAs, nucleotide IQ-TREE outputs, and MACSE nucleotide alignments.

Downstream workflows can use cleaned trees with:

```yaml
use_cleaned_gene_trees_for_alerax: true
use_cleaned_gene_trees_for_redip: true
```

If these are false, downstream workflows use the standard branch-selection logic instead.

### Standalone tree cleaning

Use:

```yaml
start_mode: "treeclean"
run_tree_cleaning: true
```

and provide existing inputs under `treeclean` as described in the [`treeclean`](#treeclean) start-mode section.

---

## Example parameter file

```yaml
genomes_tsv: "path/to/genomes.tsv"

gff_dir: "path/to/gff"
protein_dir: "path/to/protein"
cds_dir: "path/to/cds"
chr_dict_dir: "path/to/chr_dict"
fasta_dir: "path/to/fasta"

working_dir: "workingDir"
postdir: "post_genespace"

start_mode: "full"
existing_genespace_wd: null
existing_orthofinder_dir: null

run_external_orthofinder: true
orthofinder_in_blk: true
orthofinder_analysis_threads: 4
force_orthofinder: false

species_tree: "path/to/species_tree.nwk"
use_species_tree_for_orthofinder: true
use_species_tree_for_alerax: true

require_outgroup_og: true
collapse_tandems: true
tandem_max_ord_gap: 2

alignment_method: "macse_nt"
mafft_bin: "mafft"
mafft_opts: "--auto --amino"
iqtree_nt_model: "MFP"
iqtree_aa_model: "MFP"

run_tree_cleaning: false
use_cleaned_gene_trees_for_alerax: false
use_cleaned_gene_trees_for_redip: false

run_alerax: true
run_rediploidisation: false

annotation:
  run: false
  tool: "annevo"
```

---

## Optional ANNEVO annotation

ANNEVO can be run before the standard preprocessing steps. ANNEVO is disabled by default.

Enable ANNEVO:

```yaml
annotation:
  run: true
  tool: "annevo"
```

Example ANNEVO configuration:

```yaml
fasta_dir: "path/to/fasta"

annevo:
  model_path: "/path/to/ANNEVO_Magnoliopsida.pt"
  batch_size: null
  extra_args: ""
  cleanup_tmp: true
```

When ANNEVO is enabled, the pipeline uses genome FASTA files from `fasta_dir`, runs ANNEVO, converts the resulting GFF to FASTA-compatible outputs, and passes the ANNEVO-derived GFF, peptide, and CDS outputs into the rest of the pipeline.

If ANNEVO is disabled, the pipeline uses the provided GFF, protein, and CDS files directly, with Tiberius normalisation applied when appropriate.

---

## AleRax

AleRax is optional and disabled by default.

Enable AleRax:

```yaml
run_alerax: true
```

Disable AleRax:

```yaml
run_alerax: false
```

AleRax can run one or more reconciliation models. If no model list is provided, the pipeline falls back to the single-model settings:

```yaml
alerax:
  rec_model: "UndatedDL"
  model_parametrization: "PER-SPECIES"
  gene_tree_samples: 100
```

For multi-model runs:

```yaml
alerax:
  cleanup_output: true
  models:
    - model_id: "DL_global"
      rec_model: "UndatedDL"
      model_parametrization: "GLOBAL"
      gene_tree_samples: 100

    - model_id: "DL_perSpecies"
      rec_model: "UndatedDL"
      model_parametrization: "PER-SPECIES"
      gene_tree_samples: 100

    - model_id: "DTL_global"
      rec_model: "UndatedDTL"
      model_parametrization: "GLOBAL"
      gene_tree_samples: 100

    - model_id: "DTL_perSpecies"
      rec_model: "UndatedDTL"
      model_parametrization: "PER-SPECIES"
      gene_tree_samples: 100
```

Notes:

- Each model runs independently.
- Outputs are separated by `model_id`.
- `cleanup_output: true` removes large intermediate AleRax output directories that are not usually needed for downstream interpretation.
- In a standard run, AleRax uses cleaned trees when `run_tree_cleaning: true` and `use_cleaned_gene_trees_for_alerax: true`.
- Otherwise, AleRax uses the amino-acid IQ-TREE branch when `alignment_method: "mafft_aa"` or `alignment_method: "both"`; if no amino-acid branch is available, it uses the nucleotide IQ-TREE branch.
- In `start_mode: "alerax"`, AleRax requires existing IQ-TREE gene trees and matching nucleotide alignments.

---

## Rediploidisation analysis

Rediploidisation analysis can be run as part of the main workflow or as a standalone `redip` or `redip_rooted` start mode.

Enable rediploidisation during a full or downstream run:

```yaml
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"
```

Default configuration:

```yaml
rediploidisation:
  positions_source: "bed"
  positions: ""
  position_format: "auto"

  position_has_header: true
  position_key_column: "id"
  position_chr_column: "chr"
  position_start_column: "start"
  position_end_column: "end"
  position_species_column: "genome"

  gene_trees_dir: ""
  rooted_gene_trees_dir: ""
  genespace_wd: ""

  species_tree_format: 1
  gene_tree_format: 1
  tip_separator: "|"
  label_format: "species_chr_gene"

  copy_mode: "target_exactly_n"
  required_copies: 2
  min_tips: 1
  position_key_type: "gene"

  classify_mode: "recurrent"
  recent_grouping: "auto"
  min_recent_groups: 2
  min_ancestral_target_copies: null
  write_lossy_singletons: true

  species_tree_annotation: ""
  branch_colors: ""
  species_tree_plot_width: 8
  species_tree_plot_height: 6
  species_tree_plot_dpi: 300
  species_tree_circle_size: 7
  species_tree_circle_stroke: 0.5
  species_tree_branch_number_size: 4
  species_tree_tip_label_size: 5
  species_tree_tip_label_offset: 1
  species_tree_branch_label_x_offset: 0.2
  species_tree_branch_label_y_offset: 0.1
  species_tree_prune_to_branch_species: true

  cleanup_tmp: true
```

### Position sources

The rediploidisation workflow can use positions from GENESPACE-derived BED files or from a user-supplied positions file.

If using GENESPACE-derived positions, use:

```yaml
rediploidisation:
  positions_source: "bed"
  genespace_wd: "path/to/workingDir"
```

During a standard `full`, `parsed`, or `genespace` run, `genespace_wd` can be omitted and the active GENESPACE working directory is used automatically.

If using a standalone positions file:

```yaml
rediploidisation:
  positions_source: "positions"
  positions: "path/to/positions.tsv"
```

The positions file should include columns for gene ID, chromosome, start, end, and optionally species/genome.

Default column names are:

```yaml
position_key_column: "id"
position_chr_column: "chr"
position_start_column: "start"
position_end_column: "end"
position_species_column: "genome"
```

### Gene-tree rooting and classification

When input trees are unrooted, rediploidisation gene trees are rooted using outgroup information from `genomes.tsv`.

- If a single outgroup tip from the preferred available tier is present, the tree is rooted on that tip.
- If multiple outgroup tips from the preferred available tier are present and form a clade, the tree is rooted on their MRCA.
- If the MRCA of the chosen outgroup tier is already the entire tree, the tree is recorded as skipped rather than rewritten unchanged.
- Only trees successfully recorded as rooted are passed into classification.

The rooting report includes all attempted trees, including both rooted and skipped cases.

Classification is controlled by parameters including:

```yaml
copy_mode: "target_exactly_n"
required_copies: 2
classify_mode: "recurrent"
recent_grouping: "auto"
min_recent_groups: 2
min_ancestral_target_copies: null
write_lossy_singletons: true
```

### Already-rooted gene trees

Use `start_mode: "redip_rooted"` when trees have already been rooted externally or by a previous pipeline run.

```yaml
start_mode: "redip_rooted"
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"

rediploidisation:
  rooted_gene_trees_dir: "path/to/rooted_trees"
  genespace_wd: "path/to/workingDir"
```

Tree filenames must end in `.rooted.treefile`.

### Rediploidisation species-tree plots

The rediploidisation workflow can generate species-tree plots using branch definitions and optional annotation/colour files.

Relevant plot parameters include:

```yaml
rediploidisation:
  species_tree_annotation: "path/to/annotation.tsv"
  branch_colors: "path/to/branch_colors.tsv"
  species_tree_plot_width: 8
  species_tree_plot_height: 6
  species_tree_plot_dpi: 300
  species_tree_circle_size: 7
  species_tree_circle_stroke: 0.5
  species_tree_branch_number_size: 4
  species_tree_tip_label_size: 5
  species_tree_tip_label_offset: 1
  species_tree_branch_label_x_offset: 0.2
  species_tree_branch_label_y_offset: 0.1
  species_tree_prune_to_branch_species: true
```

---

## Modular execution

| Mode | Input required | Typical use |
|---|---|---|
| `full` | Raw genome annotations, sequences, and chromosome dictionaries | Complete analysis |
| `parsed` | Parsed GENESPACE working directory | Restart after parsing |
| `genespace` | Completed GENESPACE working directory | Rerun downstream analyses |
| `alerax` | Existing IQ-TREE directories and nucleotide alignments | AleRax-only analyses |
| `redip` | Existing unrooted IQ-TREE directories, species tree, and positions | Rediploidisation-only analyses with internal rooting |
| `redip_rooted` | Existing rooted gene trees, species tree, and positions | Rediploidisation-only analyses without rerooting |
| `treeclean` | Existing OG FASTAs, IQ-TREE directories, and NT alignments | Tree-cleaning-only analyses, optionally followed by AleRax and/or rediploidisation |

This allows users to:

- Avoid rerunning expensive upstream steps
- Test different alignment and tree-building strategies
- Test tree cleaning before reconciliation or rediploidisation
- Test different AleRax model configurations
- Run rediploidisation on existing unrooted or already-rooted gene trees
- Reuse completed GENESPACE analyses
- Iterate on downstream classification and plotting parameters

---

## Outputs

### Genome repository output

When `start_mode: "full"`, the staged genome repository is published as:

```text
genomeRepo/
```

Typical contents include standardised GFF, peptide, and chromosome dictionary files used by GENESPACE.

### GENESPACE output

Published as:

```text
workingDir/
```

Typical contents include:

```text
workingDir/
  bed/
  dotplots/
  orthofinder/
  pangenes/
  peptide/
  results/
  riparian/
  syntenicHits/
  tmp/
```

### Post-processing output

Published as:

```text
post_genespace/
```

Typical contents include:

```text
post_genespace/
  pangenes/
    pangenes_PASS.tsv
    og_list_min4species.txt
    pangenes_PASS.collapsed.tsv
    og_list_min4species.collapsed.txt
  tandem_collapse/
    tandem_report.tsv
  og_fasta/
    og_fastas.done
    write_og_fastas.log
  og_fasta_cds/
    og_<orthogroup>.fasta
  og_fasta_aa/
    og_<orthogroup>.fasta
```

If the MACSE nucleotide branch is run:

```text
post_genespace/
  macse/
    og_<orthogroup>_AA.fasta
    og_<orthogroup>_NT.fasta
    og_<orthogroup>.status
    og_<orthogroup>.log
    macse_report.tsv
    macse_ok_og_list.txt
  iqtree_nt/
    og_<orthogroup>/
      og_<orthogroup>_iqtree.treefile
      og_<orthogroup>_iqtree.ufboot
      og_<orthogroup>.iqtree.status
    iqtree_nt_report.tsv
    iqtree_nt_ok_og_list.txt
```

If the MAFFT amino-acid branch is run:

```text
post_genespace/
  mafft_aa/
    og_<orthogroup>_AA.fasta
    og_<orthogroup>.status
    og_<orthogroup>.log
    mafft_report.tsv
    mafft_ok_og_list.txt
  iqtree_aa/
    og_<orthogroup>/
      og_<orthogroup>_iqtree.treefile
      og_<orthogroup>_iqtree.ufboot
      og_<orthogroup>.iqtree.status
    iqtree_aa_report.tsv
    iqtree_aa_ok_og_list.txt
```

### Tree-cleaning output

If tree cleaning is enabled, typical outputs include:

```text
post_genespace/
  treeclean/
    treeclean_report.tsv
    og_<orthogroup>/
      cleaned outputs for downstream AleRax and rediploidisation
```

Exact contents depend on the local `TREECLEAN` module implementation, but the workflow publishes both a report and cleaned per-orthogroup directories.

### AleRax output

If AleRax is enabled:

```text
post_genespace/
  alerax/
    families.txt
    manifest.tsv
    <model_id>/
      output/
      alerax.log
      alerax.done
  alerax_report.tsv
```

For multi-model runs:

```text
post_genespace/
  alerax/
    DL_global/
    DL_perSpecies/
    DTL_global/
    DTL_perSpecies/
  alerax_report.tsv
```

### Rediploidisation output

If rediploidisation is enabled, outputs are published under `post_genespace/`.

Typical outputs include:

```text
post_genespace/
  rediploidisation/
    rooted_trees/
    rooting_summaries/
    branch_definitions/
    classifications/
    circos_links/
    circos_inputs/
    circos_plots/
    species_trees/
    report/
      rooting_summary.tsv
      classification_summary.tsv
      circos_links_summary.tsv
```

`rooted_trees/` contains successfully rooted trees. `rooting_summaries/` records attempted trees and rooting outcomes. In `redip_rooted` mode, the supplied rooted trees are passed through with rooting skipped.

---

## Common issues

### `species_tree` is missing

If any of the following are true, a species tree must be provided:

```yaml
use_species_tree_for_orthofinder: true
use_species_tree_for_alerax: true
run_rediploidisation: true
start_mode: "redip"
start_mode: "redip_rooted"
```

Set:

```yaml
species_tree: "path/to/species_tree.nwk"
```

### OrthoFinder error: species tree is not rooted

Ensure that the species tree is rooted before using it with OrthoFinder.

### GENESPACE reruns OrthoFinder unexpectedly

Check:

```yaml
existing_orthofinder_dir: "path/to/orthofinder"
run_external_orthofinder: true
```

Also confirm that the expected OrthoFinder results are correctly staged into the GENESPACE working directory.

If `use_species_tree_for_orthofinder: true`, the pipeline treats external OrthoFinder as required because a species tree is being supplied to OrthoFinder.

### Missing input files

Check that file names match `genome_id` exactly.

For example, if `genomes.tsv` contains:

```text
Eric_Calluna_vul
```

then default expected files are:

```text
gff_dir/Eric_Calluna_vul.gff3
protein_dir/Eric_Calluna_vul.pep
cds_dir/Eric_Calluna_vul.cds
chr_dict_dir/Eric_Calluna_vul.tsv
```

or, for a chromosome-length BED file:

```text
chr_dict_dir/Eric_Calluna_vul_chr_lengths.bed
```

If ANNEVO is enabled, the pipeline also expects a genome FASTA such as:

```text
fasta_dir/Eric_Calluna_vul.fa
```

or:

```text
fasta_dir/Eric_Calluna_vul.fa.gz
```

### Unsupported `alignment_method`

Allowed values are:

```text
macse_nt
mafft_aa
both
```

### Tree cleaning requires the MACSE nucleotide branch

Inline tree cleaning requires:

```yaml
run_tree_cleaning: true
alignment_method: "macse_nt"
```

or:

```yaml
run_tree_cleaning: true
alignment_method: "both"
```

This is because tree cleaning currently uses CDS FASTAs, nucleotide IQ-TREE results, and MACSE nucleotide alignments.

### Standalone treeclean mode is missing inputs

When using:

```yaml
start_mode: "treeclean"
```

provide:

```yaml
run_tree_cleaning: true

treeclean:
  og_fasta_dir: "path/to/og_fasta_cds"
  gene_trees_dir: "path/to/iqtree_nt"
  nt_alignments_dir: "path/to/macse"
```

### MAFFT rejects uncommon amino-acid symbols

Some protein FASTAs may contain symbols that MAFFT does not accept under the default amino-acid alphabet, for example `U`.

If this occurs, either clean or standardise the protein FASTAs before alignment, or adjust the MAFFT options deliberately if unusual symbols should be retained.

### Rediploidisation-only mode requires IQ-TREE directories

When using:

```yaml
start_mode: "redip"
```

provide:

```yaml
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"

rediploidisation:
  gene_trees_dir: "path/to/iqtree_aa_or_iqtree_nt"
```

The tree directory should contain IQ-TREE output directories named like:

```text
og_<orthogroup>/og_<orthogroup>_iqtree.treefile
og_<orthogroup>/og_<orthogroup>.iqtree.status
```

Only IQ-TREE outputs with status `OK` are used.

### `redip_rooted` mode requires `.rooted.treefile` names

When using:

```yaml
start_mode: "redip_rooted"
```

provide:

```yaml
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"

rediploidisation:
  rooted_gene_trees_dir: "path/to/rooted_trees"
```

Each tree file must end with:

```text
.rooted.treefile
```

### Rediploidisation positions are missing

If using:

```yaml
rediploidisation:
  positions_source: "positions"
```

then a positions file must be provided:

```yaml
rediploidisation:
  positions: "path/to/positions.tsv"
```

If using GENESPACE-derived positions in standalone `redip`, `redip_rooted`, or `treeclean` mode, provide:

```yaml
rediploidisation:
  genespace_wd: "path/to/workingDir"
```

### Rediploidisation produces zero classifications unexpectedly

Check that:

1. Rooted trees were successfully produced or supplied
2. Tree-tip labels match the configured `tip_separator` and `label_format`
3. Chromosome identifiers in tree tips match the chromosome identifiers in the positional inputs used downstream
4. The chosen `copy_mode`, `required_copies`, and recurrent-classification settings match the biological pattern you want to detect
5. The `redip` column in `genomes.tsv` marks the intended focal species

---

## Typical workflows

### Full WGD and gene family evolution analysis using nucleotide trees

```yaml
start_mode: "full"
alignment_method: "macse_nt"
run_alerax: true
run_rediploidisation: false
```

### Full run using amino-acid trees for deeply diverged groups

```yaml
start_mode: "full"
alignment_method: "mafft_aa"
run_alerax: true
run_rediploidisation: true
```

### Full run with both alignment/tree-building branches

```yaml
start_mode: "full"
alignment_method: "both"
run_alerax: true
run_rediploidisation: true
```

### Full run with ANNEVO annotation

```yaml
start_mode: "full"

annotation:
  run: true
  tool: "annevo"

annevo:
  model_path: "/path/to/ANNEVO_Magnoliopsida.pt"
```

### Rerun downstream steps from GENESPACE

```yaml
start_mode: "genespace"
existing_genespace_wd: "path/to/workingDir"

alignment_method: "mafft_aa"
run_alerax: true
run_rediploidisation: true
```

### Run tree cleaning inline and use cleaned trees for downstream analyses

```yaml
start_mode: "genespace"
existing_genespace_wd: "path/to/workingDir"

alignment_method: "macse_nt"
run_tree_cleaning: true
use_cleaned_gene_trees_for_alerax: true
use_cleaned_gene_trees_for_redip: true

run_alerax: true
run_rediploidisation: true
```

### Tree-cleaning-only analysis from existing outputs

```yaml
start_mode: "treeclean"
run_tree_cleaning: true

treeclean:
  og_fasta_dir: "path/to/og_fasta_cds"
  gene_trees_dir: "path/to/iqtree_nt"
  nt_alignments_dir: "path/to/macse"
```

### Tree cleaning followed by rediploidisation

```yaml
start_mode: "treeclean"
run_tree_cleaning: true
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"

treeclean:
  og_fasta_dir: "path/to/og_fasta_cds"
  gene_trees_dir: "path/to/iqtree_nt"
  nt_alignments_dir: "path/to/macse"

rediploidisation:
  genespace_wd: "path/to/workingDir"
```

### AleRax-only analysis

```yaml
start_mode: "alerax"
run_alerax: true
species_tree: "path/to/species_tree.nwk"

alerax:
  gene_trees_dir: "path/to/iqtree_nt"
  nt_alignments_dir: "path/to/macse"
```

### Rediploidisation-only analysis from unrooted IQ-TREE outputs

```yaml
start_mode: "redip"
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"
chr_dict_dir: "path/to/chr_dict"

rediploidisation:
  gene_trees_dir: "path/to/iqtree_aa"
  genespace_wd: "path/to/workingDir"
```

### Rediploidisation-only analysis from already-rooted gene trees

```yaml
start_mode: "redip_rooted"
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"
chr_dict_dir: "path/to/chr_dict"

rediploidisation:
  rooted_gene_trees_dir: "path/to/rooted_trees"
  genespace_wd: "path/to/workingDir"
```

### Rediploidisation using a standalone positions file

```yaml
start_mode: "redip_rooted"
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"

rediploidisation:
  rooted_gene_trees_dir: "path/to/rooted_trees"
  positions_source: "positions"
  positions: "path/to/positions.tsv"
```

---

## Citation

Please cite the tools used in your analysis, including:

- Nextflow
- GENESPACE
- OrthoFinder
- MACSE, if the nucleotide alignment branch was run
- MAFFT, if the amino-acid alignment branch was run
- IQ-TREE
- AleRax, if reconciliation was run
- ANNEVO, if annotation was run
- Circos, if rediploidisation plots were generated

