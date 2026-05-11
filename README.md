# wgd-nextflow

A Nextflow DSL2 pipeline for whole genome duplication (WGD), synteny-aware orthogroup construction, gene family evolution, and rediploidisation analyses (developed by Drew A. Larson).

The pipeline integrates:

- Optional ab initio annotation with ANNEVO
- Genome preprocessing and annotation parsing
- OrthoFinder
- GENESPACE
- Alternative alignment branches:
  - MACSE coding-sequence alignment followed by nucleotide gene-tree inference
  - MAFFT amino-acid alignment followed by amino-acid gene-tree inference
- IQ-TREE
- AleRax
- Rediploidisation classification and Circos visualisation

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
- [Alignment and gene-tree branches](#alignment-and-gene-tree-branches)
  - [MACSE nucleotide branch](#macse-nucleotide-branch)
  - [MAFFT amino-acid branch](#mafft-amino-acid-branch)
- [Example parameter file](#example-parameter-file)
- [Optional ANNEVO annotation](#optional-annevo-annotation)
- [AleRax](#alerax-1)
- [Rediploidisation analysis](#rediploidisation-analysis)
  - [Position sources](#position-sources)
  - [Gene-tree rooting and classification](#gene-tree-rooting-and-classification)
  - [Gene tree input](#gene-tree-input)
- [Modular execution](#modular-execution)
- [Outputs](#outputs)
  - [GENESPACE output](#genespace-output)
  - [Post-processing output](#post-processing-output)
  - [AleRax output](#alerax-output)
  - [Rediploidisation output](#rediploidisation-output)
- [Common issues](#common-issues)
  - [`species_tree` is missing](#species_tree-is-missing)
  - [OrthoFinder error: species tree is not rooted](#orthofinder-error-species-tree-is-not-rooted)
  - [GENESPACE reruns OrthoFinder unexpectedly](#genespace-reruns-orthofinder-unexpectedly)
  - [Missing input files](#missing-input-files)
  - [MAFFT rejects uncommon amino-acid symbols](#mafft-rejects-uncommon-amino-acid-symbols)
  - [Rediploidisation-only mode requires IQ-TREE directories](#rediploidisation-only-mode-requires-iq-tree-directories)
  - [Rediploidisation positions are missing](#rediploidisation-positions-are-missing)
  - [Rediploidisation produces zero classifications unexpectedly](#rediploidisation-produces-zero-classifications-unexpectedly)
- [Typical workflows](#typical-workflows)
- [Citation](#citation)

---

## Overview

This pipeline can perform:

1. Optional genome annotation with ANNEVO
2. Genome preprocessing and primary transcript selection
3. Annotation parsing for GENESPACE
4. Orthogroup inference with OrthoFinder
5. Synteny-aware orthogroup construction with GENESPACE
6. Orthogroup filtering
7. Optional tandem duplicate collapse
8. Orthogroup FASTA writing for both coding sequences and proteins
9. One or both downstream alignment/tree-building branches:
   - coding-sequence alignment with MACSE followed by nucleotide IQ-TREE analyses
   - amino-acid alignment with MAFFT followed by amino-acid IQ-TREE analyses
10. Optional gene tree - species tree reconciliation with AleRax
11. Optional rediploidisation event classification
12. Optional Circos link and plot generation for rediploidisation results

The pipeline supports several entry points, so users can run the full workflow from raw inputs or restart from existing GENESPACE, IQ-TREE, or AleRax inputs.

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

For example, the pipeline will look for local SIFs matching patterns such as:

```text
apptainer/genespace_*.sif
apptainer/macse_*.sif
apptainer/mafft_*.sif
apptainer/iqtree_*.sif
apptainer/alerax_*.sif
apptainer/annevo_*.sif
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
- `genome_source`: one of `ncbi`, `ensembl`, `phytozome`, or `helixer`
- `ploidy`: integer, for example `2`

Optional columns used by downstream analyses:

- `outgroup`: may be boolean-like (`yes`, `no`, `true`, `false`) or integer-tiered; lower positive integers are treated as more basal outgroup tiers during rediploidisation rooting
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

When integer outgroup tiers are used, the lowest tier present in a gene tree is tried first during rooting. This allows a more basal outgroup to be preferred over a secondary outgroup when both are available.

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

A species tree in Newick format is required when using:

```yaml
use_species_tree_for_orthofinder: true
```

or:

```yaml
use_species_tree_for_alerax: true
```

or:

```yaml
run_rediploidisation: true
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

The pipeline supports five start modes:

| Mode | Description |
|---|---|
| `full` | Run the full pipeline from genome inputs |
| `parsed` | Start from parsed GENESPACE inputs |
| `genespace` | Start from a completed GENESPACE working directory |
| `alerax` | Run only the AleRax workflow from existing IQ-TREE and alignment outputs |
| `redip` | Run only the rediploidisation workflow from existing IQ-TREE outputs and position data |

---

### `full`

Runs the complete workflow from genome annotations and sequence files.

This mode can optionally run ANNEVO before preprocessing if:

```yaml
annotation:
  run: true
  tool: "annevo"
```

---

### `parsed`

Starts from an existing parsed GENESPACE working directory.

Requires:

```yaml
start_mode: "parsed"
existing_genespace_wd: "path/to/workingDir"
```

The pipeline creates a parse-completion marker and validates that the expected parsed inputs exist.

---

### `genespace`

Starts from a completed GENESPACE working directory.

Requires:

```yaml
start_mode: "genespace"
existing_genespace_wd: "path/to/workingDir"
```

The pipeline validates that the directory contains completed GENESPACE output and then reruns downstream analyses only.

This is useful for rerunning alignment, tree-building, AleRax, or rediploidisation analyses without rerunning GENESPACE.

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

---

### `redip`

Runs only the rediploidisation workflow.

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

## Alignment and gene-tree branches

The pipeline can run one or both downstream alignment/tree-building branches. Select the branch with:

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

1. takes per-orthogroup CDS FASTAs
2. writes amino-acid and nucleotide MACSE alignments
3. infers nucleotide gene trees with IQ-TREE

Outputs are written under:

```text
post_genespace/macse/
post_genespace/iqtree_nt/
```

### MAFFT amino-acid branch

The MAFFT branch:

1. takes per-orthogroup protein FASTAs
2. writes amino-acid MAFFT alignments
3. infers amino-acid gene trees with IQ-TREE

Outputs are written under:

```text
post_genespace/mafft_aa/
post_genespace/iqtree_aa/
```

When both branches are run, AleRax and rediploidisation currently prefer the amino-acid IQ-TREE branch unless an explicit tree directory is supplied in the relevant parameters.

---

## Example parameter file

```yaml
genomes_tsv: "path/to/genomes.tsv"

gff_dir: "path/to/gff"
protein_dir: "path/to/protein"
cds_dir: "path/to/cds"
chr_dict_dir: "path/to/chr_dict"

working_dir: "workingDir"
postdir: "post_genespace"

start_mode: "full"
existing_genespace_wd: null
existing_orthofinder_dir: null

run_external_orthofinder: true
orthofinder_analysis_threads: 4

species_tree: "path/to/species_tree.nwk"
use_species_tree_for_orthofinder: true
use_species_tree_for_alerax: true

require_outgroup_og: true
collapse_tandems: true

alignment_method: "macse_nt"
mafft_bin: "mafft"
mafft_opts: "--auto --amino"
iqtree_nt_model: "MFP"
iqtree_aa_model: "MFP"

run_alerax: true
run_rediploidisation: false
```

---

## Optional ANNEVO annotation

ANNEVO can be run before the standard preprocessing steps. ANNEVO is disabled by default.

ANNEVO is resource heavy and may require GPU access depending on the selected profile and configuration.

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

When ANNEVO is enabled, the pipeline uses genome FASTA files from `fasta_dir`, runs ANNEVO, converts the resulting GFF to FASTA-compatible files, and passes the ANNEVO-derived GFF, peptide, and CDS outputs into the rest of the pipeline.

If ANNEVO is disabled, the pipeline uses the provided GFF, protein, and CDS files directly.

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
- In a standard full or downstream run, AleRax uses the amino-acid IQ-TREE branch when `alignment_method: "mafft_aa"` or `alignment_method: "both"`; otherwise it uses the nucleotide IQ-TREE branch.

---

## Rediploidisation analysis

Rediploidisation analysis can be run as part of the main workflow or as a standalone `redip` start mode.

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
  genespace_wd: ""

  species_tree_format: 1
  gene_tree_format: 1
  tip_separator: "|"
  label_format: "species_chr_gene"

  copy_mode: "target_exactly_n"
  required_copies: 2
  min_tips: 1
  position_key_type: "gene"

  cleanup_tmp: true
```

### Position sources

The rediploidisation workflow can use positions from GENESPACE-derived files or from a user-supplied positions file.

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

Rediploidisation gene trees are rooted using outgroup information from `genomes.tsv`.

- If a single outgroup tip from the most basal available tier is present, the tree is rooted on that tip.
- If multiple outgroup tips from the most basal available tier are present and form a clade, the tree is rooted on their MRCA.
- If the MRCA of the chosen outgroup tier is already the entire tree, the tree is recorded as `SKIPPED` rather than rewritten unchanged.
- Only trees successfully recorded as `ROOTED` are passed into classification.

The rooting report includes all attempted trees, including both `ROOTED` and `SKIPPED` cases.

### Gene tree input

When running rediploidisation as part of the full pipeline, only IQ-TREE directories with internal IQ-TREE status `OK` are passed into the rediploidisation workflow.

When running `start_mode: "redip"`, provide a directory containing IQ-TREE outputs named like:

```text
og_<orthogroup>/og_<orthogroup>_iqtree.treefile
og_<orthogroup>/og_<orthogroup>.iqtree.status
```

Example:

```yaml
start_mode: "redip"
run_rediploidisation: true

rediploidisation:
  gene_trees_dir: "path/to/iqtree_aa"
  genespace_wd: "path/to/workingDir"
```

Only directories whose `og_<orthogroup>.iqtree.status` file contains `OK` are used.

---

## Modular execution

| Mode | Input required | Typical use |
|---|---|---|
| `full` | Raw genome annotations and sequences | Complete analysis |
| `parsed` | Parsed GENESPACE input directory | Restart after parsing |
| `genespace` | Completed GENESPACE working directory | Rerun downstream analyses |
| `alerax` | Existing IQ-TREE directories and nucleotide alignments | AleRax-only analyses |
| `redip` | Existing IQ-TREE directories, species tree, and positions | Rediploidisation-only analyses |

This allows users to:

- Avoid rerunning expensive upstream steps
- Test different alignment and tree-building strategies
- Test different AleRax model configurations
- Run rediploidisation on existing IQ-TREE outputs
- Reuse completed GENESPACE analyses
- Iterate on downstream classification and plotting parameters

---

## Outputs

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

---

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

---

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

---

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
    report/
      rooting_summary.tsv
      classification_summary.tsv
      circos_links_summary.tsv
```

`rooted_trees/` contains only trees that were successfully rooted. `rooting_summaries/` records all attempted trees, including skipped cases.

---

## Common issues

### `species_tree` is missing

If any of the following are true, a species tree must be provided:

```yaml
use_species_tree_for_orthofinder: true
use_species_tree_for_alerax: true
run_rediploidisation: true
```

Set:

```yaml
species_tree: "path/to/species_tree.nwk"
```

---

### OrthoFinder error: species tree is not rooted

Ensure that the species tree is rooted before using it with OrthoFinder.

---

### GENESPACE reruns OrthoFinder unexpectedly

Check:

```yaml
existing_orthofinder_dir: "path/to/orthofinder"
run_external_orthofinder: true
```

Also confirm that the expected OrthoFinder results are correctly staged into the GENESPACE working directory.

---

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

---

### MAFFT rejects uncommon amino-acid symbols

Some protein FASTAs may contain symbols that MAFFT does not accept under the default amino-acid alphabet, for example `U`.

If this occurs, either clean or standardise the protein FASTAs before alignment, or adjust the MAFFT options deliberately if unusual symbols should be retained.

---

### Rediploidisation-only mode requires IQ-TREE directories

When using:

```yaml
start_mode: "redip"
```

you must provide:

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

---

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

If using GENESPACE-derived positions, provide:

```yaml
rediploidisation:
  genespace_wd: "path/to/workingDir"
```

This is required in `start_mode: "redip"` unless `positions_source` is set to `positions`.

---

### Rediploidisation produces zero classifications unexpectedly

Check that:

1. rooted trees were successfully produced and passed into classification
2. tree-tip labels match the configured `tip_separator` and `label_format`
3. chromosome identifiers in tree tips match the chromosome identifiers in the positional inputs used downstream
4. the chosen `copy_mode` and `required_copies` match the biological pattern you want to detect

---

## Typical workflows

### Full WGD and gene family evolution analysis using nucleotide trees

```yaml
start_mode: "full"
alignment_method: "macse_nt"
run_alerax: true
run_rediploidisation: false
```

---

### Full run using amino-acid trees for deeply diverged groups

```yaml
start_mode: "full"
alignment_method: "mafft_aa"
run_alerax: true
run_rediploidisation: true
```

---

### Full run with both alignment/tree-building branches

```yaml
start_mode: "full"
alignment_method: "both"
run_alerax: true
run_rediploidisation: true
```

---

### Full run with ANNEVO annotation

```yaml
start_mode: "full"

annotation:
  run: true
  tool: "annevo"

annevo:
  model_path: "/path/to/ANNEVO_Magnoliopsida.pt"
```

---

### Rerun downstream steps from GENESPACE

```yaml
start_mode: "genespace"
existing_genespace_wd: "path/to/workingDir"

alignment_method: "mafft_aa"
run_alerax: true
run_rediploidisation: true
```

---

### AleRax-only analysis

```yaml
start_mode: "alerax"
run_alerax: true
species_tree: "path/to/species_tree.nwk"

alerax:
  gene_trees_dir: "path/to/iqtree_nt"
  nt_alignments_dir: "path/to/macse"
```

---

### Rediploidisation-only analysis

```yaml
start_mode: "redip"
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"
chr_dict_dir: "path/to/chr_dict"

rediploidisation:
  gene_trees_dir: "path/to/iqtree_aa"
  genespace_wd: "path/to/workingDir"
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

