# wgd-nextflow

A Nextflow DSL2 pipeline for whole genome duplication (WGD), synteny-aware orthogroup construction, gene family evolution, and rediploidisation analyses (developed by Drew A. Larson).

The pipeline integrates:

- Optional ab initio annotation with ANNEVO
- Genome preprocessing and annotation parsing
- OrthoFinder
- GENESPACE
- MACSE
- IQ-TREE
- AleRax
- Rediploidisation classification and Circos visualisation

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
8. Per-orthogroup coding sequence alignment with MACSE
9. Gene tree inference with IQ-TREE
10. Optional gene tree - species tree reconciliation with AleRax
11. Optional rediploidisation event classification
12. Optional Circos link and plot generation for rediploidisation results

The pipeline supports several entry points, so users can run the full workflow from raw inputs or restart from existing GENESPACE or gene-tree outputs.

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
genome_id    genome_source    ploidy    outgroup  rediploidisation
```

Required columns:

- `genome_id`: genome identifier; must match input file names
- `genome_source`: one of `ncbi`, `ensembl`, `phytozome`, or `helixer`
- `ploidy`: integer, for example `2`

Optional column:

- `outgroup`: `yes` or `no`; defaults to `no` if absent
- `rediploidisation`: `yes` or `no`

Example:

```text
genome_id	genome_source	ploidy	outgroup rediploidisation
Eric_Calluna_vul	ensembl 2 no  yes
Sarr_Sarracenia_pur	ncbi	2	no  yes
Eben_Diospyros_lotus  ncbi  2 no  yes
Cornales_Cornus_flo	ncbi	2	no  no
Gunnerales_Myrothamnus_fla  phytozome 1 yes no
```

Providing an outgroup column allows you to require an outgroup gene copy in your gene families when filtering post-GENESPACE. Outgroup information is also used for rooting gene trees prior to rediploidisation analysis. The outgroup column is in regard to phylogenetic outgroup, not the WGD. Only species within the most basal clade should have `yes` for `outgroup`

Providing a rediploidisation column is required when running the rediploidisation subworkflow. Each species with `yes` in the rediploidisation column will be treated as an in-group to the WGD and will be treated as a focal species in the rediploidisation analysis.

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

When ANNEVO is enabled, the pipeline searches for genome FASTA files in `fasta`.

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

### 4. Chromosome dictionary

For each genome, provide either:

```text
chr_dict/<genome_id>.tsv
```

or:

```text
chr_dict/<genome_id>_chr_lengths.bed
```

These files are used to standardise or simplify chromosome names for downstream analyses and plotting. If running rediploidisation analysis, Style B is required, as chromosome lengths are used to build circos plots. 

#### Style A: two-column mapping

```text
old_seqid	new_seqid
NC_001234.1	chr1
NC_001235.1	chr2
```

#### Style B: BED-like chromosome length file

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

The pipeline supports four start modes:

| Mode | Description |
|---|---|
| `full` | Run the full pipeline from genome inputs |
| `parsed` | Start from parsed GENESPACE inputs |
| `genespace` | Start from a completed GENESPACE working directory |
| `redip` | Run only the rediploidisation workflow from existing gene trees and position data |

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

The pipeline validates that the directory contains:

```text
results/
results/combBed.txt
pangenes/
```

This is useful for rerunning downstream analyses without rerunning GENESPACE.

---

### `redip`

Runs only the rediploidisation workflow.

Requires:

```yaml
start_mode: "redip"
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"
```

Also requires existing gene trees:

```yaml
rediploidisation:
  gene_trees_dir: "path/to/iqtree_dir"
```

If positions are taken from GENESPACE or BED-like GENESPACE-derived files, provide:

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

run_alerax: true
run_rediploidisation: false
```

---

## Optional ANNEVO annotation

ANNEVO can be run before the standard preprocessing steps. ANNEVO is disabled by default. 

ANNEVO is resource heavy and requires at least 1 GPU. 

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

AleRax is optional, and disabled by default.

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

### Gene tree input

When running rediploidisation as part of the full pipeline, IQ-TREE outputs are passed directly into the rediploidisation workflow.

When running `start_mode: "redip"`, provide a directory containing IQ-TREE outputs named like:

```text
og_<orthogroup>_iqtree.treefile
```

Example:

```yaml
start_mode: "redip"
run_rediploidisation: true

rediploidisation:
  gene_trees_dir: "path/to/iqtree"
  genespace_wd: "path/to/workingDir"
```

---

## Modular execution

| Mode | Input required | Typical use |
|---|---|---|
| `full` | Raw genome annotations and sequences | Complete analysis |
| `parsed` | Parsed GENESPACE input directory | Restart after parsing |
| `genespace` | Completed GENESPACE working directory | Rerun downstream analyses |
| `redip` | Existing gene trees, species tree, and positions | Rediploidisation-only analyses |

This allows users to:

- Avoid rerunning expensive upstream steps
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
  iqtree/
  macse/
  pangenes_PASS.tsv
  og_list_min4species.txt
  og_fasta/
  macse_report.tsv
  iqtree_report.tsv
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
    branch_definitions/
    classifications/
    circos_links/
    circos_plots/
    rediploidisation_report.tsv
```

Exact file names may vary depending on the number of focal species, gene trees, and plotting configuration.

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

If ANNEVO is enabled, the pipeline also expects a genome FASTA such as:

```text
fasta_dir/Eric_Calluna_vul.fa
```

or:

```text
fasta_dir/Eric_Calluna_vul.fa.gz
```

---

### Rediploidisation-only mode requires gene trees

When using:

```yaml
start_mode: "redip"
```

you must provide:

```yaml
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"

rediploidisation:
  gene_trees_dir: "path/to/iqtree"
```

The gene tree directory should contain files named like:

```text
og_<orthogroup>_iqtree.treefile
```

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

## Typical workflows

### Full WGD and gene family evolution analysis

```yaml
start_mode: "full"
run_alerax: true
run_rediploidisation: false
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

run_alerax: true
run_rediploidisation: true
```

---

### Rediploidisation-only analysis

```yaml
start_mode: "redip"
run_rediploidisation: true
species_tree: "path/to/species_tree.nwk"

rediploidisation:
  gene_trees_dir: "path/to/iqtree"
  genespace_wd: "path/to/workingDir"
```

---

## Citation

Please cite the tools used in your analysis, including:

- Nextflow
- GENESPACE
- OrthoFinder
- MACSE
- IQ-TREE
- AleRax
- ANNEVO, if annotation was run
- Circos, if rediploidisation plots were generated
