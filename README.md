# wgd-nextflow

A Nextflow pipeline for detecting whole genome duplication (WGD) and downstream gene family evolution using GENESPACE, OrthoFinder, MACSE, IQ-TREE, and AleRax.

---

## Overview

This pipeline performs:

1. Genome preprocessing and annotation parsing  
2. Orthogroup inference (OrthoFinder)  
3. Synteny-aware orthogroup construction (GENESPACE)  
4. Orthogroup filtering  
5. Optional tandem duplicate collapse  
6. Per-orthogroup alignment (MACSE)  
7. Gene tree inference (IQ-TREE)  
8. Gene tree - species tree reconciliation (AleRax)  

---

## Requirements

- Nextflow (≥ 23.10)
- Apptainer / Singularity or Docker
- HPC recommended for large datasets

---

## Input files

### 1. `genomes.tsv`

Tab-separated file with required columns:

```
genome_id    genome_source    ploidy    outgroup
```

- `genome_id`: must match file names
- `genome_source`: one of `ncbi`, `ensembl`, `phytozome`, `helixer`
- `ploidy`: integer (e.g. 2)
- `outgroup`: `yes` or `no` (optional; defaults to `no`)

Example:

```
Eric_Calluna_vul    ncbi       2    no
Sarr_Sarracenia_pur ncbi       2    yes
Thea_Camellia_fas   ensembl    2    no
Eben_Diospyros_lot  ensembl    2    no
```

---

### 2. Genome files

For each genome:

```
gff_dir/<genome_id>.gff3
protein_dir/<genome_id>.pep
cds_dir/<genome_id>.cds
```

---

### 3. Chromosome dictionary

For each genome:

```
chr_dict_dir/<genome_id>.tsv
```
 
Suports two styles: 
Style A (2+ columns):
old_seqid <TAB> new_seqid [<TAB> ...]
(header optional)

Style B (chr_lengths.bed / chrom sizes BED-ish):
#chrom <TAB> chromStart <TAB> chromEnd <TAB> name
chr1   0   88077321   1
NC_... 0   9953      NC_...

For Style B, interpret mapping as:
old_seqid = column 4 (name)
new_seqid = column 1 (#chrom)

---

### 4. (Optional) species tree

Newick format:

```
species_tree.nwk
```

- Required if using species tree with OrthoFinder
- Must be **rooted** for OrthoFinder

---

## Running the pipeline

### Basic run

```
nextflow run main.nf \
  -profile eddie \
  -params-file params.yaml \
  -with-apptainer
```

---

### Resume a failed run

```
nextflow run main.nf -resume
```

---

## Key parameters

Example `params.yaml`:

```yaml
genomes_tsv: "path/to/genomes.tsv"

gff_dir: "path/to/gff"
protein_dir: "path/to/protein"
cds_dir: "path/to/cds"
chr_dict_dir: "path/to/chr_dict"

working_dir: "workingDir"
postdir: "post_genespace"

start_mode: "full"
existing_genespace_wd: ""

run_external_orthofinder: true
existing_orthofinder_dir: null

species_tree: "path/to/species_tree.nwk"
use_species_tree_orthofinder: true
use_species_tree_alerax: true

require_outgroup_og: true
collapse_tandems: true
tandem_max_ord_gap: 2
```

---

## Pipeline options

### Start modes

- `full` : run entire pipeline from raw inputs  
- `parsed` : start from an existing GENESPACE working directory  

---

### OrthoFinder

- `run_external_orthofinder: true`  
  Runs OrthoFinder before GENESPACE  

- `existing_orthofinder_dir`  
  Use precomputed OrthoFinder results  

Both are automatically integrated into GENESPACE.

---

### Species tree

Single parameter:

```
species_tree: path/to/tree.nwk
```

Controlled by:

- `use_species_tree_orthofinder` → passed to OrthoFinder  
- `use_species_tree_alerax` → passed to AleRax  

If not used for AleRax, a random tree is generated.

---

### Filtering

- `require_outgroup_og: true`  
  Keep only orthogroups with ≥1 outgroup gene  

- `collapse_tandems: true`  
  Collapse tandem duplicates prior to alignment  

---

## Outputs

### GENESPACE output

```
workingDir/
  ├── pangenes/
  ├── results/
  ├── riparian/
```

---

### Post-processing output

```
post_genespace/
  ├── pangenes_PASS.tsv
  ├── og_list_min4species.txt
  ├── og_fasta/
  ├── macse_report.tsv
  ├── iqtree_report.tsv
  ├── alerax/
```

---

## Common issues

### OrthoFinder error: "Species tree is not rooted"

Ensure your tree is rooted before using it with OrthoFinder.

---

### GENESPACE re-running OrthoFinder

GENESPACE only detects existing results if they are located at:

```
workingDir/orthofinder/
```

This pipeline automatically stages results correctly.

---

### Missing files

Check that:

- genome IDs match across all inputs
- all required files are present
- file extensions match configuration

---

## Typical workflow

1. Prepare input files  
2. Set parameters  
3. Run pipeline  
4. Resume if needed  
5. Inspect outputs in `post_genespace/`  

---

## Citation

If you use this pipeline, please cite:

- GENESPACE  
- OrthoFinder  
- MACSE  
- IQ-TREE  
- AleRax  

(see respective tool documentation)
