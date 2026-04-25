# wgd-nextflow

A Nextflow pipeline for detecting whole genome duplication (WGD) and downstream gene family evolution using GENESPACE, OrthoFinder, MACSE, IQ-TREE, and AleRax.

---

## Overview

This pipeline performs:

1. Genome preprocessing and annotation parsing  
2. Orthogroup inference (OrthoFinder; optional)  
3. Synteny-aware orthogroup construction (GENESPACE)  
4. Orthogroup filtering  
5. Optional tandem duplicate collapse  
6. Per-orthogroup alignment (MACSE)  
7. Gene tree inference (IQ-TREE)  
8. Gene tree - species tree reconciliation (AleRax; optional, multi-model)  

---

## Requirements

- Nextflow (≥ 23.10)
- Apptainer / Singularity or Docker
- HPC recommended for large datasets

---

## Optional: Pre-pull containers (recommended for HPC)

To avoid repeated downloads and filesystem latency issues, you can pre-pull all containers:

```
  bash utils/pull_containers.sh
```

Then point the pipeline to these local images usingsif_paths within your params file:

```
  sif_paths:
    genespace: "apptainer/genespace_1.1.sif"
    macse: "apptainer/macse_1.0.sif"
    iqtree: "apptainer/iqtree_1.0.sif"
    alerax: "apptainer/alerax_1.0.sif"
    annevo: "apptainer/annevo_1.0.sif"
    gffutils: "apptainer/gffutils_1.0.sif"
```

---

## Input files

### 1. genomes.tsv

Tab-separated file with required columns:

```
genome_id    genome_source    ploidy    outgroup
```

- genome_id: must match file names
- genome_source: one of ncbi, ensembl, phytozome, helixer
- ploidy: integer (e.g. 2)
- outgroup: yes or no (optional; defaults to no)

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

Supports two styles:

Style A:  
```
old_seqid <TAB> new_seqid
```

Style B (BED-like):  
```
#chrom <TAB> start <TAB> end <TAB> name

#chrom <TAB> chromStart <TAB> chromEnd <TAB> name chr1 0 88077321 1 NC_... 0 9953 NC_...
```
For Style B, interpret mapping as: 
old_seqid = column 4 (name) 
new_seqid = column 1 (#chrom)

---

### 4. (Optional) species tree

Newick format:

```
species_tree.nwk
```

Required if using species tree with OrthoFinder or AleRax.

---

## Running the pipeline

Basic run:

```bash
nextflow run main.nf -profile <profile_name> -params-file params.yaml -with-apptainer
```

Resume:

```bash
nextflow run main.nf -resume
```

---

## Start modes

- full: run entire pipeline from raw inputs  
- parsed: start from parsed GENESPACE inputs  
- genespace: start from completed GENESPACE working directory  

genespace mode requires:

existing_genespace_wd: "path/to/workingDir"

The pipeline validates that this directory contains:
- results/
- results/combBed.txt
- pangenes/

---

## Key parameters

Example params.yaml:

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
use_species_tree_for_orthofinder: true
use_species_tree_for_alerax: true

require_outgroup_og: true
collapse_tandems: true

run_alerax: true
```

---

## AleRax (optional, multi-model)

Enable:

```yaml
run_alerax: true
```

Example configuration:

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
- Each model runs independently
- Outputs are separated by model_id
- cleanup_output removes large intermediate directories

Disable AleRax:

```
run_alerax: false
````

---

## Modular execution

The pipeline supports starting from intermediate stages:

| Mode       | Input required                        |
|-----------|--------------------------------------|
| full      | Raw genome inputs                    |
| parsed    | Parsed GENESPACE inputs              |
| genespace | Completed GENESPACE working directory|

This allows:
- rapid re-running of downstream analyses
- testing different AleRax models
- skipping expensive upstream steps

---

## Outputs

GENESPACE output:

```
workingDir/
  ├── bed/
  ├── dotplots/
  ├── orthofinder/
  ├── pangenes/
  ├── peptide/
  ├── results/
  ├── riparian/
  ├── syntenicHits/
  ├── tmp/
```

Post-processing output:

```
post_genespace/
  ├── pangenes_PASS.tsv
  ├── og_list_min4species.txt
  ├── og_fasta/
  ├── macse_report.tsv
  ├── iqtree_report.tsv
```

AleRax output (if enabled):

```
post_genespace/
  ├── alerax/
  │   ├── DL_global/
  │   │   ├── output/
  │   │   ├── alerax.log
  │   │   └── alerax.done
  │   ├── DL_perSpecies/
  │   ├── DTL_global/
  │   └── DTL_perSpecies/
  ├── alerax_report.tsv
```

---

## Common issues

OrthoFinder error: species tree not rooted  
→ Ensure tree is rooted

GENESPACE rerunning OrthoFinder  
→ Ensure results are correctly staged into the working directory

Missing files  
→ Check filenames and extensions match genome_id

---

## Typical workflow

1. Run full pipeline  
2. Restart from GENESPACE if needed  
3. Iterate downstream analyses (e.g. AleRax models)  

---

## Citation

Please cite the following tools:

- GENESPACE  
- OrthoFinder  
- MACSE  
- IQ-TREE  
- AleRax
