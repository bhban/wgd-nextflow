nextflow.enable.dsl=2

process EXTRACT_REDIP_SPECIES {
    tag "extract_redip_species"

    input:
    path genomes_tsv

    output:
    path "redip_species.txt"

    script:
    """
    python - <<'PY'
from pathlib import Path
import sys

sys.path.insert(0, "${projectDir}/scripts/rediploidisation")
from redip_utils import read_redip_species_from_genomes_tsv

species = read_redip_species_from_genomes_tsv("${genomes_tsv}")

with open("redip_species.txt", "w") as out:
    for sp in species:
        out.write(sp + "\\n")
PY
    """
}


process ROOT_GENE_TREE {
    tag { "og_${og}" }

    input:
    tuple val(og), path(treefile)
    path genomes_tsv

    output:
    tuple val(og),
          path("og_${og}.rooted.treefile"),
          path("og_${og}.rooting_summary.tsv")

    script:
    """
    python ${projectDir}/scripts/rediploidisation/root_gene_tree.py \\
        --tree ${treefile} \\
        --genomes-tsv ${genomes_tsv} \\
        --output-tree og_${og}.rooted.treefile \\
        --summary-tsv og_${og}.rooting_summary.tsv \\
        --tip-separator '${params.rediploidisation.tip_separator}' \\
        --tree-format ${params.rediploidisation.gene_tree_format}
    """
}


process WRITE_BRANCH_DEFS {
    tag "write_branch_defs"

    input:
    path species_tree
    path genomes_tsv

    output:
    path "branch_definitions"

    script:
    """
    mkdir -p branch_definitions

    python ${projectDir}/scripts/rediploidisation/write_branch_defs.py \\
        --tree ${species_tree} \\
        --genomes-tsv ${genomes_tsv} \\
        --output-dir branch_definitions \\
        --tree-format ${params.rediploidisation.species_tree_format}
    """
}


process CLASSIFY_REDIP_EVENTS {
    tag { species }

    input:
    tuple val(species), path(rooted_trees)
    path genomes_tsv

    output:
    tuple val(species), path("${species}.redip_classification.tsv")

    script:
    def tree_list = rooted_trees.collect { it.toString() }.join(' ')

    """
    : > ${species}.rooted_gene_trees.nwk

    for tree in ${tree_list}; do
        cat "\$tree" >> ${species}.rooted_gene_trees.nwk
        echo >> ${species}.rooted_gene_trees.nwk
    done

    python ${projectDir}/scripts/rediploidisation/classify.py \\
        --target-species ${species} \\
        --treefile ${species}.rooted_gene_trees.nwk \\
        --genomes-tsv ${genomes_tsv} \\
        --output ${species}.redip_classification.tsv \\
        --tip-separator '${params.rediploidisation.tip_separator}' \\
        --label-format ${params.rediploidisation.label_format} \\
        --copy-mode ${params.rediploidisation.copy_mode} \\
        --required-copies ${params.rediploidisation.required_copies} \\
        --min-tips ${params.rediploidisation.min_tips}
    """
}


process MAKE_REDIP_LINKS {
    tag { species }

    input:
    tuple val(species), path(classification)
    path branch_definitions
    path genespace_wd

    output:
    tuple val(species), path("${species}.circos_links.tsv")

    script:
    def source = params.rediploidisation.positions_source ?: 'bed'

    def position_arg
    if (source == 'pangenes') {
        position_arg = "--pangenes-dir ${genespace_wd}/pangenes"
    } else if (source == 'positions') {
        if (!params.rediploidisation.positions?.toString()?.trim()) {
            throw new IllegalArgumentException(
                "params.rediploidisation.positions must be provided when positions_source = 'positions'"
            )
        }

        def species_col_arg = params.rediploidisation.position_species_column?.toString()?.trim()
            ? "--position-species-column ${params.rediploidisation.position_species_column}"
            : ""

        position_arg = """
            --positions ${params.rediploidisation.positions}
            --position-format ${params.rediploidisation.position_format}
            ${params.rediploidisation.position_has_header ? '--position-has-header' : '--no-position-has-header'}
            --position-key-column ${params.rediploidisation.position_key_column}
            --position-chr-column ${params.rediploidisation.position_chr_column}
            --position-start-column ${params.rediploidisation.position_start_column}
            --position-end-column ${params.rediploidisation.position_end_column}
            ${species_col_arg}
        """
    } else {
        position_arg = "--bed-dir ${genespace_wd}/bed"
    }

    """
    python ${projectDir}/scripts/rediploidisation/make_links.py \\
        --species ${species} \\
        --classification-tsv ${classification} \\
        ${position_arg} \\
        --branch-definitions ${branch_definitions}/${species}.branch_definitions.tsv \\
        --output ${species}.circos_links.tsv \\
        --tip-separator '${params.rediploidisation.tip_separator}' \\
        --label-format ${params.rediploidisation.label_format} \\
        --position-key-type ${params.rediploidisation.position_key_type} \\
        --write-header \\
        --include-metadata \\
        --on-exists overwrite \\
        --log-level INFO
    """
}


process PREP_REDIP_CIRCOS {
    tag { species }

    input:
    tuple val(species), path(circos_links)

    output:
    tuple val(species), path("${species}")

    script:
    """
    python ${projectDir}/scripts/rediploidisation/prep_circos.py \\
        --species ${species} \\
        --circos-links ${circos_links} \\
        --chr-bed ${params.chr_dict_dir}/${species}_chr_lengths.bed \\
        --output-dir ${species}
    """
}


process PLOT_REDIP_CIRCOS {
    tag { species }

    input:
    tuple val(species), path(species_dir)

    output:
    tuple val(species), path("${species}")

    script:
    """
    cp -r ${species_dir} ${species}

    cd ${species}

    ${params.circos_bin} -conf circos.conf
    """
}


process REDIP_REPORT {
    tag "redip_report"

    input:
    path rooting_summaries
    path classifications
    path circos_links

    output:
    path "redip_report"

    script:
    def root_list = rooting_summaries.collect { it.toString() }.join(' ')
    def class_list = classifications.collect { it.toString() }.join(' ')
    def link_list = circos_links.collect { it.toString() }.join(' ')

    """
    mkdir -p redip_report

    {
        first=1
        for f in ${root_list}; do
            if [[ "\$first" -eq 1 ]]; then
                cat "\$f"
                first=0
            else
                tail -n +2 "\$f"
            fi
        done
    } > redip_report/rooting_summary.tsv

    {
        echo -e "species\\tclassification_rows"
        for f in ${class_list}; do
            species=\$(basename "\$f" .redip_classification.tsv)
            n=\$(awk 'NR > 1 { count++ } END { print count + 0 }' "\$f")
            echo -e "\${species}\\t\${n}"
        done
    } > redip_report/classification_summary.tsv

    {
        echo -e "species\\tcircos_link_rows"
        for f in ${link_list}; do
            species=\$(basename "\$f" .circos_links.tsv)
            n=\$(awk 'NR > 1 { count++ } END { print count + 0 }' "\$f")
            echo -e "\${species}\\t\${n}"
        done
    } > redip_report/circos_links_summary.tsv
    """
}


workflow REDIPLOIDISATION {
    take:
    genomes_tsv
    species_tree
    iqtree_results
    genespace_wd

    main:
    EXTRACT_REDIP_SPECIES(genomes_tsv)

    redip_species_ch = EXTRACT_REDIP_SPECIES.out
        .splitText()
        .map { it.trim() }
        .filter { it }

    iqtree_tree_ch = iqtree_results.map { og, treefile, ufboot, status, nt ->
        tuple(og, treefile)
    }

    ROOT_GENE_TREE(iqtree_tree_ch, genomes_tsv)

    rooted_trees_ch = ROOT_GENE_TREE.out
        .map { og, rooted_tree, summary -> rooted_tree }
        .collect()

    rooting_summaries_ch = ROOT_GENE_TREE.out
        .map { og, rooted_tree, summary -> summary }
        .collect()

    WRITE_BRANCH_DEFS(species_tree, genomes_tsv)

    classify_input_ch = redip_species_ch.combine(rooted_trees_ch)

    CLASSIFY_REDIP_EVENTS(classify_input_ch, genomes_tsv)

    MAKE_REDIP_LINKS(
        CLASSIFY_REDIP_EVENTS.out,
        WRITE_BRANCH_DEFS.out,
        genespace_wd
    )

    PREP_REDIP_CIRCOS(MAKE_REDIP_LINKS.out)

    PLOT_REDIP_CIRCOS(PREP_REDIP_CIRCOS.out)

    classifications_ch = CLASSIFY_REDIP_EVENTS.out
        .map { species, classification -> classification }
        .collect()

    circos_links_ch = MAKE_REDIP_LINKS.out
        .map { species, links -> links }
        .collect()

    REDIP_REPORT(rooting_summaries_ch, classifications_ch, circos_links_ch)

    emit:
    redip_species = EXTRACT_REDIP_SPECIES.out
    rooted_trees = ROOT_GENE_TREE.out
    branch_definitions = WRITE_BRANCH_DEFS.out
    classifications = CLASSIFY_REDIP_EVENTS.out
    circos_links = MAKE_REDIP_LINKS.out
    circos_inputs = PREP_REDIP_CIRCOS.out
    circos_plots = PLOT_REDIP_CIRCOS.out
    report = REDIP_REPORT.out
}
