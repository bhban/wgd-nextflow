nextflow.enable.dsl=2

process EXTRACT_REDIP_SPECIES {
    tag "extract_redip_species"

    input:
    path genomes_tsv
    path redip_utils

    output:
    path "rediploidisation/redip_species.txt"

    script:
    """
    mkdir -p rediploidisation

    python - <<'PY'
import sys
from pathlib import Path

sys.path.insert(0, str(Path("${redip_utils}").parent))

from redip_utils import read_redip_species_from_genomes_tsv

species = read_redip_species_from_genomes_tsv("${genomes_tsv}")

with open("rediploidisation/redip_species.txt", "w") as out:
    for sp in species:
        out.write(sp + "\\n")
PY
    """
}


process ROOT_GENE_TREE {
    tag { og }
    array (params.array_size as int)

    input:
    tuple val(og), path(iqtree_dir)
    path genomes_tsv
    path root_script
    path redip_utils
    val redip_params

    output:
    tuple val(og),
          path("rediploidisation/rooted_trees/${og}.rooted.treefile"),
          optional: true,
          emit: rooted_trees

    tuple val(og),
          path("rediploidisation/rooting_summaries/${og}.rooting_summary.tsv"),
          emit: summaries

    script:
    """
    mkdir -p rediploidisation/rooted_trees rediploidisation/rooting_summaries

    test -s ${iqtree_dir}/${og}_iqtree.treefile

    python ${root_script} \\
        --tree ${iqtree_dir}/${og}_iqtree.treefile \\
        --genomes-tsv ${genomes_tsv} \\
        --output-tree rediploidisation/rooted_trees/${og}.rooted.treefile \\
        --summary-tsv rediploidisation/rooting_summaries/${og}.rooting_summary.tsv \\
        --tip-separator '${redip_params.tip_separator}' \\
        --tree-format ${redip_params.gene_tree_format}
    """
}


process WRITE_BRANCH_DEFS {
    tag "write_branch_defs"

    input:
    path species_tree
    path genomes_tsv
    path branch_script
    path redip_utils
    val redip_params

    output:
    path "rediploidisation/branch_definitions"

    script:
    """
    mkdir -p rediploidisation/branch_definitions

    python ${branch_script} \\
        --tree ${species_tree} \\
        --genomes-tsv ${genomes_tsv} \\
        --output-dir rediploidisation/branch_definitions \\
        --tree-format ${redip_params.species_tree_format}
    """
}


process CLASSIFY_REDIP_EVENTS {
    tag { species }
    array (params.array_size as int)

    input:
    tuple val(species), path(rooted_trees)
    path genomes_tsv
    path classify_script
    path redip_utils
    val redip_params

    output:
    tuple val(species), path("rediploidisation/classifications/${species}.redip_classification.tsv")

    script:
    def tree_list = rooted_trees.collect { it.toString() }.join(' ')

    """
    mkdir -p rediploidisation/classifications

    : > ${species}.rooted_gene_trees.nwk
    
    for tree in ${tree_list}; do
        tree_base=\$(basename "\$tree")
    
        while IFS= read -r newick || [ -n "\$newick" ]; do
            [ -z "\$newick" ] && continue
            printf "%s\\t%s\\n" "\$tree_base" "\$newick" >> ${species}.rooted_gene_trees.nwk
        done < "\$tree"
    done

    python ${classify_script} \\
        --target-species ${species} \\
        --treefile ${species}.rooted_gene_trees.nwk \\
        --genomes-tsv ${genomes_tsv} \\
        --output rediploidisation/classifications/${species}.redip_classification.tsv \\
        --tip-separator '${redip_params.tip_separator}' \\
        --label-format ${redip_params.label_format} \\
        --copy-mode ${redip_params.copy_mode} \\
        --required-copies ${redip_params.required_copies} \\
        --min-tips ${redip_params.min_tips}
    """
}


process MAKE_REDIP_LINKS {
    tag { species }
    array (params.array_size as int)

    input:
    tuple val(species), path(classification)
    path branch_definitions
    path genespace_wd
    path links_script
    path redip_utils
    val redip_params

    output:
    tuple val(species), path("rediploidisation/circos_links/${species}.circos_links.tsv")

    script:
    def source = redip_params.positions_source ?: 'bed'

    def position_arg
    if (source == 'pangenes') {
        position_arg = "--pangenes-dir ${genespace_wd}/pangenes"
    } else if (source == 'positions') {
        if (!redip_params.positions?.toString()?.trim()) {
            throw new IllegalArgumentException(
                "rediploidisation.positions must be provided when positions_source = 'positions'"
            )
        }

        def species_col_arg = redip_params.position_species_column?.toString()?.trim()
            ? "--position-species-column ${redip_params.position_species_column}"
            : ""

        position_arg = """
            --positions ${redip_params.positions}
            --position-format ${redip_params.position_format}
            ${redip_params.position_has_header ? '--position-has-header' : '--no-position-has-header'}
            --position-key-column ${redip_params.position_key_column}
            --position-chr-column ${redip_params.position_chr_column}
            --position-start-column ${redip_params.position_start_column}
            --position-end-column ${redip_params.position_end_column}
            ${species_col_arg}
        """
    } else {
        position_arg = "--bed-dir ${genespace_wd}/bed"
    }

    """
    mkdir -p rediploidisation/circos_links

    python ${links_script} \\
        --species ${species} \\
        --classification-tsv ${classification} \\
        ${position_arg} \\
        --branch-definitions ${branch_definitions}/${species}.branch_definitions.tsv \\
        --output rediploidisation/circos_links/${species}.circos_links.tsv \\
        --tip-separator '${redip_params.tip_separator}' \\
        --label-format ${redip_params.label_format} \\
        --position-key-type ${redip_params.position_key_type} \\
        --write-header \\
        --include-metadata \\
        --on-exists overwrite \\
        --log-level INFO
    """
}


process PREP_REDIP_CIRCOS {
    tag { species }
    array (params.array_size as int)

    input:
    tuple val(species), path(circos_links), path(chr_bed)
    path prep_script
    path redip_utils

    output:
    tuple val(species), path("rediploidisation/circos_inputs/${species}")

    script:
    """
    mkdir -p rediploidisation/circos_inputs

    test -s ${circos_links}
    test -s ${chr_bed}

    python ${prep_script} \\
        --species ${species} \\
        --circos-links ${circos_links} \\
        --chr-bed ${chr_bed} \\
        --output-dir rediploidisation/circos_inputs/${species}
    """
}


process PLOT_REDIP_CIRCOS {
    tag { species }
    array (params.array_size as int)

    input:
    tuple val(species), path(species_dir)

    output:
    tuple val(species), path("rediploidisation/circos_plots/${species}")

    script:
    """
    mkdir -p rediploidisation/circos_plots

    cp -r ${species_dir} rediploidisation/circos_plots/${species}

    cd rediploidisation/circos_plots/${species}

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
    path "rediploidisation/report"

    script:
    def root_list = rooting_summaries.collect { it.toString() }.join(' ')
    def class_list = classifications.collect { it.toString() }.join(' ')
    def link_list = circos_links.collect { it.toString() }.join(' ')

    """
    mkdir -p rediploidisation/report

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
    } > rediploidisation/report/rooting_summary.tsv

    {
        echo -e "species\\tclassification_rows"
        for f in ${class_list}; do
            species=\$(basename "\$f" .redip_classification.tsv)
            n=\$(awk 'NR > 1 { count++ } END { print count + 0 }' "\$f")
            echo -e "\${species}\\t\${n}"
        done
    } > rediploidisation/report/classification_summary.tsv

    {
        echo -e "species\\tcircos_link_rows"
        for f in ${link_list}; do
            species=\$(basename "\$f" .circos_links.tsv)
            n=\$(awk 'NR > 1 { count++ } END { print count + 0 }' "\$f")
            echo -e "\${species}\\t\${n}"
        done
    } > rediploidisation/report/circos_links_summary.tsv
    """
}


workflow REDIPLOIDISATION {
    take:
    genomes_tsv
    species_tree
    iqtree_results
    genespace_wd
    redip_params

    main:
    redip_utils_ch = Channel.value(file('scripts/rediploidisation/redip_utils.py'))
    root_script_ch = Channel.value(file('scripts/rediploidisation/root_gene_tree.py'))
    branch_script_ch = Channel.value(file('scripts/rediploidisation/write_branch_defs.py'))
    classify_script_ch = Channel.value(file('scripts/rediploidisation/classify.py'))
    links_script_ch = Channel.value(file('scripts/rediploidisation/make_links.py'))
    prep_script_ch = Channel.value(file('scripts/rediploidisation/prep_circos.py'))

    EXTRACT_REDIP_SPECIES(genomes_tsv, redip_utils_ch)

    redip_species_ch = EXTRACT_REDIP_SPECIES.out
        .splitText()
        .map { it.trim() }
        .filter { it }

    iqtree_tree_ch = iqtree_results.map { og, iqtree_dir ->
        tuple(og, iqtree_dir)
    }

    ROOT_GENE_TREE(
        iqtree_tree_ch,
        genomes_tsv,
        root_script_ch,
        redip_utils_ch,
        redip_params
    )

    rooted_trees_ch = ROOT_GENE_TREE.out.rooted_trees
        .map { og, rooted_tree -> rooted_tree }
        .collect()

    rooting_summaries_ch = ROOT_GENE_TREE.out.summaries
        .map { og, summary -> summary }
        .collect()

    WRITE_BRANCH_DEFS(
        species_tree,
        genomes_tsv,
        branch_script_ch,
        redip_utils_ch,
        redip_params
    )

    classify_input_ch = redip_species_ch
        .combine(rooted_trees_ch)
        .map { row ->
            def species = row[0]
            def rooted_trees = (row.size() == 2 && row[1] instanceof List)
                ? row[1]
                : row[1..-1]

            tuple(species, rooted_trees)
        }

    CLASSIFY_REDIP_EVENTS(
        classify_input_ch,
        genomes_tsv,
        classify_script_ch,
        redip_utils_ch,
        redip_params
    )

    MAKE_REDIP_LINKS(
        CLASSIFY_REDIP_EVENTS.out,
        WRITE_BRANCH_DEFS.out,
        genespace_wd,
        links_script_ch,
        redip_utils_ch,
        redip_params
    )

    chr_beds_ch = redip_species_ch.map { species ->
        tuple(
            species,
            file("${params.chr_dict_dir}/${species}_chr_lengths.bed", checkIfExists: true)
        )
    }

    circos_links_with_chr_beds_ch = MAKE_REDIP_LINKS.out
        .join(chr_beds_ch)

    PREP_REDIP_CIRCOS(
        circos_links_with_chr_beds_ch,
        prep_script_ch,
        redip_utils_ch
    )

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
    rooted_trees = ROOT_GENE_TREE.out.rooted_trees
    rooting_summaries = ROOT_GENE_TREE.out.summaries
    branch_definitions = WRITE_BRANCH_DEFS.out
    classifications = CLASSIFY_REDIP_EVENTS.out
    circos_links = MAKE_REDIP_LINKS.out
    circos_inputs = PREP_REDIP_CIRCOS.out
    circos_plots = PLOT_REDIP_CIRCOS.out
    report = REDIP_REPORT.out
}
