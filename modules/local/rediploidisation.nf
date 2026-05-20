nextflow.enable.dsl=2

/*
 * Rediploidisation workflow.
 *
 * This version supports:
 *   - default recurrent-safe classification mode
 *   - per-species text modes in genomes.tsv
 *   - backward-compatible 1/0 rediploidisation values
 *   - standard/recent/ancestral Circos layers
 *   - complete and branch-specific Circos plots for each layer
 *   - optional skip_rooting mode for already-rooted *.rooted.treefile inputs
 */

process EXTRACT_REDIP_SPECIES {
    tag "extract_redip_species"

    input:
    path genomes_tsv
    path redip_utils
    val redip_params

    output:
    path "rediploidisation/redip_species_modes.tsv", emit: modes
    path "rediploidisation/redip_species.txt", emit: species

    script:
    def classify_mode = redip_params.classify_mode ?: 'recurrent'

    """
    mkdir -p rediploidisation

    python - <<'PY'
import sys
from pathlib import Path

sys.path.insert(0, str(Path("${redip_utils}").parent))

from redip_utils import read_redip_species_modes_from_genomes_tsv

rows = read_redip_species_modes_from_genomes_tsv(
    "${genomes_tsv}",
    default_mode="${classify_mode}",
)

with open("rediploidisation/redip_species_modes.tsv", "w") as out:
    out.write("species\\tredip_mode\\n")
    for species, mode in rows:
        out.write(f"{species}\\t{mode}\\n")

# Kept for backward compatibility with earlier workflow code and quick manual checks.
with open("rediploidisation/redip_species.txt", "w") as out:
    for species, _mode in rows:
        out.write(species + "\\n")
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

    python ${root_script} \
        --tree ${iqtree_dir}/${og}_iqtree.treefile \
        --genomes-tsv ${genomes_tsv} \
        --output-tree rediploidisation/rooted_trees/${og}.rooted.treefile \
        --summary-tsv rediploidisation/rooting_summaries/${og}.rooting_summary.tsv \
        --tip-separator '${redip_params.tip_separator}' \
        --tree-format ${redip_params.gene_tree_format}
    """
}


process MARK_ROOTED_GENE_TREE {
    tag { tree_id }
    array (params.array_size as int)

    input:
    tuple val(tree_id), path(rooted_tree)

    output:
    tuple val(tree_id),
          path("rediploidisation/rooted_trees/${tree_id}.rooted.treefile"),
          emit: rooted_trees

    tuple val(tree_id),
          path("rediploidisation/rooting_summaries/${tree_id}.rooting_summary.tsv"),
          emit: summaries

    script:
    """
    mkdir -p rediploidisation/rooted_trees rediploidisation/rooting_summaries

    test -s ${rooted_tree}

    # redip_rooted mode expects already-rooted trees whose original filenames
    # end in .rooted.treefile. main.nf should pass tree_id as the full prefix
    # before that suffix, for example:
    #   og_8819_subtree_001.rooted.treefile -> tree_id=og_8819_subtree_001
    # This preserves subtree/sample identifiers instead of shortening them.
    case "\$(basename ${rooted_tree})" in
        *.rooted.treefile) ;;
        *)
            echo "ERROR: redip_rooted input does not end with .rooted.treefile: ${rooted_tree}" >&2
            exit 1
            ;;
    esac

    cp ${rooted_tree} rediploidisation/rooted_trees/${tree_id}.rooted.treefile

    {
        echo -e "tree_file\tstatus\tmessage"
        echo -e "${rooted_tree}\tALREADY_ROOTED\tRooting skipped because redip_params.skip_rooting was true"
    } > rediploidisation/rooting_summaries/${tree_id}.rooting_summary.tsv
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

    python ${branch_script} \
        --tree ${species_tree} \
        --genomes-tsv ${genomes_tsv} \
        --output-dir rediploidisation/branch_definitions \
        --tree-format ${redip_params.species_tree_format}
    """
}


process CLASSIFY_REDIP_EVENTS {
    tag { "${species}.${redip_mode}" }
    array (params.array_size as int)

    input:
    tuple val(species), val(redip_mode), path(rooted_trees)
    path genomes_tsv
    path classify_script
    path redip_utils
    val redip_params

    output:
    tuple val(species),
          val(redip_mode),
          path("rediploidisation/classifications/${species}.redip_classification.tsv")

    script:
    def tree_list = rooted_trees.collect { it.toString() }.join(' ')

    def recent_grouping = redip_params.recent_grouping ?: 'auto'
    def min_recent_groups = redip_params.min_recent_groups ?: 2
    def write_lossy_singletons = redip_params.containsKey('write_lossy_singletons') ? redip_params.write_lossy_singletons : true
    def min_ancestral_target_copies = redip_params.min_ancestral_target_copies

    def mode_args = ""
    if (redip_mode == 'standard') {
        mode_args = ""
    } else if (redip_mode == 'recurrent') {
        // Keep recurrent flags on one physical shell line so they remain part of
        // the python command instead of being executed as separate commands.
        mode_args = "--recurrent-wgd --recent-grouping ${recent_grouping} --min-recent-groups ${min_recent_groups}"
    } else if (redip_mode == 'recurrent_lossy') {
        def singleton_arg = write_lossy_singletons ? "--write-lossy-singletons" : ""
        def min_anc_arg = min_ancestral_target_copies ? "--min-ancestral-target-copies ${min_ancestral_target_copies}" : ""
        // Kept as one line for the same reason as recurrent mode.
        mode_args = "--recurrent-wgd --allow-ancestral-loss --recent-grouping ${recent_grouping} --min-recent-groups ${min_recent_groups} ${min_anc_arg} ${singleton_arg}"
    } else {
        throw new IllegalArgumentException("Unknown redip_mode: ${redip_mode}")
    }

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

    python ${classify_script} \
        --target-species ${species} \
        --treefile ${species}.rooted_gene_trees.nwk \
        --genomes-tsv ${genomes_tsv} \
        --output rediploidisation/classifications/${species}.redip_classification.tsv \
        --tip-separator '${redip_params.tip_separator}' \
        --label-format ${redip_params.label_format} \
        --copy-mode ${redip_params.copy_mode} \
        --required-copies ${redip_params.required_copies} \
        --min-tips ${redip_params.min_tips} ${mode_args}
    """
}


process MAKE_REDIP_LINKS {
    tag { "${species}.${plot_level}" }
    array (params.array_size as int)

    input:
    tuple val(species), val(redip_mode), val(plot_level), path(classification)
    path branch_definitions
    path genespace_wd
    path links_script
    path redip_utils
    val redip_params

    output:
    tuple val(species),
          val(plot_level),
          path("rediploidisation/circos_links/${plot_level}/complete/${species}.${plot_level}.circos_links.tsv")

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
    mkdir -p rediploidisation/circos_links/${plot_level}/complete

    python ${links_script} \
        --species ${species} \
        --classification-tsv ${classification} \
        ${position_arg} \
        --branch-definitions ${branch_definitions}/${species}.branch_definitions.tsv \
        --output rediploidisation/circos_links/${plot_level}/complete/${species}.${plot_level}.circos_links.tsv \
        --plot-level ${plot_level} \
        --tip-separator '${redip_params.tip_separator}' \
        --label-format ${redip_params.label_format} \
        --position-key-type ${redip_params.position_key_type} \
        --write-header \
        --include-metadata \
        --on-exists overwrite \
        --log-level INFO
    """
}


process PREP_REDIP_CIRCOS {
    tag { "${species}.${plot_level}" }
    array (params.array_size as int)

    input:
    tuple val(species), val(plot_level), path(circos_links), path(chr_bed)
    path prep_script
    path redip_utils

    output:
    tuple val(species),
          val(plot_level),
          path("rediploidisation/circos_inputs/${plot_level}/${species}"),
          emit: inputs

    tuple val(species),
          val(plot_level),
          path("rediploidisation/circos_input_summaries/${species}.${plot_level}.branch_summary.tsv"),
          emit: summaries

    script:
    """
    mkdir -p rediploidisation/circos_inputs/${plot_level}
    mkdir -p rediploidisation/circos_input_summaries

    test -s ${circos_links}
    test -s ${chr_bed}

    python ${prep_script} \
        --species ${species}.${plot_level} \
        --circos-links ${circos_links} \
        --chr-bed ${chr_bed} \
        --output-dir rediploidisation/circos_inputs/${plot_level}/${species}

    summary="rediploidisation/circos_inputs/${plot_level}/${species}/branch_summary.tsv"

    if [[ -s "\$summary" ]]; then
        cp "\$summary" rediploidisation/circos_input_summaries/${species}.${plot_level}.branch_summary.tsv
    else
        {
            echo -e "branch_id\tlink_rows"
            echo -e "complete\t0"
        } > rediploidisation/circos_input_summaries/${species}.${plot_level}.branch_summary.tsv
    fi
    """
}


process PLOT_REDIP_CIRCOS {
    tag { "${species}.${plot_level}" }
    array (params.array_size as int)

    input:
    tuple val(species), val(plot_level), path(prep_dir)

    output:
    tuple val(species),
          val(plot_level),
          path("rediploidisation/circos_plots/${plot_level}/${species}")

    script:
    """
    mkdir -p rediploidisation/circos_plots/${plot_level}

    cp -r ${prep_dir} rediploidisation/circos_plots/${plot_level}/${species}

    find rediploidisation/circos_plots/${plot_level}/${species} \
        -name circos.conf \
        -print0 | while IFS= read -r -d '' conf; do
            plot_dir=\$(dirname "\$conf")
            (
                cd "\$plot_dir"
                ${params.circos_bin} -conf circos.conf
            )
        done
    """
}


process PLOT_REDIP_SPECIES_TREE {
    tag { "${species}.${plot_level}" }

    input:
    tuple val(species), val(plot_level), path(circos_links)
    path species_tree
    path branch_definitions
    path plot_species_tree_script
    val redip_params

    output:
    tuple val(species),
          val(plot_level),
          path("rediploidisation/species_tree_plots/${plot_level}/${species}.${plot_level}.redip_species_tree.png")

    script:
    def annotation_arg = redip_params.species_tree_annotation?.toString()?.trim()
        ? "--annotation ${redip_params.species_tree_annotation}"
        : ""

    def colors_arg = redip_params.branch_colors?.toString()?.trim()
        ? "--colors ${redip_params.branch_colors}"
        : ""

    def width = redip_params.species_tree_plot_width ?: 8
    def height = redip_params.species_tree_plot_height ?: 6
    def dpi = redip_params.species_tree_plot_dpi ?: 300
    def circle_size = redip_params.species_tree_circle_size ?: 7
    def circle_stroke = redip_params.species_tree_circle_stroke ?: 0.5
    def branch_number_size = redip_params.species_tree_branch_number_size ?: 4
    def tip_label_size = redip_params.species_tree_tip_label_size ?: 5
    def tip_label_offset = redip_params.species_tree_tip_label_offset ?: 1
    def branch_label_x_offset = redip_params.species_tree_branch_label_x_offset ?: 0.15
    def branch_label_y_offset = redip_params.species_tree_branch_label_y_offset ?: 0.1
    def prune_to_branch_species = redip_params.containsKey('species_tree_prune_to_branch_species')
        ? redip_params.species_tree_prune_to_branch_species
        : true

    """
    mkdir -p rediploidisation/species_tree_plots/${plot_level}

    Rscript ${plot_species_tree_script} \
        --species-tree ${species_tree} \
        --branch-definitions ${branch_definitions}/${species}.branch_definitions.tsv \
        --circos-links ${circos_links} \
        --output rediploidisation/species_tree_plots/${plot_level}/${species}.${plot_level}.redip_species_tree.png \
        --width ${width} \
        --height ${height} \
        --dpi ${dpi} \
        --circle-size ${circle_size} \
        --circle-stroke ${circle_stroke} \
        --branch-number-size ${branch_number_size} \
        --tip-label-size ${tip_label_size} \
        --tip-label-offset ${tip_label_offset} \
        --branch-label-x-offset ${branch_label_x_offset} \
        --branch-label-y-offset ${branch_label_y_offset} \
        --prune-to-branch-species ${prune_to_branch_species} \
        ${annotation_arg} \
        ${colors_arg}
    """
}


process REDIP_REPORT {
    tag "redip_report"

    input:
    path rooting_summaries
    path classifications
    path circos_links
    path circos_input_summaries

    output:
    path "rediploidisation/report"

    script:
    def root_list = rooting_summaries.collect { it.toString() }.join(' ')
    def class_list = classifications.collect { it.toString() }.join(' ')
    def link_list = circos_links.collect { it.toString() }.join(' ')
    def summary_list = circos_input_summaries.collect { it.toString() }.join(' ')

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

    python - <<'PY'
import csv
from collections import Counter
from pathlib import Path

class_files = '${class_list}'.split()
link_files = '${link_list}'.split()
summary_files = '${summary_list}'.split()

report_dir = Path("rediploidisation/report")

with open(report_dir / "classification_summary.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["species", "event_level", "event_type", "loss_status", "rows"])
    for path in class_files:
        species = Path(path).name.replace(".redip_classification.tsv", "")
        counts = Counter()
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                counts[(
                    row.get("event_level", "standard"),
                    row.get("event_type", "standard_exact_n"),
                    row.get("loss_status", "complete"),
                )] += 1
        if not counts:
            writer.writerow([species, "none", "none", "none", 0])
        else:
            for key, n in sorted(counts.items()):
                writer.writerow([species, *key, n])

with open(report_dir / "circos_links_summary.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["species", "plot_level", "circos_link_rows"])
    for path in link_files:
        name = Path(path).name.replace(".circos_links.tsv", "")
        if "." in name:
            species, plot_level = name.rsplit(".", 1)
        else:
            species, plot_level = name, "unknown"
        with open(path, newline="") as handle:
            n = max(sum(1 for _ in handle) - 1, 0)
        writer.writerow([species, plot_level, n])

with open(report_dir / "circos_branch_summary.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["species", "plot_level", "branch_id", "link_rows"])
    for summary_file in summary_files:
        p = Path(summary_file)
        name = p.name.replace(".branch_summary.tsv", "")

        if "." in name:
            species, plot_level = name.rsplit(".", 1)
        else:
            species, plot_level = name, "unknown"

        if not p.exists():
            continue

        with open(p, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                writer.writerow([
                    species,
                    plot_level,
                    row.get("branch_id", ""),
                    row.get("link_rows", ""),
                ])
PY
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
    plot_species_tree_script_ch = Channel.value(file('scripts/rediploidisation/plot_redip_species_tree.R'))

    def skip_rooting = redip_params.containsKey('skip_rooting')
        ? redip_params.skip_rooting
        : false

    def rooted_trees_emit_ch = Channel.empty()
    def rooting_summaries_emit_ch = Channel.empty()
    def rooted_trees_ch
    def rooting_summaries_ch

    EXTRACT_REDIP_SPECIES(genomes_tsv, redip_utils_ch, redip_params)

    /*
     * Parse per-species redip modes from the helper output.
     *
     * Important:
     * Do not collect and then index species strings with item[0].
     * That caused species like Eric_Calluna_vul to become E.
     *
     * These channels are parsed separately from the same small TSV so each
     * downstream branch has its own queue channel.
     */
    redip_species_modes_for_species_ch = EXTRACT_REDIP_SPECIES.out.modes
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.species as String, row.redip_mode as String)
        }

    redip_species_modes_for_classify_ch = EXTRACT_REDIP_SPECIES.out.modes
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.species as String, row.redip_mode as String)
        }

    redip_species_modes_for_plot_ch = EXTRACT_REDIP_SPECIES.out.modes
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.species as String, row.redip_mode as String)
        }

    redip_species_ch = redip_species_modes_for_species_ch
        .map { species_name, redip_mode ->
            species_name
        }

    if (skip_rooting) {
        /*
         * redip_rooted mode:
         * iqtree_results is expected to be tuple(tree_id, rooted_tree), where
         * rooted_tree must end in .rooted.treefile and tree_id is the full prefix
         * before that suffix.
         */
        MARK_ROOTED_GENE_TREE(iqtree_results)

        rooted_trees_emit_ch = MARK_ROOTED_GENE_TREE.out.rooted_trees
        rooting_summaries_emit_ch = MARK_ROOTED_GENE_TREE.out.summaries

        rooted_trees_ch = MARK_ROOTED_GENE_TREE.out.rooted_trees
            .map { tree_id, rooted_tree -> rooted_tree }
            .collect()

        rooting_summaries_ch = MARK_ROOTED_GENE_TREE.out.summaries
            .map { tree_id, summary -> summary }
            .collect()

    } else {
        /*
         * Standard redip mode:
         * iqtree_results is expected to be tuple(og, iqtree_dir).
         */
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

        rooted_trees_emit_ch = ROOT_GENE_TREE.out.rooted_trees
        rooting_summaries_emit_ch = ROOT_GENE_TREE.out.summaries

        rooted_trees_ch = ROOT_GENE_TREE.out.rooted_trees
            .map { og, rooted_tree -> rooted_tree }
            .collect()

        rooting_summaries_ch = ROOT_GENE_TREE.out.summaries
            .map { og, summary -> summary }
            .collect()
    }

    WRITE_BRANCH_DEFS(
        species_tree,
        genomes_tsv,
        branch_script_ch,
        redip_utils_ch,
        redip_params
    )

    /*
     * Pair each focal species/mode with the one collected list of rooted trees.
     *
     * combine() can emit either:
     *   [species, mode, rooted_trees]
     * or:
     *   [[species, mode], rooted_trees]
     *
     * This row-unpacking form handles both and avoids the LinkedList closure
     * error seen when the collected rooted-tree list is interpreted as
     * multiple closure arguments.
     */
    classify_input_ch = redip_species_modes_for_classify_ch
        .combine(rooted_trees_ch)
        .map { row ->
            def species_name
            def redip_mode
            def rooted_trees

            if (row.size() == 2 && row[0] instanceof List) {
                species_name = row[0][0]
                redip_mode = row[0][1]
                rooted_trees = row[1]
            } else {
                species_name = row[0]
                redip_mode = row[1]
                // If combine() flattens the collected rooted-tree list, keep
                // everything from index 2 onward as the tree list. If it does
                // not flatten, row[2] is already the collected tree list.
                rooted_trees = (row.size() == 3 && row[2] instanceof List)
                    ? row[2]
                    : row[2..-1]
            }

            tuple(species_name as String, redip_mode as String, rooted_trees)
        }

    CLASSIFY_REDIP_EVENTS(
        classify_input_ch,
        genomes_tsv,
        classify_script_ch,
        redip_utils_ch,
        redip_params
    )

    /*
     * Create Circos plot layers.
     *
     * standard species:
     *   standard only
     *
     * recurrent / recurrent_lossy species:
     *   standard, recent, ancestral
     */
    plot_levels_ch = redip_species_modes_for_plot_ch
        .flatMap { species_name, redip_mode ->
            if (redip_mode == 'standard') {
                return [
                    tuple(species_name, redip_mode, 'standard')
                ]
            }

            return [
                tuple(species_name, redip_mode, 'standard'),
                tuple(species_name, redip_mode, 'recent'),
                tuple(species_name, redip_mode, 'ancestral')
            ]
        }

    classifications_for_join_ch = CLASSIFY_REDIP_EVENTS.out
        .map { species_name, redip_mode, classification ->
            tuple(species_name, classification)
        }
    
    links_input_ch = plot_levels_ch
        .combine(classifications_for_join_ch, by: 0)
        .map { species_name, redip_mode, plot_level, classification ->
            tuple(species_name, redip_mode, plot_level, classification)
        }

    MAKE_REDIP_LINKS(
        links_input_ch,
        WRITE_BRANCH_DEFS.out,
        genespace_wd,
        links_script_ch,
        redip_utils_ch,
        redip_params
    )

    /*
     * Pair each species/layer Circos links file with the chromosome BED for
     * that focal species.
     */
    chr_beds_ch = redip_species_ch.map { species_name ->
        tuple(
            species_name,
            file("${params.chr_dict_dir}/${species_name}_chr_lengths.bed", checkIfExists: true)
        )
    }

    circos_links_with_chr_beds_ch = MAKE_REDIP_LINKS.out
        .combine(chr_beds_ch, by: 0)
        .map { species_name, plot_level, links, chr_bed ->
            tuple(species_name, plot_level, links, chr_bed)
        }

    PREP_REDIP_CIRCOS(
        circos_links_with_chr_beds_ch,
        prep_script_ch,
        redip_utils_ch
    )

    PLOT_REDIP_CIRCOS(PREP_REDIP_CIRCOS.out.inputs)

    PLOT_REDIP_SPECIES_TREE(
        MAKE_REDIP_LINKS.out,
        species_tree,
        WRITE_BRANCH_DEFS.out,
        plot_species_tree_script_ch,
        redip_params
    )

    classifications_ch = CLASSIFY_REDIP_EVENTS.out
        .map { species_name, redip_mode, classification -> classification }
        .collect()

    circos_links_ch = MAKE_REDIP_LINKS.out
        .map { species_name, plot_level, links -> links }
        .collect()

    circos_input_summaries_ch = PREP_REDIP_CIRCOS.out.summaries
        .map { species_name, plot_level, summary -> summary }
        .collect()

    REDIP_REPORT(
        rooting_summaries_ch,
        classifications_ch,
        circos_links_ch,
        circos_input_summaries_ch
    )

    emit:
    redip_species_modes = EXTRACT_REDIP_SPECIES.out.modes
    redip_species = EXTRACT_REDIP_SPECIES.out.species
    rooted_trees = rooted_trees_emit_ch
    rooting_summaries = rooting_summaries_emit_ch
    branch_definitions = WRITE_BRANCH_DEFS.out
    classifications = CLASSIFY_REDIP_EVENTS.out
    circos_links = MAKE_REDIP_LINKS.out
    circos_inputs = PREP_REDIP_CIRCOS.out.inputs
    circos_input_summaries = PREP_REDIP_CIRCOS.out.summaries
    circos_plots = PLOT_REDIP_CIRCOS.out
    species_tree_plots = PLOT_REDIP_SPECIES_TREE.out
    report = REDIP_REPORT.out
}
