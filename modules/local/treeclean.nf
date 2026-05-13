/*
 * Tree-cleaning subworkflow for wgd-nextflow.
 *
 * Expected inputs:
 *
 * ch_og_raw_fasta:
 *   tuple val(og), path(raw_cds_fasta)
 *
 * ch_og_iqtree_dir:
 *   tuple val(og), path(iqtree_dir)
 *
 * genomes_tsv:
 *   path genomes.tsv
 */

process DECOMPOSE_LONG_BRANCH_TREES {
    tag { og }

    input:
    tuple val(og), path(raw_fasta), path(iqtree_dir)
    path decompose_script

    output:
    tuple val(og),
          path("${og}.subtree_membership.tsv"),
          path("${og}.decomposition_report.tsv")

    script:
    """
    set -euo pipefail

    treefile="${iqtree_dir}/og_${og}_iqtree.treefile"

    if [ ! -s "\$treefile" ]; then
      treefile="${iqtree_dir}/${og}_iqtree.treefile"
    fi

    if [ ! -s "\$treefile" ]; then
      echo "Could not find IQ-TREE treefile for OG ${og} in ${iqtree_dir}" >&2
      echo "Expected one of:" >&2
      echo "  ${iqtree_dir}/og_${og}_iqtree.treefile" >&2
      echo "  ${iqtree_dir}/${og}_iqtree.treefile" >&2
      exit 1
    fi

    python3 ${decompose_script} \
      --og "og_${og}" \
      --tree "\$treefile" \
      --min-branch ${params.treeclean_min_branch} \
      --min-tips ${params.treeclean_min_tips} \
      --out-membership "${og}.subtree_membership.tsv" \
      --out-report "${og}.decomposition_report.tsv"
    """
}


process FILTER_DECOMPOSED_SUBTREES {
    tag { og }

    input:
    tuple val(og), path(membership), path(decomp_report)
    path genomes_tsv
    path filter_script

    output:
    tuple val(og),
          path("${og}.passing_subtree_membership.tsv"),
          path("${og}.subtree_filter_report.tsv")

    script:
    def require_outgroup_arg = params.treeclean_require_outgroup ? "--require-outgroup" : ""

    """
    set -euo pipefail

    python3 ${filter_script} \
      --membership "${membership}" \
      --genomes-tsv "${genomes_tsv}" \
      --min-species ${params.treeclean_min_species} \
      ${require_outgroup_arg} \
      --tip-separator "${params.treeclean_tip_separator}" \
      --tip-label-format "${params.treeclean_tip_label_format}" \
      --out-membership "${og}.passing_subtree_membership.tsv" \
      --out-report "${og}.subtree_filter_report.tsv"
    """
}


process WRITE_DECOMPOSED_FASTAS {
    tag { og }

    input:
    tuple val(og), path(raw_fasta), path(passing_membership)
    path write_fastas_script

    output:
    path "decomposed_fasta/*.fasta", emit: fastas, optional: true
    path "${og}.decomposed_fasta_manifest.tsv", emit: manifest

    script:
    """
    set -euo pipefail

    mkdir -p decomposed_fasta

    python3 ${write_fastas_script} \
      --og "og_${og}" \
      --raw-fasta "${raw_fasta}" \
      --membership "${passing_membership}" \
      --out-dir decomposed_fasta \
      --out-manifest "${og}.decomposed_fasta_manifest.tsv"
    """
}


process MACSE_ALIGN_DECOMPOSED {
    tag { unit }

    input:
    tuple val(unit), path(raw_fasta)

    output:
    tuple val(unit),
          path("${unit}_AA.fasta"),
          path("${unit}_NT.fasta"),
          path("${unit}.macse.status"),
          path("${unit}.macse.log")

    script:
    """
    set -euo pipefail

    rm -f \
      "${unit}.macse.status" \
      "${unit}_AA.fasta" \
      "${unit}_NT.fasta" \
      "${unit}.macse.log"

    set +e
    ${params.macse_bin} \
      -prog alignSequences \
      -seq "${raw_fasta}" \
      -out_AA "${unit}_AA.fasta" \
      -out_NT "${unit}_NT.fasta" \
      > "${unit}.macse.log" 2>&1
    exit_code=\$?
    set -e

    if [ "\$exit_code" -eq 0 ] && [ -s "${unit}_NT.fasta" ]; then
      echo "OK" > "${unit}.macse.status"
    else
      echo "FAIL" > "${unit}.macse.status"
      exit "\$exit_code"
    fi
    """
}


process IQTREE_DECOMPOSED {
    tag { unit }

    input:
    tuple val(unit), path(aa_aln), path(nt_aln), path(macse_status), path(macse_log)

    output:
    tuple val(unit),
          path("iqtree_decomposed/${unit}"),
          path(nt_aln)

    script:
    """
    set -euo pipefail

    mkdir -p "iqtree_decomposed/${unit}"

    rm -f iqtree_decomposed/${unit}/${unit}_iqtree.*

    set +e
    ${params.iqtree_bin} \
      -s "${nt_aln}" \
      -nt ${task.cpus} \
      -m MFP \
      -bb 1000 \
      -wbtl \
      -redo \
      -pre "iqtree_decomposed/${unit}/${unit}_iqtree" \
      > "iqtree_decomposed/${unit}/${unit}_iqtree.log" 2>&1
    exit_code=\$?
    set -e

    if [ "\$exit_code" -eq 0 ] && [ -s "iqtree_decomposed/${unit}/${unit}_iqtree.treefile" ]; then
      echo "OK" > "iqtree_decomposed/${unit}/${unit}_iqtree.status"
    else
      echo "FAIL" > "iqtree_decomposed/${unit}/${unit}_iqtree.status"
      exit "\$exit_code"
    fi
    """
}


process MAKE_TREESHRINK_INPUTS {
    tag "treeshrink_inputs"

    input:
    path iqtree_dirs
    path make_treeshrink_inputs_script

    output:
    path "treeshrink_input.trees"
    path "treeshrink_tree_order.tsv"

    script:
    """
    set -euo pipefail

    python3 ${make_treeshrink_inputs_script} \
      --iqtree-dirs ${iqtree_dirs} \
      --out-trees treeshrink_input.trees \
      --out-order treeshrink_tree_order.tsv
    """
}


process WRITE_TREESHRINK_PROTECTED_SPECIES {
    tag "protected_species"

    input:
    path genomes_tsv
    path write_protected_outgroups_script

    output:
    path "treeshrink_protected_species.txt"

    script:
    """
    set -euo pipefail

    python3 ${write_protected_outgroups_script} \
      --genomes-tsv "${genomes_tsv}" \
      --out treeshrink_protected_species.txt
    """
}


process RUN_TREESHRINK_DECOMPOSED {
    tag "treeshrink"

    input:
    path trees
    path tree_order
    path protected_species
    path parse_treeshrink_removals_script

    output:
    path "treeshrink_out", emit: outdir
    path "treeshrink.removals.tsv", emit: removals
    path "treeshrink.log", emit: log

    script:
    """
    set -euo pipefail

    mkdir -p treeshrink_out

    echo "TreeShrink mode: ${params.treeshrink_mode}" > treeshrink.log
    echo "TreeShrink alpha: ${params.treeshrink_alpha}" >> treeshrink.log
    echo "Tree count: \$(grep -c . ${trees})" >> treeshrink.log
    echo "Protected species will be handled during pruning, not by TreeShrink -x." >> treeshrink.log
    echo "Protected species:" >> treeshrink.log
    cat "${protected_species}" >> treeshrink.log || true
    echo "" >> treeshrink.log

    run_treeshrink.py \\
      -t "${trees}" \\
      -m ${params.treeshrink_mode} \\
      -q ${params.treeshrink_alpha} \\
      -o treeshrink_out \\
      -O treeshrink \\
      >> treeshrink.log 2>&1

    python3 ${parse_treeshrink_removals_script} \\
      --treeshrink-dir treeshrink_out \\
      --prefix treeshrink \\
      --alpha ${params.treeshrink_alpha} \\
      --tree-order "${tree_order}" \\
      --out treeshrink.removals.tsv
    """
}


process WRITE_TREESHRINK_PRUNED_FASTAS {
    tag "pruned_fastas"

    input:
    path decomposed_fastas
    path removals
    path genomes_tsv
    path prune_fastas_script

    output:
    path "pruned_fasta/*.fasta", emit: fastas, optional: true
    path "post_treeshrink_membership.tsv", emit: membership
    path "treeshrink_pruning_report.tsv", emit: report

    script:
    def protect_outgroups_arg = params.treeshrink_protect_outgroups ? "--protect-outgroups" : ""

    """
    set -euo pipefail

    mkdir -p pruned_fasta

    python3 ${prune_fastas_script} \\
      --fastas ${decomposed_fastas} \\
      --removals "${removals}" \\
      --genomes-tsv "${genomes_tsv}" \\
      --tip-separator "${params.treeclean_tip_separator}" \\
      --tip-label-format "${params.treeclean_tip_label_format}" \\
      ${protect_outgroups_arg} \\
      --out-dir pruned_fasta \\
      --out-membership post_treeshrink_membership.tsv \\
      --out-report treeshrink_pruning_report.tsv
    """
}


process FILTER_POST_TREESHRINK_SUBTREES {
    tag "post_treeshrink_filter"

    input:
    path membership
    path genomes_tsv
    path filter_script

    output:
    path "final_passing_membership.tsv"
    path "post_treeshrink_filter_report.tsv"

    script:
    def require_outgroup_arg = params.treeclean_require_outgroup ? "--require-outgroup" : ""

    """
    set -euo pipefail

    python3 ${filter_script} \
      --membership "${membership}" \
      --genomes-tsv "${genomes_tsv}" \
      --min-species ${params.treeclean_min_species} \
      ${require_outgroup_arg} \
      --tip-separator "${params.treeclean_tip_separator}" \
      --tip-label-format "${params.treeclean_tip_label_format}" \
      --out-membership final_passing_membership.tsv \
      --out-report post_treeshrink_filter_report.tsv
    """
}


process SELECT_FINAL_FASTAS {
    tag "select_final_fastas"

    input:
    path pruned_fastas
    path final_passing_membership
    path select_fastas_script

    output:
    path "final_raw_fasta/*.fasta", emit: fastas, optional: true
    path "final_raw_fasta_manifest.tsv", emit: manifest

    script:
    """
    set -euo pipefail

    mkdir -p final_raw_fasta

    python3 ${select_fastas_script} \
      --fastas ${pruned_fastas} \
      --membership "${final_passing_membership}" \
      --out-dir final_raw_fasta \
      --out-manifest final_raw_fasta_manifest.tsv
    """
}


process MACSE_ALIGN_TREESHRINK_PRUNED {
    tag { unit }

    input:
    tuple val(unit), path(raw_fasta)

    output:
    tuple val(unit),
          path("${unit}.cleaned_AA.fasta"),
          path("${unit}.cleaned_NT.fasta"),
          path("${unit}.cleaned_macse.status"),
          path("${unit}.cleaned_macse.log")

    script:
    """
    set -euo pipefail

    rm -f \
      "${unit}.cleaned_macse.status" \
      "${unit}.cleaned_AA.fasta" \
      "${unit}.cleaned_NT.fasta" \
      "${unit}.cleaned_macse.log"

    set +e
    ${params.macse_bin} \
      -prog alignSequences \
      -seq "${raw_fasta}" \
      -out_AA "${unit}.cleaned_AA.fasta" \
      -out_NT "${unit}.cleaned_NT.fasta" \
      > "${unit}.cleaned_macse.log" 2>&1
    exit_code=\$?
    set -e

    if [ "\$exit_code" -eq 0 ] && [ -s "${unit}.cleaned_NT.fasta" ]; then
      echo "OK" > "${unit}.cleaned_macse.status"
    else
      echo "FAIL" > "${unit}.cleaned_macse.status"
      exit "\$exit_code"
    fi
    """
}


process IQTREE_CLEANED_FINAL {
    tag { unit }

    input:
    tuple val(unit), path(aa_aln), path(nt_aln), path(macse_status), path(macse_log)

    output:
    tuple val(unit),
          path("iqtree_cleaned/${unit}"),
          path(nt_aln)

    script:
    """
    set -euo pipefail

    mkdir -p "iqtree_cleaned/${unit}"

    rm -f iqtree_cleaned/${unit}/${unit}_iqtree.*

    set +e
    ${params.iqtree_bin} \
      -s "${nt_aln}" \
      -nt ${task.cpus} \
      -m MFP \
      -bb 1000 \
      -wbtl \
      -redo \
      -pre "iqtree_cleaned/${unit}/${unit}_iqtree" \
      > "iqtree_cleaned/${unit}/${unit}_iqtree.log" 2>&1
    exit_code=\$?
    set -e

    if [ "\$exit_code" -eq 0 ] && [ -s "iqtree_cleaned/${unit}/${unit}_iqtree.treefile" ]; then
      echo "OK" > "iqtree_cleaned/${unit}/${unit}_iqtree.status"
    else
      echo "FAIL" > "iqtree_cleaned/${unit}/${unit}_iqtree.status"
      exit "\$exit_code"
    fi
    """
}


process TREE_CLEANING_REPORT {
    tag "tree_cleaning_report"

    input:
    path decomp_reports
    path subtree_filter_reports
    path treeshrink_pruning_report
    path post_treeshrink_filter_report
    path combine_reports_script

    output:
    path "tree_cleaning_report.tsv"

    script:
    """
    set -euo pipefail

    python3 ${combine_reports_script} \
      --decomposition-reports ${decomp_reports} \
      --subtree-filter-reports ${subtree_filter_reports} \
      --treeshrink-pruning-report "${treeshrink_pruning_report}" \
      --post-treeshrink-filter-report "${post_treeshrink_filter_report}" \
      --out tree_cleaning_report.tsv
    """
}


workflow TREECLEAN {
    take:
    ch_og_raw_fasta
    ch_og_iqtree_dir
    ch_og_nt_alignment
    genomes_tsv

    main:
    decompose_script_ch = Channel.value(file("${projectDir}/scripts/treeclean/decompose_tree_by_long_branches.py", checkIfExists: true))
    filter_script_ch = Channel.value(file("${projectDir}/scripts/treeclean/filter_cleaning_units.py", checkIfExists: true))
    write_fastas_script_ch = Channel.value(file("${projectDir}/scripts/treeclean/write_fastas_from_membership.py", checkIfExists: true))
    make_treeshrink_inputs_script_ch = Channel.value(file("${projectDir}/scripts/treeclean/make_treeshrink_inputs.py", checkIfExists: true))
    write_protected_outgroups_script_ch = Channel.value(file("${projectDir}/scripts/treeclean/write_protected_outgroups.py", checkIfExists: true))
    parse_treeshrink_removals_script_ch = Channel.value(file("${projectDir}/scripts/treeclean/parse_treeshrink_removals.py", checkIfExists: true))
    prune_fastas_script_ch = Channel.value(file("${projectDir}/scripts/treeclean/prune_fastas_from_treeshrink.py", checkIfExists: true))
    select_fastas_script_ch = Channel.value(file("${projectDir}/scripts/treeclean/select_fastas_from_membership.py", checkIfExists: true))
    combine_reports_script_ch = Channel.value(file("${projectDir}/scripts/treeclean/combine_treeclean_reports.py", checkIfExists: true))

    ch_og_raw_tree = ch_og_raw_fasta
        .join(ch_og_iqtree_dir)
        .map { og, raw_fasta, iqtree_dir ->
            tuple(og, raw_fasta, iqtree_dir)
        }

    DECOMPOSE_LONG_BRANCH_TREES(
        ch_og_raw_tree,
        decompose_script_ch
    )

    FILTER_DECOMPOSED_SUBTREES(
        DECOMPOSE_LONG_BRANCH_TREES.out,
        genomes_tsv,
        filter_script_ch
    )

    ch_raw_plus_passing_membership = ch_og_raw_fasta
        .join(
            FILTER_DECOMPOSED_SUBTREES.out.map { og, passing_membership, report ->
                tuple(og, passing_membership)
            }
        )
        .map { og, raw_fasta, passing_membership ->
            tuple(og, raw_fasta, passing_membership)
        }

    WRITE_DECOMPOSED_FASTAS(
        ch_raw_plus_passing_membership,
        write_fastas_script_ch
    )

    /*
     * Candidate units after decomposition:
     *
     * unchanged units:
     *   og_29
     *
     * true decomposed units:
     *   og_123_subtree_001
     */
    ch_candidate_fastas = WRITE_DECOMPOSED_FASTAS.out.fastas
        .flatten()
        .map { fasta ->
            tuple(fasta.baseName, fasta)
        }

    ch_true_subtree_fastas = ch_candidate_fastas.filter { unit, fasta ->
        unit.contains('_subtree_')
    }

    ch_unchanged_units = ch_candidate_fastas.filter { unit, fasta ->
        !unit.contains('_subtree_')
    }

    /*
     * Only true decomposed subtrees are realigned/rebuilt before TreeShrink.
     */
    MACSE_ALIGN_DECOMPOSED(ch_true_subtree_fastas)

    IQTREE_DECOMPOSED(MACSE_ALIGN_DECOMPOSED.out)

    ch_rebuilt_decomposed_iqtree_dirs = IQTREE_DECOMPOSED.out
        .map { unit, iqtree_dir, nt_aln ->
            tuple(unit, iqtree_dir)
        }

    ch_original_iqtree_dirs_for_unchanged = ch_unchanged_units
        .map { unit, fasta ->
            def og = unit.replaceFirst(/^og_/, '')
            tuple(og, unit)
        }
        .join(ch_og_iqtree_dir)
        .map { og, unit, iqtree_dir ->
            tuple(unit, iqtree_dir)
        }

    ch_treeclean_iqtree_dirs_for_treeshrink = ch_rebuilt_decomposed_iqtree_dirs
        .mix(ch_original_iqtree_dirs_for_unchanged)

    ch_iqtree_dirs_for_treeshrink = ch_treeclean_iqtree_dirs_for_treeshrink
        .map { unit, iqtree_dir -> iqtree_dir }
        .collect()

    MAKE_TREESHRINK_INPUTS(
        ch_iqtree_dirs_for_treeshrink,
        make_treeshrink_inputs_script_ch
    )

    WRITE_TREESHRINK_PROTECTED_SPECIES(
        genomes_tsv,
        write_protected_outgroups_script_ch
    )

    RUN_TREESHRINK_DECOMPOSED(
        MAKE_TREESHRINK_INPUTS.out[0],
        MAKE_TREESHRINK_INPUTS.out[1],
        WRITE_TREESHRINK_PROTECTED_SPECIES.out,
        parse_treeshrink_removals_script_ch
    )

    ch_all_candidate_fastas = WRITE_DECOMPOSED_FASTAS.out.fastas
        .flatten()
        .collect()

    WRITE_TREESHRINK_PRUNED_FASTAS(
        ch_all_candidate_fastas,
        RUN_TREESHRINK_DECOMPOSED.out.removals,
        genomes_tsv,
        prune_fastas_script_ch
    )

    FILTER_POST_TREESHRINK_SUBTREES(
        WRITE_TREESHRINK_PRUNED_FASTAS.out.membership,
        genomes_tsv,
        filter_script_ch
    )

    SELECT_FINAL_FASTAS(
        WRITE_TREESHRINK_PRUNED_FASTAS.out.fastas.flatten().collect(),
        FILTER_POST_TREESHRINK_SUBTREES.out[0],
        select_fastas_script_ch
    )

    /*
     * Decide which units need a final rebuild.
     *
     * Rebuild if:
     *   - the unit was decomposed, or
     *   - TreeShrink pruned at least one tip
     *
     * Reuse original/rebuilt input if:
     *   - the unit was not decomposed and TreeShrink pruned nothing
     */
    ch_pruning_status = WRITE_TREESHRINK_PRUNED_FASTAS.out.report
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.cleaning_unit_id,
                row.was_decomposed,
                row.was_pruned
            )
        }

    ch_final_raw_fastas_all = SELECT_FINAL_FASTAS.out.fastas
        .flatten()
        .map { fasta ->
            def unit = fasta.baseName.replaceFirst(/\\.pruned$/, '')
            tuple(unit, fasta)
        }

    ch_final_rebuild_fastas = ch_final_raw_fastas_all
        .join(ch_pruning_status)
        .filter { unit, fasta, was_decomposed, was_pruned ->
            was_decomposed == 'true' || was_pruned == 'true'
        }
        .map { unit, fasta, was_decomposed, was_pruned ->
            tuple(unit, fasta)
        }

    ch_final_passthrough_units = ch_final_raw_fastas_all
        .join(ch_pruning_status)
        .filter { unit, fasta, was_decomposed, was_pruned ->
            was_decomposed != 'true' && was_pruned != 'true'
        }
        .map { unit, fasta, was_decomposed, was_pruned ->
            unit
        }

    /*
     * Rebuild units that were decomposed or TreeShrink-pruned.
     */
    MACSE_ALIGN_TREESHRINK_PRUNED(ch_final_rebuild_fastas)

    IQTREE_CLEANED_FINAL(MACSE_ALIGN_TREESHRINK_PRUNED.out)

    ch_rebuilt_final_for_alerax = IQTREE_CLEANED_FINAL.out.map { unit, iqtree_dir, nt_aln ->
        tuple(unit, iqtree_dir, nt_aln)
    }

    ch_rebuilt_final_for_redip = IQTREE_CLEANED_FINAL.out.map { unit, iqtree_dir, nt_aln ->
        tuple(unit, iqtree_dir)
    }

    ch_rebuilt_final_dirs = IQTREE_CLEANED_FINAL.out.map { unit, iqtree_dir, nt_aln ->
        tuple(unit, iqtree_dir)
    }

    /*
     * Passthrough units reuse the original tree and original NT alignment.
     * These are only unchanged and unpruned original OGs.
     */
    ch_passthrough_keyed = ch_final_passthrough_units
        .map { unit ->
            def unit_id = unit instanceof List ? unit[0].toString() : unit.toString()
            def og = unit_id.replaceFirst(/^og_/, '')
            tuple(og, unit_id)
        }

    ch_passthrough_with_tree = ch_passthrough_keyed
        .join(ch_og_iqtree_dir)
        .map { og, unit, iqtree_dir ->
            tuple(og, unit, iqtree_dir)
        }

    ch_passthrough_for_alerax = ch_passthrough_with_tree
        .join(ch_og_nt_alignment)
        .map { og, unit, iqtree_dir, nt_aln ->
            tuple(unit, iqtree_dir, nt_aln)
        }

    ch_passthrough_for_redip = ch_passthrough_with_tree
        .map { og, unit, iqtree_dir ->
            tuple(unit, iqtree_dir)
        }

    ch_passthrough_dirs = ch_passthrough_with_tree
        .map { og, unit, iqtree_dir ->
            tuple(unit, iqtree_dir)
        }

    TREE_CLEANING_REPORT(
        DECOMPOSE_LONG_BRANCH_TREES.out.map { og, membership, report -> report }.collect(),
        FILTER_DECOMPOSED_SUBTREES.out.map { og, membership, report -> report }.collect(),
        WRITE_TREESHRINK_PRUNED_FASTAS.out.report,
        FILTER_POST_TREESHRINK_SUBTREES.out[1],
        combine_reports_script_ch
    )

    emit:
    cleaned_for_alerax = ch_rebuilt_final_for_alerax.mix(ch_passthrough_for_alerax)

    cleaned_for_redip = ch_rebuilt_final_for_redip.mix(ch_passthrough_for_redip)

    cleaned_dirs = ch_rebuilt_final_dirs.mix(ch_passthrough_dirs)

    report = TREE_CLEANING_REPORT.out
}
