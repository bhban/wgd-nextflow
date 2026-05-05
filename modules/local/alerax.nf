nextflow.enable.dsl=2

process WRITE_ALERAX_MAPPING {
    tag { "og_${og}" }
    array (params.array_size as int)

    input:
    tuple val(og), path(iqtree_dir), path(nt)

    output:
    path("alerax_payload_og_${og}")

    script:
    """
    mkdir -p alerax_payload_og_${og}

    treefile="${iqtree_dir}/og_${og}_iqtree.treefile"
    ufboot="${iqtree_dir}/og_${og}_iqtree.ufboot"
    status="${iqtree_dir}/og_${og}.iqtree.status"

    if [ ! -s ${nt} ]; then
        echo "Missing NT fasta for og_${og}" >&2
        exit 1
    fi

    if [ ! -s "\$treefile" ]; then
        echo "Missing IQ-TREE treefile for og_${og}: \$treefile" >&2
        exit 1
    fi

    if [ ! -s "\$status" ] || [ "\$(tr -d '\\r\\n' < "\$status")" != "OK" ]; then
        echo "IQ-TREE status is not OK for og_${og}: \$status" >&2
        exit 1
    fi

    if [ ! -s "\$ufboot" ]; then
        echo "Missing ufboot for og_${og}: \$ufboot" >&2
        exit 1
    fi

    cp "\$ufboot" alerax_payload_og_${og}/og_${og}_iqtree.ufboot

    awk '
      /^>/{
        h=substr(\$1,2)
        split(h,a,/\\|/)
        if(length(a[1]) == 0) {
          exit 1
        }
        if(!(h in seen)) {
          print h "\\t" a[1]
          seen[h]=1
        }
      }
    ' ${nt} > alerax_payload_og_${og}/og_${og}.mapping.tsv

    test -s alerax_payload_og_${og}/og_${og}.mapping.tsv
    test -s alerax_payload_og_${og}/og_${og}_iqtree.ufboot
    """
}


process WRITE_ALERAX_FAMILIES {
    tag "write_alerax_families"

    input:
    path mapping_payload_dirs

    output:
    path("${params.postdir}/alerax")

    script:
    """
    mkdir -p ${params.postdir}/alerax
    mkdir -p ${params.postdir}/alerax/gene_trees
    mkdir -p ${params.postdir}/alerax/gene_to_species_mapping

    echo "Staged files in WRITE_ALERAX_FAMILIES:" >&2
    find . -maxdepth 4 -type f | head -100 >&2

    while IFS= read -r f; do
        cp "\$f" ${params.postdir}/alerax/gene_to_species_mapping/
    done < <(find . -path './alerax_payload_og_*/og_*.mapping.tsv' -type f | sort)

    while IFS= read -r f; do
        cp "\$f" ${params.postdir}/alerax/gene_trees/
    done < <(find . -path './alerax_payload_og_*/og_*_iqtree.ufboot' -type f | sort)

    {
      echo "[FAMILIES]"

      for uf in \$(find ${params.postdir}/alerax/gene_trees -name 'og_*_iqtree.ufboot' -type f | sort); do
        [ -s "\$uf" ] || continue

        og=\$(basename "\$uf" _iqtree.ufboot | sed 's/^og_//')
        mp=${params.postdir}/alerax/gene_to_species_mapping/og_\${og}.mapping.tsv

        [ -s "\$mp" ] || continue

        echo "- family_\${og}"
        echo "gene_tree = alerax_input/gene_trees/og_\${og}_iqtree.ufboot"
        echo "mapping = alerax_input/gene_to_species_mapping/og_\${og}.mapping.tsv"
      done
    } > ${params.postdir}/alerax/families.txt

    test -s ${params.postdir}/alerax/families.txt
    test -d ${params.postdir}/alerax/gene_trees
    test -d ${params.postdir}/alerax/gene_to_species_mapping

    n_gene_trees=\$(find ${params.postdir}/alerax/gene_trees -name 'og_*_iqtree.ufboot' -type f | wc -l)
    n_mappings=\$(find ${params.postdir}/alerax/gene_to_species_mapping -name 'og_*.mapping.tsv' -type f | wc -l)
    n_families=\$(grep -c '^- family_' ${params.postdir}/alerax/families.txt || true)

    echo "AleRax family input summary:" >&2
    echo "  ufboot files: \$n_gene_trees" >&2
    echo "  mapping files: \$n_mappings" >&2
    echo "  families: \$n_families" >&2

    if [ "\$n_gene_trees" -eq 0 ]; then
        echo "No ufboot files found for AleRax families" >&2
        echo "Contents of work directory:" >&2
        find . -maxdepth 5 -type f | head -200 >&2
        exit 1
    fi

    if [ "\$n_mappings" -eq 0 ]; then
        echo "No mapping files found for AleRax families" >&2
        echo "Contents of work directory:" >&2
        find . -maxdepth 5 -type f | head -200 >&2
        exit 1
    fi

    if [ "\$n_families" -eq 0 ]; then
        echo "No valid AleRax families written" >&2
        exit 1
    fi
    """
}


process WRITE_ALERAX_MANIFEST {
    tag "write_alerax_manifest"

    input:
    val models

    output:
    path("${params.postdir}/alerax/model_manifest.tsv")

    script:
    def lines = models.collect { m ->
        "${m['model_id']}\t${m['rec_model']}\t${m['model_parametrization']}\t${m['gene_tree_samples']}"
    }.join('\n')

    """
    mkdir -p ${params.postdir}/alerax

    printf '%b\\n' 'model_id\\trec_model\\tmodel_parametrization\\tgene_tree_samples' > ${params.postdir}/alerax/model_manifest.tsv
    printf '%b\\n' '${lines}' >> ${params.postdir}/alerax/model_manifest.tsv

    test -s ${params.postdir}/alerax/model_manifest.tsv
    """
}


process RUN_ALERAX {
    tag { "alerax_${model['model_id']}" }

    input:
    tuple path(alerax_input_dir, stageAs: "alerax_input"), path(species_tree), val(model)

    output:
    tuple val(model['model_id']), path("alerax/${model['model_id']}")

    script:
    def model_id = model['model_id']
    def rec_model = model['rec_model']
    def model_param = model['model_parametrization']
    def gene_tree_samples = model['gene_tree_samples']
    def cleanup = params.alerax.cleanup_output ? "true" : "false"

    """
    mkdir -p alerax/${model_id}/output
    mkdir -p mpi_tmp

    test -s alerax_input/families.txt
    test -d alerax_input/gene_trees
    test -d alerax_input/gene_to_species_mapping
    test -s ${species_tree}

    export TMPDIR="\$PWD/mpi_tmp"
    export TEMP="\$PWD/mpi_tmp"
    export TMP="\$PWD/mpi_tmp"
    export OMPI_MCA_orte_tmpdir_base="\$PWD/mpi_tmp"

    mpiexec --mca orte_tmpdir_base "\$PWD/mpi_tmp" -np ${task.cpus} ${params.alerax_bin} \\
      -f alerax_input/families.txt \\
      -s ${species_tree} \\
      -p alerax/${model_id}/output \\
      -r ${rec_model} \\
      --model-parametrization ${model_param} \\
      --gene-tree-samples ${gene_tree_samples} \\
      > alerax/${model_id}/alerax.log 2>&1

    test -d alerax/${model_id}/output
    test -s alerax/${model_id}/alerax.log

    if [ "${cleanup}" = "true" ]; then
      if [ -s alerax/${model_id}/output/reconciliations/totalSpeciesEventCounts.txt ]; then
        rm -rf alerax/${model_id}/output/ccps
        rm -rf alerax/${model_id}/output/reconciliations/all
      fi
    fi

    touch alerax/${model_id}/alerax.done
    """
}


process RUN_ALERAX_RANDOM {
    tag { "alerax_${model['model_id']}" }

    input:
    tuple path(alerax_input_dir, stageAs: "alerax_input"), val(model)

    output:
    tuple val(model['model_id']), path("alerax/${model['model_id']}")

    script:
    def model_id = model['model_id']
    def rec_model = model['rec_model']
    def model_param = model['model_parametrization']
    def gene_tree_samples = model['gene_tree_samples']
    def cleanup = params.alerax.cleanup_output ? "true" : "false"

    """
    mkdir -p alerax/${model_id}/output
    mkdir -p mpi_tmp

    test -s alerax_input/families.txt
    test -d alerax_input/gene_trees
    test -d alerax_input/gene_to_species_mapping

    export TMPDIR="\$PWD/mpi_tmp"
    export TEMP="\$PWD/mpi_tmp"
    export TMP="\$PWD/mpi_tmp"
    export OMPI_MCA_orte_tmpdir_base="\$PWD/mpi_tmp"

    mpiexec --mca orte_tmpdir_base "\$PWD/mpi_tmp" -np ${task.cpus} ${params.alerax_bin} \\
      -f alerax_input/families.txt \\
      -s random \\
      -p alerax/${model_id}/output \\
      -r ${rec_model} \\
      --model-parametrization ${model_param} \\
      --gene-tree-samples ${gene_tree_samples} \\
      > alerax/${model_id}/alerax.log 2>&1

    test -d alerax/${model_id}/output
    test -s alerax/${model_id}/alerax.log

    if [ "${cleanup}" = "true" ]; then
      if [ -s alerax/${model_id}/output/reconciliations/totalSpeciesEventCounts.txt ]; then
        rm -rf alerax/${model_id}/output/ccps
        rm -rf alerax/${model_id}/output/reconciliations/all
      fi
    fi

    touch alerax/${model_id}/alerax.done
    """
}


process ALERAX_REPORT {
    tag "alerax_report"

    input:
    path alerax_results
    path model_manifest

    output:
    path("${params.postdir}/alerax/alerax_report.tsv")
    path("${params.postdir}/alerax/alerax.done")

    script:
    """
    mkdir -p ${params.postdir}/alerax

    {
      echo -e "model_id\\trec_model\\tmodel_parametrization\\tgene_tree_samples\\tstatus\\tresult_dir\\tlog_file"

      for d in */; do
        [ -d "\$d" ] || continue

        model_id=\$(basename "\$d")

        status="FAIL"
        [ -e "\$d/alerax.done" ] && status="OK"

        rec_model=\$(awk -F'\\t' -v id="\$model_id" 'NR>1 && \$1==id {print \$2}' ${model_manifest})
        model_param=\$(awk -F'\\t' -v id="\$model_id" 'NR>1 && \$1==id {print \$3}' ${model_manifest})
        gts=\$(awk -F'\\t' -v id="\$model_id" 'NR>1 && \$1==id {print \$4}' ${model_manifest})

        result_dir=""
        [ -d "\$d/output" ] && result_dir=\$(realpath "\$d/output")

        log_file=""
        [ -s "\$d/alerax.log" ] && log_file=\$(realpath "\$d/alerax.log")

        echo -e "\${model_id}\\t\${rec_model}\\t\${model_param}\\t\${gts}\\t\${status}\\t\${result_dir}\\t\${log_file}"
      done | sort -k1,1
    } > ${params.postdir}/alerax/alerax_report.tsv

    touch ${params.postdir}/alerax/alerax.done
    """
}


workflow ALERAX_WORKFLOW {

    take:
    iqtree_results
    species_tree
    models

    main:

    alerax_map_out = WRITE_ALERAX_MAPPING(iqtree_results)

    families_out = WRITE_ALERAX_FAMILIES(
        alerax_map_out.collect()
    )

    families_dir_ch = families_out
    families_publish_ch = families_dir_ch

    manifest_out = WRITE_ALERAX_MANIFEST(models.collect())
    manifest_file_ch = manifest_out[0]

    def alerax_results_out

    if (params.use_species_tree_for_alerax) {
        alerax_in = models
            .combine(families_dir_ch)
            .combine(species_tree)
            .map { row ->
                def model = row[0][0]
                def families_dir = row[0][1]
                def tree = row[1]

                tuple(families_dir, tree, model)
            }

        alerax_results_out = RUN_ALERAX(alerax_in)

    } else {
        alerax_in = models
            .combine(families_dir_ch)
            .map { model, families_dir ->
                tuple(families_dir, model)
            }

        alerax_results_out = RUN_ALERAX_RANDOM(alerax_in)
    }

    alerax_report_out = ALERAX_REPORT(
        alerax_results_out
            .map { model_id, model_dir ->
                model_dir
            }
            .collect(),
        manifest_file_ch
    )

    alerax_report_tsv_ch = alerax_report_out[0]
    alerax_done_ch = alerax_report_out[1]
    alerax_report_publish_ch = alerax_report_tsv_ch.mix(alerax_done_ch)

    emit:
    families = families_publish_ch
    manifest = manifest_file_ch
    results = alerax_results_out
    report = alerax_report_publish_ch
}
