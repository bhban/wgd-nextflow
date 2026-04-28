nextflow.enable.dsl=2

process WRITE_ALERAX_MAPPING {
    tag { "og_${og}" }
    array (params.array_size as int)

    input:
    tuple val(og), path(treefile), path(ufboot), path(status), path(nt), path(iqtree_log), path(iqtree_all)

    output:
    tuple val(og),
          path("alerax/gene_to_species_mapping/og_${og}.mapping.tsv"),
          path(ufboot)

    script:
    """
    mkdir -p alerax/gene_to_species_mapping

    if [ ! -s ${nt} ]; then
        echo "Missing NT fasta for og_${og}" > alerax/gene_to_species_mapping/og_${og}.mapping.log
        exit 1
    fi

    if [ ! -s ${status} ] || [ "\$(tr -d '\\r' < ${status})" != "OK" ]; then
        echo "IQ-TREE status is not OK for og_${og}" > alerax/gene_to_species_mapping/og_${og}.mapping.log
        exit 1
    fi

    if [ ! -s ${ufboot} ]; then
        echo "Missing ufboot for og_${og}" > alerax/gene_to_species_mapping/og_${og}.mapping.log
        exit 1
    fi

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
    ' ${nt} > alerax/gene_to_species_mapping/og_${og}.mapping.tsv

    test -s alerax/gene_to_species_mapping/og_${og}.mapping.tsv
    test -s ${ufboot}
    """
}

process WRITE_ALERAX_FAMILIES {
    tag "write_alerax_families"

    input:
    path mapping_results

    output:
    path("${params.postdir}/alerax/families.txt")
    path("${params.postdir}/alerax/gene_to_species_mapping")

    script:
    """
    mkdir -p ${params.postdir}/alerax
    mkdir -p ${params.postdir}/alerax/gene_to_species_mapping

    find . -name 'og_*.mapping.tsv' -type f -exec cp {} ${params.postdir}/alerax/gene_to_species_mapping/ \\;

    {
      echo "[FAMILIES]"

      for uf in \$(find . -name 'og_*_iqtree.ufboot' -type f | sort); do
        [ -s "\$uf" ] || continue

        og=\$(basename "\$uf" _iqtree.ufboot | sed 's/^og_//')
        mp=${params.postdir}/alerax/gene_to_species_mapping/og_\${og}.mapping.tsv

        [ -s "\$mp" ] || continue

        echo "- family_\${og}"
        echo "gene_tree = \$(realpath "\$uf")"
        echo "mapping = \$(realpath "\$mp")"
      done
    } > ${params.postdir}/alerax/families.txt

    test -s ${params.postdir}/alerax/families.txt
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
        "${m.model_id}\t${m.rec_model}\t${m.model_parametrization}\t${m.gene_tree_samples}"
    }.join('\n')

    """
    mkdir -p ${params.postdir}/alerax

    {
      echo -e "model_id\trec_model\tmodel_parametrization\tgene_tree_samples"
      cat <<'EOF'
${lines}
EOF
    } > ${params.postdir}/alerax/model_manifest.tsv

    test -s ${params.postdir}/alerax/model_manifest.tsv
    """
}

process RUN_ALERAX {
    tag { "alerax_${model.model_id}" }

    input:
    tuple path(families), path(species_tree), val(model)

    output:
    tuple val(model.model_id), path("alerax/${model.model_id}")

    script:
    def model_id = model.model_id
    def rec_model = model.rec_model
    def model_param = model.model_parametrization
    def gene_tree_samples = model.gene_tree_samples
    def cleanup = params.alerax.cleanup_output ? "true" : "false"

    """
    mkdir -p alerax/${model_id}/output
    mkdir -p mpi_tmp

    export TMPDIR="\$PWD/mpi_tmp"
    export TEMP="\$PWD/mpi_tmp"
    export TMP="\$PWD/mpi_tmp"
    export OMPI_MCA_orte_tmpdir_base="\$PWD/mpi_tmp"

    mpiexec --mca orte_tmpdir_base "\$PWD/mpi_tmp" -np ${task.cpus} ${params.alerax_bin} \\
      -f ${families} \\
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
    tag { "alerax_${model.model_id}" }

    input:
    tuple path(families), val(model)

    output:
    tuple val(model.model_id), path("alerax/${model.model_id}")

    script:
    def model_id = model.model_id
    def rec_model = model.rec_model
    def model_param = model.model_parametrization
    def gene_tree_samples = model.gene_tree_samples
    def cleanup = params.alerax.cleanup_output ? "true" : "false"

    """
    mkdir -p alerax/${model_id}/output
    mkdir -p mpi_tmp

    export TMPDIR="\$PWD/mpi_tmp"
    export TEMP="\$PWD/mpi_tmp"
    export TMP="\$PWD/mpi_tmp"
    export OMPI_MCA_orte_tmpdir_base="\$PWD/mpi_tmp"

    mpiexec --mca orte_tmpdir_base "\$PWD/mpi_tmp" -np ${task.cpus} ${params.alerax_bin} \\
      -f ${families} \\
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
        alerax_map_out
            .flatMap { og, mapping, ufboot -> [mapping, ufboot] }
            .collect()
    )

    manifest_out = WRITE_ALERAX_MANIFEST(models)

    def alerax_results_out

    if (species_tree) {
        alerax_in = models
            .combine(families_out)
            .combine(species_tree)
            .map { model, fam, tree ->
                tuple(fam, tree, model)
            }

        alerax_results_out = RUN_ALERAX(alerax_in)

    } else {
        alerax_in = models
            .combine(families_out)
            .map { model, fam ->
                tuple(fam, model)
            }

        alerax_results_out = RUN_ALERAX_RANDOM(alerax_in)
    }

    alerax_report_out = ALERAX_REPORT(
        alerax_results_out
            .map { model_id, model_dir -> model_dir }
            .collect(),
        manifest_out
    )

    emit:
    families = families_out
    manifest = manifest_out
    results = alerax_results_out
    report = alerax_report_out
}
