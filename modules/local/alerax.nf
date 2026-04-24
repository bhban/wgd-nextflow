process WRITE_ALERAX_MAPPING {
    tag { "og_${og}" }
    array (params.array_size as int)

    input:
    tuple val(og), path(treefile), path(ufboot), path(status), path(nt)

    output:
    tuple val(og), path("og_${og}.mapping.tsv"), path("og_${og}_iqtree.ufboot")

    script:
    """
    if [ ! -s ${nt} ]; then
        echo "Missing NT fasta for og_${og}" > og_${og}.mapping.log
        exit 1
    fi

    if [ ! -s ${status} ] || [ "\$(tr -d '\\r' < ${status})" != "OK" ]; then
        echo "IQ-TREE status is not OK for og_${og}" > og_${og}.mapping.log
        exit 1
    fi

    if [ ! -s ${ufboot} ]; then
        echo "Missing ufboot for og_${og}" > og_${og}.mapping.log
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
    ' ${nt} > og_${og}.mapping.tsv

    test -s og_${og}.mapping.tsv
    test -s ${ufboot}
    """
}

process WRITE_ALERAX_FAMILIES {
    tag "write_alerax_families"

    input:
    path mapping_results

    output:
    path("${params.postdir}/alerax/families.txt")

    script:
    """
    mkdir -p ${params.postdir}/alerax
    mkdir -p ${params.postdir}/iqtree
    mkdir -p ${params.postdir}/alerax/gene_to_species_mapping

    for f in og_*_iqtree.ufboot; do
      [ -e "\$f" ] || continue
      cp "\$f" ${params.postdir}/iqtree/
    done

    for f in og_*.mapping.tsv; do
      [ -e "\$f" ] || continue
      cp "\$f" ${params.postdir}/alerax/gene_to_species_mapping/
    done

    {
      echo "[FAMILIES]"
      for uf in ${params.postdir}/iqtree/og_*_iqtree.ufboot; do
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
    path("model_manifest.tsv")

    script:
    def lines = models.collect { m ->
        "${m.model_id}\t${m.rec_model}\t${m.model_parametrization}\t${m.gene_tree_samples}"
    }.join('\n')

    """
    {
      echo -e "model_id\trec_model\tmodel_parametrization\tgene_tree_samples"
      cat <<'EOF'
${lines}
EOF
    } > model_manifest.tsv

    test -s model_manifest.tsv
    """
}

process RUN_ALERAX {
    cache false
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

    mpiexec --mca orte_tmpdir_base "\$PWD/mpi_tmp" -np ${task.cpus} ${params.alerax_bin} \
      -f ${families} \
      -s ${species_tree} \
      -p alerax/${model_id}/output \
      -r ${rec_model} \
      --model-parametrization ${model_param} \
      --gene-tree-samples ${gene_tree_samples} \
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

    mpiexec --mca orte_tmpdir_base "\$PWD/mpi_tmp" -np ${task.cpus} ${params.alerax_bin} \
      -f ${families} \
      -s random \
      -p alerax/${model_id}/output \
      -r ${rec_model} \
      --model-parametrization ${model_param} \
      --gene-tree-samples ${gene_tree_samples} \
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
    path("${params.postdir}/alerax_report.tsv")
    path("${params.postdir}/alerax.done")

    script:
    """
    mkdir -p ${params.postdir}

    {
      echo -e "model_id\trec_model\tmodel_parametrization\tgene_tree_samples\tstatus\tresult_dir\tlog_file"
      for d in alerax/*; do
        [ -d "\$d" ] || continue
        model_id=\$(basename "\$d")

        status="FAIL"
        [ -s "\$d/alerax.done" ] && status="OK"

        rec_model=\$(awk -F'\\t' -v id="\$model_id" 'NR>1 && \$1==id {print \$2}' model_manifest.tsv)
        model_param=\$(awk -F'\\t' -v id="\$model_id" 'NR>1 && \$1==id {print \$3}' model_manifest.tsv)
        gts=\$(awk -F'\\t' -v id="\$model_id" 'NR>1 && \$1==id {print \$4}' model_manifest.tsv)

        result_dir=""
        [ -d "\$d/output" ] && result_dir=\$(realpath "\$d/output")

        log_file=""
        [ -s "\$d/alerax.log" ] && log_file=\$(realpath "\$d/alerax.log")

        echo -e "\${model_id}\t\${rec_model}\t\${model_param}\t\${gts}\t\${status}\t\${result_dir}\t\${log_file}"
      done | sort -k1,1
    } > ${params.postdir}/alerax_report.tsv

    touch ${params.postdir}/alerax.done
    """
}
