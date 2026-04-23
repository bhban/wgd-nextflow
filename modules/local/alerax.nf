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

process RUN_ALERAX {
    tag "alerax"

    input:
    path families
    path species_tree

    output:
    path("${params.postdir}/alerax/output")
    path("${params.postdir}/alerax/alerax.done")
    path("${params.postdir}/post_genespace.done")

    script:
    """
    mkdir -p ${params.postdir}/alerax/output
    mkdir -p mpi_tmp

    export TMPDIR="\$PWD/mpi_tmp"
    export TEMP="\$PWD/mpi_tmp"
    export TMP="\$PWD/mpi_tmp"
    export OMPI_MCA_orte_tmpdir_base="\$PWD/mpi_tmp"

    mpiexec --mca orte_tmpdir_base "\$PWD/mpi_tmp" -np ${task.cpus} ${params.alerax_bin} \
      -f ${families} \
      -s ${species_tree} \
      -p ${params.postdir}/alerax/output \
      -r ${params.alerax.rec_model} \
      --model-parametrization ${params.alerax.model_parametrization} \
      --gene-tree-samples ${params.alerax.gene_tree_samples} \
      > alerax.log 2>&1

    test -d ${params.postdir}/alerax/output
    touch ${params.postdir}/alerax/alerax.done
    touch ${params.postdir}/post_genespace.done
    """
}

process RUN_ALERAX_RANDOM {
    tag "alerax"

    input:
    path families

    output:
    path("${params.postdir}/alerax/output")
    path("${params.postdir}/alerax/alerax.done")
    path("${params.postdir}/post_genespace.done")

    script:
    """
    mkdir -p ${params.postdir}/alerax/output
    mkdir -p mpi_tmp

    export TMPDIR="\$PWD/mpi_tmp"
    export TEMP="\$PWD/mpi_tmp"
    export TMP="\$PWD/mpi_tmp"
    export OMPI_MCA_orte_tmpdir_base="\$PWD/mpi_tmp"

    mpiexec --mca orte_tmpdir_base "\$PWD/mpi_tmp" -np ${task.cpus} ${params.alerax_bin} \
      -f ${families} \
      -s random \
      -p ${params.postdir}/alerax/output \
      -r ${params.alerax.rec_model} \
      --model-parametrization ${params.alerax.model_parametrization} \
      --gene-tree-samples ${params.alerax.gene_tree_samples} \
      > alerax.log 2>&1

    test -d ${params.postdir}/alerax/output
    touch ${params.postdir}/alerax/alerax.done
    touch ${params.postdir}/post_genespace.done
    """
}
