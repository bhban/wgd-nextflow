process PANGENES_PASS_FILTER {
    tag "pangenes_pass_filter"

    input:
    path genespace_wd
    path genespace_done

    output:
    path("${params.postdir}/pangenes_PASS.tsv")
    path("${params.postdir}/og_list_min4species.txt")

    script:
    """
    mkdir -p ${params.postdir}

    python ${projectDir}/scripts/pangenes_pass_filter.py \
      --genespace-wd genespace/${params.working_dir} \
      --out-tsv ${params.postdir}/pangenes_PASS.tsv \
      --out-og-list ${params.postdir}/og_list_min4species.txt \
      > pangenes_pass_filter.log 2>&1

    test -s ${params.postdir}/pangenes_PASS.tsv
    test -s ${params.postdir}/og_list_min4species.txt
    """
}

process WRITE_OG_FASTAS {
    tag "write_og_fastas"

    input:
    path pass_tsv
    path og_list
    path genomes_tsv
    path cds_files

    output:
    path("${params.postdir}/og_fasta")
    path("${params.postdir}/og_list_min4species.txt")
    path("${params.postdir}/og_fastas.done")

    script:
    """
    mkdir -p ${params.postdir}/og_fasta
    mkdir -p cds
    cp ${cds_files.join(' ')} cds/

    python ${projectDir}/scripts/write_og_fastas.py \
      --pangenes-pass ${pass_tsv} \
      --genomes-tsv ${genomes_tsv} \
      --cds-dir cds \
      --outdir ${params.postdir}/og_fasta \
      --og-list ${params.postdir}/og_list_min4species.txt \
      > write_og_fastas.log 2>&1

    find ${params.postdir}/og_fasta -maxdepth 1 -name 'og_*.fasta' | grep -q .
    touch ${params.postdir}/og_fastas.done
    """
}

process MACSE_ALIGN_OG {
    tag { "og_${og}" }
    array (params.array_size as int)

    input:
    tuple val(og), path(fasta)

    output:
    tuple val(og),
          path("og_${og}_AA.fasta"),
          path("og_${og}_NT.fasta"),
          path("og_${og}.status")

    script:
    """
    rm -f og_${og}.status og_${og}_AA.fasta og_${og}_NT.fasta

    echo "STARTED" > og_${og}.status

    on_exit() {
        rc=\$?
        if [ ! -s og_${og}.status ] || grep -qx 'STARTED' og_${og}.status; then
            rm -f og_${og}_AA.fasta og_${og}_NT.fasta
            : > og_${og}_AA.fasta
            : > og_${og}_NT.fasta
            echo "FAIL" > og_${og}.status
        fi
        exit \$rc
    }
    trap on_exit EXIT TERM INT

    set +e
    ${params.macse_bin} -prog alignSequences \
      -seq ${fasta} \
      -out_AA og_${og}_AA.fasta \
      -out_NT og_${og}_NT.fasta \
      > og_${og}.log 2>&1
    rc=\$?
    set -e

    if [ \$rc -eq 0 ] && [ -s og_${og}_AA.fasta ] && [ -s og_${og}_NT.fasta ]; then
        echo "OK" > og_${og}.status
    else
        rm -f og_${og}_AA.fasta og_${og}_NT.fasta
        : > og_${og}_AA.fasta
        : > og_${og}_NT.fasta
        echo "FAIL" > og_${og}.status
    fi

    trap - EXIT TERM INT
    """
}

process MACSE_REPORT {
    tag "macse_report"

    input:
    path macse_results

    output:
    path("${params.postdir}/macse_report.tsv")
    path("${params.postdir}/macse_ok_og_list.txt")

    script:
    """
    mkdir -p ${params.postdir}

    ls -1 *.status >/dev/null 2>&1

    {
      echo -e "og\tstatus"
      for f in *.status; do
        [ -e "\$f" ] || continue
        og=\$(basename "\$f" .status | sed 's/^og_//')
        st=\$(tr -d '\\r' < "\$f")
        [ -n "\$st" ] || st="FAIL"
        echo -e "\${og}\t\${st}"
      done | sort -k1,1n
    } > ${params.postdir}/macse_report.tsv

    awk 'NR>1 && \$2=="OK"{print \$1}' ${params.postdir}/macse_report.tsv \
      > ${params.postdir}/macse_ok_og_list.txt
    """
}

process IQTREE_OG {
    tag { "og_${og}" }
    array (params.array_size as int)

    input:
    tuple val(og), path(aa), path(nt), path(status)

    output:
    tuple val(og),
          path("og_${og}_iqtree.treefile"),
          path("og_${og}_iqtree.ufboot"),
          path("og_${og}.iqtree.status"),
          path("og_${og}_NT.fasta")

    script:
    """
    rm -f og_${og}.iqtree.status og_${og}_iqtree.treefile og_${og}_iqtree.ufboot
    cp ${nt} og_${og}_NT.fasta

    echo "STARTED" > og_${og}.iqtree.status

    on_exit() {
        rc=\$?
        if [ ! -s og_${og}.iqtree.status ] || grep -qx 'STARTED' og_${og}.iqtree.status; then
            : > og_${og}_iqtree.treefile
            : > og_${og}_iqtree.ufboot
            echo "FAIL" > og_${og}.iqtree.status
        fi
        exit \$rc
    }
    trap on_exit EXIT TERM INT

    if [ ! -s ${nt} ]; then
        echo "FAIL" > og_${og}.iqtree.status
        : > og_${og}_iqtree.treefile
        : > og_${og}_iqtree.ufboot
        trap - EXIT TERM INT
        exit 0
    fi

    if [ "\$(tr -d '\\r' < ${status})" != "OK" ]; then
        echo "FAIL" > og_${og}.iqtree.status
        : > og_${og}_iqtree.treefile
        : > og_${og}_iqtree.ufboot
        trap - EXIT TERM INT
        exit 0
    fi

    set +e
    ${params.iqtree_bin} \
      -s ${nt} \
      -nt ${task.cpus} \
      -m MFP \
      -bb 1000 \
      -wbtl \
      -redo \
      -pre og_${og}_iqtree \
      > og_${og}.log 2>&1
    rc=\$?
    set -e

    if [ \$rc -eq 0 ] && [ -s og_${og}_iqtree.treefile ] && [ -s og_${og}_iqtree.ufboot ]; then
        echo "OK" > og_${og}.iqtree.status
    else
        : > og_${og}_iqtree.treefile
        : > og_${og}_iqtree.ufboot
        echo "FAIL" > og_${og}.iqtree.status
    fi

    trap - EXIT TERM INT
    """
}

process IQTREE_REPORT {
    tag "iqtree_report"

    input:
    path iqtree_results

    output:
    path("${params.postdir}/iqtree_report.tsv")
    path("${params.postdir}/iqtree_ok_og_list.txt")

    script:
    """
    mkdir -p ${params.postdir}

    ls -1 *.iqtree.status >/dev/null 2>&1

    {
      echo -e "og\tstatus"
      for f in *.iqtree.status; do
        [ -e "\$f" ] || continue
        og=\$(basename "\$f" .iqtree.status | sed 's/^og_//')
        st=\$(tr -d '\\r' < "\$f")
        [ -n "\$st" ] || st="FAIL"
        echo -e "\${og}\t\${st}"
      done | sort -k1,1n
    } > ${params.postdir}/iqtree_report.tsv

    awk 'NR>1 && \$2=="OK"{print \$1}' ${params.postdir}/iqtree_report.tsv \
      > ${params.postdir}/iqtree_ok_og_list.txt
    """
}

process WRITE_ALERAX_MAPPING {
    tag { "og_${og}" }

    input:
    tuple val(og), path(treefile), path(ufboot), path(status), path(nt)

    output:
    tuple val(og), path("og_${og}.mapping.tsv"), path("og_${og}_iqtree.ufboot")

    script:
    """
    cp ${ufboot} og_${og}_iqtree.ufboot

    if [ ! -s ${nt} ]; then
        echo "Missing NT fasta for og_${og}" > og_${og}.mapping.log
        exit 1
    fi

    if [ ! -s ${status} ] || [ "\$(tr -d '\\r' < ${status})" != "OK" ]; then
        echo "IQ-TREE status is not OK for og_${og}" > og_${og}.mapping.log
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

    mpiexec -np ${task.cpus} ${params.alerax_bin} \
      -f ${families} \
      -s ${species_tree} \
      -p ${params.postdir}/alerax/output \
      -r ${params.alerax.rec_model} \
      --model-parametrization ${params.alerax.model_parametrization} \
      --gene-tree-samples ${params.alerax.gene_tree_samples} \
      > alerax.log 2>&1

    touch ${params.postdir}/alerax/alerax.done
    touch ${params.postdir}/post_genespace.done
    """
}
