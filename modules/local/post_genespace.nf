process PANGENES_PASS_FILTER {
    tag "pangenes_pass_filter"

    input:
    path genespace_wd
    path genespace_done
    path genomes_tsv
    path pangenes_pass_filter_script

    output:
    path("${params.postdir}/pangenes_PASS.tsv")
    path("${params.postdir}/og_list_min4species.txt")

    script:
    def requireOutgroupArg = params.require_outgroup_og ? "--require-outgroup" : ""
    """
    mkdir -p ${params.postdir}

    python ${pangenes_pass_filter_script} \
      --genespace-wd ${genespace_wd} \
      --genomes-tsv ${genomes_tsv} \
      --out-tsv ${params.postdir}/pangenes_PASS.tsv \
      --out-og-list ${params.postdir}/og_list_min4species.txt \
      ${requireOutgroupArg} \
      > pangenes_pass_filter.log 2>&1

    test -s ${params.postdir}/pangenes_PASS.tsv
    test -s ${params.postdir}/og_list_min4species.txt
    """
}

process COLLAPSE_TANDEMS {
    tag "collapse_tandems"

    input:
    path pass_tsv
    path genespace_wd
    path genespace_done
    path genomes_tsv
    path collapse_tandems_script

    output:
    path("${params.postdir}/pangenes_PASS.collapsed.tsv")
    path("${params.postdir}/tandem_report.tsv")
    path("${params.postdir}/og_list_min4species.collapsed.txt")

    script:
    def requireOutgroupArg = params.require_outgroup_og ? "--require-outgroup" : ""
    """
    mkdir -p ${params.postdir}

    test -s ${genespace_wd}/results/combBed.txt

    python ${collapse_tandems_script} \
      --infile ${pass_tsv} \
      --combBed ${genespace_wd}/results/combBed.txt \
      --genomes-tsv ${genomes_tsv} \
      --outfile_filtered ${params.postdir}/pangenes_PASS.collapsed.tsv \
      --outfile_tandems ${params.postdir}/tandem_report.tsv \
      --outfile_og_list ${params.postdir}/og_list_min4species.collapsed.txt \
      ${requireOutgroupArg} \
      > collapse_tandems.log 2>&1

    test -s ${params.postdir}/pangenes_PASS.collapsed.tsv
    test -s ${params.postdir}/tandem_report.tsv
    test -s ${params.postdir}/og_list_min4species.collapsed.txt
    """
}

process WRITE_OG_FASTAS {
    tag "write_og_fastas"

    input:
    path pass_tsv
    path og_list
    path genomes_tsv
    path cds_files
    path write_og_fastas_script

    output:
    path("${params.postdir}/og_fasta")
    path("${params.postdir}/*min4species*.txt")
    path("${params.postdir}/og_fastas.done")

    script:
    def cdsList = cds_files.collect { "\"${it}\"" }.join(' ')
    def ogListOut = "${params.postdir}/${og_list.getFileName()}"
    """
    mkdir -p ${params.postdir}/og_fasta
    mkdir -p cds

    for f in ${cdsList}; do
      cp "\$f" cds/
    done

    python ${write_og_fastas_script} \
      --pangenes-pass ${pass_tsv} \
      --genomes-tsv ${genomes_tsv} \
      --cds-dir cds \
      --outdir ${params.postdir}/og_fasta \
      --og-list ${ogListOut} \
      > write_og_fastas.log 2>&1

    compgen -G "${params.postdir}/og_fasta/og_*.fasta" > /dev/null
    test -s ${ogListOut}
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

    test -s ${params.postdir}/macse_report.tsv
    test -s ${params.postdir}/macse_ok_og_list.txt
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

    if [ ! -s ${status} ] || [ "\$(tr -d '\\r' < ${status})" != "OK" ]; then
        echo "FAIL" > og_${og}.iqtree.status
        : > og_${og}_iqtree.treefile
        : > og_${og}_iqtree.ufboot
        trap - EXIT TERM INT
        exit 0
    fi

    set +e
    ${params.iqtree_bin} \
      -s ${nt} \
      -T ${task.cpus} \
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

    test -s ${params.postdir}/iqtree_report.tsv
    test -s ${params.postdir}/iqtree_ok_og_list.txt
    """
}
