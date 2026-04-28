process PANGENES_PASS_FILTER {
    tag "pangenes_pass_filter"

    input:
    path genespace_wd
    path genespace_done
    path genomes_tsv
    path pangenes_pass_filter_script

    output:
    path("${params.postdir}/pangenes/pangenes_PASS.tsv")
    path("${params.postdir}/pangenes/og_list_min4species.txt")
    path("${params.postdir}/pangenes/pangenes_pass_filter.log")

    script:
    def requireOutgroupArg = params.require_outgroup_og ? "--require-outgroup" : ""
    """
    mkdir -p ${params.postdir}/pangenes

    python ${pangenes_pass_filter_script} \\
      --genespace-wd ${genespace_wd} \\
      --genomes-tsv ${genomes_tsv} \\
      --out-tsv ${params.postdir}/pangenes/pangenes_PASS.tsv \\
      --out-og-list ${params.postdir}/pangenes/og_list_min4species.txt \\
      ${requireOutgroupArg} \\
      > ${params.postdir}/pangenes/pangenes_pass_filter.log 2>&1

    test -s ${params.postdir}/pangenes/pangenes_PASS.tsv
    test -s ${params.postdir}/pangenes/og_list_min4species.txt
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
    path("${params.postdir}/pangenes/pangenes_PASS.collapsed.tsv")
    path("${params.postdir}/pangenes/og_list_min4species.collapsed.txt")
    path("${params.postdir}/tandem_collapse/tandem_report.tsv")
    path("${params.postdir}/tandem_collapse/collapse_tandems.log")

    script:
    def requireOutgroupArg = params.require_outgroup_og ? "--require-outgroup" : ""
    """
    mkdir -p ${params.postdir}/pangenes
    mkdir -p ${params.postdir}/tandem_collapse

    test -s ${genespace_wd}/results/combBed.txt

    python ${collapse_tandems_script} \\
      --infile ${pass_tsv} \\
      --combBed ${genespace_wd}/results/combBed.txt \\
      --genomes-tsv ${genomes_tsv} \\
      --outfile_filtered ${params.postdir}/pangenes/pangenes_PASS.collapsed.tsv \\
      --outfile_tandems ${params.postdir}/tandem_collapse/tandem_report.tsv \\
      --outfile_og_list ${params.postdir}/pangenes/og_list_min4species.collapsed.txt \\
      ${requireOutgroupArg} \\
      > ${params.postdir}/tandem_collapse/collapse_tandems.log 2>&1

    test -s ${params.postdir}/pangenes/pangenes_PASS.collapsed.tsv
    test -s ${params.postdir}/tandem_collapse/tandem_report.tsv
    test -s ${params.postdir}/pangenes/og_list_min4species.collapsed.txt
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
    path("${params.postdir}/pangenes/*min4species*.txt")
    path("${params.postdir}/og_fasta/og_fastas.done")
    path("${params.postdir}/og_fasta/write_og_fastas.log")

    script:
    def cdsList = cds_files.collect { "\"${it}\"" }.join(' ')
    def ogListOut = "${params.postdir}/pangenes/${og_list.getFileName()}"
    """
    mkdir -p ${params.postdir}/og_fasta
    mkdir -p ${params.postdir}/pangenes
    mkdir -p cds

    for f in ${cdsList}; do
      cp "\$f" cds/
    done

    python ${write_og_fastas_script} \\
      --pangenes-pass ${pass_tsv} \\
      --genomes-tsv ${genomes_tsv} \\
      --cds-dir cds \\
      --outdir ${params.postdir}/og_fasta \\
      --og-list ${ogListOut} \\
      > ${params.postdir}/og_fasta/write_og_fastas.log 2>&1

    compgen -G "${params.postdir}/og_fasta/og_*.fasta" > /dev/null
    test -s ${ogListOut}
    touch ${params.postdir}/og_fasta/og_fastas.done
    """
}

process MACSE_ALIGN_OG {
    tag { "og_${og}" }
    array (params.array_size as int)

    input:
    tuple val(og), path(fasta)

    output:
    tuple val(og),
          path("macse/og_${og}_AA.fasta"),
          path("macse/og_${og}_NT.fasta"),
          path("macse/og_${og}.status"),
          path("macse/og_${og}.log")

    script:
    """
    mkdir -p macse

    rm -f macse/og_${og}.status macse/og_${og}_AA.fasta macse/og_${og}_NT.fasta macse/og_${og}.log

    echo "STARTED" > macse/og_${og}.status

    on_exit() {
        rc=\$?
        if [ ! -s macse/og_${og}.status ] || grep -qx 'STARTED' macse/og_${og}.status; then
            rm -f macse/og_${og}_AA.fasta macse/og_${og}_NT.fasta
            : > macse/og_${og}_AA.fasta
            : > macse/og_${og}_NT.fasta
            echo "FAIL" > macse/og_${og}.status
        fi
        exit \$rc
    }
    trap on_exit EXIT TERM INT

    set +e
    ${params.macse_bin} -prog alignSequences \\
      -seq ${fasta} \\
      -out_AA macse/og_${og}_AA.fasta \\
      -out_NT macse/og_${og}_NT.fasta \\
      > macse/og_${og}.log 2>&1
    rc=\$?
    set -e

    if [ \$rc -eq 0 ] && [ -s macse/og_${og}_AA.fasta ] && [ -s macse/og_${og}_NT.fasta ]; then
        echo "OK" > macse/og_${og}.status
    else
        rm -f macse/og_${og}_AA.fasta macse/og_${og}_NT.fasta
        : > macse/og_${og}_AA.fasta
        : > macse/og_${og}_NT.fasta
        echo "FAIL" > macse/og_${og}.status
    fi

    trap - EXIT TERM INT
    """
}

process MACSE_REPORT {
    tag "macse_report"

    input:
    path macse_results

    output:
    path("${params.postdir}/macse/macse_report.tsv")
    path("${params.postdir}/macse/macse_ok_og_list.txt")

    script:
    """
    mkdir -p ${params.postdir}/macse

    find . -name '*.status' -type f | sort > macse_status_files.txt
    test -s macse_status_files.txt

    {
      echo -e "og\tstatus"
      while read -r f; do
        og=\$(basename "\$f" .status | sed 's/^og_//')
        st=\$(tr -d '\\r' < "\$f")
        [ -n "\$st" ] || st="FAIL"
        echo -e "\${og}\t\${st}"
      done < macse_status_files.txt | sort -k1,1n
    } > ${params.postdir}/macse/macse_report.tsv

    awk 'NR>1 && \$2=="OK"{print \$1}' ${params.postdir}/macse/macse_report.tsv \\
      > ${params.postdir}/macse/macse_ok_og_list.txt

    test -s ${params.postdir}/macse/macse_report.tsv
    test -s ${params.postdir}/macse/macse_ok_og_list.txt
    """
}

process IQTREE_OG {
    tag { "og_${og}" }
    array (params.array_size as int)

    input:
    tuple val(og), path(aa), path(nt), path(status), path(macse_log)

    output:
    tuple val(og),
          path("iqtree/og_${og}_iqtree.treefile"),
          path("iqtree/og_${og}_iqtree.ufboot"),
          path("iqtree/og_${og}.iqtree.status"),
          path("iqtree/og_${og}_NT.fasta"),
          path("iqtree/og_${og}.log"),
          path("iqtree/og_${og}_iqtree.*")

    script:
    """
    mkdir -p iqtree

    rm -f iqtree/og_${og}.iqtree.status iqtree/og_${og}_iqtree.treefile iqtree/og_${og}_iqtree.ufboot iqtree/og_${og}.log

    cp ${nt} iqtree/og_${og}_NT.fasta

    echo "STARTED" > iqtree/og_${og}.iqtree.status

    on_exit() {
        rc=\$?
        if [ ! -s iqtree/og_${og}.iqtree.status ] || grep -qx 'STARTED' iqtree/og_${og}.iqtree.status; then
            : > iqtree/og_${og}_iqtree.treefile
            : > iqtree/og_${og}_iqtree.ufboot
            echo "FAIL" > iqtree/og_${og}.iqtree.status
        fi
        exit \$rc
    }
    trap on_exit EXIT TERM INT

    if [ ! -s ${nt} ]; then
        echo "FAIL" > iqtree/og_${og}.iqtree.status
        : > iqtree/og_${og}_iqtree.treefile
        : > iqtree/og_${og}_iqtree.ufboot
        trap - EXIT TERM INT
        exit 0
    fi

    if [ ! -s ${status} ] || [ "\$(tr -d '\\r' < ${status})" != "OK" ]; then
        echo "FAIL" > iqtree/og_${og}.iqtree.status
        : > iqtree/og_${og}_iqtree.treefile
        : > iqtree/og_${og}_iqtree.ufboot
        trap - EXIT TERM INT
        exit 0
    fi

    set +e
    ${params.iqtree_bin} \\
      -s ${nt} \\
      -T ${task.cpus} \\
      -m MFP \\
      -bb 1000 \\
      -wbtl \\
      -redo \\
      -pre iqtree/og_${og}_iqtree \\
      > iqtree/og_${og}.log 2>&1
    rc=\$?
    set -e

    if [ \$rc -eq 0 ] && [ -s iqtree/og_${og}_iqtree.treefile ] && [ -s iqtree/og_${og}_iqtree.ufboot ]; then
        echo "OK" > iqtree/og_${og}.iqtree.status
    else
        : > iqtree/og_${og}_iqtree.treefile
        : > iqtree/og_${og}_iqtree.ufboot
        echo "FAIL" > iqtree/og_${og}.iqtree.status
    fi

    trap - EXIT TERM INT
    """
}

process IQTREE_REPORT {
    tag "iqtree_report"

    input:
    path iqtree_results

    output:
    path("${params.postdir}/iqtree/iqtree_report.tsv")
    path("${params.postdir}/iqtree/iqtree_ok_og_list.txt")

    script:
    """
    mkdir -p ${params.postdir}/iqtree

    find . -name '*.iqtree.status' -type f | sort > iqtree_status_files.txt
    test -s iqtree_status_files.txt

    {
      echo -e "og\tstatus"
      while read -r f; do
        og=\$(basename "\$f" .iqtree.status | sed 's/^og_//')
        st=\$(tr -d '\\r' < "\$f")
        [ -n "\$st" ] || st="FAIL"
        echo -e "\${og}\t\${st}"
      done < iqtree_status_files.txt | sort -k1,1n
    } > ${params.postdir}/iqtree/iqtree_report.tsv

    awk 'NR>1 && \$2=="OK"{print \$1}' ${params.postdir}/iqtree/iqtree_report.tsv \\
      > ${params.postdir}/iqtree/iqtree_ok_og_list.txt

    test -s ${params.postdir}/iqtree/iqtree_report.tsv
    test -s ${params.postdir}/iqtree/iqtree_ok_og_list.txt
    """
}
