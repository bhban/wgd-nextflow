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

    python ${pangenes_pass_filter_script} \
      --genespace-wd ${genespace_wd} \
      --genomes-tsv ${genomes_tsv} \
      --out-tsv ${params.postdir}/pangenes/pangenes_PASS.tsv \
      --out-og-list ${params.postdir}/pangenes/og_list_min4species.txt \
      ${requireOutgroupArg} \
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

    python ${collapse_tandems_script} \
      --infile ${pass_tsv} \
      --combBed ${genespace_wd}/results/combBed.txt \
      --genomes-tsv ${genomes_tsv} \
      --outfile_filtered ${params.postdir}/pangenes/pangenes_PASS.collapsed.tsv \
      --outfile_tandems ${params.postdir}/tandem_collapse/tandem_report.tsv \
      --outfile_og_list ${params.postdir}/pangenes/og_list_min4species.collapsed.txt \
      ${requireOutgroupArg} \
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
    path protein_files
    path write_og_fastas_script

    output:
    path("${params.postdir}/og_fasta_cds")
    path("${params.postdir}/og_fasta_aa")
    path("${params.postdir}/pangenes/*min4species*.txt")
    path("${params.postdir}/og_fasta/og_fastas.done")
    path("${params.postdir}/og_fasta/write_og_fastas.log")

    script:
    def cdsList = cds_files.collect { "\"${it}\"" }.join(' ')
    def proteinList = protein_files.collect { "\"${it}\"" }.join(' ')
    def ogListOut = "${params.postdir}/pangenes/${og_list.getFileName()}"
    """
    mkdir -p ${params.postdir}/og_fasta
    mkdir -p ${params.postdir}/og_fasta_cds
    mkdir -p ${params.postdir}/og_fasta_aa
    mkdir -p ${params.postdir}/pangenes
    mkdir -p cds
    mkdir -p protein

    for f in ${cdsList}; do
      cp "\$f" cds/
    done

    for f in ${proteinList}; do
      cp "\$f" protein/
    done

    python ${write_og_fastas_script} \
      --pangenes-pass ${pass_tsv} \
      --genomes-tsv ${genomes_tsv} \
      --cds-dir cds \
      --protein-dir protein \
      --outdir-cds ${params.postdir}/og_fasta_cds \
      --outdir-aa ${params.postdir}/og_fasta_aa \
      --og-list ${ogListOut} \
      > ${params.postdir}/og_fasta/write_og_fastas.log 2>&1

    compgen -G "${params.postdir}/og_fasta_cds/og_*.fasta" > /dev/null
    compgen -G "${params.postdir}/og_fasta_aa/og_*.fasta" > /dev/null
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
          path("${params.postdir}/macse/og_${og}_AA.fasta"),
          path("${params.postdir}/macse/og_${og}_NT.fasta"),
          path("${params.postdir}/macse/og_${og}.status"),
          path("${params.postdir}/macse/og_${og}.log")

    script:
    """
    mkdir -p ${params.postdir}/macse

    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export VECLIB_MAXIMUM_THREADS=${task.cpus}
    export NUMEXPR_NUM_THREADS=${task.cpus}
    export MALLOC_ARENA_MAX=2

    aa_out="${params.postdir}/macse/og_${og}_AA.fasta"
    nt_out="${params.postdir}/macse/og_${og}_NT.fasta"
    status_out="${params.postdir}/macse/og_${og}.status"
    log_out="${params.postdir}/macse/og_${og}.log"

    rm -f "\$status_out" "\$aa_out" "\$nt_out" "\$log_out"

    echo "STARTED" > "\$status_out"

    fail_outputs() {
        rm -f "\$aa_out" "\$nt_out"
        : > "\$aa_out"
        : > "\$nt_out"
        echo "FAIL" > "\$status_out"
    }

    on_exit() {
        rc=\$?

        if [ -s "\$status_out" ] && grep -qx "OK" "\$status_out"; then
            echo 0 > .exitcode
            exit 0
        fi

        if [ ! -s "\$status_out" ] || grep -qx "STARTED" "\$status_out"; then
            fail_outputs
        fi

        echo 0 > .exitcode
        exit 0
    }

    trap on_exit EXIT TERM INT USR1 USR2

    set +e
    /opt/conda/envs/macse/bin/java \
      -Xmx6g \
      -XX:ActiveProcessorCount=${task.cpus} \
      -jar /opt/conda/envs/macse/share/macse-2.07-0/macse_v2.07.jar \
      -prog alignSequences \
      -seq ${fasta} \
      -out_AA "\$aa_out" \
      -out_NT "\$nt_out" \
      > "\$log_out" 2>&1
    rc=\$?
    set -e

    if [ \$rc -eq 0 ] && [ -s "\$aa_out" ] && [ -s "\$nt_out" ]; then
        echo "OK" > "\$status_out"
    else
        fail_outputs
    fi

    echo 0 > .exitcode
    trap - EXIT TERM INT USR1 USR2
    exit 0
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

    echo -e "og\tstatus" > ${params.postdir}/macse/macse_report.tsv

    find . -name 'og_*.status' -type f | sort | while read -r f; do
      og=\$(basename "\$f" .status | sed 's/^og_//')
      st=\$(tr -d '\r\n' < "\$f")
      [ -n "\$st" ] || st="FAIL"
      echo -e "\${og}\t\${st}" >> ${params.postdir}/macse/macse_report.tsv
    done

    awk -F'\t' 'NR > 1 && \$2 == "OK" { print \$1 }' \
      ${params.postdir}/macse/macse_report.tsv \
      > ${params.postdir}/macse/macse_ok_og_list.txt

    test -s ${params.postdir}/macse/macse_report.tsv
    test -f ${params.postdir}/macse/macse_ok_og_list.txt
    """
}

process MAFFT_ALIGN_AA {
    tag { "og_${og}" }
    array (params.array_size as int)

    input:
    tuple val(og), path(fasta)

    output:
    tuple val(og),
          path("${params.postdir}/mafft_aa/og_${og}_AA.fasta"),
          path("${params.postdir}/mafft_aa/og_${og}.status"),
          path("${params.postdir}/mafft_aa/og_${og}.log")

    script:
    """
    mkdir -p ${params.postdir}/mafft_aa
    mkdir -p tmp

    export TMPDIR="\$PWD/tmp"
    export TMP="\$TMPDIR"
    export TEMP="\$TMPDIR"
    export MAFFT_TMPDIR="\$TMPDIR"

    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export VECLIB_MAXIMUM_THREADS=${task.cpus}
    export NUMEXPR_NUM_THREADS=${task.cpus}
    export MALLOC_ARENA_MAX=2

    aa_out="${params.postdir}/mafft_aa/og_${og}_AA.fasta"
    status_out="${params.postdir}/mafft_aa/og_${og}.status"
    log_out="${params.postdir}/mafft_aa/og_${og}.log"

    rm -f "\$status_out" "\$aa_out" "\$log_out"

    echo "STARTED" > "\$status_out"

    fail_outputs() {
        rm -f "\$aa_out"
        : > "\$aa_out"
        echo "FAIL" > "\$status_out"
    }

    on_exit() {
        rc=\$?

        if [ -s "\$status_out" ] && grep -qx "OK" "\$status_out"; then
            echo 0 > .exitcode
            exit 0
        fi

        if [ ! -s "\$status_out" ] || grep -qx "STARTED" "\$status_out"; then
            fail_outputs
        fi

        echo 0 > .exitcode
        exit 0
    }

    trap on_exit EXIT TERM INT USR1 USR2

    set +e
    ${params.mafft_bin} \
      ${params.mafft_opts} \
      --thread ${task.cpus} \
      ${fasta} \
      > "\$aa_out" \
      2> "\$log_out"
    rc=\$?
    set -e

    if [ \$rc -eq 0 ] && [ -s "\$aa_out" ]; then
        echo "OK" > "\$status_out"
    else
        fail_outputs
    fi

    echo 0 > .exitcode
    trap - EXIT TERM INT USR1 USR2
    exit 0
    """
}

process MAFFT_REPORT {
    tag "mafft_report"

    input:
    path mafft_results

    output:
    path("${params.postdir}/mafft_aa/mafft_report.tsv")
    path("${params.postdir}/mafft_aa/mafft_ok_og_list.txt")

    script:
    """
    mkdir -p ${params.postdir}/mafft_aa

    echo -e "og\tstatus" > ${params.postdir}/mafft_aa/mafft_report.tsv

    find . -name 'og_*.status' -type f | sort | while read -r f; do
      og=\$(basename "\$f" .status | sed 's/^og_//')
      st=\$(tr -d '\r\n' < "\$f")
      [ -n "\$st" ] || st="FAIL"
      echo -e "\${og}\t\${st}" >> ${params.postdir}/mafft_aa/mafft_report.tsv
    done

    awk -F'\t' 'NR > 1 && \$2 == "OK" { print \$1 }' \
      ${params.postdir}/mafft_aa/mafft_report.tsv \
      > ${params.postdir}/mafft_aa/mafft_ok_og_list.txt

    test -s ${params.postdir}/mafft_aa/mafft_report.tsv
    test -f ${params.postdir}/mafft_aa/mafft_ok_og_list.txt
    """
}

process IQTREE_NT_OG {
    tag { "og_${og}" }
    array (params.array_size as int)

    input:
    tuple val(og), path(aa), path(nt), path(status), path(macse_log)

    output:
    tuple val(og), path("${params.postdir}/iqtree_nt/og_${og}")

    script:
    """
    outdir="${params.postdir}/iqtree_nt/og_${og}"
    mkdir -p "\$outdir"

    status_out="\$outdir/og_${og}.iqtree.status"
    tree_out="\$outdir/og_${og}_iqtree.treefile"
    ufboot_out="\$outdir/og_${og}_iqtree.ufboot"
    log_out="\$outdir/og_${og}.log"

    rm -f "\$status_out" "\$tree_out" "\$ufboot_out" "\$log_out"

    echo "STARTED" > "\$status_out"

    fail_outputs() {
        : > "\$tree_out"
        : > "\$ufboot_out"
        echo "FAIL" > "\$status_out"
    }

    on_exit() {
        rc=\$?

        if [ -s "\$tree_out" ] && [ -s "\$ufboot_out" ]; then
            echo "OK" > "\$status_out"
            echo 0 > .exitcode
            exit 0
        fi

        if [ ! -s "\$status_out" ] || grep -qx "STARTED" "\$status_out"; then
            fail_outputs
        fi

        echo 0 > .exitcode
        exit 0
    }

    trap on_exit EXIT TERM INT USR1 USR2

    if [ ! -s ${nt} ] || [ "\$(tr -d '\r\n' < ${status})" != "OK" ]; then
        fail_outputs
        echo "Skipped: missing NT alignment or upstream MACSE status was not OK" > "\$log_out"
        echo 0 > .exitcode
        trap - EXIT TERM INT USR1 USR2
        exit 0
    fi

    gap_fraction=\$(awk '
      /^>/ {next}
      {
        line = \$0
        total += length(line)
        gaps += gsub(/[-XxNn?]/, "", line)
      }
      END {
        if (total == 0) print 1
        else print gaps / total
      }
    ' ${nt})

    if awk -v gf="\$gap_fraction" 'BEGIN {exit !(gf > 0.50)}'; then
        fail_outputs
        echo "Skipped: gap_fraction=\$gap_fraction" > "\$log_out"
        echo 0 > .exitcode
        trap - EXIT TERM INT USR1 USR2
        exit 0
    fi

    set +e
    ${params.iqtree_bin} \
      -s ${nt} \
      -st DNA \
      -T ${task.cpus} \
      -m ${params.iqtree_nt_model} \
      -bb 1000 \
      -wbtl \
      -redo \
      -pre "\$outdir/og_${og}_iqtree" \
      > "\$log_out" 2>&1
    rc=\$?
    set -e

    if [ \$rc -eq 0 ] && [ -s "\$tree_out" ] && [ -s "\$ufboot_out" ]; then
        echo "OK" > "\$status_out"
    else
        fail_outputs
    fi

    echo 0 > .exitcode
    trap - EXIT TERM INT USR1 USR2
    exit 0
    """
}

process IQTREE_NT_REPORT {
    tag "iqtree_nt_report"

    input:
    path iqtree_dirs

    output:
    path("${params.postdir}/iqtree_nt/iqtree_nt_report.tsv")
    path("${params.postdir}/iqtree_nt/iqtree_nt_ok_og_list.txt")

    script:
    """
    mkdir -p ${params.postdir}/iqtree_nt

    echo -e "og\tstatus" > ${params.postdir}/iqtree_nt/iqtree_nt_report.tsv

    find . -name 'og_*.iqtree.status' -type f | sort | while read -r f; do
      og=\$(basename "\$f" .iqtree.status | sed 's/^og_//')
      st=\$(tr -d '\r\n' < "\$f")
      [ -n "\$st" ] || st="FAIL"
      echo -e "\${og}\t\${st}" >> ${params.postdir}/iqtree_nt/iqtree_nt_report.tsv
    done

    awk -F'\t' 'NR > 1 && \$2 == "OK" { print \$1 }' \
      ${params.postdir}/iqtree_nt/iqtree_nt_report.tsv \
      > ${params.postdir}/iqtree_nt/iqtree_nt_ok_og_list.txt

    test -s ${params.postdir}/iqtree_nt/iqtree_nt_report.tsv
    test -f ${params.postdir}/iqtree_nt/iqtree_nt_ok_og_list.txt
    """
}

process IQTREE_AA_OG {
    tag { "og_${og}" }
    array (params.array_size as int)

    input:
    tuple val(og), path(aa), path(status), path(mafft_log)

    output:
    tuple val(og), path("${params.postdir}/iqtree_aa/og_${og}")

    script:
    """
    outdir="${params.postdir}/iqtree_aa/og_${og}"
    mkdir -p "\$outdir"

    status_out="\$outdir/og_${og}.iqtree.status"
    tree_out="\$outdir/og_${og}_iqtree.treefile"
    ufboot_out="\$outdir/og_${og}_iqtree.ufboot"
    log_out="\$outdir/og_${og}.log"

    rm -f "\$status_out" "\$tree_out" "\$ufboot_out" "\$log_out"

    echo "STARTED" > "\$status_out"

    fail_outputs() {
        : > "\$tree_out"
        : > "\$ufboot_out"
        echo "FAIL" > "\$status_out"
    }

    on_exit() {
        rc=\$?

        if [ -s "\$tree_out" ] && [ -s "\$ufboot_out" ]; then
            echo "OK" > "\$status_out"
            echo 0 > .exitcode
            exit 0
        fi

        if [ ! -s "\$status_out" ] || grep -qx "STARTED" "\$status_out"; then
            fail_outputs
        fi

        echo 0 > .exitcode
        exit 0
    }

    trap on_exit EXIT TERM INT USR1 USR2

    if [ ! -s ${aa} ] || [ "\$(tr -d '\r\n' < ${status})" != "OK" ]; then
        fail_outputs
        echo "Skipped: missing AA alignment or upstream MAFFT status was not OK" > "\$log_out"
        echo 0 > .exitcode
        trap - EXIT TERM INT USR1 USR2
        exit 0
    fi

    gap_fraction=\$(awk '
      /^>/ {next}
      {
        line = \$0
        total += length(line)
        gaps += gsub(/[-Xx?*]/, "", line)
      }
      END {
        if (total == 0) print 1
        else print gaps / total
      }
    ' ${aa})

    if awk -v gf="\$gap_fraction" 'BEGIN {exit !(gf > 0.50)}'; then
        fail_outputs
        echo "Skipped: gap_fraction=\$gap_fraction" > "\$log_out"
        echo 0 > .exitcode
        trap - EXIT TERM INT USR1 USR2
        exit 0
    fi

    set +e
    ${params.iqtree_bin} \
      -s ${aa} \
      -st AA \
      -T ${task.cpus} \
      -m ${params.iqtree_aa_model} \
      -bb 1000 \
      -wbtl \
      -redo \
      -pre "\$outdir/og_${og}_iqtree" \
      > "\$log_out" 2>&1
    rc=\$?
    set -e

    if [ \$rc -eq 0 ] && [ -s "\$tree_out" ] && [ -s "\$ufboot_out" ]; then
        echo "OK" > "\$status_out"
    else
        fail_outputs
    fi

    echo 0 > .exitcode
    trap - EXIT TERM INT USR1 USR2
    exit 0
    """
}

process IQTREE_AA_REPORT {
    tag "iqtree_aa_report"

    input:
    path iqtree_dirs

    output:
    path("${params.postdir}/iqtree_aa/iqtree_aa_report.tsv")
    path("${params.postdir}/iqtree_aa/iqtree_aa_ok_og_list.txt")

    script:
    """
    mkdir -p ${params.postdir}/iqtree_aa

    echo -e "og\tstatus" > ${params.postdir}/iqtree_aa/iqtree_aa_report.tsv

    find . -name 'og_*.iqtree.status' -type f | sort | while read -r f; do
      og=\$(basename "\$f" .iqtree.status | sed 's/^og_//')
      st=\$(tr -d '\r\n' < "\$f")
      [ -n "\$st" ] || st="FAIL"
      echo -e "\${og}\t\${st}" >> ${params.postdir}/iqtree_aa/iqtree_aa_report.tsv
    done

    awk -F'\t' 'NR > 1 && \$2 == "OK" { print \$1 }' \
      ${params.postdir}/iqtree_aa/iqtree_aa_report.tsv \
      > ${params.postdir}/iqtree_aa/iqtree_aa_ok_og_list.txt

    test -s ${params.postdir}/iqtree_aa/iqtree_aa_report.tsv
    test -f ${params.postdir}/iqtree_aa/iqtree_aa_ok_og_list.txt
    """
}
