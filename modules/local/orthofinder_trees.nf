/*
 * Raw OrthoFinder trees: all OG members, no synteny filtering or tandem collapse.
 * Expected per-OG application failures become OK/FAIL statuses rather than
 * failing the entire Nextflow run. Infrastructure/staging errors still fail.
 */

process ORTHOFINDER_RAW_FASTAS {
    tag 'orthofinder_raw_fastas'
    input:
    path orthofinder_dir
    path genomes_tsv
    path cds_files
    path protein_files
    path writer
    output:
    path 'orthofinder_raw/og_fasta_aa'
    path 'orthofinder_raw/og_fasta_cds'
    path 'orthofinder_raw/og_report.tsv'
    script:
    def cdsList = cds_files.collect { "\"${it}\"" }.join(' ')
    def proteinList = protein_files.collect { "\"${it}\"" }.join(' ')
    """
    mkdir -p input_cds input_protein
    for f in ${cdsList}; do ln -s "\$(realpath "\$f")" input_cds/; done
    for f in ${proteinList}; do ln -s "\$(realpath "\$f")" input_protein/; done
    # Writer must recognise Species.final.pep as well as Species.pep.
    python ${writer} \
      --orthofinder-dir ${orthofinder_dir} \
      --genomes-tsv ${genomes_tsv} \
      --cds-dir input_cds --protein-dir input_protein \
      --cds-ext '${params.ext.cds}' --protein-ext '${params.ext.pep}' \
      --min-species ${params.orthofinder_trees_min_species ?: 4} \
      --outdir orthofinder_raw
    """
}

process ORTHOFINDER_MACSE_NT {
    tag { og }
    input:
    tuple val(og), path(fasta)
    output:
    tuple val(og), path('alignment_NT.fasta'), path('alignment.status')
    // Log and AA output are kept in this task's work directory, if generated.
    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export VECLIB_MAXIMUM_THREADS=${task.cpus}
    export NUMEXPR_NUM_THREADS=${task.cpus}
    export MALLOC_ARENA_MAX=2

    echo STARTED > alignment.status
    set +e
    /opt/conda/envs/macse/bin/java \
      -XX:-UsePerfData -XX:ActiveProcessorCount=${task.cpus} -Xmx6g \
      -jar /opt/conda/envs/macse/share/macse-2.07-0/macse_v2.07.jar \
      -prog alignSequences -seq ${fasta} \
      -out_AA alignment_AA.fasta -out_NT alignment_NT.fasta \
      > macse.log 2>&1
    rc=\$?
    set -e

    if [ "\$rc" -eq 0 ] && [ -s alignment_NT.fasta ] && [ -s alignment_AA.fasta ]; then
      echo OK > alignment.status
    else
      printf 'MACSE failed or produced incomplete output (exit status %s)\\n' "\$rc" >> macse.log
      rm -f alignment_NT.fasta alignment_AA.fasta
      : > alignment_NT.fasta
      echo FAIL > alignment.status
    fi
    """
}

process ORTHOFINDER_MAFFT_AA {
    tag { og }
    input:
    tuple val(og), path(fasta)
    output:
    tuple val(og), path('alignment_AA.fasta'), path('alignment.status')
    script:
    """
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

    echo STARTED > alignment.status
    set +e
    ${params.mafft_bin} ${params.mafft_opts ?: ''} --thread ${task.cpus} \
      ${fasta} > alignment_AA.fasta 2> mafft.log
    rc=\$?
    set -e

    if [ "\$rc" -eq 0 ] && [ -s alignment_AA.fasta ]; then
      echo OK > alignment.status
    else
      printf 'MAFFT failed or produced no alignment (exit status %s)\\n' "\$rc" >> mafft.log
      : > alignment_AA.fasta
      echo FAIL > alignment.status
    fi
    """
}

process ORTHOFINDER_IQTREE_NT {
    tag { og }
    input:
    tuple val(og), path(aln), path(status)
    output:
    tuple val(og), path("orthofinder_nt/${og}")
    script:
    """
    outdir="orthofinder_nt/${og}"
    mkdir -p "\$outdir"
    tree="\$outdir/${og}.treefile"
    ufboot="\$outdir/${og}.ufboot"
    logfile="\$outdir/iqtree.log"
    tree_status="\$outdir/${og}.iqtree.status"
    cp ${aln} "\$outdir/alignment_NT.fasta"
    echo STARTED > "\$tree_status"

    fail_tree() {
      : > "\$tree"
      : > "\$ufboot"
      echo FAIL > "\$tree_status"
    }

    if [ ! -s "\$outdir/alignment_NT.fasta" ] || [ "\$(tr -d '\\r\\n' < ${status})" != OK ]; then
      fail_tree
      echo 'Skipped: missing NT alignment or upstream MACSE status not OK' > "\$logfile"
      exit 0
    fi

    gap_fraction=\$(awk '
      /^>/ {next}
      { total += length(\$0); line=\$0; gaps += gsub(/[-XxNn?]/, "", line) }
      END {if (total == 0) print 1; else print gaps / total}
    ' "\$outdir/alignment_NT.fasta")
    if awk -v gf="\$gap_fraction" 'BEGIN {exit !(gf > 0.50)}'; then
      fail_tree
      echo "Skipped: gap_fraction=\$gap_fraction (>0.50)" > "\$logfile"
      exit 0
    fi

    set +e
    ${params.iqtree_bin} -s "\$outdir/alignment_NT.fasta" -st DNA \
       -m '${params.iqtree_nt_model ?: 'GTR+F+G4'}' -T ${task.cpus} \
       -bb 1000 -wbtl -pre "\$outdir/${og}" \
       > "\$logfile" 2>&1
    rc=\$?
    set -e
    if [ "\$rc" -eq 0 ] && [ -s "\$tree" ] && [ -s "\$ufboot" ]; then
      echo OK > "\$tree_status"
    else
      printf 'IQ-TREE failed or missing outputs (exit status %s)\\n' "\$rc" >> "\$logfile"
      fail_tree
    fi
    """
}

process ORTHOFINDER_IQTREE_AA {
    tag { og }
    input:
    tuple val(og), path(aln), path(status)
    output:
    tuple val(og), path("orthofinder_aa/${og}")
    script:
    """
    outdir="orthofinder_aa/${og}"
    mkdir -p "\$outdir"
    tree="\$outdir/${og}.treefile"
    ufboot="\$outdir/${og}.ufboot"
    logfile="\$outdir/iqtree.log"
    tree_status="\$outdir/${og}.iqtree.status"
    cp ${aln} "\$outdir/alignment_AA.fasta"
    echo STARTED > "\$tree_status"

    fail_tree() {
      : > "\$tree"
      : > "\$ufboot"
      echo FAIL > "\$tree_status"
    }

    if [ ! -s "\$outdir/alignment_AA.fasta" ] || [ "\$(tr -d '\\r\\n' < ${status})" != OK ]; then
      fail_tree
      echo 'Skipped: missing AA alignment or upstream MAFFT status not OK' > "\$logfile"
      exit 0
    fi

    gap_fraction=\$(awk '
      /^>/ {next}
      { total += length(\$0); line=\$0; gaps += gsub(/[-Xx?*]/, "", line) }
      END {if (total == 0) print 1; else print gaps / total}
    ' "\$outdir/alignment_AA.fasta")
    if awk -v gf="\$gap_fraction" 'BEGIN {exit !(gf > 0.50)}'; then
      fail_tree
      echo "Skipped: gap_fraction=\$gap_fraction (>0.50)" > "\$logfile"
      exit 0
    fi

    set +e
    ${params.iqtree_bin} -s "\$outdir/alignment_AA.fasta" -st AA \
       -m '${params.iqtree_aa_model ?: 'LG+G4'}' -T ${task.cpus} \
       -bb 1000 -wbtl -pre "\$outdir/${og}" \
       > "\$logfile" 2>&1
    rc=\$?
    set -e
    if [ "\$rc" -eq 0 ] && [ -s "\$tree" ] && [ -s "\$ufboot" ]; then
      echo OK > "\$tree_status"
    else
      printf 'IQ-TREE failed or missing outputs (exit status %s)\\n' "\$rc" >> "\$logfile"
      fail_tree
    fi
    """
}
