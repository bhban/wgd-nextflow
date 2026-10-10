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
    for f in ${cdsList}; do ln -s "\$(realpath \"\$f\")" input_cds/; done
    for f in ${proteinList}; do ln -s "\$(realpath \"\$f\")" input_protein/; done
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
    script:
    """
    /opt/conda/envs/macse/bin/java -Xmx6g \
       -jar /opt/conda/envs/macse/share/macse-2.07-0/macse_v2.07.jar \
       -prog alignSequences -seq ${fasta} \
       -out_AA alignment_AA.fasta -out_NT alignment_NT.fasta \
       > macse.log 2>&1
    test -s alignment_NT.fasta
    echo OK > alignment.status
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

    ${params.mafft_bin} \
        ${params.mafft_opts ?: ''} \
        --thread ${task.cpus} \
        ${fasta} \
        > alignment_AA.fasta \
        2> mafft.log

    test -s alignment_AA.fasta

    echo OK > alignment.status
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
    mkdir -p orthofinder_nt/${og}
    cp ${aln} orthofinder_nt/${og}/alignment_NT.fasta
    cd orthofinder_nt/${og}
    ${params.iqtree_bin} -s alignment_NT.fasta -st DNA \
       -m '${params.iqtree_nt_model ?: 'GTR+F+G4'}' -T ${task.cpus} \
       -bb 1000 -wbtl -pre og_tree > iqtree.log 2>&1
    test -s og_tree.treefile
    test -s og_tree.ufboot
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
    mkdir -p orthofinder_aa/${og}
    cp ${aln} orthofinder_aa/${og}/alignment_AA.fasta
    cd orthofinder_aa/${og}
    ${params.iqtree_bin} -s alignment_AA.fasta -st AA \
       -m '${params.iqtree_aa_model ?: 'LG+G4'}' -T ${task.cpus} \
       -bb 1000 -wbtl -pre og_tree > iqtree.log 2>&1
    test -s og_tree.treefile
    test -s og_tree.ufboot
    """
}
