process NORMALISE_TIBERIUS_ANNOTATION {
    tag { genome }

    input:
    tuple val(genome), val(source), val(ploidy), path(gff), path(cds), path(pep), path(chr)
    path normalise_tiberius_script

    output:
    tuple val(genome), val(source), val(ploidy),
          path("${genome}.tiberius.gff3"),
          path("${genome}.tiberius.cds"),
          path("${genome}.tiberius.pep"),
          path(chr)

    script:
    """
    python ${normalise_tiberius_script} \
        --gff ${gff} \
        --cds ${cds} \
        --pep ${pep} \
        --out-gff ${genome}.tiberius.gff3 \
        --out-cds ${genome}.tiberius.cds \
        --out-pep ${genome}.tiberius.pep

    test -s ${genome}.tiberius.gff3
    test -s ${genome}.tiberius.cds
    test -s ${genome}.tiberius.pep
    """
}

process PRIMARY_TRANSCRIPT {
    tag { genome }

    input:
    tuple val(genome), val(source), val(ploidy), path(gff), path(pep), path(chr)
    path primary_transcript_script

    output:
    tuple val(genome), val(source), val(ploidy), path("${genome}.primary.pep"), path(gff), path(chr)

    script:
    def source_lc = source.toString().toLowerCase()
    def skip_primary_filter_sources = ['helixer', 'tiberius']

    def cmd

    if (skip_primary_filter_sources.contains(source_lc)) {
        cmd = """
        cp ${pep} ${genome}.primary.pep
        """
    } else if (source_lc == 'phytozome') {
        cmd = """
        python ${primary_transcript_script} ${pep} \
          --mode phytozome \
          --phytozome-gff ${gff} \
          > ${genome}.log 2>&1

        test -s primary_transcripts/${pep.getName()}
        cp primary_transcripts/${pep.getName()} ${genome}.primary.pep
        """
    } else {
        cmd = """
        python ${primary_transcript_script} ${pep} \
          > ${genome}.log 2>&1

        test -s primary_transcripts/${pep.getName()}
        cp primary_transcripts/${pep.getName()} ${genome}.primary.pep
        """
    }

    """
    ${cmd}
    test -s ${genome}.primary.pep
    """
}

process FINALIZE_REPO_IDS {
    tag { genome }

    input:
    tuple val(genome), val(source), val(ploidy), path(primary_pep), path(gff), path(chr)
    path apply_chr_dict_script
    path finalize_repo_ids_script

    output:
    tuple val(genome), val(source), val(ploidy),
          path("${genome}.final.gff3"),
          path("${genome}.final.pep"),
          path("${genome}.chr.tsv")

    script:
    def source_lc = source.toString().toLowerCase()

    def finalize_cmd

    if (source_lc == 'tiberius') {
        finalize_cmd = """
        cp ${genome}.chr.gff3 ${genome}.final.gff3
        cp ${primary_pep} ${genome}.final.pep
        """
    } else {
        finalize_cmd = """
        python ${finalize_repo_ids_script} \
          --source ${source_lc} \
          --in-gff ${genome}.chr.gff3 \
          --out-gff ${genome}.final.gff3 \
          --in-pep ${primary_pep} \
          --out-pep ${genome}.final.pep \
          > ${genome}.finalize.log 2>&1
        """
    }

    """
    cp ${chr} ${genome}.chr.tsv

    python ${apply_chr_dict_script} \
      --in-gff ${gff} \
      --out-gff ${genome}.chr.gff3 \
      --chr-dict ${genome}.chr.tsv \
      > ${genome}.chrify.log 2>&1

    ${finalize_cmd}

    test -s ${genome}.final.gff3
    test -s ${genome}.final.pep
    """
}
