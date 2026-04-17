process PRIMARY_TRANSCRIPT {
    tag { genome }

    input:
    tuple val(genome), val(source), val(ploidy), path(gff), path(pep), path(chr)
    path primary_transcript_script

    output:
    tuple val(genome), val(source), val(ploidy), path("${genome}.primary.pep"), path(gff), path(chr)

    script:
    def cmd =
        source == 'helixer'
            ? """
              cp ${pep} ${genome}.primary.pep
              """
        : source == 'phytozome'
            ? """
              python ${primary_transcript_script} ${pep} \
                --mode phytozome \
                --phytozome-gff ${gff} \
                > ${genome}.log 2>&1

              test -s primary_transcripts/${pep.getName()}
              cp primary_transcripts/${pep.getName()} ${genome}.primary.pep
              """
        : """
          python ${primary_transcript_script} ${pep} \
            > ${genome}.log 2>&1

          test -s primary_transcripts/${pep.getName()}
          cp primary_transcripts/${pep.getName()} ${genome}.primary.pep
          """

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
    """
    cp ${chr} ${genome}.chr.tsv

    python ${apply_chr_dict_script} \
      --in-gff ${gff} \
      --out-gff ${genome}.chr.gff3 \
      --chr-dict ${genome}.chr.tsv \
      > ${genome}.chrify.log 2>&1

    python ${finalize_repo_ids_script} \
      --source ${source} \
      --in-gff ${genome}.chr.gff3 \
      --out-gff ${genome}.final.gff3 \
      --in-pep ${primary_pep} \
      --out-pep ${genome}.final.pep \
      > ${genome}.finalize.log 2>&1

    test -s ${genome}.final.gff3
    test -s ${genome}.final.pep
    """
}
