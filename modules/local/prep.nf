process PRIMARY_TRANSCRIPT {
    tag { genome }

    input:
    tuple val(genome), val(source), val(ploidy), path(gff), path(pep), path(chr)

    output:
    tuple val(genome), val(source), val(ploidy), path("${genome}.primary.pep"), path(gff), path(chr)

    script:
    def cmd = source == 'phytozome'
        ? """
          python ${projectDir}/scripts/primary_transcript.py ${pep} \
            --mode phytozome \
            --phytozome-gff ${gff} \
            > ${genome}.primary.pep 2> ${genome}.log
          """
        : """
          python ${projectDir}/scripts/primary_transcript.py ${pep} \
            > ${genome}.primary.pep 2> ${genome}.log
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

    output:
    tuple val(genome), val(source), val(ploidy),
          path("${genome}.final.gff3"),
          path("${genome}.final.pep"),
          path("${genome}.chr.tsv")

    script:
    """
    cp ${chr} ${genome}.chr.tsv
    cp ${gff} ${genome}.gff3
    cp ${primary_pep} ${genome}.pep

    python ${projectDir}/scripts/apply_chr_dict_to_gff.py \
      --in-gff ${genome}.gff3 \
      --out-gff ${genome}.chr.gff3 \
      --chr-dict ${genome}.chr.tsv \
      > ${genome}.chrify.log 2>&1

    python ${projectDir}/scripts/finalize_repo_ids.py \
      --source ${source} \
      --in-gff ${genome}.chr.gff3 \
      --out-gff ${genome}.final.gff3 \
      --in-pep ${genome}.pep \
      --out-pep ${genome}.final.pep \
      > ${genome}.finalize.log 2>&1

    test -s ${genome}.final.gff3
    test -s ${genome}.final.pep
    """
}
