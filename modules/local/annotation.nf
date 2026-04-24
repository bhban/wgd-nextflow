process RUN_ANNEVO {
    tag { genome }

    input:
    tuple val(genome), val(source), val(ploidy), path(genome_fasta), path(chr)

    output:
    tuple val(genome), val('annevo'), val(ploidy),
          path("${genome}.annevo.gff3"),
          path(genome_fasta),
          path(chr)

    script:
    def extra = params.annevo.extra_args ?: ""
    def batch = params.annevo.batch_size ? "--batch_size ${params.annevo.batch_size}" : ""

    """
    python annotation.py \\
      --genome ${genome_fasta} \\
      --model_path ${params.annevo.model_path} \\
      --output ${genome}.annevo.gff3 \\
      --threads ${task.cpus} \\
      ${batch} \\
      ${extra}

    test -s ${genome}.annevo.gff3
    """
}

process ANNEVO_GFF_TO_FASTA {
    tag { genome }

    input:
    tuple val(genome), val(source), val(ploidy), path(gff), path(genome_fasta), path(chr)

    output:
    tuple val(genome), val(source), val(ploidy),
          path("${genome}.annevo.gff3"),
          path("${genome}.annevo.pep"),
          path(chr),
          path("${genome}.annevo.cds")

    script:
    """
    cp ${gff} ${genome}.annevo.gff3

    gffread ${genome}.annevo.gff3 \\
      -g ${genome_fasta} \\
      -x ${genome}.annevo.cds \\
      -y ${genome}.annevo.pep

    test -s ${genome}.annevo.gff3
    test -s ${genome}.annevo.pep
    test -s ${genome}.annevo.cds
    """
}
