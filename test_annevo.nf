nextflow.enable.dsl = 2

include { RUN_ANNEVO; ANNEVO_GFF_TO_FASTA } from './modules/local/annotation'


def resolveChrDict(genome) {
    def tsv = file("${params.chr_dict_dir}/${genome}.tsv")
    def bed = file("${params.chr_dict_dir}/${genome}_chr_lengths.bed")

    if (tsv.exists()) return tsv
    if (bed.exists()) return bed

    throw new IllegalArgumentException(
        "No chr_dict found for ${genome}: expected ${tsv} or ${bed}"
    )
}


def resolveGenomeFasta(genome) {
    def exts = [
        params.ext.fasta,
        "${params.ext.fasta}.gz",
        "fasta",
        "fasta.gz",
        "fna",
        "fna.gz",
        "fa",
        "fa.gz"
    ]

    def found_ext = exts.find { ext ->
        file("${params.fasta_dir}/${genome}.${ext}").exists()
    }

    if (found_ext) {
        return file("${params.fasta_dir}/${genome}.${found_ext}")
    }

    throw new IllegalArgumentException(
        "No genome FASTA found for ${genome}: expected one of " +
        exts.collect { "${params.fasta_dir}/${genome}.${it}" }.join(", ")
    )
}


def readGenomesTable(tsvPath) {
    def lines = file(tsvPath).readLines().findAll { it.trim() }

    if (!lines) {
        throw new IllegalArgumentException("Empty genomes TSV: ${tsvPath}")
    }

    def header = lines[0].split('\t', -1) as List
    def idxGenome = header.indexOf('genome_id')
    def idxSource = header.indexOf('genome_source')
    def idxPloidy = header.indexOf('ploidy')

    if (idxGenome < 0 || idxSource < 0 || idxPloidy < 0) {
        throw new IllegalArgumentException(
            "genomes TSV must contain columns: genome_id, genome_source, ploidy"
        )
    }

    def rows = []

    lines.drop(1).each { line ->
        def cols = line.split('\t', -1) as List

        rows << [
            genome : cols[idxGenome].trim(),
            source : cols[idxSource].trim().toLowerCase(),
            ploidy : cols[idxPloidy].trim() as Integer
        ]
    }

    rows
}


workflow {

    genomes_rows = readGenomesTable(params.genomes_tsv)

    if (params.test_genome?.toString()?.trim()) {
        genomes_rows = genomes_rows.findAll { row ->
            row.genome == params.test_genome
        }

        if (genomes_rows.isEmpty()) {
            throw new IllegalArgumentException(
                "No row in ${params.genomes_tsv} matched --test_genome ${params.test_genome}"
            )
        }
    }

    annevo_input_ch = Channel
        .fromList(genomes_rows)
        .map { row ->
            tuple(
                row.genome,
                row.source,
                row.ploidy,
                resolveGenomeFasta(row.genome),
                resolveChrDict(row.genome)
            )
        }

    annevo_gff_out = RUN_ANNEVO(annevo_input_ch)

    annevo_fasta_out = ANNEVO_GFF_TO_FASTA(annevo_gff_out)

    annevo_fasta_out.view { genome, source, ploidy, gff, pep, chr, cds ->
        "ANNOTATION_TEST_OK\t${genome}\t${source}\t${ploidy}\t${gff}\t${pep}\t${cds}\t${chr}"
    }
}
