nextflow.enable.dsl=2

include { PRIMARY_TRANSCRIPT; FINALIZE_REPO_IDS } from './modules/local/prep'
include { PARSE_ANNOTATIONS; VALIDATE_PARSE_OUTPUTS; ORTHOFINDER_OR_SKIP; RUN_GENESPACE } from './modules/local/genespace'
include { PANGENES_PASS_FILTER; WRITE_OG_FASTAS; MACSE_ALIGN_OG; MACSE_REPORT; IQTREE_OG; IQTREE_REPORT; WRITE_ALERAX_MAPPING; WRITE_ALERAX_FAMILIES; RUN_ALERAX } from './modules/local/post_genespace'

def resolveChrDict(genome) {
    def tsv = file("${params.chr_dict_dir}/${genome}.tsv")
    def bed = file("${params.chr_dict_dir}/${genome}_chr_lengths.bed")
    if (tsv.exists()) return tsv
    if (bed.exists()) return bed
    throw new IllegalArgumentException("No chr_dict found for ${genome}: expected ${tsv} or ${bed}")
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
        throw new IllegalArgumentException("genomes TSV must contain columns: genome_id, genome_source, ploidy")
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

    genomes_ch = Channel
        .fromList(genomes_rows)
        .map { row ->
            tuple(
                row.genome,
                row.source,
                row.ploidy,
                file("${params.gff_dir}/${row.genome}.${params.ext.gff}"),
                file("${params.protein_dir}/${row.genome}.${params.ext.pep}"),
                resolveChrDict(row.genome)
            )
        }

    genomes_tsv_ch  = Channel.value(file(params.genomes_tsv))
    species_tree_ch = Channel.value(file(params.alerax.species_tree))
    cds_files_ch    = Channel.fromPath("${params.cds_dir}/*.cds").collect()

    primary_out   = PRIMARY_TRANSCRIPT(genomes_ch)
    finalized_out = FINALIZE_REPO_IDS(primary_out)

    parse_out       = PARSE_ANNOTATIONS(finalized_out.collect(), genomes_tsv_ch)
    validate_out    = VALIDATE_PARSE_OUTPUTS(parse_out)
    orthofinder_out = ORTHOFINDER_OR_SKIP(validate_out, genomes_tsv_ch)
    genespace_out   = RUN_GENESPACE(validate_out, orthofinder_out, genomes_tsv_ch)

    pass_out      = PANGENES_PASS_FILTER(genespace_out)
    og_fastas_out = WRITE_OG_FASTAS(pass_out[0], pass_out[1], genomes_tsv_ch, cds_files_ch)

    og_fasta_ch = Channel
        .fromPath("${params.postdir}/og_fasta/og_*.fasta")
        .map { fasta ->
            def m = (fasta.baseName =~ /^og_(.+)$/)
            if (!m) {
                throw new IllegalArgumentException("Could not parse OG from filename: ${fasta}")
            }
            tuple(m[0][1], fasta)
        }

    macse_out        = MACSE_ALIGN_OG(og_fasta_ch)
    macse_report_out = MACSE_REPORT(macse_out.collect())

    iqtree_in = macse_out.filter { og, aa, nt, status ->
        status.text.trim() == 'OK'
    }

    iqtree_out        = IQTREE_OG(iqtree_in)
    iqtree_report_out = IQTREE_REPORT(iqtree_out.collect())

    alerax_map_out = WRITE_ALERAX_MAPPING(iqtree_out)
    families_out   = WRITE_ALERAX_FAMILIES(alerax_map_out.collect())
    RUN_ALERAX(families_out, species_tree_ch)
}
