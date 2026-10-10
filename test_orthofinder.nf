nextflow.enable.dsl=2

include { ORTHOFINDER_RAW_FASTAS; ORTHOFINDER_MACSE_NT; ORTHOFINDER_MAFFT_AA; ORTHOFINDER_IQTREE_NT; ORTHOFINDER_IQTREE_AA } from './modules/local/orthofinder_trees'

def rawFastaChannel(dirChannel, limit) {
    dirChannel.flatMap { dir ->
        dir.listFiles()
            .findAll { it.name.endsWith('.fasta') }
            .sort { a, b -> a.name <=> b.name }
            .take(limit)
            .collect { fasta -> tuple(fasta.baseName, fasta) }
    }
}

workflow {
    def n = (params.test_og_count ?: 5) as int
    if (n < 1) error 'test_og_count must be positive'
    if (!params.existing_orthofinder_dir) error 'Set existing_orthofinder_dir'
    if (!params.genomes_tsv) error 'Set genomes_tsv'
    if (!params.cds_dir || !params.protein_dir) error 'Set cds_dir and protein_dir'

    def fastaOut = ORTHOFINDER_RAW_FASTAS(
        Channel.value(file(params.existing_orthofinder_dir, checkIfExists: true)),
        Channel.value(file(params.genomes_tsv, checkIfExists: true)),
        Channel.fromPath("${params.cds_dir}/*.${params.ext.cds}", checkIfExists: true).collect(),
        Channel.fromPath("${params.protein_dir}/*.${params.ext.pep}", checkIfExists: true).collect(),
        Channel.value(file('scripts/write_orthofinder_og_fastas.py', checkIfExists: true))
    )

    def cds = rawFastaChannel(fastaOut[1], n)
    def aa = rawFastaChannel(fastaOut[0], n)

    def macse = ORTHOFINDER_MACSE_NT(cds)
    def mafft = ORTHOFINDER_MAFFT_AA(aa)
    ORTHOFINDER_IQTREE_NT(macse)
    ORTHOFINDER_IQTREE_AA(mafft)
}
