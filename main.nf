nextflow.enable.dsl=2

include { PRIMARY_TRANSCRIPT; FINALIZE_REPO_IDS } from './modules/local/prep'
include { STAGE_GENOMEREPO; PARSE_ANNOTATIONS_BY_SOURCE; VALIDATE_PARSE_OUTPUTS; ORTHOFINDER_OR_SKIP; RUN_GENESPACE } from './modules/local/genespace'
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

    primary_transcript_script_ch     = Channel.value(file('scripts/primary_transcript.py'))
    apply_chr_dict_script_ch         = Channel.value(file('scripts/apply_chr_dict_to_gff.py'))
    finalize_repo_ids_script_ch      = Channel.value(file('scripts/finalize_repo_ids.py'))
    parse_annotations_script_ch      = Channel.value(file('scripts/run_parse_annotations_by_source.R'))
    validate_parse_outputs_script_ch = Channel.value(file('scripts/validate_parse_outputs.py'))
    orthofinder_or_skip_script_ch    = Channel.value(file('scripts/orthofinder_or_skip.py'))
    run_genespace_script_ch          = Channel.value(file('scripts/run_genespace.R'))
    pangenes_pass_filter_script_ch   = Channel.value(file('scripts/pangenes_pass_filter.py'))
    write_og_fastas_script_ch        = Channel.value(file('scripts/write_og_fastas.py'))

    genome_ids_ch = Channel.value(genomes_rows.collect { it.genome })

    primary_out = PRIMARY_TRANSCRIPT(
        genomes_ch,
        primary_transcript_script_ch
    )

    finalized_out = FINALIZE_REPO_IDS(
        primary_out,
        apply_chr_dict_script_ch,
        finalize_repo_ids_script_ch
    )

    staged_repo_out = STAGE_GENOMEREPO(
        finalized_out
            .flatMap { genome, source, ploidy, gff, pep, chr ->
                [gff, pep, chr]
            }
            .collect()
    )

    parsed_out = PARSE_ANNOTATIONS_BY_SOURCE(
        staged_repo_out,
        genomes_tsv_ch,
        parse_annotations_script_ch
    )

    validated_out = VALIDATE_PARSE_OUTPUTS(
        parsed_out,
        validate_parse_outputs_script_ch,
        genome_ids_ch
    )

    def orthofinder_dir_arg
    def orthofinder_out = null

    if (params.run_external_orthofinder) {
        orthofinder_out = ORTHOFINDER_OR_SKIP(
            validated_out,
            genome_ids_ch,
            orthofinder_or_skip_script_ch
        )
        orthofinder_dir_arg = orthofinder_out[0]
    } else {
        orthofinder_dir_arg = params.orthofinder_dir_override ? params.orthofinder_dir_override : ''
    }

    genespace_out = RUN_GENESPACE(
        validated_out,
        orthofinder_dir_arg,
        genomes_tsv_ch,
        run_genespace_script_ch
    )

    pass_out = PANGENES_PASS_FILTER(
        genespace_out,
        pangenes_pass_filter_script_ch
    )

    og_fastas_out = WRITE_OG_FASTAS(
        pass_out[0],
        pass_out[1],
        genomes_tsv_ch,
        cds_files_ch,
        write_og_fastas_script_ch
    )

    og_fasta_ch = og_fastas_out[0]
        .flatMap { og_dir ->
            og_dir.listFiles()
                .findAll { it.name ==~ /og_.+\.fasta/ }
                .sort { a, b -> a.name <=> b.name }
                .collect { fasta ->
                    def m = (fasta.baseName =~ /^og_(.+)$/)
                    if (!m) {
                        throw new IllegalArgumentException("Could not parse OG from filename: ${fasta}")
                    }
                    tuple(m[0][1], fasta)
                }
        }

    macse_out = MACSE_ALIGN_OG(og_fasta_ch)

    macse_report_out = MACSE_REPORT(
        macse_out
            .map { og, aa, nt, status -> status }
            .collect()
    )
    
    iqtree_in = macse_out.filter { og, aa, nt, status ->
        status.text.trim() == 'OK'
    }
    
    iqtree_out = IQTREE_OG(iqtree_in)
    
    iqtree_report_out = IQTREE_REPORT(
        iqtree_out
            .map { og, treefile, ufboot, status, nt -> status }
            .collect()
    )
    
    alerax_map_out = WRITE_ALERAX_MAPPING(iqtree_out)
    
    families_out = WRITE_ALERAX_FAMILIES(
        alerax_map_out
            .flatMap { og, mapping, ufboot -> [mapping, ufboot] }
            .collect()
    )
    
    RUN_ALERAX(families_out, species_tree_ch)
}
