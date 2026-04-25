nextflow.enable.dsl=2

include { RUN_ANNEVO; ANNEVO_GFF_TO_FASTA } from './modules/local/annotation'
include { PRIMARY_TRANSCRIPT; FINALIZE_REPO_IDS } from './modules/local/prep'
include { STAGE_GENOMEREPO; PARSE_ANNOTATIONS_BY_SOURCE; MAKE_PARSE_DONE; VALIDATE_PARSE_OUTPUTS; VALIDATE_GENESPACE_RESULTS; ORTHOFINDER_OR_SKIP; RUN_GENESPACE } from './modules/local/genespace'
include { PANGENES_PASS_FILTER; COLLAPSE_TANDEMS; WRITE_OG_FASTAS; MACSE_ALIGN_OG; MACSE_REPORT; IQTREE_OG; IQTREE_REPORT } from './modules/local/post_genespace'
include { WRITE_ALERAX_MAPPING; WRITE_ALERAX_FAMILIES; WRITE_ALERAX_MANIFEST; RUN_ALERAX; RUN_ALERAX_RANDOM; ALERAX_REPORT } from './modules/local/alerax'

// Helper functions
def resolveChrDict(genome) {
    def tsv = file("${params.chr_dict_dir}/${genome}.tsv")
    def bed = file("${params.chr_dict_dir}/${genome}_chr_lengths.bed")
    if (tsv.exists()) return tsv
    if (bed.exists()) return bed
    throw new IllegalArgumentException("No chr_dict found for ${genome}: expected ${tsv} or ${bed}")
}

def resolveGenomeFasta(genome) {
    def fasta = file("${params.fasta_dir}/${genome}.${params.ext.fasta}")
    if (fasta.exists()) return fasta
    throw new IllegalArgumentException("No genome FASTA found for ${genome}: expected ${fasta}")
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
    def idxOutgroup = header.indexOf('outgroup')

    if (idxGenome < 0 || idxSource < 0 || idxPloidy < 0) {
        throw new IllegalArgumentException("genomes TSV must contain columns: genome_id, genome_source, ploidy")
    }

    def rows = []
    lines.drop(1).each { line ->
        def cols = line.split('\t', -1) as List
        rows << [
            genome   : cols[idxGenome].trim(),
            source   : cols[idxSource].trim().toLowerCase(),
            ploidy   : cols[idxPloidy].trim() as Integer,
            outgroup : idxOutgroup >= 0 ? cols[idxOutgroup].trim().toLowerCase() : "no"
        ]
    }

    rows
}

def validateAleraxModels(models) {
    if (!(models instanceof List) || models.isEmpty()) {
        throw new IllegalArgumentException("AleRax models must be a non-empty list")
    }

    def seen = [] as Set

    models.each { m ->
        if (!m.model_id) {
            throw new IllegalArgumentException("Each AleRax model must define model_id")
        }
        if (!m.rec_model) {
            throw new IllegalArgumentException("AleRax model '${m.model_id}' is missing rec_model")
        }
        if (!m.model_parametrization) {
            throw new IllegalArgumentException("AleRax model '${m.model_id}' is missing model_parametrization")
        }
        if (m.gene_tree_samples == null) {
            throw new IllegalArgumentException("AleRax model '${m.model_id}' is missing gene_tree_samples")
        }
        if (seen.contains(m.model_id)) {
            throw new IllegalArgumentException("Duplicate AleRax model_id: ${m.model_id}")
        }
        seen << m.model_id
    }

    models
}

def resolveAleraxModels() {
    def models = params.alerax?.models

    if (models) {
        return validateAleraxModels(models)
    }

    return validateAleraxModels([
        [
            model_id: 'default',
            rec_model: params.alerax.rec_model,
            model_parametrization: params.alerax.model_parametrization,
            gene_tree_samples: params.alerax.gene_tree_samples
        ]
    ])
}

// Workflow
workflow {
    main:
    genome_repo_publish_ch = Channel.empty()
    genespace_publish_ch = Channel.empty()
    post_outputs_ch = Channel.empty()

    genomes_rows = readGenomesTable(params.genomes_tsv)

    def species_tree_path = params.species_tree?.toString()?.trim()

    if (params.use_species_tree_for_orthofinder && !species_tree_path) {
        error "--species_tree must be provided when --use_species_tree_for_orthofinder is true"
    }

    if (params.use_species_tree_for_alerax && !species_tree_path) {
        error "--species_tree must be provided when --use_species_tree_for_alerax is true"
    }

    genomes_ch = Channel
        .fromList(genomes_rows)
        .map { row ->
            tuple(
                row.genome,
                row.source,
                row.ploidy,
                file("${params.gff_dir}/${row.genome}.${params.ext.gff}"),
                file("${params.protein_dir}/${row.genome}.${params.ext.pep}"),
                resolveChrDict(row.genome),
                resolveGenomeFasta(row.genome)
            )
        }

    genomes_tsv_ch = Channel.value(file(params.genomes_tsv))
    cds_files_ch   = Channel.fromPath("${params.cds_dir}/*.cds").collect()

    def orthofinder_species_tree_arg = (
        params.use_species_tree_for_orthofinder && species_tree_path
    ) ? species_tree_path : ""

    primary_transcript_script_ch     = Channel.value(file('scripts/primary_transcript.py'))
    apply_chr_dict_script_ch         = Channel.value(file('scripts/apply_chr_dict_to_gff.py'))
    finalize_repo_ids_script_ch      = Channel.value(file('scripts/finalize_repo_ids.py'))
    parse_annotations_script_ch      = Channel.value(file('scripts/run_parse_annotations_by_source.R'))
    validate_parse_outputs_script_ch = Channel.value(file('scripts/validate_parse_outputs.py'))
    orthofinder_or_skip_script_ch    = Channel.value(file('scripts/orthofinder_or_skip.py'))
    run_genespace_script_ch          = Channel.value(file('scripts/run_genespace.R'))
    pangenes_pass_filter_script_ch   = Channel.value(file('scripts/pangenes_pass_filter.py'))
    collapse_tandems_script_ch       = Channel.value(file('scripts/collapse_tandems.py'))
    write_og_fastas_script_ch        = Channel.value(file('scripts/write_og_fastas.py'))

    genome_ids_ch = Channel.value(genomes_rows.collect { it.genome })

    def alerax_models = resolveAleraxModels()
    alerax_models_ch = Channel.fromList(alerax_models)
    alerax_models_list_ch = Channel.value(alerax_models)

    def validated_out = null
    def genespace_ready_out

    if (params.start_mode == 'full') {
        def prep_input_ch

        if (params.annotation?.run) {
            if (params.annotation.tool != 'annevo') {
                error "Unsupported annotation tool: ${params.annotation.tool}. Currently supported: annevo"
            }

            annevo_input_ch = genomes_ch.map { genome, source, ploidy, gff, pep, chr, fasta ->
                tuple(genome, source, ploidy, fasta, chr)
            }

            annevo_gff_out = RUN_ANNEVO(annevo_input_ch)

            annevo_fasta_out = ANNEVO_GFF_TO_FASTA(annevo_gff_out)

            prep_input_ch = annevo_fasta_out.map { genome, source, ploidy, gff, pep, chr, cds ->
                tuple(genome, source, ploidy, gff, pep, chr)
            }
            
            cds_files_ch = annevo_fasta_out
                .map { genome, source, ploidy, gff, pep, chr, cds -> cds }
                .collect()

        } else {
            prep_input_ch = genomes_ch.map { genome, source, ploidy, gff, pep, chr, fasta ->
                tuple(genome, source, ploidy, gff, pep, chr)
            }
        }

        primary_out = PRIMARY_TRANSCRIPT(
            prep_input_ch,
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

        genome_repo_publish_ch = staged_repo_out

        parsed_out = PARSE_ANNOTATIONS_BY_SOURCE(
            staged_repo_out,
            genomes_tsv_ch,
            parse_annotations_script_ch
        )

        validated_out = VALIDATE_PARSE_OUTPUTS(
            parsed_out[0],
            parsed_out[1],
            validate_parse_outputs_script_ch,
            genome_ids_ch
        )

    } else if (params.start_mode == 'parsed') {
        if (!params.existing_genespace_wd) {
            error "When --start_mode parsed is used, --existing_genespace_wd must be provided"
        }

        existing_wd_ch = Channel.value(file(params.existing_genespace_wd))
        parse_done_ch  = MAKE_PARSE_DONE()

        validated_out = VALIDATE_PARSE_OUTPUTS(
            existing_wd_ch,
            parse_done_ch,
            validate_parse_outputs_script_ch,
            genome_ids_ch
        )

    } else if (params.start_mode == 'genespace') {
        if (!params.existing_genespace_wd) {
            error "When --start_mode genespace is used, --existing_genespace_wd must be provided"
        }

        existing_wd_ch = Channel.value(file(params.existing_genespace_wd))
        genespace_ready_out = VALIDATE_GENESPACE_RESULTS(existing_wd_ch)

    } else {
        error "Unsupported start_mode: ${params.start_mode}. Use 'full', 'parsed', or 'genespace'."
    }

    if (params.start_mode != 'genespace') {
        def orthofinder_dir_arg
        def orthofinder_out = null

        def use_external_orthofinder = (
            params.run_external_orthofinder ||
            params.use_species_tree_for_orthofinder
        )

        if (params.existing_orthofinder_dir) {
            orthofinder_dir_arg = params.existing_orthofinder_dir

        } else if (use_external_orthofinder) {
            orthofinder_out = ORTHOFINDER_OR_SKIP(
                validated_out,
                genome_ids_ch,
                orthofinder_or_skip_script_ch,
                orthofinder_species_tree_arg
            )
            orthofinder_dir_arg = orthofinder_out[0]

        } else {
            orthofinder_dir_arg = ''
        }

        genespace_ready_out = RUN_GENESPACE(
            validated_out,
            orthofinder_dir_arg,
            genomes_tsv_ch,
            run_genespace_script_ch
        )
    }

    genespace_publish_ch = genespace_ready_out[0]

    pass_out = PANGENES_PASS_FILTER(
        genespace_ready_out,
        genomes_tsv_ch,
        pangenes_pass_filter_script_ch
    )

    post_outputs_ch = post_outputs_ch.mix(pass_out)

    pass_tsv_for_og = pass_out[0]
    og_list_for_og  = pass_out[1]

    if (params.collapse_tandems) {
        tandem_out = COLLAPSE_TANDEMS(
            pass_tsv_for_og,
            genespace_ready_out[0],
            genespace_ready_out[1],
            genomes_tsv_ch,
            collapse_tandems_script_ch
        )

        post_outputs_ch = post_outputs_ch.mix(tandem_out)

        pass_tsv_for_og = tandem_out[0]
        og_list_for_og  = tandem_out[2]
    }

    og_fastas_out = WRITE_OG_FASTAS(
        pass_tsv_for_og,
        og_list_for_og,
        genomes_tsv_ch,
        cds_files_ch,
        write_og_fastas_script_ch
    )

    post_outputs_ch = post_outputs_ch.mix(og_fastas_out)

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

    post_outputs_ch = post_outputs_ch.mix(macse_report_out)

    iqtree_in = macse_out.filter { og, aa, nt, status ->
        status.text.trim() == 'OK'
    }

    iqtree_out = IQTREE_OG(iqtree_in)

    iqtree_report_out = IQTREE_REPORT(
        iqtree_out
            .map { og, treefile, ufboot, status, nt -> status }
            .collect()
    )

    post_outputs_ch = post_outputs_ch.mix(iqtree_report_out)

    if (params.run_alerax) {
        alerax_map_out = WRITE_ALERAX_MAPPING(iqtree_out)

        families_out = WRITE_ALERAX_FAMILIES(
            alerax_map_out
                .flatMap { og, mapping, ufboot -> [mapping, ufboot] }
                .collect()
        )

        post_outputs_ch = post_outputs_ch.mix(families_out)

        manifest_out = WRITE_ALERAX_MANIFEST(alerax_models_list_ch)

        post_outputs_ch = post_outputs_ch.mix(manifest_out)

        if (params.use_species_tree_for_alerax) {
            species_tree_ch = Channel.value(file(species_tree_path))

            alerax_in = alerax_models_ch
                .combine(families_out)
                .combine(species_tree_ch)
                .map { model, fam, tree ->
                    tuple(fam, tree, model)
                }

            alerax_results_out = RUN_ALERAX(alerax_in)
        } else {
            alerax_in = alerax_models_ch
                .combine(families_out)
                .map { model, fam ->
                    tuple(fam, model)
                }

            alerax_results_out = RUN_ALERAX_RANDOM(alerax_in)
        }

        post_outputs_ch = post_outputs_ch.mix(alerax_results_out)

        alerax_report_inputs = alerax_results_out
            .map { model_id, model_dir -> model_dir }
            .collect()

        alerax_report_out = ALERAX_REPORT(
            alerax_report_inputs,
            manifest_out
        )

        post_outputs_ch = post_outputs_ch.mix(alerax_report_out)
    }

    publish:
    genomeRepo = genome_repo_publish_ch
    genespace_wd = genespace_publish_ch
    post_genespace = post_outputs_ch
}

output {
    genomeRepo {
        path '.'
    }

    genespace_wd {
        path { wd -> wd >> "${params.working_dir}/" }
    }

    post_genespace {
        path '.'
    }
}
