nextflow.enable.dsl=2

include { RUN_ANNEVO; ANNEVO_GFF_TO_FASTA } from './modules/local/annotation'
include { PRIMARY_TRANSCRIPT; FINALIZE_REPO_IDS } from './modules/local/prep'
include { STAGE_GENOMEREPO; PARSE_ANNOTATIONS_BY_SOURCE; MAKE_PARSE_DONE; VALIDATE_PARSE_OUTPUTS; VALIDATE_GENESPACE_RESULTS; ORTHOFINDER_OR_SKIP; RUN_GENESPACE } from './modules/local/genespace'
include { PANGENES_PASS_FILTER; COLLAPSE_TANDEMS; WRITE_OG_FASTAS; MACSE_ALIGN_OG; MACSE_REPORT; IQTREE_OG; IQTREE_REPORT } from './modules/local/post_genespace'
include { ALERAX_WORKFLOW } from './modules/local/alerax'
include { REDIPLOIDISATION } from './modules/local/rediploidisation'


// Helper functions
def resolveChrDict(genome) {
    def tsv = file("${params.chr_dict_dir}/${genome}.tsv")
    def bed = file("${params.chr_dict_dir}/${genome}_chr_lengths.bed")
    if (tsv.exists()) return tsv
    if (bed.exists()) return bed
    throw new IllegalArgumentException("No chr_dict found for ${genome}: expected ${tsv} or ${bed}")
}

def resolveGenomeFasta(genome) {
    def exts = [
        params.ext.fasta,
        "${params.ext.fasta}.gz",
        "fasta",
        "fasta.gz",
        "fna",
        "fna.gz"
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

def makeIqtreeDirChannelFromDir(treeDir) {
    Channel
        .fromPath("${treeDir}/og_*/og_*_iqtree.treefile", checkIfExists: true)
        .map { treefile ->
            def m = (treefile.baseName =~ /^og_(.+)_iqtree$/)
            if (!m) {
                throw new IllegalArgumentException("Could not parse OG from IQ-TREE filename: ${treefile}")
            }

            def og = m[0][1]
            tuple(og, treefile.parent)
        }
}

def makeAleraxInputChannelFromDirs(treeDir, ntDir) {
    Channel
        .fromPath("${treeDir}/og_*/og_*_iqtree.treefile", checkIfExists: true)
        .map { treefile ->
            def m = (treefile.baseName =~ /^og_(.+)_iqtree$/)
            if (!m) {
                throw new IllegalArgumentException("Could not parse OG from IQ-TREE filename: ${treefile}")
            }

            def og = m[0][1]
            def iqtree_dir = treefile.parent
            def nt = file("${ntDir}/og_${og}_NT.fasta", checkIfExists: true)

            tuple(og, iqtree_dir, nt)
        }
}

// Workflow
workflow {
    main:
    genome_repo_publish_ch = Channel.empty()
    genespace_publish_ch = Channel.empty()
    post_outputs_ch = Channel.empty()

    def species_tree_path = params.species_tree?.toString()?.trim()

    def redipDefaults = [
        positions_source: 'bed',
        positions: '',
        position_format: 'auto',

        position_has_header: true,
        position_key_column: 'id',
        position_chr_column: 'chr',
        position_start_column: 'start',
        position_end_column: 'end',
        position_species_column: 'genome',

        gene_trees_dir: '',
        genespace_wd: '',

        species_tree_format: 1,
        gene_tree_format: 1,
        tip_separator: '|',
        label_format: 'species_chr_gene',
        copy_mode: 'target_exactly_n',
        required_copies: 2,
        min_tips: 1,
        position_key_type: 'gene',
        cleanup_tmp: true
    ]

    def redipParams = redipDefaults + (params.rediploidisation ?: [:])

    /*
     * =========================
     * ALERAX-ONLY MODE
     * =========================
     */
    
    if (params.start_mode == 'alerax') {
        if (!params.run_alerax) {
            error "--run_alerax must be true when --start_mode alerax"
        }
    
        if (params.use_species_tree_for_alerax && !species_tree_path) {
            error "--species_tree must be provided when --start_mode alerax and --use_species_tree_for_alerax is true"
        }
    
        if (!params.alerax?.gene_trees_dir?.toString()?.trim()) {
            error "--alerax.gene_trees_dir must be provided when --start_mode alerax"
        }
    
        if (!params.alerax?.nt_alignments_dir?.toString()?.trim()) {
            error "--alerax.nt_alignments_dir must be provided when --start_mode alerax"
        }
    
        def alerax_models = resolveAleraxModels()
        def alerax_models_ch = Channel.fromList(alerax_models)
    
        def species_tree_ch = params.use_species_tree_for_alerax
            ? Channel.value(file(species_tree_path, checkIfExists: true))
            : Channel.empty()
    
        def alerax_input_ch = makeAleraxInputChannelFromDirs(
            params.alerax.gene_trees_dir,
            params.alerax.nt_alignments_dir
        )
    
        def alerax_out = ALERAX_WORKFLOW(
            alerax_input_ch,
            species_tree_ch,
            alerax_models_ch
        )
    
        post_outputs_ch = post_outputs_ch.mix(alerax_out.families)
        post_outputs_ch = post_outputs_ch.mix(alerax_out.manifest)
        post_outputs_ch = post_outputs_ch.mix(
            alerax_out.results.map { model_id, dir -> dir }.collect()
        )
        post_outputs_ch = post_outputs_ch.mix(alerax_out.report)

    /*
     * =========================
     * REDIP-ONLY MODE
     * =========================
     */

    } else if (params.start_mode == 'redip') {
    
        genomes_rows = readGenomesTable(params.genomes_tsv)
        genomes_tsv_ch = Channel.value(file(params.genomes_tsv, checkIfExists: true))
    
        if (!params.run_rediploidisation) {
            error "--run_rediploidisation must be true when --start_mode redip"
        }

        if (!species_tree_path) {
            error "--species_tree must be provided when --start_mode redip"
        }

        if (!redipParams.gene_trees_dir?.toString()?.trim()) {
            error "--rediploidisation.gene_trees_dir must be provided when --start_mode redip"
        }
        
        if (
            redipParams.positions_source != 'positions' &&
            !redipParams.genespace_wd?.toString()?.trim()
        ) {
            error "--rediploidisation.genespace_wd must be provided unless positions_source = positions"
        }
        
        if (
            redipParams.positions_source == 'positions' &&
            !redipParams.positions?.toString()?.trim()
        ) {
            error "--rediploidisation.positions must be provided when positions_source = positions"
        }
        
        def redip_genespace_wd_ch = redipParams.genespace_wd?.toString()?.trim()
            ? Channel.value(file(redipParams.genespace_wd))
            : Channel.value(file('.'))
        
        def redip_out = REDIPLOIDISATION(
            genomes_tsv_ch,
            Channel.value(file(species_tree_path)),
            makeIqtreeDirChannelFromDir(redipParams.gene_trees_dir),
            redip_genespace_wd_ch,
            redipParams
        )

        post_outputs_ch = post_outputs_ch.mix(
            redip_out.rooted_trees.flatMap { og, tree, summary -> [tree, summary] }.collect()
        )
        post_outputs_ch = post_outputs_ch.mix(redip_out.branch_definitions)
        post_outputs_ch = post_outputs_ch.mix(
            redip_out.classifications.map { species, file -> file }.collect()
        )
        post_outputs_ch = post_outputs_ch.mix(
            redip_out.circos_links.map { species, file -> file }.collect()
        )
        post_outputs_ch = post_outputs_ch.mix(
            redip_out.circos_plots.map { species, dir -> dir }.collect()
        )
        post_outputs_ch = post_outputs_ch.mix(redip_out.report)

    } else {
    
        /*
         * =========================
         * FULL / PARSED / GENESPACE PIPELINE
         * =========================
         */
    
        genomes_rows = readGenomesTable(params.genomes_tsv)
        genomes_tsv_ch = Channel.value(file(params.genomes_tsv, checkIfExists: true))
    
        if (params.use_species_tree_for_orthofinder && !species_tree_path) {
            error "--species_tree must be provided when --use_species_tree_for_orthofinder is true"
        }
    
        if (params.use_species_tree_for_alerax && params.run_alerax && !species_tree_path) {
            error "--species_tree must be provided when --run_alerax is true and --use_species_tree_for_alerax is true"
        }

        genomes_ch = Channel
            .fromList(genomes_rows)
            .map { row ->
                def fasta = params.annotation?.run
                    ? resolveGenomeFasta(row.genome)
                    : null
        
                tuple(
                    row.genome,
                    row.source,
                    row.ploidy,
                    file("${params.gff_dir}/${row.genome}.${params.ext.gff}"),
                    file("${params.protein_dir}/${row.genome}.${params.ext.pep}"),
                    resolveChrDict(row.genome),
                    fasta
                )
            }

        cds_files_ch = Channel.fromPath("${params.cds_dir}/*.cds").collect()

        def orthofinder_species_tree_ch = (
            params.use_species_tree_for_orthofinder && species_tree_path
        ) ? Channel.value(file(species_tree_path, checkIfExists: true)) : Channel.value([])

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
            error "Unsupported start_mode: ${params.start_mode}. Use 'full', 'parsed', 'genespace', 'alerax', or 'redip'."
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
                    validated_out[0],
                    validated_out[1],
                    genome_ids_ch,
                    orthofinder_or_skip_script_ch,
                    orthofinder_species_tree_ch
                )
                orthofinder_dir_arg = orthofinder_out[0]

            } else {
                orthofinder_dir_arg = ''
            }

            genespace_ready_out = RUN_GENESPACE(
                validated_out[0],
                validated_out[1],
                orthofinder_dir_arg,
                genomes_tsv_ch,
                run_genespace_script_ch
            )
        }

        genespace_publish_ch = genespace_ready_out[0]

        pass_out = PANGENES_PASS_FILTER(
            genespace_ready_out[0],
            genespace_ready_out[1],
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
            og_list_for_og  = tandem_out[1]
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

        macse_all_outputs_ch = macse_out
            .flatMap { og, aa, nt, status, log -> [aa, nt, status, log] }
            .collect()

        macse_report_out = MACSE_REPORT(macse_all_outputs_ch)

        post_outputs_ch = post_outputs_ch.mix(macse_all_outputs_ch)
        post_outputs_ch = post_outputs_ch.mix(macse_report_out)

        iqtree_in = macse_out.filter { og, aa, nt, status, log ->
            status.text.trim() == 'OK'
        }
        
        iqtree_out = IQTREE_OG(iqtree_in)
        
        iqtree_dirs_ch = iqtree_out
            .map { og, iqtree_dir -> iqtree_dir }
            .collect()
        
        iqtree_report_out = IQTREE_REPORT(iqtree_dirs_ch)
        
        post_outputs_ch = post_outputs_ch.mix(iqtree_dirs_ch)
        post_outputs_ch = post_outputs_ch.mix(iqtree_report_out)
        
        iqtree_for_alerax_ch = iqtree_out
            .join(
                macse_out
                    .filter { og, aa, nt, status, log ->
                        status.text.trim() == 'OK'
                    }
                    .map { og, aa, nt, status, log ->
                        tuple(og, nt)
                    }
            )
            .map { og, iqtree_dir, nt ->
                tuple(og, iqtree_dir, nt)
            }
        
        iqtree_for_redip_ch = iqtree_out
            .map { og, iqtree_dir ->
                tuple(og, iqtree_dir)
            }
        
        if (params.run_alerax) {
            def alerax_models = resolveAleraxModels()
            def alerax_models_ch = Channel.fromList(alerax_models)
        
            def species_tree_ch = params.use_species_tree_for_alerax
                ? Channel.value(file(species_tree_path))
                : Channel.empty()
        
            def alerax_out = ALERAX_WORKFLOW(
                iqtree_for_alerax_ch,
                species_tree_ch,
                alerax_models_ch
            )
        
            post_outputs_ch = post_outputs_ch.mix(alerax_out.families)
            post_outputs_ch = post_outputs_ch.mix(alerax_out.manifest)
            post_outputs_ch = post_outputs_ch.mix(
                alerax_out.results.map { model_id, dir -> dir }.collect()
            )
            post_outputs_ch = post_outputs_ch.mix(alerax_out.report)
        }

        if (params.run_rediploidisation) {
            if (!species_tree_path) {
                error "--species_tree must be provided when --run_rediploidisation is true"
            }
        
            if (
                redipParams.positions_source == 'positions' &&
                !redipParams.positions?.toString()?.trim()
            ) {
                error "--rediploidisation.positions must be provided when positions_source = positions"
            }
        
            def redip_gene_trees_dir = redipParams.gene_trees_dir?.toString()?.trim()
        
            def redip_iqtree_ch = redip_gene_trees_dir
                ? makeIqtreeDirChannelFromDir(redip_gene_trees_dir)
                : iqtree_for_redip_ch
        
            def redip_genespace_wd_ch = redipParams.genespace_wd?.toString()?.trim()
                ? Channel.value(file(redipParams.genespace_wd))
                : genespace_ready_out[0]
        
            def redip_out = REDIPLOIDISATION(
                genomes_tsv_ch,
                Channel.value(file(species_tree_path)),
                redip_iqtree_ch,
                redip_genespace_wd_ch,
                redipParams
            )
        
            post_outputs_ch = post_outputs_ch.mix(
                redip_out.rooted_trees.flatMap { og, tree, summary -> [tree, summary] }.collect()
            )
            post_outputs_ch = post_outputs_ch.mix(redip_out.branch_definitions)
            post_outputs_ch = post_outputs_ch.mix(
                redip_out.classifications.map { species, file -> file }.collect()
            )
            post_outputs_ch = post_outputs_ch.mix(
                redip_out.circos_links.map { species, file -> file }.collect()
            )
            post_outputs_ch = post_outputs_ch.mix(
                redip_out.circos_plots.map { species, dir -> dir }.collect()
            )
            post_outputs_ch = post_outputs_ch.mix(redip_out.report)
        }
    }
    /*
     * =========================
     * SINGLE PUBLISH BLOCK
     * =========================
     */

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
