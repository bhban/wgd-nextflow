process STAGE_GENOMEREPO {
    tag "stage_genomerepo"

    input:
    path finalized_files

    output:
    path("genomeRepo")

    script:
    def stageLines = finalized_files.collect { f ->
        def name = f.getName()

        if (name.endsWith('.final.gff3')) {
            def genome = name.replaceFirst(/\.final\.gff3$/, '')
            return "mkdir -p genomeRepo/${genome}\ncp ${f} genomeRepo/${genome}/${genome}.final.gff3"
        }
        if (name.endsWith('.final.pep')) {
            def genome = name.replaceFirst(/\.final\.pep$/, '')
            return "mkdir -p genomeRepo/${genome}\ncp ${f} genomeRepo/${genome}/${genome}.final.pep"
        }
        if (name.endsWith('.chr.tsv')) {
            def genome = name.replaceFirst(/\.chr\.tsv$/, '')
            return "mkdir -p genomeRepo/chr_dict\ncp ${f} genomeRepo/chr_dict/${genome}.tsv"
        }

        return null
    }.findAll { it != null }.join('\n')

    """
    rm -rf genomeRepo
    mkdir -p genomeRepo

    ${stageLines}

    test -d genomeRepo
    """
}

process PARSE_ANNOTATIONS_BY_SOURCE {
    tag "parse_annotations_by_source"

    input:
    path genomeRepo
    path genomes_tsv
    path parse_annotations_script

    output:
    path("${params.working_dir}")
    path("parse_annotations.done")

    script:
    """
    mkdir -p ${params.working_dir}

    Rscript --vanilla ${parse_annotations_script} \
      --genomes-tsv ${genomes_tsv} \
      --raw-genomerepo ${genomeRepo} \
      --genespace-wd ${params.working_dir} \
      > parse_annotations_by_source.log 2>&1

    touch parse_annotations.done
    """
}

process MAKE_PARSE_DONE {
    tag "make_parse_done"

    output:
    path("parse_annotations.done")

    script:
    """
    touch parse_annotations.done
    """
}

process VALIDATE_PARSE_OUTPUTS {
    tag "validate_parse_outputs"

    input:
    path genespace_wd
    path parse_done
    path validate_parse_outputs_script
    val genomes

    output:
    path("${genespace_wd.getFileName()}")
    path("parse_outputs.ok")

    script:
    def genomes_arg = genomes.join(' ')

    """
    python ${validate_parse_outputs_script} \
      --genespace-wd ${genespace_wd} \
      --genomes ${genomes_arg} \
      > validate_parse_outputs.log 2>&1

    touch parse_outputs.ok
    """
}

process ORTHOFINDER_OR_SKIP {
    tag "orthofinder"

    input:
    path genespace_wd
    path parse_ok
    val genomes
    path orthofinder_or_skip_script
    val species_tree_arg

    output:
    path("orthofinder")
    path("orthofinder.done")

    script:
    def genomes_arg = genomes.join(' ')
    def of_species_tree_arg = species_tree_arg ? "--species-tree ${species_tree_arg}" : ""

    """
    rm -rf orthofinder

    python ${orthofinder_or_skip_script} \
      --threads ${task.cpus} \
      --analysis-threads ${params.orthofinder_analysis_threads} \
      --peptide-dir ${genespace_wd}/peptide \
      --orthofinder-dir orthofinder \
      --orthofinder-bin ${params.orthofinder_bin} \
      --genomes ${genomes_arg} \
      --force ${params.force_orthofinder} \
      ${of_species_tree_arg} \
      > orthofinder.log 2>&1

    touch orthofinder.done
    """
}

process RUN_GENESPACE {
    tag "genespace"

    input:
    path genespace_wd
    path parse_ok
    val orthofinder_dir_arg
    path genomes_tsv
    path run_genespace_script

    output:
    path("${genespace_wd.getFileName()}")
    path("genespace.done")

    script:
    def raw_of_arg = params.genespace.rawOrthofinderDir ? "--raw-orthofinder-dir ${params.genespace.rawOrthofinderDir}" : ""
    def stage_external_of = orthofinder_dir_arg ? """
    rm -rf ${genespace_wd}/orthofinder
    cp -r ${orthofinder_dir_arg} ${genespace_wd}/orthofinder
    """ : ""

    """
    ${stage_external_of}

    Rscript --vanilla ${run_genespace_script} \
      --genespace-wd ${genespace_wd} \
      ${raw_of_arg} \
      --genomes-tsv ${genomes_tsv} \
      --blk-size ${params.genespace.blkSize} \
      --orthofinder-in-blk ${params.genespace.orthofinder_in_blk} \
      > genespace.log 2>&1

    touch genespace.done
    """
}
