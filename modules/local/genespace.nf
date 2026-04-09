process STAGE_GENOMEREPO {
    tag "stage_genomerepo"

    input:
    path finalized_files

    output:
    path("genespace/genomeRepo")

    script:
    def stageLines = finalized_files.collect { f ->
        def name = f.getName()

        if (name.endsWith('.final.gff3')) {
            def genome = name.replaceFirst(/\.final\.gff3$/, '')
            return "mkdir -p genespace/genomeRepo/${genome}\ncp ${f} genespace/genomeRepo/${genome}/${genome}.final.gff3"
        }
        if (name.endsWith('.final.pep')) {
            def genome = name.replaceFirst(/\.final\.pep$/, '')
            return "mkdir -p genespace/genomeRepo/${genome}\ncp ${f} genespace/genomeRepo/${genome}/${genome}.final.pep"
        }
        if (name.endsWith('.chr.tsv')) {
            def genome = name.replaceFirst(/\.chr\.tsv$/, '')
            return "mkdir -p genespace/genomeRepo/chr_dict\ncp ${f} genespace/genomeRepo/chr_dict/${genome}.tsv"
        }

        return null
    }.findAll { it != null }.join('\n')

    """
    rm -rf genespace/genomeRepo
    mkdir -p genespace/genomeRepo

    ${stageLines}

    test -d genespace/genomeRepo
    """
}

process PARSE_ANNOTATIONS_BY_SOURCE {
    tag "parse_annotations_by_source"

    input:
    path genomeRepo
    path genomes_tsv
    path parse_annotations_script

    output:
    path("genespace/${params.working_dir}")
    path("genespace/${params.working_dir}/parse_annotations.done")

    script:
    """
    mkdir -p genespace/${params.working_dir}

    Rscript --vanilla ${parse_annotations_script} \
      --genomes-tsv ${genomes_tsv} \
      --raw-genomerepo ${genomeRepo} \
      --genespace-wd genespace/${params.working_dir} \
      > parse_annotations_by_source.log 2>&1

    touch genespace/${params.working_dir}/parse_annotations.done
    """
}

process VALIDATE_PARSE_OUTPUTS {
    tag "validate_parse_outputs"

    input:
    path genespace_wd
    path parse_done
    path validate_parse_outputs_script

    output:
    path("${genespace_wd.getFileName()}")
    path("parse_outputs.ok")

    script:
    def genomes = file(params.genomes_tsv)
        .readLines()
        .drop(1)
        .findAll { it.trim() }
        .collect { it.split('\t', -1)[0].trim() }
        .join(' ')

    """
    python ${validate_parse_outputs_script} \
      --genespace-wd ${genespace_wd} \
      --genomes ${genomes} \
      > validate_parse_outputs.log 2>&1

    touch parse_outputs.ok
    """
}

process ORTHOFINDER_OR_SKIP {
    tag "orthofinder"

    input:
    path genespace_wd
    path parse_ok
    path genomes_tsv
    path orthofinder_or_skip_script

    output:
    path("genespace/orthofinder")
    path("genespace/orthofinder.done")

    script:
    def genomes = genomes_tsv
        .readLines()
        .drop(1)
        .findAll { it.trim() }
        .collect { it.split('\t', -1)[0].trim() }
        .join(' ')

    """
    mkdir -p genespace/orthofinder

    python ${orthofinder_or_skip_script} \
      --threads ${task.cpus} \
      --peptide-dir genespace/${params.working_dir}/peptide \
      --orthofinder-dir genespace/orthofinder \
      --orthofinder-bin ${params.orthofinder_bin} \
      --genomes ${genomes} \
      --force ${params.force_orthofinder} \
      > orthofinder.log 2>&1

    touch genespace/orthofinder.done
    """
}

process RUN_GENESPACE {
    tag "genespace"

    input:
    path genespace_wd
    path parse_ok
    path orthofinder_dir
    path orthofinder_done
    path genomes_tsv
    path run_genespace_script

    output:
    path("genespace/${params.working_dir}")
    path("genespace/${params.working_dir}/genespace.done")

    script:
    """
    Rscript --vanilla ${run_genespace_script} \
      --config ${params.config_for_r ?: 'config.yaml'} \
      --genespace-wd genespace/${params.working_dir} \
      --orthofinder-dir genespace/orthofinder \
      --genomes-tsv ${genomes_tsv} \
      > genespace.log 2>&1

    touch genespace/${params.working_dir}/genespace.done
    """
}
