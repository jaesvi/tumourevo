process CONIPHER_CLUSTERING {
    tag "$meta.id"
    label "process_high"
    label "error_ignore"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/conipher:2.2.0--r40hdfd78af_0' :
        'quay.io/biocontainers/conipher:2.2.0--r40hdfd78af_0' }"

    input:
    tuple val(meta), path(input_tsv), val(tumour_samples)

    output:
    tuple val(meta), path("*.SCoutput.CLEAN.tsv"), emit: cluster_tsv
    tuple val(meta), path("*.SCoutput.FULL.tsv") , emit: cluster_full , optional: true
    tuple val(meta), path("*.SCoutput.DIRTY.tsv"), emit: cluster_dirty, optional: true
    path "versions.yml"                          , emit: versions    , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def case_id = meta.patient ?: prefix
    """
    touch ${case_id}.SCoutput.CLEAN.tsv
    touch ${case_id}.SCoutput.FULL.tsv
    touch ${case_id}.SCoutput.DIRTY.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CONIPHER: \$(Rscript -e "library(CONIPHER); cat(as.character(packageVersion('CONIPHER')))")
    END_VERSIONS
    """
}
