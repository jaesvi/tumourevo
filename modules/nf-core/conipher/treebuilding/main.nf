process CONIPHER_TREEBUILDING {
    tag "$meta.id"
    label "process_medium"
    label "error_ignore"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/conipher:2.2.0--r40hdfd78af_0' :
        'quay.io/biocontainers/conipher:2.2.0--r40hdfd78af_0' }"

    input:
    tuple val(meta), path(cluster_tsv)

    output:
    tuple val(meta), path("trees/*.tree.RDS")             , emit: tree_rds , optional: true
    tuple val(meta), path("trees/pytree_and_bar.pdf")     , emit: tree_pdf , optional: true
    tuple val(meta), path("trees/pytree_multipletrees.pdf"), emit: trees_pdf, optional: true
    path "versions.yml"                               , emit: versions  , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def case_id = meta.patient ?: prefix
    """
    mkdir -p trees
    touch trees/${case_id}.tree.RDS
    touch trees/pytree_and_bar.pdf
    touch trees/pytree_multipletrees.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CONIPHER: \$(Rscript -e "library(CONIPHER); cat(as.character(packageVersion('CONIPHER')))")
    END_VERSIONS
    """
}
