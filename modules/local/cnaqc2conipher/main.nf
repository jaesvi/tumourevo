process CNAQC2CONIPHER {
    tag "$meta.id"
    label "process_single"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/63f20f0f8ea78a7aecb1baa9aa35c03cdff9e18d19181d34561ed6d7e6376a0c/data' :
        'community.wave.seqera.io/library/r-cnaqc_r-cli_r-dplyr_r-readr_pruned:aebdc556849a14c7' }"

    input:
    tuple val(meta), path(rds_join), val(tumour_samples)

    output:
    tuple val(meta), path("*_conipher_input.tsv"), val(tumour_samples), emit: tsv
    path "versions.yml",                                               emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_conipher_input.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CNAqc: \$(Rscript -e "library(CNAqc); cat(as.character(packageVersion('CNAqc')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
