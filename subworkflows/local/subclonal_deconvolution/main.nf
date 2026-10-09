// SUBCLONAL DECONVOLUTION WORKFLOW

include { MOBSTER } from "../../../modules/nf-core/mobster/main"
include { VIBER } from "../../../modules/nf-core/viber/main"
include { PYCLONEVI } from "../../../modules/nf-core/pyclonevi/main"
include { FORMATTER } from "../../../subworkflows/local/formatter/main"
include { CNAQC2CONIPHER } from "../../../modules/local/cnaqc2conipher/main"
include { CONIPHER_CLUSTERING } from "../../../modules/nf-core/conipher/clustering/main"
include { CONIPHER_TREEBUILDING } from "../../../modules/nf-core/conipher/treebuilding/main"
include { CTREE as CTREE_MOBSTER } from "../../../modules/nf-core/ctree/main"
include { CTREE as CTREE_PYCLONEVI } from "../../../modules/nf-core/ctree/main"
include { CTREE as CTREE_VIBER } from "../../../modules/nf-core/ctree/main"

workflow SUBCLONAL_DECONVOLUTION {
    take:
    rds_join // tuple val(meta), path("*.rds"), val(tumour_samples), emit: rds

    main:
    ch_versions = Channel.empty()
    mobster_pdf = null
    ctree_mobster_pdf = null
    viber_pdf = null
    ctree_viber_pdf = null
    pyclone_fits = null
    pyclone_best = null
    pyclone_table = null
    ctree_pyclone_pdf = null
    conipher_cluster_tsv = null
    conipher_tree_pdf = null

    if (params.tools && params.tools.split(",").contains("mobster")) {
        joinCNAqc = rds_join.transpose().map{ meta, rds, sample ->
            meta = meta + ["tumour_sample": sample, "id":"${meta.dataset}_${meta.patient}_${sample}"]
            [meta, rds]}

        MOBSTER(joinCNAqc)
        ch_versions = ch_versions.mix(MOBSTER.out.versions_mobster)

        CTREE_MOBSTER(MOBSTER.out.mobster_best_rds)
        ch_versions = ch_versions.mix(CTREE_MOBSTER.out.versions)

        mobster_pdf = MOBSTER.out.mobster_report_pdf
        ctree_mobster_pdf = CTREE_MOBSTER.out.ctree_report_pdf

    }

    if (params.tools && params.tools.split(",").contains("viber")) {
        VIBER(rds_join)
        ch_versions = ch_versions.mix(VIBER.out.versions_viber)

        CTREE_VIBER(VIBER.out.viber_rds)
        ch_versions = ch_versions.mix(CTREE_VIBER.out.versions)

        viber_pdf = VIBER.out.viber_report_pdf
        ctree_viber_pdf = CTREE_VIBER.out.ctree_report_pdf
    }

    if (params.tools && params.tools.split(",").contains("pyclone-vi")) {
        FORMATTER(rds_join, "rds")
        PYCLONEVI(FORMATTER.out.out_data)
        CTREE_PYCLONEVI(PYCLONEVI.out.ctree_input)

        ch_versions = ch_versions.mix(FORMATTER.out.versions)
        ch_versions = ch_versions.mix(PYCLONEVI.out.versions_pyclonevi)
        ch_versions = ch_versions.mix(CTREE_PYCLONEVI.out.versions)

        pyclone_fits = PYCLONEVI.out.pyclone_all_fits
        pyclone_best = PYCLONEVI.out.pyclone_best_fit
        pyclone_table = FORMATTER.out.out_data
        ctree_pyclone_pdf = CTREE_PYCLONEVI.out.ctree_report_pdf
    }

    if (params.tools && params.tools.split(",").contains("conipher")) {
        CNAQC2CONIPHER(rds_join)
        CONIPHER_CLUSTERING(CNAQC2CONIPHER.out.tsv)
        CONIPHER_TREEBUILDING(CONIPHER_CLUSTERING.out.cluster_tsv)

        ch_versions = ch_versions.mix(CNAQC2CONIPHER.out.versions)
        ch_versions = ch_versions.mix(CONIPHER_CLUSTERING.out.versions)
        ch_versions = ch_versions.mix(CONIPHER_TREEBUILDING.out.versions)

        conipher_cluster_tsv = CONIPHER_CLUSTERING.out.cluster_tsv
        conipher_tree_pdf = CONIPHER_TREEBUILDING.out.tree_pdf
    }

    emit:
    pyclone_fits
    pyclone_best
    ctree_pyclone_pdf
    viber_pdf
    ctree_viber_pdf
    mobster_pdf
    ctree_mobster_pdf
    pyclone_table
    conipher_cluster_tsv
    conipher_tree_pdf
    versions = ch_versions
}
