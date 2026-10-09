#!/usr/bin/env Rscript
parse_args = function(x) {
  x = gsub("\\\\[","",x)
  x = gsub("\\\\]","",x)
  args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
  args_vals = lapply(args_list, function(x) {
    x_splt = strsplit(x, split=":")[[1]]
    c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
  })
  args_vals = lapply(args_vals, function(z){ length(z) = 2; z})
  parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
  parsed_args[! is.na(parsed_args)]
}

opt = list(
    prefix = ifelse("$task.ext.prefix" == "null", "$meta.id", "$task.ext.prefix"),
    qc_chr = "FALSE"
)
args_opt = parse_args("$task.ext.args")
for (ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]

library(dplyr)
library(tidyr)
library(CNAqc)
library(cli)

requested_samples = gsub("\\\\[|\\\\]", "", "$tumour_samples")
requested_samples = trimws(unlist(strsplit(requested_samples, ",")))
requested_samples = requested_samples[nzchar(requested_samples)]

extract_sample_table = function(cnaqc_sample, sample_name) {
    mutations = cnaqc_sample[["mutations"]]
    if (is.null(mutations) || nrow(mutations) == 0) {
        return(NULL)
    }
    if ("blacklisted" %in% colnames(mutations)) {
        mutations = mutations %>% dplyr::filter(is.na(blacklisted) | blacklisted == FALSE)
    }
    if (!all(c("chr", "from", "ref", "alt") %in% colnames(mutations))) {
        cli::cli_alert_warning("Sample {sample_name} is missing required mutation columns")
        return(NULL)
    }
    if (!"karyotype" %in% colnames(mutations)) {
        cli::cli_alert_warning("Sample {sample_name} has no karyotype column; skipping")
        return(NULL)
    }

    nr = if ("NR" %in% colnames(mutations)) mutations[["NR"]] else NA_real_
    nv = if ("NV" %in% colnames(mutations)) mutations[["NV"]] else NA_real_
    dp = if ("DP" %in% colnames(mutations)) mutations[["DP"]] else nr + nv

    ploidy = cnaqc_sample[["ploidy"]]
    if (is.null(ploidy) || length(ploidy) == 0 || is.na(ploidy[1])) {
        cli::cli_alert_warning("Sample {sample_name} has no ploidy; defaulting to 2")
        ploidy = 2
    }

    mutations %>%
        dplyr::mutate(
            CASE_ID = "$meta.patient",
            SAMPLE = gsub("-", ".", sample_name),
            CHR = gsub("^chr", "", chr, ignore.case = TRUE),
            POS = as.integer(from),
            REF = ref,
            ALT = alt,
            REF_COUNT = as.numeric(nr),
            VAR_COUNT = as.numeric(nv),
            DEPTH = as.numeric(dp),
            ACF = as.numeric(cnaqc_sample[["purity"]]),
            PLOIDY = as.numeric(ploidy[1]),
            karyotype = as.character(karyotype)
        ) %>%
        tidyr::separate(karyotype, into = c("COPY_NUMBER_A", "COPY_NUMBER_B"), sep = ":", remove = TRUE, convert = TRUE, extra = "drop", fill = "right") %>%
        dplyr::filter(
            CHR %in% as.character(1:22),
            !is.na(REF_COUNT),
            !is.na(VAR_COUNT),
            !is.na(DEPTH),
            !is.na(COPY_NUMBER_A),
            !is.na(COPY_NUMBER_B),
            DEPTH > 0
        ) %>%
        dplyr::select(CASE_ID, SAMPLE, CHR, POS, REF, ALT, REF_COUNT, VAR_COUNT, DEPTH, COPY_NUMBER_A, COPY_NUMBER_B, ACF, PLOIDY)
}

tables_from_mcnaqc = function(obj) {
    sample_names = get_sample_name(obj)
    if (length(requested_samples) > 0) {
        sample_names = intersect(sample_names, requested_samples)
    }
    if (length(sample_names) == 0) {
        return(NULL)
    }
    samples = get_sample(m_cnaqc_obj = obj, sample = sample_names, which_obj = "original")
    lapply(sample_names, function(s) extract_sample_table(samples[[s]], s)) %>%
        dplyr::bind_rows()
}

multi_cnaqc = readRDS(file = "$rds_join")
if (as.logical(opt[["qc_chr"]])) {
    parts = lapply(names(multi_cnaqc), function(c) tables_from_mcnaqc(multi_cnaqc[[c]]))
} else {
    parts = list(tables_from_mcnaqc(multi_cnaqc))
}
parts = parts[!sapply(parts, is.null)]
joint_table = if (length(parts) == 0) {
    data.frame()
} else {
    dplyr::bind_rows(parts)
}

if (nrow(joint_table) == 0) {
    cli::cli_abort("No autosomal mutations with copy-number and read counts could be formatted for CONIPHER")
}

options(scipen = 999)
write.table(
    joint_table,
    file = paste0(opt[["prefix"]], "_conipher_input.tsv"),
    append = FALSE,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
)

f = file("versions.yml", "w")
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    CNAqc:", sessionInfo()\$otherPkgs\$CNAqc\$Version), f)
writeLines(paste("    dplyr:", sessionInfo()\$otherPkgs\$dplyr\$Version), f)
close(f)
