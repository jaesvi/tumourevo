#!/usr/bin/env Rscript
parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
    args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
    args_vals = lapply(args_list, function(x) {
        x_splt = strsplit(x, split=":")[[1]]
        c(x_splt[1], paste(x_splt[2:length(x_splt)], collapse=":"))
    })
    args_vals = lapply(args_vals, function(z){ length(z) = 2; z})
    parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

as_bool = function(x) {
    if (is.logical(x) && length(x) == 1) return(x)
    toupper(as.character(x)) %in% c("TRUE", "T", "1")
}

opt = list(
    prefix = "$meta.patient",
    ccf_buffer = "10",
    pval_cutoff = "0.01",
    use_boot = "TRUE",
    merge_clusters = "TRUE",
    min_ccf = "0.01",
    min_cluster_size = "5",
    multi_trees = "TRUE"
)
args_opt = parse_args("$task.ext.args")
for (ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]

library(CONIPHER)

out_dir = "trees"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

conipher_treebuilding(
    input_tsv_loc = "$cluster_tsv",
    out_dir = out_dir,
    prefix = opt[["prefix"]],
    ccf_buffer = as.numeric(opt[["ccf_buffer"]]),
    pval_cutoff = as.numeric(opt[["pval_cutoff"]]),
    use_boot = as_bool(opt[["use_boot"]]),
    merge_clusters = as_bool(opt[["merge_clusters"]]),
    min_ccf = as.numeric(opt[["min_ccf"]]),
    min_cluster_size = as.integer(opt[["min_cluster_size"]]),
    multi_trees = as_bool(opt[["multi_trees"]])
)

f = file("versions.yml", "w")
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    CONIPHER:", as.character(packageVersion("CONIPHER"))), f)
close(f)
