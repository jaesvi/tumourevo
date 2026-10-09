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
    min_cluster_size = "5",
    clonal_cutOff = "0.9",
    propClonal_threshold = "0.25",
    burn_in = "1000",
    seed = "1024",
    subclonal_copy_correction = "TRUE",
    clean_clusters = "TRUE"
)
args_opt = parse_args("$task.ext.args")
for (ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]

# PyClone 0.13.1 uses @jit(cache=True). Numba then tries to write a cache next to
# site-packages, which is read-only in Singularity/Apptainer. Prefer a work-dir
# cache; disable JIT only if that import still fails.
numba_dir = file.path(getwd(), "numba_cache")
dir.create(numba_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(NUMBA_CACHE_DIR = normalizePath(numba_dir))
pyclone_import = system("python -c \"from pyclone.math_utils import log_beta_binomial_likelihood\"")
if (!identical(as.integer(pyclone_import), 0L)) {
    Sys.setenv(NUMBA_DISABLE_JIT = "1")
}

library(CONIPHER)

out_dir = "clustering"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

conipher_clustering(
    case_id = "$meta.patient",
    out_dir = out_dir,
    input_tsv_loc = "$input_tsv",
    min_cluster_size = as.integer(opt[["min_cluster_size"]]),
    clonal_cutOff = as.numeric(opt[["clonal_cutOff"]]),
    propClonal_threshold = as.numeric(opt[["propClonal_threshold"]]),
    burn_in = as.integer(opt[["burn_in"]]),
    seed = as.integer(opt[["seed"]]),
    nProcs = as.integer("$task.cpus"),
    subclonal_copy_correction = as_bool(opt[["subclonal_copy_correction"]]),
    clean_clusters = as_bool(opt[["clean_clusters"]])
)

outputs = list.files(out_dir, pattern = "SCoutput", full.names = TRUE)
outputs = outputs[grepl("tsv\$", outputs)]
if (length(outputs) == 0) {
    stop("CONIPHER clustering produced no SCoutput TSV files")
}
file.copy(outputs, ".", overwrite = TRUE)

f = file("versions.yml", "w")
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    CONIPHER:", as.character(packageVersion("CONIPHER"))), f)
close(f)
