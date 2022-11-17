#!/usr/bin/env Rscript

library(jsonlite)
library(sitePath)
library(optparse)

parser <- OptionParser()
parser <- add_option(
    parser,
    "--treeFile",
    help = "The input tree file"
)
parser <- add_option(
    parser,
    "--msaFile",
    help = "The input multiple sequence alignment file"
)
parser <- add_option(
    parser,
    "--Nmin",
    type = "double",
    default = NULL,
    help = paste(
        "The parameter to resolve phylogenetic lineages",
        "(default is recommended)"
    ),
    metavar = "number"
)
parser <- add_option(
    parser,
    "--outFile",
    default = "fixed_and_parallel.json",
    help = paste(
        "Path of the output file",
        "(default is 'fixed_and_parallel.json')"
    )
)
args <- parse_args(parser)

tree_file <- args$treeFile
msa_file <- args$msaFile
n_min <- args$Nmin
out_file <- args$outFile

if (is.null(tree_file)) {
    stop("treeFile is missing")
}
if (is.null(msa_file)) {
    stop("msaFile is missing")
}

tree <- read.tree(file = tree_file)
paths <- addMSA(tree = tree, msaPath = msa_file)
if (!is.null(n_min)) {
    paths <- lineagePath(tree = paths, similarity = n_min)
}
min_entropy <- sitesMinEntropy(x = paths)

fixed_sites <- fixationSites(paths = min_entropy)
para_sites <- parallelSites(x = min_entropy)

fixed_and_parallel <- list(
    "fixed" = allSitesName(fixed_sites),
    "parallel" = allSitesName(para_sites)
)

write_json(fixed_and_parallel, out_file)
