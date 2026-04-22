#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

# Read source and libraries
library(SAIGE)
library(data.table)
library(Matrix)
library(optparse)

option_list <- list(
    make_option(c("--rdaFile"), type="character", default="",
        help="Path to rda file from SAIGE step 1"),
    make_option(c("--chrom"), type="character", default="",
        help="Chromosome"),
    make_option(c("--geneName"), type="character", default="",
        help="Gene name to analyze"),
    make_option(c("--groupFile"), type="character", default="",
        help="Path to group file (containing functional annotation of variants)"),
    make_option(c("--traitType"), type="character", default="",
        help="Trait type (quantitative or binary)"),
    make_option(c("--bedFile"), type="character", default="",
        help="Path to bed file containing rare variants"),
    make_option(c("--bimFile"), type="character", default="",
        help="Path to bim file containing rare variants"),
    make_option(c("--famFile"), type="character", default="",
        help="Path to fam file containing rare variants"),
    make_option(c("--macThreshold"), type="integer", default=10,
        help="MAC threshold for ultra-rare variant collapsing (default: 10)"),
    make_option(c("--collapseGroups"), type="character", default="",
        help="Comma-separated list of group names to collapse into super-variants (e.g. 'lof,missense'). Default: none collapsed."),
    make_option(c("--apply_AR"), type="logical", default=FALSE,
        help="Apply adaptive ridge (approximated L0-regularization) when estimating effect size (default: FALSE)"),
    make_option(c("--outputPrefix"), type="character", default="",
        help="Path to save output (without extension and gene name)"),
    make_option(c("--groups"), type="character", default="lof,missense,synonymous",
        help="Comma-separated list of annotation group labels to analyze (default: 'lof,missense,synonymous')")
)

## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = FALSE)
print(args)

# Parse collapseGroups into a named logical vector
collapseGroups <- NULL
if (nchar(args$collapseGroups) > 0) {
    collapse_names <- trimws(strsplit(args$collapseGroups, ",")[[1]])
    collapseGroups <- setNames(rep(TRUE, length(collapse_names)), collapse_names)
}

run_RareEffect(
    rdaFile = args$rdaFile,
    chrom = args$chrom,
    geneName = args$geneName,
    groupFile = args$groupFile,
    traitType = args$traitType,
    bedFile = args$bedFile,
    bimFile = args$bimFile,
    famFile = args$famFile,
    macThreshold = args$macThreshold,
    collapseGroups = collapseGroups,
    apply_AR = args$apply_AR,
    outputPrefix = args$outputPrefix,
    groups = args$groups
)
