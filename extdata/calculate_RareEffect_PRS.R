#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

# Calculate RareEffect PRS for each individual

## Load packages and sources

library(SAIGE)
library(data.table)
library(Matrix)
library(optparse)

option_list <- list(
    make_option(c("--effectFile"), type="character", default="",
        help="Path to effect file from RareEffect main function"),
    make_option(c("--bedFile"), type="character", default="",
        help="Path to bed file containing rare variants"),
    make_option(c("--bimFile"), type="character", default="",
        help="Path to bim file containing rare variants"),
    make_option(c("--famFile"), type="character", default="",
        help="Path to fam file containing rare variants"),
    make_option(c("--groupFile"), type="character", default="",
        help="Path to group file (containing functional annotation of variants)"),
    make_option(c("--geneName"), type="character", default="",
        help="Gene name to analyze"),
    make_option(c("--variantListFile"), type="character", default="",
        help="Path to file containing list of variants to exclude (non-ultra-rare variants)"),
    make_option(c("--outputPrefix"), type="character", default="",
        help="Path to save output (without extension and gene name)"),
    make_option(c("--groups"), type="character", default="lof,missense,synonymous",
        help="Comma-separated list of annotation group labels to analyze (default: 'lof,missense,synonymous')")
)

## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = FALSE)
print(args)

## Parse groups
group_names <- trimws(strsplit(args$groups, ",")[[1]])
n_groups <- length(group_names)

## Assign varibles
bedFile <- args$bedFile
bimFile <- args$bimFile
famFile <- args$famFile
effectFile <- args$effectFile
groupFile <- args$groupFile
geneName <- args$geneName

## Set PLINK object
bim <- fread(bimFile)
fam <- fread(famFile)
sampleID <- as.character(fam$V2)
modglmm <- list()
modglmm$sampleID <- sampleID
n_samples <- length(sampleID)
var_by_func_anno <- read_groupfile(groupFile, geneName, group_names)

objGeno <- SAIGE::setGenoInput(bgenFile = "",
        bgenFileIndex = "",
        vcfFile = "",
        vcfFileIndex = "",
        vcfField = "",
        savFile = "",
        savFileIndex = "",
        sampleFile = "",
        bedFile=bedFile,
        bimFile=bimFile,
        famFile=famFile,
        idstoIncludeFile = "",
        rangestoIncludeFile = "",
        chrom = "",
        AlleleOrder = "alt-first",
        sampleInModel = sampleID
    )

effect <- fread(effectFile)

# Print which groups have no variants
for (g in group_names) {
    if (length(var_by_func_anno[[g]]) == 0) {
        print(paste0(g, " variant does not exist."))
    }
}

# Remove variants not in plink file
for (g in group_names) {
    var_by_func_anno[[g]] <- var_by_func_anno[[g]][which(var_by_func_anno[[g]] %in% objGeno$markerInfo$ID)]
}
print(str(var_by_func_anno))

# Read genotype matrix for each group
mat_collapsed <- list()
mat_flipped <- list()
for (g in group_names) {
    tmp <- collapse_matrix(objGeno, var_by_func_anno[[g]], sampleID, modglmm, macThreshold = 0)
    mat_collapsed[[g]] <- tmp[[1]]
    mat_flipped[[g]] <- tmp[[2]]
}
print("Genotype matrix is loaded.")

# Remove common variants
if (args$variantListFile == "") {
    print("No variant list file is provided.")
} else {
    var_list <- fread(args$variantListFile)
    common_var <- var_list$V1
    for (g in group_names) {
        exclude_idx <- which(colnames(mat_collapsed[[g]]) %in% common_var)
        if (length(exclude_idx) > 0) {
            mat_collapsed[[g]] <- mat_collapsed[[g]][, -exclude_idx, drop = FALSE]
        }
    }
    print("Common variants are removed.")
}

# Read effect size
effect <- as.data.frame(effect)

# Build effect vectors for each group
effect_per_group <- list()
for (g in group_names) {
    ur_label <- paste0(g, "_UR")
    if (length(which(effect$variant == ur_label)) == 0) {
        effect_UR <- 0
    } else {
        effect_UR <- effect[which(effect$variant == ur_label), ]$effect
    }

    n_col <- ncol(mat_collapsed[[g]])
    effect_g <- rep(effect_UR, n_col)

    if (n_col > 1) {
        for (i in 1:(n_col - 1)) {
            match_idx <- which(effect$variant == colnames(mat_collapsed[[g]])[i])
            if (length(match_idx) != 0) {
                effect_g[i] <- effect[match_idx, ]$effect
            }
        }
    }

    effect_per_group[[g]] <- effect_g
}

print("Effect size is loaded.")

## Calculate PRS
total_prs <- Matrix::Matrix(0, nrow = n_samples, ncol = 1, sparse = TRUE)
rownames(total_prs) <- sampleID

for (g in group_names) {
    if (ncol(mat_collapsed[[g]]) > 0 && length(effect_per_group[[g]]) > 0) {
        g_prs <- mat_collapsed[[g]] %*% effect_per_group[[g]]
        # Align rows by sample ID
        total_prs[match(rownames(g_prs), sampleID), 1] <- total_prs[match(rownames(g_prs), sampleID), 1] + g_prs[,1]
    }
}

out <- data.frame(IID = sampleID, PRS = as.numeric(total_prs[,1]))

## Save PRS
write.table(out, paste0(args$outputPrefix, "_", geneName, "_score.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
write.table(out, "test.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
