#!/usr/bin/env Rscript

# Written by Harshil Patel and Gavin Kelly and released under the MIT license.

################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

## PCA, HEATMAP AND SCATTERPLOTS FOR SAMPLES IN COUNTS FILE
## - SAMPLE NAMES HAVE TO END IN e.g. "_R1" REPRESENTING REPLICATE ID. LAST 3 CHARACTERS OF SAMPLE NAME WILL BE TRIMMED TO OBTAIN GROUP ID FOR DESEQ2 COMPARISONS.
## - PACKAGES BELOW NEED TO BE AVAILABLE TO LOAD WHEN RUNNING R

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(
    make_option(c("-i", "--count_file"    ), type="character", default=NULL    , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."                        ),
    make_option(c("-f", "--count_col"     ), type="integer"  , default=3       , metavar="integer", help="First column containing sample count data."                                             ),
    make_option(c("-d", "--id_col"        ), type="integer"  , default=1       , metavar="integer", help="Column containing identifiers to be used."                                              ),
    make_option(c("-r", "--sample_suffix" ), type="character", default=''      , metavar="string" , help="Suffix to remove after sample name in columns e.g. '.rmDup.bam' if 'DRUG_R1.rmDup.bam'."),
    make_option(c("-m", "--metadata"      ), type="character", default=NULL    , metavar="path"   , help="Optional metadata file with sample information. Must contain a column matching sample names."),
    make_option(c("--pca_color"           ), type="character", default=NULL    , metavar="string" , help="Optional metadata column name to use for PCA plot colors."                             ),
    make_option(c("--pca_shape"           ), type="character", default=NULL    , metavar="string" , help="Optional metadata column name to use for PCA plot shapes."                             ),
    make_option(c("--deg_formula"         ), type="character", default=NULL    , metavar="string" , help="Optional formula string for differential expression analysis (e.g., '~condition' or '~condition + batch'). Must reference columns from metadata."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='deseq2', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-v", "--vst"           ), type="logical"  , default=FALSE   , metavar="boolean", help="Run vst transform instead of rlog."                                                     ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$count_file)){
    print_help(opt_parser)
    stop("Please provide a counts file.", call.=FALSE)
}

################################################
################################################
## READ IN COUNTS FILE                        ##
################################################
################################################

count.table           <- read.delim(file=opt$count_file,header=TRUE, row.names=NULL)
rownames(count.table) <- count.table[,opt$id_col]
count.table           <- count.table[,opt$count_col:ncol(count.table),drop=FALSE]
colnames(count.table) <- gsub(opt$sample_suffix,"",colnames(count.table))
colnames(count.table) <- gsub(pattern='\\.$', replacement='', colnames(count.table))

################################################
################################################
## RUN DESEQ2                                 ##
################################################
################################################

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir, recursive=TRUE)
}

# Store original working directory for metadata file path resolution
original_wd <- getwd()
setwd(opt$outdir)

samples.vec     <- colnames(count.table)
name_components <- strsplit(samples.vec, "_")
n_components    <- length(name_components[[1]])
decompose       <- n_components!=1 && all(sapply(name_components, length)==n_components)
coldata         <- data.frame(samples.vec, sample=samples.vec, row.names=1)
if (decompose) {
    groupings        <- as.data.frame(lapply(1:n_components, function(i) sapply(name_components, "[[", i)))
    n_distinct       <- sapply(groupings, function(grp) length(unique(grp)))
    groupings        <- groupings[n_distinct!=1 & n_distinct!=length(samples.vec)]
    if (ncol(groupings)!=0) {
        names(groupings) <- paste0("Group", 1:ncol(groupings))
        coldata <- cbind(coldata, groupings)
    } else {
        decompose <- FALSE
    }
}

################################################
################################################
## READ AND MERGE METADATA                    ##
################################################
################################################

if (!is.null(opt$metadata)) {
    # Check metadata file path (try relative to original working directory first, then relative to output directory)
    metadata_path <- opt$metadata
    if (!file.exists(metadata_path)) {
        metadata_path <- file.path(original_wd, opt$metadata)
    }
    if (!file.exists(metadata_path)) {
        stop(paste("Metadata file not found:", opt$metadata), call.=FALSE)
    }
    metadata <- read.delim(file=metadata_path, header=TRUE, row.names=NULL, check.names=FALSE)
    # Try to find sample column (could be "sample", "Sample", or first column)
    sample_col <- NULL
    if ("sample" %in% colnames(metadata)) {
        sample_col <- "sample"
    } else if ("Sample" %in% colnames(metadata)) {
        sample_col <- "Sample"
    } else {
        # Use first column as sample identifier
        sample_col <- colnames(metadata)[1]
    }

    # Merge metadata with coldata based on sample names
    if (sample_col %in% colnames(metadata)) {
        # Create a copy of metadata excluding the sample column to avoid duplication
        metadata_cols_to_add <- setdiff(colnames(metadata), sample_col)
        metadata_to_merge <- metadata[, c(sample_col, metadata_cols_to_add), drop=FALSE]

        # Add row names to coldata for matching
        coldata$rowname_backup <- rownames(coldata)
        metadata_to_merge$rowname_backup <- metadata_to_merge[,sample_col]

        # Merge, keeping all rows from coldata
        coldata <- merge(coldata, metadata_to_merge[,c("rowname_backup", metadata_cols_to_add), drop=FALSE],
                        by.x="rowname_backup", by.y="rowname_backup", all.x=TRUE, sort=FALSE)

        # Restore row names and remove backup column
        rownames(coldata) <- coldata$rowname_backup
        coldata$rowname_backup <- NULL

        # Convert character columns to factors (needed for DESeq2 formulas)
        if (!is.null(opt$deg_formula)) {
            for (col in colnames(coldata)) {
                if (col != "sample" && is.character(coldata[,col])) {
                    coldata[,col] <- as.factor(coldata[,col])
                }
            }
        }
    } else {
        warning(paste("Could not find sample column in metadata file. Expected 'sample' or 'Sample' or using first column."), call.=FALSE)
    }

    # Validate PCA color and shape columns if specified
    if (!is.null(opt$pca_color)) {
        if (!opt$pca_color %in% colnames(coldata)) {
            stop(paste("PCA color column '", opt$pca_color, "' not found in metadata. Available columns:", paste(colnames(coldata), collapse=", ")), call.=FALSE)
        }
    }
    if (!is.null(opt$pca_shape)) {
        if (!opt$pca_shape %in% colnames(coldata)) {
            stop(paste("PCA shape column '", opt$pca_shape, "' not found in metadata. Available columns:", paste(colnames(coldata), collapse=", ")), call.=FALSE)
        }
    }
}

DDSFile <- paste(opt$outprefix,".dds.RData",sep="")

counts  <- count.table[,samples.vec,drop=FALSE]

# Determine design formula
if (!is.null(opt$deg_formula)) {
    # Parse the formula string provided by user
    design_formula <- as.formula(opt$deg_formula)
    # Validate that all variables in the formula exist in coldata
    formula_vars <- all.vars(design_formula)
    missing_vars <- formula_vars[!formula_vars %in% colnames(coldata)]
    if (length(missing_vars) > 0) {
        stop(paste("Variables in DEG formula not found in metadata:", paste(missing_vars, collapse=", "),
                   ". Available columns:", paste(colnames(coldata), collapse=", ")), call.=FALSE)
    }
    dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design=design_formula)
} else {
    # `design=~1` creates intercept-only model, equivalent to setting `blind=TRUE` for transformation.
    dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design=~1)
}

# If DEG formula is provided, run full DESeq analysis; otherwise just estimate size factors
if (!is.null(opt$deg_formula)) {
    dds <- DESeq(dds)
} else {
    dds <- estimateSizeFactors(dds)
}
if (min(dim(count.table))<=1)  { # No point if only one sample, or one gene
    save(dds,file=DDSFile)
    saveRDS(dds, file=sub("\\.dds\\.RData$", ".rds", DDSFile))
    warning("Not enough samples or genes in counts file for PCA.", call.=FALSE)
    quit(save = "no", status = 0, runLast = FALSE)
}
if (!opt$vst) {
    vst_name <- "rlog"
    rld      <- rlog(dds)
} else {
    vst_name <- "vst"
    rld      <- varianceStabilizingTransformation(dds)
}

assay(dds, vst_name) <- assay(rld)
save(dds,file=DDSFile)
saveRDS(dds, file=sub("\\.dds\\.RData$", ".rds", DDSFile))

################################################
################################################
## EXTRACT DEG RESULTS                        ##
################################################
################################################

if (!is.null(opt$deg_formula)) {
    # Extract results (default comparison: last level vs first level)
    res <- results(dds)

    # Convert to data frame and add gene IDs
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)

    # Reorder columns to have gene first
    res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

    # Sort by adjusted p-value
    res_df <- res_df[order(res_df$padj, na.last=TRUE), ]

    # Write DEG table
    DEGFile <- paste(opt$outprefix, ".deg_results.txt", sep="")
    write.table(res_df, file=DEGFile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

    # Write summary statistics
    DEGSummaryFile <- paste(opt$outprefix, ".deg_summary.txt", sep="")
    sink(DEGSummaryFile)
    cat("DESeq2 Differential Expression Summary\n")
    cat("=======================================\n\n")
    cat("Formula:", opt$deg_formula, "\n\n")
    cat("Results summary:\n")
    print(summary(res))
    sink()
}

################################################
################################################
## PLOT QC                                    ##
################################################
################################################

##' PCA pre-processeor
##'
##' Generate all the necessary information to plot PCA from a DESeq2 object
##' in which an assay containing a variance-stabilised matrix of counts is
##' stored. Copied from DESeq2::plotPCA, but with additional ability to
##' say which assay to run the PCA on.
##'
##' @param object The DESeq2DataSet object.
##' @param ntop number of top genes to use for principla components, selected by highest row variance.
##' @param assay the name or index of the assay that stores the variance-stabilised data.
##' @return A data.frame containing the projected data alongside the grouping columns.
##' A 'percentVar' attribute is set which includes the percentage of variation each PC explains,
##' and additionally how much the variation within that PC is explained by the grouping variable.
##' @author Gavin Kelly
plotPCA_vst <- function (object,  ntop = 500, assay=length(assays(object))) {
    rv         <- rowVars(assay(object, assay))
    select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca        <- prcomp(t(assay(object, assay)[select, ]), center=TRUE, scale=FALSE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    df         <- cbind( as.data.frame(colData(object)), pca$x)
    #Order points so extreme samples are more likely to get label
    ord        <- order(abs(rank(df$PC1)-median(df$PC1)), abs(rank(df$PC2)-median(df$PC2)))
    df         <- df[ord,]
    attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar)
    return(df)
}

PlotFile <- paste(opt$outprefix,".plots.pdf",sep="")

pdf(file=PlotFile, onefile=TRUE, width=7, height=7)
## PCA
ntop <- c(500, Inf)
for (n_top_var in ntop) {
    pca.data      <- plotPCA_vst(dds, assay=vst_name, ntop=n_top_var)
    percentVar    <- round(attr(pca.data, "percentVar")$percentVar)
    plot_subtitle <- ifelse(n_top_var==Inf, "All genes", paste("Top", n_top_var, "genes"))

    # Build aesthetics mapping based on optional color and shape columns
    # Build base mapping
    base_mapping <- aes(x=PC1, y=PC2, label=paste0(" ", sample, " "))

    # Add color if specified
    if (!is.null(opt$pca_color) && opt$pca_color %in% colnames(pca.data)) {
        base_mapping$colour <- as.name(opt$pca_color)
    }
    # Add shape if specified
    if (!is.null(opt$pca_shape) && opt$pca_shape %in% colnames(pca.data)) {
        base_mapping$shape <- as.name(opt$pca_shape)
    }

    pl <- ggplot(pca.data, base_mapping) +
        geom_point() +
        geom_text(check_overlap=TRUE, vjust=0.5, hjust="inward", show.legend=FALSE) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        labs(title = paste0("First PCs on ", vst_name, "-transformed data"), subtitle = plot_subtitle) +
        theme(legend.position="top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
    print(pl)

    if (decompose) {
        pc_names <- paste0("PC", attr(pca.data, "percentVar")$PC)
        long_pc <- reshape(pca.data, varying=pc_names, direction="long", sep="", timevar="component", idvar="pcrow")
        long_pc <- subset(long_pc, component<=5)
        long_pc_grp <- reshape(long_pc, varying=names(groupings), direction="long", sep="", timevar="grouper")
        long_pc_grp <- subset(long_pc_grp, grouper<=5)
        long_pc_grp$component <- paste("PC", long_pc_grp$component)
        long_pc_grp$grouper <- paste0(long_pc_grp$grouper, c("st","nd","rd","th","th")[long_pc_grp$grouper], " prefix")
        pl <- ggplot(long_pc_grp, aes(x=Group, y=PC)) +
            geom_point() +
            stat_summary(fun=mean, geom="line", aes(group = 1)) +
            labs(x=NULL, y=NULL, subtitle = plot_subtitle, title="PCs split by sample-name prefixes") +
            facet_grid(component~grouper, scales="free_x") +
            scale_x_discrete(guide = guide_axis(n.dodge = 3))
        print(pl)
    }
} # at end of loop, we'll be using the user-defined ntop if any, else all genes

## WRITE PC1 vs PC2 VALUES TO FILE
pca.vals           <- pca.data[,c("PC1","PC2")]
colnames(pca.vals) <- paste0(colnames(pca.vals), ": ", percentVar[1:2], '% variance')
pca.vals           <- cbind(sample = rownames(pca.vals), pca.vals)
write.table(pca.vals, file = paste(opt$outprefix, ".pca.vals.txt", sep=""),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = TRUE)

## SAMPLE CORRELATION HEATMAP
sampleDists      <- dist(t(assay(dds, vst_name)))
sampleDistMatrix <- as.matrix(sampleDists)
colors           <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(
    sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors,
    main=paste("Euclidean distance between", vst_name, "of samples")
)

## WRITE SAMPLE DISTANCES TO FILE
write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),file=paste(opt$outprefix, ".sample.dists.txt", sep=""),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
dev.off()

################################################
################################################
## SAVE SIZE FACTORS                          ##
################################################
################################################

SizeFactorsDir <- "size_factors/"
if (file.exists(SizeFactorsDir) == FALSE) {
    dir.create(SizeFactorsDir, recursive=TRUE)
}

NormFactorsFile <- paste(SizeFactorsDir,opt$outprefix, ".size_factors.RData", sep="")

normFactors <- sizeFactors(dds)
save(normFactors, file=NormFactorsFile)

for (name in names(sizeFactors(dds))) {
    sizeFactorFile <- paste(SizeFactorsDir,name, ".txt", sep="")
    write(as.numeric(sizeFactors(dds)[name]), file=sizeFactorFile)
}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

RLogFile <- "R_sessionInfo.log"

sink(RLogFile)
a <- sessionInfo()
print(a)
sink()

################################################
################################################
################################################
################################################
