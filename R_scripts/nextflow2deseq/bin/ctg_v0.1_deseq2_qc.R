#!/usr/bin/env Rscript

# Main R script for DESeq2 analysis in bulkRNA seq analyses 
# based on nfcore pipeline 3.0 deseq2_qc.r
#

# added functionality 

# Core R functions are defied in source script. ctg_r_functions.R


################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

# NFCORE 
##  PCA, HEATMAP AND SCATTERPLOTS FOR SAMPLES IN COUNTS FILE
## - SAMPLE NAMES HAVE TO END IN e.g. "_R1" REPRESENTING REPLICATE ID. LAST 3 CHARACTERS OF SAMPLE NAME WILL BE TRIMMED TO OBTAIN GROUP ID FOR DESEQ2 COMPARISONS.
## - PACKAGES BELOW NEED TO BE AVAILABLE TO LOAD WHEN RUNNING R


## ========================================== ## 
## NF CORE vs CTG
## ========================================== ## 
# Diferences bethee

## ========================================== ## 
## LOAD PACKAGES                             ##
## ========================================== ## 
library(optparse)
library(DESeq2)
library(BiocParallel)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)


## ========================================== ## 
## PARSE COMMAND-LINE PARAMETERS              ##
## ========================================== ## 
# sample_sheet added 
# should be add gene/tx info reference file?

option_list <- list(
  make_option(c("-i", "--count_file"    ), type="character", default=NULL    , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."                        ),
  make_option(c("-i", "--sample_sheet"  ), type="character", default=NULL    , metavar="path"   , help="Sample Sheet, colData, where samples are rows and columns are sample annotations."                        ), ## added 
  # make_option(c("-i", "--gene_info"     ), type="character", default=NULL    , metavar="path"   , help="Transcript info file where genes/transcripts are rows and columns are sample annotations."        ), ## added 
  
  make_option(c("-f", "--count_col"     ), type="integer"  , default=2       , metavar="integer", help="First column containing sample count data."                                             ),
  make_option(c("-d", "--id_col"        ), type="integer"  , default=1       , metavar="integer", help="Column containing identifiers to be used."                                              ),
  
  make_option(c("-r", "--sample_suffix" ), type="character", default=''      , metavar="string" , help="Suffix to remove after sample name in columns e.g. '.rmDup.bam' if 'DRUG_R1.rmDup.bam'."),
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

if (is.null(opt$count_file)){
  print_help(opt_parser)
  stop("Please provide a counts file.", call.=FALSE)
}


## convert featureCounts output to nfore style count matrix file

## Mach up counts_file with 



## ========================================== ## 
## PROCESS COUNTS FILE & SAMPLE SHEET         ##
## ========================================== ## 

count.table           <- read.delim(file=opt$count_file,header=TRUE)
rownames(count.table) <- count.table[,opt$id_col]
count.table           <- count.table[,opt$count_col:ncol(count.table),drop=FALSE]
colnames(count.table) <- gsub(opt$sample_suffix,"",colnames(count.table))
colnames(count.table) <- gsub(pattern='\\.$', replacement='', colnames(count.table))



## ========================================== ## 
## RUN DESEQ2                                 ##
## ========================================== ## 

if (file.exists(opt$outdir) == FALSE) {
  dir.create(opt$outdir,recursive=TRUE)
}
setwd(opt$outdir)

samples.vec <- sort(colnames(count.table))
groups      <- sub("_[^_]+$", "", samples.vec)
if (length(unique(groups)) == 1 || length(unique(groups)) == length(samples.vec)) {
  quit(save = "no", status = 0, runLast = FALSE)
}

DDSFile <- paste(opt$outprefix,".dds.RData",sep="")
if (file.exists(DDSFile) == FALSE) {
  counts  <- count.table[,samples.vec,drop=FALSE]
  coldata <- data.frame(row.names=colnames(counts), condition=groups)
  dds     <- DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design=~ condition)
  dds     <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(opt$cores))
  if (!opt$vst) {
    vst_name <- "rlog"
    rld      <- rlog(dds)
  } else {
    vst_name <- "vst"
    rld      <- varianceStabilizingTransformation(dds)
  }
  assay(dds, vst_name) <- assay(rld)
  save(dds,file=DDSFile)
} else {
  load(DDSFile)
  vst_name <- intersect(assayNames(dds), c("vst", "rlog"))
  if (length(vst_name)==0) { # legacy might mean vst was saved as a separate object called rld
    vst_name <- "loaded_rld"
    assay(dds, vst_name) <- assay(rld)
  } else {
    vst_name==vst_name[1]
  }
}
