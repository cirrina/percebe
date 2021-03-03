#!/usr/bin/env Rscript

## ========================================== ## 
## LOAD PACKAGES                             ##
## ========================================== ## 
library(optparse)
#library(DESeq2)
#library(BiocParallel)
#library(ggplot2)
#library(RColorBrewer)
#library(pheatmap)


## ========================================== ## 
## PARSE COMMAND-LINE PARAMETERS              ##
## ========================================== ## 
# sample_sheet added 
# should be add gene/tx info reference file?

option_list <- list(
  make_option(c("-i", "--sample_sheet"  ), type="character", default=NULL   , metavar="path"   , help="Sample Sheet, colData, where samples are rows and columns are sample annotations."    ), 
  make_option(c("-b", "--iem_format"    ), type="logical"  , default=TRUE   , metavar="boolean", help="If Sample sheet is in illumina IEM format."                                            ),
  make_option(c("-s", "--sep"     ), type="character", default=','    , metavar="string" , help="Column separator, defaults to comma, ','"                                                     ),
  make_option(c("-a", "--require_adapter"), type="logical"  , default=TRUE  , metavar="boolean", help="If to check for Illumina adapter. Will only check if there is a character vector here or not"),
  make_option(c("-d", "--flag_dup_Names"), type="logical"  , default=TRUE  , metavar="boolean", help="If to check for duplicate sample names."),
  make_option(c("-h", "--header"    ), type="logical"  , default=TRUE   , metavar="boolean", help="If Sample sheet data section includes a header"                                            ),
  make_option(c("-n", "--header_names"    ), type="character"  , default=c("Sample_ID","Sample_Name","Sample_Project"), metavar="string", help="string of column header names required for data section" )
  
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt        <- parse_args(opt_parser)

if (is.null(opt$sample_sheet)){
  print_help(opt_parser)
  stop("Please provide a sample_sheet file.", call.=FALSE)
}

cat("\n ... CHECKING SAMPLE SHEET\n\n")

## ========================================== ## 
## PROCESS COUNTS FILE  - IEM STYLE
## ========================================== ## 
## comment table is read and hashes are commented, hashes used when featureCounts produce output
if(opt$iem_format){
  #ssheet <- read.delim()
  # opt$sample_sheet <- "/Users/david/CTG_projects/bulkRNA/CTG_2020_177/samplesheet_2020_177.csv"
  all.lines <- scan(file = opt$sample_sheet, what = "character", sep="\n", nmax = 250)
  
  ## check nlines of document
  if(length(all.lines)>=250) stop("Number of rows in sample sheet exceeds what is read by the scan function ('nmax' set to 250)")
  
  ## check file delimiter, column separators
  # u <- grep(opt$sep, all.lines[1])
  # if(!length(u)) stop("Cant find separators in string, check if correct file delimiter format is used using --sep")
  
  ## check IEM sections
  iem.headers <- list(Header="[[]Header[]]", Reads="[[]Reads[]]",Settings="[[]Settings[]]",Data="[[]Data[]]")
  iem.index <- sapply(iem.headers, function(x) grep(x, all.lines))
  if(!any(unlist(lapply(iem.index, length)==1))){
    stop("Sample Sheet Not IEM format - check Illumina IEM sections, must include - [Header], [Reads], [Settings], [Data] ...\n ... or set --iem_format to FALSE")
  }
  iem.index <- sapply(iem.headers, function(x) grep(x, all.lines))
  
  ## Check Adapter if adapter sequence slot is NA or not
  # Note, will not check if different Adapter reads for forward or reverse
  if(opt$require_adapter){
    settings.section <- read.delim(file=opt$sample_sheet, header = F, sep = opt$sep, 
                                 skip = iem.index["Settings"], nrows = iem.index["Data"]-iem.index["Settings"]-1)
    u <- grep("Adapter", settings.section[,1])
    if(!length(u)) stop("Cant find 'Adapter' row in [Settings] section ")
    if(is.na(settings.section[u[1],1])) stop("No Adapter sequence detected. Define adapter or set --require_adapter FALSE")
  }
  
  # load data sheet 
  data.section <- read.delim(file=opt$sample_sheet, header = opt$header, sep = opt$sep, 
                                 skip = iem.index["Data"])
} # end if iem sheet format

# if not iem style. load delimited file
if(!opt$iem_format){
  data.section <- read.delim(file=opt$sample_sheet, header = opt$header, sep = opt$sep)
} # end if not iem format

# Check Data section. 
## 1. Columns must include Sample_ID,Sample_Name,Sample_Project columns
# opt$header <- T
  required.columns <- c("Sample_ID","Sample_Name","Sample_Project")
  if(!all (opt$header_names %in% colnames(data.section))){
    stop("not all required header names are present. Check --header_names or set --header to FALSE")
  }
  
  ## 2. Check for non allowed characters * \ / . + - ' ' ( ) " '
  specials <- c("[*]", "[/]", "[.]","[+]", "[-]", "[(]", "[)]", ",", " ")
  # matrix(data.section)
  special.flag <- F
  for(i in 1:length(specials)){
    for(j in 1:ncol(data.section)){
      u<-grepl(specials[i], data.section[,j]) 
      
      if(any(unlist(u))){
        cat("\n ... Column  '", colnames(data.section)[j], "'  contains: ", ifelse(specials[i]==" ", "white space", specials[i]))  
        special.flag <- T
      } 
    }
  } 
cat("\n ...")

## Check if Sampl_Name duplicates
if(opt$header && opt$flag_dup_Names){
  u <- duplicated(data.section$Sample_Name)
  if(any(u)) stop("Duplicated Sample_Name present: ", data.section$Sample_Name[u])
}


if(special.flag) stop("Sample sheet contains special characters. Fix this !!")



cat("\n ... ### --------------- ### " )
cat("\n ... ### Sample Sheet OK ### " )
cat("\n ... ### --------------- ### \n\n" )

  

