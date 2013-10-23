#!/usr/local/R.2.9.2/bin/Rscript --no-save

source("funcs.R")
args=commandArgs(TRUE)

usage="\n  Usage: ./Rae.R
    --genome=[Required: hg18|mm9]
    --project=[Required: unique project name]
    --data=[Required: absolute path of data directory]
    --qv=[Optional (default=0.1): Q-value threshold for significant alterations]
    --exclude=[Optional (default=FALSE): TRUE|FALSE, to exclude samples by pre-determined QC]
    --SNP6=[Optional (default=FALSE): TRUE|FALSE, enable special SNP6 array handling]\n
		--normOnly=[Optional: runs only the first normalization/parameterization steps, skips latter steps of pipeline]\n
"

# Initialize default parameters
p=list()
p$Qval=0.1
p$Exclude=FALSE
p$SNP6=FALSE
p$Dir=getwd()
p$NormOnly=FALSE

# Verify that at least required parameters and their values are passed
if(length(args)<3) {
	cat(usage)
	stop("Incorrect or missing required input!")
}

# Verify genome parameter
idx=grep("--genome=",args)
if(is.integer(idx)) {
	genome=gsub("--genome=","",args[idx])
	if(!genome%in%c("hg18","mm9")) {
		stop("Invalid/unsupported genome [options: hg18, mm9]")
	} else {
		p$Build=genome
	}
} else {
	stop("Missing required --genome parameter!")
}

# Verify project parameter
idx=grep("--project=",args)
if(is.integer(idx)) {
	project=gsub("--project=","",args[idx])
	if(length(grep("[!@#$%^&*|\"\':;/]",project))>0) {
		stop("Invalid --project parameter: no special characters allowed!")
	} else {
		p$Project=project
	}
} else {
	stop("Missing required --project parameter!")
}

# Verify data directory
idx=grep("--data=",args)
if(is.integer(idx)) {
	p$Data=gsub("--data=","",args[idx])
	if(nchar(p$Data)==0) {
		stop("Missing value for --data parameter!")
	}
	if(!isDirectory(p$Data)) {
		stop("Directory supplied to --data does not exist!")
	}
	if(isDirectory(p$Data) & length(list.files(p$Data))==0) {
		stop("Directory supplied to --data is empty!")
	}
} else {
	stop("Missing required --data parameter!")
}

# Custom FDR provided, reset from default
if(sum(grepl("--qv=",args))>0) p$Qval=as.numeric(gsub("--qv=","",args[grep("--qv=",args)]))

# Provided value for exclude, reset from default
if(sum(grepl("--exclude=",args))>0) p$Exclude=as.logical(gsub("--exclude=","",args[grep("--exclude=",args)]))

# Provided value for SNP6, reset from default
if(sum(grepl("--SNP6=",args))>0) p$SNP6=as.logical(gsub("--SNP6=","",args[grep("--SNP6=",args)]))

if(sum(grepl("--normOnly",args))>0) p$NormOnly=TRUE

# Load necessary packages
options(warn=-1)
required="\n  Required packages:
    DNAcopy
    Biostrings
    BSgenome
    BSgenome.Hsapiens.UCSC.hg18
    BSgenome.Mmusculus.UCSC.mm9
    IRanges\n
"
if(!suppressMessages(library(DNAcopy,logical.return=TRUE)) |
	 !suppressMessages(library(Biostrings,logical.return=TRUE)) |
	 !suppressMessages(library(BSgenome,logical.return=TRUE))) {
		cat(required)
		stop("Required packages are not installed!")
}

if(p$Build=="hg18") {
	if(!suppressMessages(library(BSgenome.Hsapiens.UCSC.hg18,logical.return=TRUE))) {
			cat(required)
			stop("Required packages are not installed!")
	}
	Genome=Hsapiens
}

if(p$Build=="mm9") {
	if(!suppressMessages(library(BSgenome.Mmusculus.UCSC.mm9,logical.return=TRUE))) {
			cat(required)
			stop("Required packages are not installed!")
	}
	Genome=Mmusculus
}
options(warn=0)

BETA.tol=0.001
N.density=2048
MAXDBL=.Machine$double.xmax
MINDBL=5e-324
started=proc.time()

cat("\n")
cat("--------------------------------------------------------------\n")
cat("Starting RAE analysis...\n")
cat("--------------------------------------------------------------\n")

if(grepl("hg",p$Build)) {
	cat(">> Step 0 of 5: Get genome features...\n")
	feat=read.table("FEAT.file",header=T,skip=1,as.is=T,colClasses=c(rep("character",2),rep("numeric",3)))
	cnv=read.table("CNV.file",header=T,skip=1,as.is=T)
	cnv=cnv[order(cnv$Chr,cnv$Start),]
}

cat(">> Step 1 of 5: Gathering project segmentation...\n")
dd=gatherProjectSegmentation(p)

cat(">> Step 2 of 5: Parameterizing discriminators per sample...\n")
pp=parameterizeMultiComponentModel(p,dd)

#save(p,dd,pp,
#     file=paste("CHECKPOINT",p$Project,gsub(" ","_",date()),".Rdata",sep="__"),
#     compress=T)

if(p$NormOnly==FALSE) {
	cat(">> Step 3 of 5: Scoring and assessing the recurrence of alterations...\n")
	dd=scoreProjectAnalysis(p,dd,pp)

	cat(">> Step 4 of 5: Selecting regions of interest...\n")
	try({getRegionsOfInterest(p,dd,pp,errorFactor=0.25)})

	cat(">> Step 5 of 5: Plotting output...\n")
	try({plotGenome(p,dd,save="pdf")})

	if(grepl("hg",p$Build)) {
		cat(">> Writing gene/tumor map...\n")
		writeGeneMap(p)
	}
} else {
	cat(">> Running normalization/parameterization-only mode, skip steps 3-5...\n")
}

ended=proc.time()

cat("--------------------------------------------------------------\n")
cat("Finished...\n")
cat("Elapsed (seconds): ",as.numeric((ended-started)[3]),"\n",sep="")
cat("Elapsed (hrs): ",as.numeric((ended-started)[3])/3600,"\n",sep="")
cat("--------------------------------------------------------------\n")

save.image(file="raeImage.Rdata")

