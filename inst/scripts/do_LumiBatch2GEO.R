################################################################################
# HOWTO create an Illumina GEO submission.
# In this example I have raw data in the form of GenomeStudio-like output, and
# a tsv, exported from R, after lumiN.
################################################################################
library(lumi)
library(microarrays)

raw.f <- "GEX_for_methyaltion_108_raw.txt" # GenomeStudio-like output
norm.f <- "GEX_methylation_108_Normalized_q_Back_affy_log2.txt" # this is a tsv, exported via lumiN
options(stringsAsFactors=F)

x.raw <- lumiR(raw.f)
LumiBatch2GEOarchive(x.raw, "raw_data.tsv")

x.norm.df <- read.delim(norm.f, check.names=F)

# print.venn(sampleNames(x.raw), colnames(x.norm.df))
# sampleNames(x.raw) == colnames(x.norm.df)
x.norm <- new("LumiBatch", exprs=as.matrix(x.norm.df), se.exprs=se.exprs(x.raw), detection=detection(x.raw))

all(featureNames(x.norm) == featureNames(x.raw))
all(sampleNames(x.norm) == sampleNames(x.raw))

#
# export
# 
LumiBatch2GEOarchive(x.raw, "Matrix non-normalized.tsv")
LumiBatch2GEOarchive(x.norm, "Matrix normalized.tsv")

#
# merge
# 
$ which tab2xls
/Library/Frameworks/R.framework/Resources/library/excelIO/bin/tab2xls
$ tab2xls "Matrix normalized.tsv" "Matrix non-normalized.tsv" GEOsub.xls

#
# insert into a GA_illumina_expression.xls sheet downloaded from GEO
#

################################################################################
################################################################################
