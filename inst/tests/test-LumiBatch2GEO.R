context("testing LumiBatch2GEOarchive")

test_that(".LumiBatch2GEOarchive_df test data exists", {
	f <- "sample_GenomeStudio_output_9x3.txt"
	expect_that(file.exists(f), is_true())
})

test_that(".LumiBatch2GEOarchive_df test data is valid Lumi data", {
	f <- "sample_GenomeStudio_output_9x3.txt"
	raw <- lumiR(f, verbose=FALSE)
	expect_that(is(raw, "LumiBatch"), is_true())
})

test_that(".LumiBatch2GEOarchive_df works", {
	f <- "sample_GenomeStudio_output_9x3.txt"
	raw <- lumiR(f, verbose=FALSE)
	res <- microarrays:::.LumiBatch2GEOarchive_df(raw)
	expect_is(res, "data.frame")
	
	expected.res <- structure(list(ID_REF = c("ILMN_1762337", "ILMN_2055271", "ILMN_1736007", 
	"ILMN_2383229", "ILMN_1806310", "ILMN_1779670", "ILMN_1653355", 
	"ILMN_1717783", "ILMN_1705025"), `5445316009_D` = c(136.0626, 
	193.0642, 134.0098, 380.6206, 642.526, 118.9719, 138.6569, 96.33086, 
	155.3983), `Detection Pval` = c(0.45584, 0.00519, 0.48831, 0, 
	0, 0.71948, 0.38961, 0.92597, 0.14545), `5445316009_E` = c(172.9379, 
	149.0012, 136.9481, 287.6883, 389.5496, 104.1193, 165.0129, 100.7001, 
	153.3313), `Detection Pval` = c(0.03636, 0.28961, 0.5039, 0, 
	0, 0.88182, 0.07662, 0.91558, 0.19481), `5445316009_G` = c(212.3289, 
	177.6925, 139.0936, 325.9922, 408.0934, 102.0919, 152.1122, 113.2898, 
	146.5808), `Detection Pval` = c(0.00519, 0.02338, 0.32857, 0, 
	0, 0.87403, 0.15195, 0.74805, 0.21558)), .Names = c("ID_REF", 
	"5445316009_D", "Detection Pval", "5445316009_E", "Detection Pval", 
	"5445316009_G", "Detection Pval"), row.names = c(NA, -9L), class = "data.frame")
	
	expect_equivalent(res, expected.res)
	
})

test_that(".LumiBatch2GEOarchive_df,ids works", {
	f <- "sample_GenomeStudio_output_9x3.txt"
	x <- lumiR(f, verbose=FALSE)
	ids <- c("ILMN_1653355", "ILMN_1705025", "ILMN_1717783", "ILMN_1736007", 
	"ILMN_1762337", "ILMN_1779670", "ILMN_1806310", "ILMN_2055271", 
	"ILMN_2383229")
	res <- microarrays:::.LumiBatch2GEOarchive_df(x, ids=ids)
	expect_equivalent(res$ID_REF, ids)
})

test_that(".LumiBatch2GEOarchive_df,ids,subset works", {
	f <- "sample_GenomeStudio_output_9x3.txt"
	x <- lumiR(f, verbose=FALSE)
	ids <- c("ILMN_1653355", "ILMN_1705025", "ILMN_1717783", 
	"ILMN_2383229")
	res <- microarrays:::.LumiBatch2GEOarchive_df(x, ids=ids)
	expect_equivalent(res$ID_REF, ids)
})

test_that(".LumiBatch2GEOarchive_df,ids,extra works", {

	f <- "sample_GenomeStudio_output_9x3.txt"
	x <- lumiR(f, verbose=FALSE)
	ids <- c("ILMN_1653355", "ILMN_1705025", "ILMN_1717783", 
	"ILMN_2383229", "ILMN_FAKE")
	res <- microarrays:::.LumiBatch2GEOarchive_df(x, ids=ids)
	expect_equivalent(res$ID_REF, ids)
	expect_equal(sum(is.na(res[5,])), 6)
})

test_that(".LumiBatch2GEOarchive_df,ids,1NA works", {
	f <- "sample_GenomeStudio_output_9x3.txt"
	x <- lumiR(f, verbose=FALSE)
	ids <- c("ILMN_1653355", "ILMN_1705025", "ILMN_1717783", 
	"ILMN_2383229", NA)
	res <- microarrays:::.LumiBatch2GEOarchive_df(x, ids=ids)
	expect_equivalent(res$ID_REF, ids)
	expect_equal(sum(is.na(res[5,])), 7)
})

test_that(".LumiBatch2GEOarchive_df,ids,with dup works", {

	f <- "sample_GenomeStudio_output_9x3.txt"
	x <- lumiR(f, verbose=FALSE)
	ids <- c("ILMN_1653355", "ILMN_1653355", "ILMN_2383229")
	expect_error(res <- microarrays:::.LumiBatch2GEOarchive_df(x, ids=ids)
				, "ids can't contain duplicates")
	
})

test_that(".LumiBatch2GEOarchive_df,ids,2NA works", {

	f <- "sample_GenomeStudio_output_9x3.txt"
	x <- lumiR(f, verbose=FALSE)
	ids <- c(NA, "ILMN_1653355", "ILMN_1705025", "ILMN_1717783", 
	"ILMN_2383229", NA)
	expect_error(res <- microarrays:::.LumiBatch2GEOarchive_df(x, ids=ids)
				, "ids can't contain duplicates")
})
