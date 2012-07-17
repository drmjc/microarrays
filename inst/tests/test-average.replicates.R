context("testing average.replicates")

################################################################################
# #
# # create a small LumiBatch suitable for testing average.replicates
# #
# load("/Volumes/GRIW/ICGCPancreas/Mark/projects/2012-05-29 Combine iDAT/Rmisc/x.norm.RDa.gz")
# # ids <- sample(featureNames(x.norm),40)
# ids <- c("ILMN_3237291", "ILMN_3214052", "ILMN_1821724", "ILMN_1891277", 
# "ILMN_1757137", "ILMN_3210171", "ILMN_1774937", "ILMN_1727880", 
# "ILMN_2326376", "ILMN_1741977", "ILMN_3243953", "ILMN_1737760", 
# "ILMN_3232701", "ILMN_1748836", "ILMN_1736056", "ILMN_1721087", 
# "ILMN_3242786", "ILMN_1681979", "ILMN_2259223", "ILMN_1691467", 
# "ILMN_3182275", "ILMN_1794818", "ILMN_1732782", "ILMN_3294263", 
# "ILMN_1746720", "ILMN_1760982", "ILMN_2262462", "ILMN_1706094", 
# "ILMN_1799137", "ILMN_2282477", "ILMN_2132309", "ILMN_1719340", 
# "ILMN_2402972", "ILMN_3225649", "ILMN_1736448", "ILMN_1651949", 
# "ILMN_2072178", "ILMN_1790688", "ILMN_1815083", "ILMN_1665873"
# )
# 
# x.LumiBatch <- x.norm
# x.LumiBatch <- x.LumiBatch[ids,]
# x.LumiBatch@controlData <- x.LumiBatch@controlData[1:35,]
# x.LumiBatch <- x.LumiBatch[,1:7]
# # bug fix in the subsetting of the controlData slot.
# x.LumiBatch@controlData <- data.frame(x.norm@controlData[1:35,1:2], x.LumiBatch@controlData)
# sampleNames(x.LumiBatch) <- c("normal.1", "normal.2", "normal.3", "mutant.1", "mutant.2", "mutant.3", "reference")
# sampleNames( x.LumiBatch )
# save(x.LumiBatch, file="~/src/R/pwbc/inst/examples/x.LumiBatch.example.RDa")
# 
# x.ExpressionSet <- as(x.LumiBatch, "ExpressionSet")
# save(x.ExpressionSet, file="~/src/R/pwbc/inst/examples/x.ExpressionSet.example.RDa")
#
# # save the expected result too.
# classes <- sub("\\.[1-9]$", "", sampleNames(x.LumiBatch)) # character
# x.LumiBatch.avg <- average.replicates(x.LumiBatch, classes)
# save(x.LumiBatch.avg, file="~/src/R/pwbc/inst/examples/x.LumiBatch.avg.RDa")
# save(x.LumiBatch.avg, file="~/src/R/quicktest/inst/examples/x.LumiBatch.avg.RDa")
# 
# x.ExpressionSet.avg <- average.replicates(x.ExpressionSet, classes)
# save(x.ExpressionSet.avg, file="~/src/R/pwbc/inst/examples/x.ExpressionSet.avg.RDa")
#
################################################################################

test_that("test that input data exists", {
	expect_true( file.exists(system.file("examples", "x.LumiBatch.example.RDa", package="pwbc")) )
	expect_true( file.exists(system.file("examples", "x.LumiBatch.avg.RDa", package="pwbc")) )
})

test_that("average.replicates on matrix", {
	require(lumi) || stop("required package 'lumi' is not installed")
	load( system.file("examples", "x.LumiBatch.example.RDa", package="pwbc") ) # load x.LumiBatch object
	classes <- as.factor(sub("\\.[1-9]$", "", sampleNames(x.LumiBatch))) # character
	res <- average.replicates(exprs(x.LumiBatch), classes)
	expected.res <- structure(c(8.15150004530165, 10.1328745467003, 7.49775545803208, 
	7.4476481696329, 7.04401424384069, 11.2706112785709, 7.2284567269886, 
	7.38306934877342, 6.69251279620315, 6.99606650081166, 6.90200128207641, 
	6.84088257649328, 7.22454179317623, 8.13375079270194, 7.07309579565193, 
	7.78473033452568, 6.68312521707511, 7.08583171554452, 7.17478966392928, 
	6.87224219059137, 7.55816547408239, 5.87133975309216, 6.38738173698081, 
	7.01746351320833, 8.13554336425842, 8.7633395728979, 8.59823938194037, 
	11.3141877007204, 7.01052330080514, 7.24637689363769, 6.96449355710694, 
	6.57244311986803, 8.12546021072627, 7.21166738792694, 6.78906158755257, 
	8.05447404273069, 7.18752414981699, 7.16502260541611, 8.90414783930211, 
	6.78241837686025, 7.97501862079713, 9.58381220345539, 7.73357125947685, 
	7.20230792005529, 6.96649702440178, 11.2626563397547, 7.02700789930857, 
	6.55471603165014, 7.06798699771778, 7.03565621421496, 6.79018594678216, 
	6.64928854710145, 7.59169971170484, 8.04900209030033, 7.18510803926087, 
	8.49105372559396, 6.81295287383096, 6.91728578158643, 7.07703512676991, 
	6.53758355899177, 6.94811938176567, 6.20150859798392, 6.69471729652901, 
	6.82188137672293, 9.53758666042858, 8.92632167154977, 8.13770103101301, 
	11.0424415047994, 6.82414075153545, 7.08124634263355, 6.89704497671146, 
	6.72792696685374, 7.84103761036026, 7.57872959530072, 7.11182560014886, 
	7.9156382547289, 9.25011075284244, 7.00720109610439, 8.45333471175059, 
	6.71520734486919, 7.88670416113808, 10.0554641736314, 7.72255348977541, 
	7.07452949857107, 6.98005762075189, 11.5314993457204, 6.60602750697996, 
	7.34434175947847, 6.92351515359686, 7.03194766178152, 6.44249547823365, 
	6.72215549098823, 7.05833441095459, 8.00123439924826, 7.08254248378704, 
	7.98160373369768, 7.05574713457588, 6.97210459411254, 7.05505325667835, 
	6.61076537039296, 6.96423140172482, 6.13166540602059, 6.35104387908652, 
	7.09954036052535, 7.83808584689372, 8.49691204366696, 8.30822350540834, 
	10.6557273216795, 6.84183546107577, 7.44228076509898, 7.15071167550413, 
	6.75540718029703, 8.76705076327991, 7.07157384567065, 6.68776764659472, 
	8.26250063654079, 8.84353211114571, 7.10453985398356, 8.96533459829041, 
	6.81786893733246), .Dim = c(40L, 3L), .Dimnames = list(c("ILMN_3237291", 
	"ILMN_3214052", "ILMN_1821724", "ILMN_1891277", "ILMN_1757137", 
	"ILMN_3210171", "ILMN_1774937", "ILMN_1727880", "ILMN_2326376", 
	"ILMN_1741977", "ILMN_3243953", "ILMN_1737760", "ILMN_3232701", 
	"ILMN_1748836", "ILMN_1736056", "ILMN_1721087", "ILMN_3242786", 
	"ILMN_1681979", "ILMN_2259223", "ILMN_1691467", "ILMN_3182275", 
	"ILMN_1794818", "ILMN_1732782", "ILMN_3294263", "ILMN_1746720", 
	"ILMN_1760982", "ILMN_2262462", "ILMN_1706094", "ILMN_1799137", 
	"ILMN_2282477", "ILMN_2132309", "ILMN_1719340", "ILMN_2402972", 
	"ILMN_3225649", "ILMN_1736448", "ILMN_1651949", "ILMN_2072178", 
	"ILMN_1790688", "ILMN_1815083", "ILMN_1665873"), c("mutant", 
	"normal", "reference")))
	expect_equivalent(res, expected.res)
	
	res <- average.replicates(exprs(x.LumiBatch), as.character(classes))
	expect_equivalent(res, expected.res[,unique(classes)])
	
})


test_that("average.replicates on data.frame", {
	require(lumi) || stop("required package 'lumi' is not installed")
	load( system.file("examples", "x.LumiBatch.example.RDa", package="pwbc") ) # load x.LumiBatch object
	classes <- as.factor(sub("\\.[1-9]$", "", sampleNames(x.LumiBatch))) # character
	res <- average.replicates(as.data.frame(exprs(x.LumiBatch)), classes)
	expected.res <- structure(list(mutant = c(8.15150004530165, 10.1328745467003, 
	7.49775545803208, 7.4476481696329, 7.04401424384069, 11.2706112785709, 
	7.2284567269886, 7.38306934877342, 6.69251279620315, 6.99606650081166, 
	6.90200128207641, 6.84088257649328, 7.22454179317623, 8.13375079270194, 
	7.07309579565193, 7.78473033452568, 6.68312521707511, 7.08583171554452, 
	7.17478966392928, 6.87224219059137, 7.55816547408239, 5.87133975309216, 
	6.38738173698081, 7.01746351320833, 8.13554336425842, 8.7633395728979, 
	8.59823938194037, 11.3141877007204, 7.01052330080514, 7.24637689363769, 
	6.96449355710694, 6.57244311986803, 8.12546021072627, 7.21166738792694, 
	6.78906158755257, 8.05447404273069, 7.18752414981699, 7.16502260541611, 
	8.90414783930211, 6.78241837686025), normal = c(7.97501862079713, 
	9.58381220345539, 7.73357125947685, 7.20230792005529, 6.96649702440178, 
	11.2626563397547, 7.02700789930857, 6.55471603165014, 7.06798699771778, 
	7.03565621421496, 6.79018594678216, 6.64928854710145, 7.59169971170484, 
	8.04900209030033, 7.18510803926087, 8.49105372559396, 6.81295287383096, 
	6.91728578158643, 7.07703512676991, 6.53758355899177, 6.94811938176567, 
	6.20150859798392, 6.69471729652901, 6.82188137672293, 9.53758666042858, 
	8.92632167154977, 8.13770103101301, 11.0424415047994, 6.82414075153545, 
	7.08124634263355, 6.89704497671146, 6.72792696685374, 7.84103761036026, 
	7.57872959530072, 7.11182560014886, 7.9156382547289, 9.25011075284244, 
	7.00720109610439, 8.45333471175059, 6.71520734486919), reference = c(7.88670416113808, 
	10.0554641736314, 7.72255348977541, 7.07452949857107, 6.98005762075189, 
	11.5314993457204, 6.60602750697996, 7.34434175947847, 6.92351515359686, 
	7.03194766178152, 6.44249547823365, 6.72215549098823, 7.05833441095459, 
	8.00123439924826, 7.08254248378704, 7.98160373369768, 7.05574713457588, 
	6.97210459411254, 7.05505325667835, 6.61076537039296, 6.96423140172482, 
	6.13166540602059, 6.35104387908652, 7.09954036052535, 7.83808584689372, 
	8.49691204366696, 8.30822350540834, 10.6557273216795, 6.84183546107577, 
	7.44228076509898, 7.15071167550413, 6.75540718029703, 8.76705076327991, 
	7.07157384567065, 6.68776764659472, 8.26250063654079, 8.84353211114571, 
	7.10453985398356, 8.96533459829041, 6.81786893733246)), .Names = c("mutant", 
	"normal", "reference"), row.names = c("ILMN_3237291", "ILMN_3214052", 
	"ILMN_1821724", "ILMN_1891277", "ILMN_1757137", "ILMN_3210171", 
	"ILMN_1774937", "ILMN_1727880", "ILMN_2326376", "ILMN_1741977", 
	"ILMN_3243953", "ILMN_1737760", "ILMN_3232701", "ILMN_1748836", 
	"ILMN_1736056", "ILMN_1721087", "ILMN_3242786", "ILMN_1681979", 
	"ILMN_2259223", "ILMN_1691467", "ILMN_3182275", "ILMN_1794818", 
	"ILMN_1732782", "ILMN_3294263", "ILMN_1746720", "ILMN_1760982", 
	"ILMN_2262462", "ILMN_1706094", "ILMN_1799137", "ILMN_2282477", 
	"ILMN_2132309", "ILMN_1719340", "ILMN_2402972", "ILMN_3225649", 
	"ILMN_1736448", "ILMN_1651949", "ILMN_2072178", "ILMN_1790688", 
	"ILMN_1815083", "ILMN_1665873"), class = "data.frame")
	expect_equivalent(res, expected.res)
	
	res <- average.replicates(as.data.frame(exprs(x.LumiBatch)), as.character(classes))
	expect_equivalent(res, expected.res[,unique(classes)])
	
})

test_that("average.replicates on AnnotatedDataFrame", {
	require(lumi) || stop("required package 'lumi' is not installed")
	load( system.file("examples", "x.LumiBatch.example.RDa", package="pwbc") ) # load x.LumiBatch object
	classes <- as.factor(sub("\\.[1-9]$", "", sampleNames(x.LumiBatch))) # character
	res <- average.replicates(phenoData(x.LumiBatch), classes)
	expected.res.df <- structure(list(sampleID = c("APGI_1839-CR", "APGI_1151-TR", "APGI_1953-TR"
	)), .Names = "sampleID", row.names = c("mutant", "normal", "reference"
	), class = "data.frame")
	
	expect_equivalent(as(res,"data.frame"), expected.res.df)
	
	res <- average.replicates(phenoData(x.LumiBatch), as.character(classes))
	expected.res.df <- structure(list(sampleID = c("APGI_1151-TR", "APGI_1839-CR", "APGI_1953-TR"
	)), .Names = "sampleID", row.names = c("normal", "mutant", "reference"
	), class = "data.frame")
	
	expect_equivalent(as(res, "data.frame"), expected.res.df)
	
	#
	# use the example from ?AnnotatedDataFrame
	# 
	df <- data.frame(x=1:6,
                     y=rep(c("Low", "High"),3),
                     z=I(LETTERS[1:6]),
                     row.names=paste("Sample", 1:6, sep="_"))
    metaData <- data.frame(labelDescription=c(
                   "Numbers",
                   "Factor levels",
                   "Characters"))
    adf <- AnnotatedDataFrame(data=df, varMetadata=metaData)
    res <- average.replicates(adf, rep(c("Low", "High"),3))

	expect_equivalent(
		pData(res), 
		structure(list(x = 1:2, y = structure(c(2L, 1L), .Label = c("High", 
	"Low"), class = "factor"), z = structure(c("A", "B"), class = "AsIs")), .Names = c("x", 
	"y", "z"), row.names = c("Low", "High"), class = "data.frame")
	)

	expect_equivalent(
		varMetadata(res),
		structure(list(labelDescription = c("Numbers", "Factor levels", 
		"Characters")), .Names = "labelDescription", row.names = c("x", 
		"y", "z"), class = "data.frame")
	)
	
	expect_equivalent(
		varLabels(res),
		c("x", "y", "z")
	)
	
	
})


test_that("average.replicates on ExpressionSet", {
	load( system.file("examples", "x.ExpressionSet.example.RDa", package="pwbc") ) # load x.ExpressionSet object
	classes <- as.factor(sub("\\.[1-9]$", "", sampleNames(x.ExpressionSet)))
	
	res <- average.replicates(x.ExpressionSet, as.character(classes))
	
	load( system.file("examples", "x.ExpressionSet.avg.RDa", package="pwbc") ) # load x.ExpressionSet.avg object
	expected.res <- x.ExpressionSet.avg
	
	expect_equivalent(dims(res), dims(expected.res))
	expect_equivalent(assayData(res), assayData(expected.res))
	expect_equivalent(pData(res), pData(expected.res))
	expect_equivalent(phenoData(res), phenoData(expected.res))
	expect_equivalent(featureData(res), featureData(expected.res))
	expect_equivalent(protocolData(res), protocolData(expected.res))
	
})


test_that("average.replicates on LumiBatch", {
	load( system.file("examples", "x.LumiBatch.example.RDa", package="pwbc") )
	classes <- as.factor(sub("\\.[1-9]$", "", sampleNames(x.LumiBatch)))
	
	res <- average.replicates(x.LumiBatch, classes)
	
	load( system.file("examples", "x.LumiBatch.avg.RDa", package="pwbc") )
	expected.res <- x.LumiBatch.avg
	
	expect_equivalent(dims(res), dims(expected.res))
	expect_equivalent(assayData(res), assayData(expected.res))
	expect_equivalent(pData(res), pData(expected.res))
	expect_equivalent(phenoData(res), phenoData(expected.res))
	expect_equivalent(featureData(res), featureData(expected.res))
	expect_equivalent(protocolData(res), protocolData(expected.res))
	# LumiBatch specific slots, except history where we'll never match the timestamps.
	expect_equivalent(res@QC$sampleSummary, expected.res@QC$sampleSummary)
	expect_equivalent(res@controlData, expected.res@controlData)
	
})
