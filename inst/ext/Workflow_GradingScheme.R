# #### OLD CODE FOR USING TO BUILD GRADING SCHEME 


#create the grading system [ideally we should pull this out from an xml file....one day]
#codes.list.NOX <- c('E','D','C','B','A')
#grep ('_score',names(df.qualityAssessment),value = T)


# conditionals.codes.list.NOX <- list (
#   E= "(.lazylogic (e = 'HumanDetection_score >= test.strictness.value')) | (.lazylogic (e = 'institutionLocality_score >= test.strictness.value')) | (.lazylogic (e = 'centroidDetection_score >= test.strictness.value')) | (.lazylogic (e = 'unknownRange_score >= test.strictness.value')) | (.lazylogic (e = 'wrongReportCtry_score >= test.strictness.value')) ",
#   D= "(.lazylogic (e = 'geoOutliers_score >= test.strictness.value') ) &  (.lazylogic ('envOutliers_score  >= test.strictness.value') )",
#   C= ".lazylogic ('envOutliers_score  >= test.strictness.value')",
#   B= ".lazylogic (e = 'geoOutliers_score >= test.strictness.value') ",
#   A=NA
#   )
# 
# 
# descriptors.codes.list.NOX <- list(
#   E=list("placeIssues",
#          "if(.lazylogic (e = 'HumanDetection_score >= test.strictness.value'))  {'HumanEnvironment'}",
#          "if(.lazylogic (e = 'institutionLocality_score >= test.strictness.value'))  {'instituionLocation'}",
#          "if(.lazylogic (e = 'centroidDetection_score >= test.strictness.value')) {'centroid'}",
#          "if(.lazylogic (e = 'unknownRange_score >= test.strictness.value')) {'unknownRange'}",
#          "if(.lazylogic (e = 'wrongReportCtry_score >= test.strictness.value')) {'wrongReportedCtry'}"
#          ),
#   D="GeoEnvOutlier",
#   C="EnvOutlier",
#   B="GeoOutlier",
#   A=NA
#   )
# 
# 
# grading.scheme.vNOX <- occTest:::.make.qgrading.system(
#                                             codes = codes.list.NOX,
#                                             conditionals.codes =
#                                               conditionals.codes.list.NOX,
#                                             descriptors.codes =
#                                               descriptors.codes.list.NOX)
# ### label the different occurrences
# dat.qgraded <- occTest:::.qgrade.data.frame(df = df.qualityAssessment,
#                                   grading.scheme = grading.scheme.vNOX)
# ### add lables to the data
# 


