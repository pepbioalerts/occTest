codes.list.NOX <- c('E','D','C','B','A')

conditionals.codes.list.NOX <- list (
  E='missingEnvironment == 1 | botgarden == 1',
  D="(lazylogic (e = 'geooutlier == 1') | lazylogic (e = 'human == 1') | centroids == 1 ) &  (lazylogic ('outlierEnvironment == 1') )",
  C="lazylogic('outlierEnvironment == 1')",
  B="(lazylogic (e = 'geooutlier == 1')  |  lazylogic (e = 'human == 1') | centroids == 1)",
  A=NA
)


descriptors.codes.list.NOX <- list(
  E=list("if(missingEnvironment == 1) {'MissingEnvironment'}",
         "if(botgarden == 1) {'LikelyBotanicalGarden'}"),
  D=list ("GeoEnvOutlier",
          "if(lazylogic (e = 'geooutlier == 1')) {'OutAlphaHull'}",
          "if(lazylogic (e = 'human == 1')) {'HighHumanEnvironment'}",
          "if(lazylogic (e=  'centroids  == 1')) {'InCentroid'}"),
  C="EnvironmentalOutlier",
  B=list ("GeoOutlier",
          "if(lazylogic (e = 'geooutlier == 1')) {'OutAlphaHull'}",
          "if(lazylogic (e = 'human == 1')) {'HighHumanEnvironment'}",
          "if(lazylogic (e=  'centroids  == 1')) {'InCentroid'}"),
  A=NA
) 








