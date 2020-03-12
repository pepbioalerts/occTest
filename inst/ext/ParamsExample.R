#params example

#load necessary data
occ.species.nogeoscrub = system.file('ext/SampleData/Sp3Occurrence_v3.csv'
                                     ,package='occProfileR')
occ.species.nogeoscrub <-  read.csv (occ.species.nogeoscrub)
ras.env = system.file('ext/AllEnv.tif',package='occProfileR')
ras.env = raster (ras.env)
ras.dem = system.file('ext/DEM.tif',package='occProfileR')
ras.dem = raster (ras.dem)
ras.humaninfluence = system.file('ext/HII.tif',package='occProfileR')
ras.humaninfluence = raster(ras.humaninfluence)
countryshapefile.dir =  system.file('ext/CountryLevel',package='occProfileR')
countries.pol <- rgdal::readOGR(dsn =countryshapefile.dir,
                                layer = 'Countries_SPDF_MED')
occ.sp <- occ.species.nogeoscrub[,c('MAPX','MAPY')]
occ.sp <- occ.sp [complete.cases(occ.sp),]
coordinates (occ.sp) <-  ~ MAPX + MAPY
occ.sp@proj4string <- sp:::CRS (projection (ras.dem))


#occProfile default parameters
output.dir; sp.table; sp.name;r.env; r.dem;
countries.shapefile;countryfield.shapefile = 'ISO';ntv.ctry; inv.ctry; ras.hii;
taxonobservation.id = 'ID_GRAFIC';x.field = 'x'; y.field = 'y'; t.field = NULL;l.field = NULL;c.field =NULL;e.field = NULL;
coordinate.decimal.precision = 4;points.proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0");
alpha.parameter = 2;th.human.influence = 45;th.perc.outenv =  0.2;elev.quality.threshold = 100;
write.simple.output=T
write.full.output=T
excludeUnknownRanges= F ; excludeNotmatchCountry= F ; do.centroidDetection=T ; centroidDetection.method='all';
do.HumanDetection=T ; do.institutionLocality=T ; institutionLocality.method='all'; do.geoOutliers=T ; geoOutliers.method='all'
output.base.filename="QAQC"


#occProfile function parameters in the example
output.dir = getwd();
sp.table = occ.species.nogeoscrub;
sp.name = "my Species";
r.env = ras.env;
r.dem = ras.dem ;
countries.shapefile = countries.pol;
countryfield.shapefile = 'ISO';
ntv.ctry = c('ESP');
inv.ctry = c('FRA') ;
ras.hii = ras.humaninfluence;
x.field = 'MAPX';
y.field = 'MAPY'


#hey


