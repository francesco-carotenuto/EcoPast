#' @title Creates model's predictors
#' @description This function creates both the response variable and the predictor variables to be used with \code{\link{minosse.target}} function.
#' @usage minosse.data(obj,species_name,domain,time.overlap=0.95,coc.by="locality",
#' min.occs=3,abiotic.covs=NULL,combine.covs=TRUE,reduce_covs_by="pca",covs_th=0.95,
#' c.size="mean",bkg.predictors="presence",min.bkg=NULL,sampling.by.distance=TRUE,
#' prediction.ground=NULL,crop.by.mcp=FALSE,constrain.predictors=FALSE,
#' temporal.tolerance=NULL,projection=NULL,lon_0=NULL,lat_0=NULL,n.clusters=NULL,seed=NULL)
#' @param obj A \eqn{n x m} dataframe where \eqn{n} are the single occurrences and \eqn{m} are the following columns: spec (the species name), x and y (longitude and latitude in decimal degrees, respectively) and loc_id (an id identifying the fossil locality).
#' @param species_name Character. The name of the species whose geographic range is to be estimated.
#' @param domain Character or \code{NULL}. Only used if no prediction ground is provided. If set as \code{"land"}, then present day mainland portions are selected according to fossil data spatial distribution, if \code{"sea"}, marine domain portion is used as prediction ground. Default \code{NULL}.
#' @param time.overlap Numeric. The proportion of temporal intersection between the target and the predictors' time span. Default is 0.95.
#' @param coc.by Character. Either \code{"locality"} or \code{"cell"} to enable cooccurence analysis. See details below.
#' @param min.occs Either numeric or numeric vector of length 2. The number occurrences below which to discard a species from being either valid predictors either a target. If ony one value is provided, the threshold is the same for both target and predicors.
#' @param abiotic.covs the raster or rasters' stack of additional environmental predictors.
#' @param combine.covs Logical. Should \code{\link{minosse.data}} collate species and abiotic predictors when performing variables' number reduction? Default \code{TRUE} See details.
#' @param reduce_covs_by Character. The method used for predictors' number reduction. Available strategies are \code{"pca"}, \code{"variance"} or \code{"corr"}. See details.
#' @param covs_th Numeric. The threshold value used for predictors' number reduction strategy. See details.
#' @param c.size Numeric.This is the (square) cell resolution in meters for spatial interpolations. Some character values are possible: "mean", "semimean" and "max" (see details for forther explanations). If \code{prediction.ground} is not null and is a raster it is possible to use the raster resolution by setting "raster" as c.size.
#' @param bkg.predictors The number of pseudo absences to be simulated for each predictor species. If \code{"presence"}, the pseudo absences number equals the presences in each species.
#' @param min.bkg Numeric. If \code{bkg.predictors} is set to \code{"presence"}, this is the minimum number of pseudo absences to simulate if a species occurrence number is below this value.
#' @param sampling.by.distance Logical. \code{TRUE} for a distace-based density pseudo absences simulation. \code{FALSE} for a pure spatial random distribution of pseudo absences.
#' @param prediction.ground Either a raster or a SpatialPolygons class object where to perform all the spatial interpolations and target species prediction.
#' @param crop.by.mcp Logical. If \code{TRUE}, interpoalations and prediction are restricted to the \code{prediction.grund} area delimited by the MCP enclosing the fossil occurrences of the whole dataset. Default \code{FALSE}.
#' @param constrain.predictors Logical. Removing from the predictors' record all the localities not complying with spatial and temporal restrictions? Default is \code{FALSE}. See details.
#' @param temporal.tolerance Numeric. If \code{constrain.predictors} is \code{TRUE} this is the maximum difference (expressed in Million years unit, i.e. 0.1 = one kylo years) allowed between target and predictors species' age estimate of the localites. See details.
#' @param projection Character. This argument works only if \code{prediction.ground} is \code{NULL}. This is the equal-area projection for spatial interpolations. A character string in the proj4 format or either \code{"moll"} (Mollweide) or \code{"laea"} (Lambert Azimuthal equal area) projections (see details).
#' @param lon_0 Numeric. Only if \code{prediction.ground} is \code{NULL}. The longitude of the projection centre used when setting either \code{"moll"} or \code{"laea"} projections.  If \code{NULL} the mean longitude of the whole fossil record is used. Default \code{NULL}.
#' @param lat_0 Numeric. Only if \code{prediction.ground} is \code{NULL}. The latitude of the projection centre used when setting \code{"laea"} projection.  If \code{NULL} the mean latitude of whole fossil record is used. Default \code{NULL}.
#' @param n.clusters Numeric or \code{NULL}. The number of cores to use during spatial interpolations. If "automatic", the number of used cores is equal to the number of predictors. If predictors' number > the avaialble cores, all cores - 1 is then used. Default is \code{NULL}.
#' @param seed Numeric. The \code{seed} number for experiment replication.
#' @importFrom stats complete.cases optim prcomp var quantile
#' @importFrom methods as
#' @importFrom utils data
#' @importFrom raster rasterFromXYZ values crop extent crs mask projectRaster extract stack rasterToPoints raster
#' @importFrom sp coordinates CRS spTransform zerodist proj4string split
#' @importFrom cooccur cooccur pair
#' @importFrom fossil create.matrix
#' @importFrom parallelDist parDist
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom plotKML vect2rast
#' @importFrom maptools unionSpatialPolygons
#' @importFrom rgeos gDifference gBuffer gUnion
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% foreach
#' @importFrom usdm exclude vifcor vifstep
#' @importFrom smoothr smooth
#' @importFrom lava wait
#' @importFrom fields rdist.earth
#' @importFrom stats dist
#' @importClassesFrom sp SpatialPoints SpatialPointsDataFrame
#' @importClassesFrom sp Spatial SpatialPixels SpatialPixelsDataFrame
#' @importClassesFrom raster Raster RasterLayer
#' @export
#' @details In \code{minosse.data} there are different strategies for predictor species (covariates) dimension reduction. The first one considers only the species that are significantly related (positively or negatively) to the target species, then discarding all the others.
#' This first stratery uses the cooccurrence analysis that can be performed either at the locality level, i.e. by seeking pattern of cooccurrence whithin the species list of any single fossil locality, or at the cell level, i.e. by considering lists of
#' unique species occurring inside the squared cell of the prediction ground. A cell based analysis is useful when having many low-richness fossil localities. If the significantly relationships is less than 4, then all the species are considered.  Other
#' strategies can be used for predictors' dimensionality reduction. These additional strategies are performed over the predictors'maps and can employ one of the following methods: Principal Component Analysis (\code{"pca"}), Variance  Inflation Factor (\code{"variance"}) and correlation (\code{"corr"}).
#' These strategies need a threshold value (\code{"covs_th"}) to be set in order to select the predictors to retain. If the strategy is \code{"pca"}, then the \code{covs_th} is the percentage (from 1 to 100) of variance to be explained by PCA axes. If the strategy is \code{"corr"}, then \code{covs_th} is any
#' number between 0 and 1 indicating the correlation between predictors below which predictor species can be retained. If the strategy is \code{"variance"}, then \code{covs_th} is any mumber higher than one indicating the higher variance inflation that can be achieved by the predictor.
#' See details of \code{vif} function in the usdm package for further explanations. For \code{c.size} some automatic values are available: by setting "mean", the algorithm uses the mean nearest neighbour distance between fossil localities as cell resolution; by setting "semimean" it uses half of the average nearest neighbour distance, whereas, by setting "max" it uses the maximum nearest neighbour distance.
#' If \code{abiotic.covs} is not \code{NULL}, the \code{combine.covs} argument indicates if performing predictors maps number reduction by including (\code{TRUE}) or excluding (\code{FALSE}) abiotic covariates. If \code{FALSE}, abiotic covariates
#' are always included in the final dataset of predictors. Spatial interpolations always need equal area coordinates reference system to be used. The user can specify its own projected CRS (in the proj4 format, see https://proj4.org/operations/projections/index.html) or can use predefined choices like \code{"laea"} (for Lambert Azimuthal equal area) or \code{"moll"} (for Mollweide) projections. When setting predefined projections,
#' the user can specify the projection centre's coordinates in decimal degrees by \code{lon_0} and \code{lat_0} arguments. If both \code{lon_0} and \code{lat_0} are \code{NULL}, the mean longitude and latitude of the whole fossil record are used. Warning: If not \code{NULL}, the \code{prediction.ground}'s coordinates reference system has the priority over all the other projection settings.
#' time.overlap indicates the percentage of target and predictors species temporal overalp. Each predictor temporally overlapping target species' time span is automatically ruled out from prediction.
#' The argument \code{constrain.predictors} enables setting spatial and temporal restriction to predictors' fossil localities in order to be considered synchronous and syntopic to the target occurrences. The spatial restriction is set as the average nearest neighbour's distance between
#' target species fossil sites and cannot be changed. The user is allowed to set a temporal restriction by the argument \code{temporal.tolerance}. By this argument it is possible to set the maximum temporal differnce between target and predictors' fossil localities estimated ages. All the
#' sites exciding this value are ruled out from any analysis. See the reference paper and related Supporting Information for further details.
#' @return A list of three objects to be used with \code{minosse.target} function. The first element of the list is the dataset of target species occurrences. The second object is the raster stack of predictor species. The third object, if present, is the result of the cooccurrence analysis.
#' @author Francesco Carotenuto, francesco.carotenuto@unina.it
#' @examples
#'   \dontrun{
#' library(raster)
#' data(lgm)
#' raster(system.file("exdata/prediction_ground.gri", package="EcoPast"))->prediction_ground
#'
#' minosse_dat<-minosse.data(obj=lgm,species_name="Mammuthus_primigenius",
#' domain=NULL,time.overlap=0.95,coc.by="locality",min.occs=3,abiotic.covs=NULL,
#' combine.covs=TRUE,reduce_covs_by=NULL,covs_th=0.95,c.size="mean",
#' bkg.predictors="presence",min.bkg=100,sampling.by.distance=TRUE,
#' prediction.ground=prediction_ground,crop.by.mcp=FALSE,constrain.predictors=FALSE,
#' temporal.tolerance=NULL,projection=NULL,lon_0=NULL,lat_0=NULL,
#' n.clusters=3,seed=625)
#'
#'   }

minosse.data<-function(obj,
                       species_name,
                       domain=NULL,
                       time.overlap=0.95,
                       coc.by="locality",
                       min.occs=3,
                       abiotic.covs=NULL,
                       combine.covs=TRUE,
                       reduce_covs_by="pca",
                       covs_th=0.95,
                       c.size="mean",
                       bkg.predictors="presence",
                       min.bkg=NULL,
                       sampling.by.distance=TRUE,
                       prediction.ground=NULL,
                       crop.by.mcp=FALSE,
                       constrain.predictors=FALSE,
                       temporal.tolerance=NULL,
                       projection=NULL,
                       lon_0=NULL,
                       lat_0=NULL,
                       n.clusters=NULL,
                       seed=NULL) {

  # library(fossil)
  # library(cooccur)
  # library(sp)
  # library(raster)
  # library(maptools)
  # library(intamap)
  # library(plotKML)
  # library(parallel)
  # library(parallelDist)
  # library(rworldmap)
  # library(maptools)
  # library(raster)
  # library(automap)
  # library(dismo)
  # library(foreach)
  # library(doSNOW)
  # library(doParallel)
  # library(parallel)
  # library(usdm)
  # library(R.utils)
  # library(rgeos)
  # library(smoothr)
  # library(lava)
  # library(spatstat)
  # library(fields)
  # library(rangeBuilder)

  if (!requireNamespace("rangeBuilder", quietly = TRUE)) {
    stop("Package \"rangeBuilder\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  length(getLoadedDLLs())
  Sys.getenv("R_MAX_NUM_DLLs", unset = NA)

  if(isTRUE(constrain.predictors)) {
    if(is.null(temporal.tolerance)) stop("Enabling the predictors constraint requires setting a temporal tolerance")
  }

  if(is.null(prediction.ground)){
  if(is.null(projection)) stop("Spatial interpolation needs projected coordinates. Please provide a CRS object or set either  'laea' (for Lambert Azimuthal Equal Area projection) or 'moll' (for Mollweide projection)")
  if(all(projection!=c("laea","moll"))) my.prj<-sp::CRS(projection)
  }
  projection.ground<-NULL
  if(!is.null(abiotic.covs)){
    rasterFromXYZ(data.frame(sp::coordinates(abiotic.covs[[1]])[!is.na(values(abiotic.covs[[1]])),],z=1))->grid_mask;
    crop(abiotic.covs,grid_mask)->abiotic.covs
  }

  if(!is.null(abiotic.covs)){
    if(combine.covs==TRUE){
      abiotic.covs2<-abiotic.covs
    }
    if(combine.covs==FALSE) {abiotic.covs2<-abiotic.covs}
  } else abiotic.covs2<-abiotic.covs

  obj<-obj[which(obj$age>=range(obj[which(obj$spec==species_name),]$age)[1] & obj$age<=range(obj[which(obj$spec==species_name),]$age)[2]),]
  if(nrow(obj)==0) stop("Sorry, no species occurred in the selected time frame")

  as.character(obj$spec)->obj$spec
  if(!is.null(time.overlap)) {

  tapply(obj$age,obj$spec,function(x)range(x))->range_list
  target_range<-range(obj[which(obj$spec==species_name),]$age)
  sapply(range_list,function(x)sum(abs(x-target_range)))->range_diff
  quantile(target_range,(1-as.numeric(time.overlap)))-target_range[1]->target_tolerance
  names(target_tolerance)<-NULL
  valid_specs<-which(sapply(range_diff,function(x)x<=target_tolerance))
  names(range_list[valid_specs])->valid_names
  if(length(valid_names[-which((valid_names==species_name))])<4) stop("The number of selected predictors is less than the minimum required (4) to perform a good prediction. Please, consider to decrease the proportion of temporal overlap between target and predictor species")
  obj[which(obj$spec%in%valid_names),]->obj
  as.character(obj$loc_id)->obj$loc_id

  } else obj<-obj

  if(constrain.predictors==TRUE) {
  target_locs<-obj[which(obj$spec==species_name),c("x","y","age","loc_id")]
  target_locs[order(target_locs$loc_id),]->target_locs
  target_locs<-target_locs[!duplicated(target_locs$loc_id),]
  rownames(target_locs)<-NULL
  predictors_locs<-obj[-which(obj$loc_id%in%target_locs$loc_id),c("x","y","age","loc_id")]
  predictors_locs<-predictors_locs[!duplicated(predictors_locs$loc_id),]
  predictors_locs[order(predictors_locs$loc_id),]->predictors_locs
  rownames(predictors_locs)<-NULL

  target_pol<-rangeBuilder::getDynamicAlphaHull(target_locs, coordHeaders=c("x","y"),
                      clipToCoast = 'no',partCount=1)[[1]]
  predictors_locs->temp_predictors_locs
  sp::coordinates(temp_predictors_locs)<-~x+y
  raster::crs(temp_predictors_locs)<-raster::crs(target_pol)
  #plot(target_pol)
  quarantine<-temp_predictors_locs[which(temp_predictors_locs$loc_id%in%temp_predictors_locs[target_pol,]$loc_id),]
  #points(quarantine,pch=19,col="red")
  #points(target_locs[,c("x","y")],pch=19,col="green")
  as.data.frame(quarantine)->quarantine
  if(length(temp_predictors_locs[-which(temp_predictors_locs$loc_id%in%temp_predictors_locs[target_pol,]$loc_id),])>0){
  healthy<-temp_predictors_locs[-which(temp_predictors_locs$loc_id%in%temp_predictors_locs[target_pol,]$loc_id),]$loc_id
  obj[which(obj$loc_id%in%healthy),]->healthy
  #points(healthy[,c("x","y")],pch=19,col="blue")
  }
  if(length(temp_predictors_locs[-which(temp_predictors_locs$loc_id%in%temp_predictors_locs[target_pol,]$loc_id),])==0) {healthy<-NULL}
  target_locs->nearest_dat
  rbind(target_locs,quarantine)->original_dist_dat
  rownames(nearest_dat)<-NULL
  rownames(original_dist_dat)<-NULL
  rdist.earth(nearest_dat[,c("x","y")])->dist_mat
  rownames(dist_mat)<-nearest_dat$loc_id
  colnames(dist_mat)<-nearest_dat$loc_id
  diag(dist_mat)<-NA
  fields::rdist.earth(original_dist_dat[,c("x","y")])->original_dist_mat
  as.data.frame(original_dist_mat)->original_dist_mat
  rownames(original_dist_mat)<-original_dist_dat$loc_id
  colnames(original_dist_mat)<-original_dist_dat$loc_id

  original_dist_mat<-original_dist_mat[1:nrow(target_locs),-(1:nrow(target_locs))]
  mean(apply(dist_mat,1,function(z)unique(z[which(z==min(z,na.rm=TRUE))])))->target_aver_dist
  apply(original_dist_mat,1,function(z){
    j<-rep(0,length(z));
    j[which(z<=(target_aver_dist/2))]<-1;
    return(j)})->dist_mat_bin
  as.data.frame(dist_mat_bin)->dist_mat_bin
  colSums(dist_mat_bin)
  as.matrix(dist(original_dist_dat$age))->age_dist_mat
  age_dist_mat[1:nrow(target_locs),-(1:nrow(target_locs))]->age_dist_mat
  apply(age_dist_mat,1,function(z)names(z[which(z<=temporal.tolerance)]))->age_nearest_mat
  matrix(data=0,nrow(age_dist_mat),ncol=ncol(age_dist_mat))->age_dist_mat_bin
  rownames(age_dist_mat_bin)<-rownames(age_dist_mat)
  colnames(age_dist_mat_bin)<-colnames(age_dist_mat)
  lapply(1:length(age_nearest_mat),function(z){age_dist_mat_bin[rownames(age_dist_mat_bin)==names(age_nearest_mat)[z],colnames(age_dist_mat_bin)%in%age_nearest_mat[z][[1]]]<-1;
  return(age_dist_mat_bin[rownames(age_dist_mat_bin)==names(age_nearest_mat)[z],])})->age_dist_mat_bin
  do.call(rbind,age_dist_mat_bin)->age_dist_mat_bin
  dist_mat_bin*age_dist_mat_bin->valid_predictors_locs
  names(which(colSums(valid_predictors_locs)>0))->valid_predictors_locs
  if(length(which(valid_predictors_locs%in%unique(target_locs$loc_id)==FALSE))) print("No locality selected outside target species' sites according to the set restrictions")
  obj[obj$loc_id%in%obj[which(obj$spec==species_name),"loc_id"],]->target_obj
  obj[-which(obj$loc_id%in%target_obj$loc_id),]->predictors_obj
  predictors_obj<-predictors_obj[which(predictors_obj$loc_id%in%valid_predictors_locs),]
  rbind(predictors_obj,healthy)->predictors_obj
  if(nrow(predictors_obj)!=0){
  rbind(target_obj,predictors_obj)->obj
  } else {obj<-obj}
  }
  if(constrain.predictors==FALSE) {obj<-obj}

  if(is.null(coc.by)) {obj<-obj; valid_pairs<-NA}
  if(!is.null(coc.by)) {
    if(coc.by=="locality"){
      create.matrix(obj,"spec","loc_id")->mat
      set.seed(seed)
      if(class(try(cooccur(mat,spp_names =TRUE)->coc_mat))=="try-error") { cat("\n","Unable to perform cooccurrence analysis");
        obj<-obj} else {
          pair(coc_mat, spp=species_name, all = FALSE)->valid_pairs
          as.character(valid_pairs$sp2)->valid_species;
          if(length(valid_species)>=4) {
            valid_species<-c(species_name,valid_species)
            obj[which(obj$spec%in%valid_species),]->obj
          } else if(length(valid_species)<4) {obj<-obj}
        }
    }
  }

  rownames(obj)<-NULL

  row.names(obj)<-NULL
  obj$spec<-as.character(obj$spec)
  projection<-projection
  obj<-obj


  rownames(obj)<-NULL
  sp::coordinates(obj)<-~x+y
  sp::proj4string(obj)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  original_obj<-obj

  if(is.null(lon_0)){
    mean(extent(obj)[1:2])->lon_0
    mean(extent(obj)[3:4])->lat_0
  }
  if(!is.null(lon_0)) {
    lon_0->lon_0
    lat_0->lat_0
  }

  if (!is.null(prediction.ground)) {
    temp_obj <- obj
    obj <- sp::spTransform(obj, CRSobj = raster::crs(prediction.ground))
    if(class(prediction.ground)%in%c("SpatialPolygonsDataFrame","SpatialPolygons"))
      obj <- obj[prediction.ground,]
    if(class(prediction.ground)=="RasterLayer") {
      obj <- obj[as(prediction.ground,"SpatialPixels"),]
    }

  }
  if(is.null(prediction.ground)) {
    temp_obj<-obj
    laea_prj<-CRS(paste("+proj=laea +lat_0=",lat_0, " +lon_0=",lon_0," +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs",sep=""))
    moll_prj<-CRS(paste("+proj=moll +lon_0=",lon_0," +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs",sep=""))
    if(projection=="laea") {obj<-spTransform(obj,CRSobj=laea_prj)}
    if(projection=="moll") {obj<-spTransform(obj,CRSobj=moll_prj)}
    if(all(projection!=c("laea","moll"))) {obj<-spTransform(obj,CRSobj=my.prj)}
    }


  if(c.size=="mean") {
    dist_res<-parDist(as.matrix(sp::coordinates(obj[-zerodist(obj)[,2],])), method = "euclidean", threads = detectCores()-1)
    dist_res<-as.matrix(dist_res)
    c.size<-mean(apply(dist_res,1,function(x) min(x[x!=0])))
  } else c.size<-c.size


  if(c.size=="max") {
    dist_res<-parDist(as.matrix(sp::coordinates(obj[-zerodist(obj)[,2],])), method = "euclidean", threads = 3)
    dist_res<-as.matrix(dist_res)
    c.size<-max(apply(dist_res,1,function(x) min(x[x!=0])))
  } else c.size<-c.size

  if(c.size=="semimean") {
    vect2rast(obj[-zerodist(obj)[,2],],fname="age")->temp_res;
    summary(temp_res)$grid$cellsize[[1]]->c.size
  }

  if(is.null(prediction.ground)) {
    rworldmap::countriesCoarseLessIslands->countriesCoarseLessIslands
    w_temp<-countriesCoarseLessIslands[grep("Africa|Antarctica|Australia|Eurasia|North America|South America", countriesCoarseLessIslands$continent),]
    rgeos::gConvexHull(temp_obj)->conv_pol
    smoothr::smooth(conv_pol,method="densify")->conv_pol
    if(domain=="land") {
      w_pol<-w_temp[which(w_temp$continent%in%unique(w_temp[temp_obj,]$continent)),]
      w_pol<-unionSpatialPolygons(w_pol, cut(sp::coordinates(w_pol)[,1], range(sp::coordinates(w_pol)[,1]), include.lowest=FALSE))
      if(isTRUE(crop.by.mcp)) rgeos::gIntersection(conv_pol,w_pol)->w_pol else w_pol->w_pol
      if(projection=="laea"){
        w_pol<-spTransform(w_pol,CRSobj=laea_prj)
        prj_obj<-spTransform(temp_obj,CRSobj=laea_prj)
      }
      if(projection=="moll") {
        w_pol<-spTransform(w_pol,CRSobj=moll_prj)
        prj_obj<-spTransform(temp_obj,CRSobj=moll_prj)
      }
      if(all(projection!=c("laea","moll"))) {
        w_pol<-spTransform(w_pol,CRSobj=my.prj)
        prj_obj<-spTransform(temp_obj,CRSobj=my.prj)
      }
    }

    if(domain=="sea") {
      if(isTRUE(crop.by.mcp)) gDifference(conv_pol,w_temp)->w_pol else w_temp->w_pol
      if(projection=="laea"){
        w_pol<-spTransform(w_pol,CRSobj=laea_prj)
        prj_obj<-spTransform(temp_obj,CRSobj=laea_prj)
      }
      if(projection=="moll") {
        w_pol<-spTransform(w_pol,CRSobj=moll_prj)
        prj_obj<-spTransform(temp_obj,CRSobj=moll_prj)
      }
      if(all(projection!=c("laea","moll"))) {
        w_pol<-spTransform(w_pol,CRSobj=my.prj)
        prj_obj<-spTransform(temp_obj,CRSobj=my.prj)
      }
      dista<-apply(rgeos::gDistance(prj_obj, w_pol,byid=TRUE),2,min)
      dista<-dista[as.numeric(rownames(prj_obj[which(is.na(over(prj_obj,w_pol))),]@data))]
      dista+(dista/3)->dista
      gBuffer(prj_obj[prj_obj[which(is.na(over(prj_obj,w_pol))),],], byid=TRUE,width=dista)->to_add
      gUnion(w_pol,to_add)->w_pol
    }

    max.range1<-function(w_pol,c.size){
      vect2rast(as(extent(w_pol),"SpatialPolygons"),cell.size=c.size)->temp_w_rast
      if(projection=="laea") sp::proj4string(temp_w_rast)<-laea_prj else if(projection=="moll") sp::proj4string(temp_w_rast)<-moll_prj
      temp_w_rast<-as(temp_w_rast,"SpatialPixels")
      land<-temp_w_rast[w_pol]
      return(land)
    }

    max.range1(w_pol,c.size)->rast
    rast->ra
    ra<-raster(ra)
    raster::values(ra)<-1
    new_ra<-mask(ra,rast)
  }
  if (!is.null(prediction.ground)) {
    max.range <- function(w_pol, c.size, crop.by.mcp = crop.by.mcp) {
      temp_w_rast <- plotKML::vect2rast(w_pol, cell.size = c.size)
      if (isTRUE(crop.by.mcp)) {
        prj_obj <- sp::spTransform(temp_obj, CRSobj = raster::crs(temp_w_rast))
        conv_pol <- rgeos::gConvexHull(prj_obj)
        temp_w_rast <- raster::mask(temp_w_rast, conv_pol)
      }
      if (projection == "laea")
        raster::crs(temp_w_rast) <- laea_prj
      else if (projection == "moll")
        raster::crs(temp_w_rast) <- moll_prj
      temp_w_rast <- as(temp_w_rast, "SpatialPixels")
      land <- temp_w_rast[w_pol]
      return(land)
    }
    if (class(prediction.ground)%in%c("SpatialPolygonsDataFrame","SpatialPolygons")) {
      rast <- plotKML::vect2rast(prediction.ground, cell.size = c.size)
      rast<-raster(rast)
      if (isTRUE(crop.by.mcp)) {
        prj_obj <- sp::spTransform(temp_obj, CRSobj = raster::crs(rast))
        conv_pol <- rgeos::gConvexHull(prj_obj)
        rast <- raster::mask(rast, conv_pol)
        if (!is.null(abiotic.covs))
          abiotic.covs <- raster::mask(abiotic.covs,
                                       conv_pol)
      }
      rast <- as(rast, "SpatialPixels")
    }
    if (class(prediction.ground) == "RasterLayer") {
      if (c.size == "raster") {
        rast <- prediction.ground
        if (isTRUE(crop.by.mcp)) {
          prj_obj <- sp::spTransform(temp_obj, CRSobj = raster::crs(rast))
          conv_pol <- rgeos::gConvexHull(prj_obj)
          rast <- raster::mask(rast, conv_pol)
          if (!is.null(abiotic.covs))
            abiotic.covs <- raster::mask(abiotic.covs,
                                         conv_pol)
        }
        rast <- as(rast, "SpatialPixels")
      }
      if (is.numeric(c.size)) {
        gr_ext <- extent(prediction.ground)
        x_gr <- seq(from = gr_ext[1], to = gr_ext[2],
                    by = c.size)
        y_gr <- seq(from = gr_ext[3], to = gr_ext[4],
                    by = c.size)
        target_ras <- expand.grid(x_gr, y_gr)
        sp::coordinates(target_ras) <- ~Var1 + Var2
        raster::crs(target_ras) <- raster::crs(prediction.ground)
        sp::gridded(target_ras) <- TRUE
        target_ras <- raster(target_ras)
        prediction.ground <- raster::resample(prediction.ground,
                                              target_ras, method = "ngb")
        if (!is.null(abiotic.covs))
          abiotic.covs <- raster::resample(abiotic.covs,
                                           target_ras)
        rast <- prediction.ground
        if (isTRUE(crop.by.mcp)) {
          prj_obj <- sp::spTransform(temp_obj, CRSobj = raster::crs(rast))
          conv_pol <- rgeos::gConvexHull(prj_obj)
          rast <- raster::mask(rast, conv_pol)
          if (!is.null(abiotic.covs))
            abiotic.covs <- raster::mask(abiotic.covs,
                                         conv_pol)
        }
        rast <- as(rast, "SpatialPixels")
      }
    }
    ra <- rast
    ra <- raster(ra)
    raster::values(ra) <- 1
    new_ra <- mask(ra, rast)
    new_ra<-raster::crop(new_ra,extent(rasterToPolygons(new_ra)))
  }

  x.min<-obj@bbox[1,1]
  x.max<-obj@bbox[1,2]
  y.min<-obj@bbox[2,1]
  y.max<-obj@bbox[2,2]

  obj_list<-split(obj,obj$spec)

  lapply(obj_list,function(x)as.data.frame(x))->obj_list

  lapply(obj_list,function(x){
    x$value<-1;
    return(x)
  })->interp_dataset

  names(interp_dataset)<-names(obj_list)

  if(is.null(coc.by)) interp_dataset<-interp_dataset
  if(!is.null(coc.by)) {
    if(coc.by=="cell"){
      dat_for_coc<-data.frame(spec=obj$spec,locality=raster::extract(ra,obj,cellnumbers=TRUE)[,1]);
      dat_for_coc<-dat_for_coc[order(dat_for_coc$locality),];
      create.matrix(dat_for_coc,"spec","locality")->dat_for_coc
      set.seed(seed)
      if(class(try(cooccur(dat_for_coc,spp_names=TRUE)->coc_res))=="try-error") { cat("\n","Unable to perform cooccurrence analysis");
        interp_dataset<-interp_dataset} else {
          pair(coc_res, spp=species_name, all = FALSE)->valid_pairs
          as.character(valid_pairs$sp2)->valid_species;
          if(length(valid_species)>=4) {
            valid_species<-c(species_name,valid_species);
            obj[which(obj$spec%in%valid_species),]->obj;
            obj_list<-split(obj,obj$spec)
            lapply(obj_list,function(x)as.data.frame(x))->obj_list
            lapply(obj_list,function(x){
              x$value<-1;
              return(x)
            })->interp_dataset
            names(interp_dataset)<-names(obj_list)
          } else if(length(valid_species)<4) {interp_dataset<-interp_dataset}
        }
    }
  }

  rownames(obj)<-NULL

  ### For each species remove replicated occurences in each cell
  for(i in 1:length(interp_dataset)){
    sp::coordinates(interp_dataset[[i]])<-~x+y;
    raster::crs(interp_dataset[[i]])<-raster::crs(ra)
    raster::extract(ra,interp_dataset[[i]],cellnumbers=TRUE)[,1]->cells;
    interp_dataset[[i]]$cells<-cells
    interp_dataset[[i]][order(interp_dataset[[i]]$cells,interp_dataset[[i]]$value),]->interp_dataset[[i]]
    rownames(interp_dataset[[i]])<-NULL
    interp_dataset[[i]]<-interp_dataset[[i]][which(duplicated(interp_dataset[[i]]$cells,fromLast=TRUE)==FALSE),]
    interp_dataset[[i]]<-interp_dataset[[i]][,-which(names(interp_dataset[[i]]) %in% c("optional","cells"))]
  }
  gc()

  if(length(min.occs)==1) {min.occsX<-min.occs;min.occsY<-min.occs}
  if(length(min.occs)==2) {min.occsX<-min.occs[1];min.occsY<-min.occs[2]}

  interp_dataset[[which(names(interp_dataset)==species_name)]]->Y
  if(length(Y)<min.occsY) stop(paste("Unable to perform predictions with less than ", min.occsY, " occurrences")) ##
  interp_dataset[-which(names(interp_dataset)==species_name)]->X
  valid_X<-sapply(X,function(x)length(x)>=min.occsX)
  if(length(which(valid_X))<2) stop("Warning! Only one predictor detected. Target species prediction could be inaccurate. MInOSSE stops")
  if(length(which(valid_X))==0) stop("No predictor satisfies the minimum number of occurrences")
  if(length(which(is.na(raster::extract(new_ra,Y))))>length(Y)/2) {print("Warning, more than 50% of target species occurrences falls out of the prediction ground domain! Consider fixing target species geographic coordinates by using fix.coastal.points function.")}
  #if(length(which(is.na(raster::extract(new_ra,Y))))>length(Y)/2) {print("Warning, more than 50% of target species occurrences falls out of the prediction ground domain! By pressing ENTER, you allow MInOSSE computing the predictors and, then, you can try to fix target species geographic coordinates by using fix.coastal.points function. By pressing ESC MInOSSE stops"); lava::wait()}



  ####### PREDICTORS' PREDICTION CORE

  X_po_dis<-list()
  zero<-list()
  for(i in 1:length(X)){
    set.seed(seed)
    X_po_dis[[i]]<-raster::distanceFromPoints(new_ra, X[[i]])
    X_po_dis[[i]]<-raster::crop(X_po_dis[[i]],new_ra)
    X_po_dis[[i]]<-raster::mask(X_po_dis[[i]],new_ra)
    if(is.numeric(bkg.predictors)) num_zero<-bkg.predictors
    if(bkg.predictors=="presence") {
      num_zero<-length(X[[i]])
      if(!is.null(min.bkg)) {
        if(length(X[[i]])<min.bkg) num_zero<-min.bkg else num_zero<-length(X[[i]])
      }
      if(is.null(min.bkg)) num_zero<-num_zero
    }

    set.seed(seed)
    dismo::randomPoints(X_po_dis[[i]],
                        n=num_zero,p=X[[i]],
                        prob=sampling.by.distance)->zero[[i]]
    data.frame(zero[[i]],value=0)->zero[[i]]
    sp::coordinates(zero[[i]])<-~x+y
    sp::proj4string(zero[[i]])<-raster::crs(X[[i]])
    X[[i]]<-rbind(X[[i]][,"value"],zero[[i]])
  }

  ############ Interpolations of the valid predictors
  cat("\n","Performing interpolations of the predictors' occurreces")
  cat(paste("\n","Interpolating",length(X),"predictors",sep=" "))

  if(is.null(n.clusters)) {
    interps<-list()
    for(i in 1:length(X)) {
      seed<-seed
      set.seed(seed)
      methodName<-"automap"
      if(class(try(interps[[i]]<-minosse_interpolate(X[[i]],
                                                as(new_ra,"SpatialPixels"),
                                                methodName=methodName,
                                                outputWhat=list(mean=TRUE,
                                                                variance=TRUE,MOK=1),
                                                optList = list(confProj = TRUE,set.seed=seed)),silent=FALSE))=="try-error") interps[[i]]<-NULL else interps[[i]]<-interps[[i]]
      print(i)
    }
  }


  if(!is.null(n.clusters)) {
  if(n.clusters=="automatic") {
    if((detectCores()-1)>length(X)) n.clusters<-length(X) else (detectCores()-1)->n.clusters
  }

  cl <- makeCluster(n.clusters, type="SOCK")
  registerDoSNOW (cl)
  interps<-list()

  system.time(foreach(i=1:length(X),.packages=c("intamap","Metrics","raster"),.export=c("minosse_interpolate")) %dopar% {
    seed<-seed
    set.seed(seed)
    methodName<-"automap"
    if(class(try(interps<-minosse_interpolate(X[[i]],
                                              as(new_ra,"SpatialPixels"),
                                              methodName=methodName,
                                              outputWhat=list(mean=TRUE,
                                                              variance=TRUE,MOK=1),
                                              optList = list(confProj = TRUE,set.seed=seed)),silent=TRUE))=="try-error") interps<-NULL else interps<-interps
    print(i)
    return(interps)}->interps)
  stopCluster(cl)
  gc()
  }

  sapply(interps,function(x)!is.null(x))->valid_interps
  interps[which(valid_interps)]->interps
  X[which(valid_interps)]->X

  MOK<-1
  if(!is.null(MOK)){
    mappa<-lapply(interps,function(x)rasterFromXYZ(x$predictions[,"MOK1"]))
  }
  if(is.null(MOK)){
    mappa<-lapply(interps,function(x)rasterFromXYZ(x$predictions[,1]))
  }

  mappa<-lapply(mappa,function(x){
    if(!is.null(prediction.ground)){
      sp::proj4string(x)<-raster::crs(prediction.ground)
    }
    if(is.null(prediction.ground)) {
      if(projection=="laea") sp::proj4string(x)<-laea_prj;
      if(projection=="moll") sp::proj4string(x)<-moll_prj
    }
    return(x)
  })
  gc()

  names(mappa)<-names(X)
  raster::stack(mappa)->X_maps

  if(!is.null(reduce_covs_by)){
    if(!is.null(abiotic.covs2)){
  if(combine.covs==TRUE){
    abiotic.covs2->abiotic.covs;
    raster::unstack(X_maps)->res_maps
    res_maps[length(res_maps)+1]<-abiotic.covs
    raster::resample(res_maps[[length(res_maps)]],res_maps[[1]])->res_maps[[length(res_maps)]]
    the_min_maps<-which.min(sapply(res_maps,function(z)length(na.omit(values(z)))))
    lapply(res_maps,function(z)raster::crop(z,res_maps[[the_min_maps]]))->res_maps
    lapply(res_maps,function(z)raster::mask(z,res_maps[[the_min_maps]]))->res_maps
    the_min_maps<-which.min(sapply(res_maps,function(z)length(na.omit(values(z)))))
    lapply(res_maps,function(z)raster::mask(z,res_maps[[the_min_maps]]))->res_maps
    raster::stack(raster::scale(raster::stack(res_maps)))->X_maps
    } else {X_maps->X_maps}}
  } else {X_maps->X_maps}
  if(length(names(X_maps))<=2) {
    reduce_covs_by<-NULL;
    cat(paste("\n","Predictors' reduction not necessary; number of predictors = ",length(names(X_maps)),sep="" ))
  }


  if(is.null(reduce_covs_by)) {
    X_maps->mappa_pix;
    X_maps->mappa_ras;
    abiotic.covs2->abiotic.covs
    if(!is.null(abiotic.covs)) {
      raster::unstack(mappa_ras)->res_ras
      res_ras[length(res_ras)+1]<-abiotic.covs
      raster::resample(res_ras[[length(res_ras)]],res_ras[[1]])->res_ras[[length(res_ras)]]
      the_min<-which.min(sapply(res_ras,function(z)length(na.omit(values(z)))))
      lapply(res_ras,function(z)raster::crop(z,res_ras[[the_min]]))->res_ras
      lapply(res_ras,function(z)raster::mask(z,res_ras[[the_min]]))->res_ras
      the_min<-which.min(sapply(res_ras,function(z)length(na.omit(values(z)))))
      lapply(res_ras,function(z)raster::mask(z,res_ras[[the_min]]))->res_ras
      raster::stack(raster::scale(raster::stack(res_ras)))->mappa_ras}
    as(mappa_pix,"SpatialPixelsDataFrame")->mappa_pix
    cat(paste("\n","Number of predictors without reduction:",length(names(mappa_ras))))
  }
  if(!is.null(reduce_covs_by)) {
    if(reduce_covs_by=="corr"){
      if(covs_th=="automatic"){
        gimme.th<-function(par,optim.maps){
          abs(length(names(exclude(optim.maps,vifcor(optim.maps,th=par))))-8)
        }
        optim.res<-optim(par=0.5,fn=gimme.th,optim.maps=X_maps,method="Brent",lower=0.1,upper=0.95)$par
        exclude(X_maps,vifcor(X_maps,th=optim.res))->mappa_ras;
        mappa_ras->mappa_pix
      }
      if(is.numeric(covs_th)) {
        if(length(names(exclude(X_maps,vifcor(X_maps,th=covs_th))))<2) {
          X_maps->mappa_pix;
          X_maps->mappa_ras
        }
        if(length(names(exclude(X_maps,vifcor(X_maps,th=covs_th))))>=2){
          set.seed(seed)
          exclude(X_maps,vifcor(X_maps,th=covs_th))->mappa_pix;
          set.seed(seed)
          exclude(X_maps,vifcor(X_maps,th=covs_th))->mappa_ras
        }
      }
      as(mappa_pix,"SpatialPixelsDataFrame")->mappa_pix;
      cat(paste("\n","Number of predictors after reduction:",ncol(mappa_pix)))
    }
    if(reduce_covs_by=="variance") {
      if(length(names(exclude(X_maps,vifstep(X_maps,th=covs_th))))<2){
        X_maps->mappa_pix;
        X_maps->mappa_ras
      }
      if(length(names(exclude(X_maps,vifstep(X_maps,th=covs_th))))>=2) {
        set.seed(seed)
        exclude(X_maps,vifstep(X_maps,th=covs_th))->mappa_pix;
        set.seed(seed)
        exclude(X_maps,vifstep(X_maps,th=covs_th))->mappa_ras
      }
      as(mappa_pix,"SpatialPixelsDataFrame")->mappa_pix;
      cat(paste("\n","Number of predictors after reduction:",ncol(mappa_pix)))
    }
    if(reduce_covs_by=="pca") {
      X_maps->original_X_maps
      as.data.frame(rasterToPoints(X_maps))->X_maps
      pca_res <- prcomp(X_maps[complete.cases(X_maps),-c(1:2)], scale = TRUE)
      vars <- apply(pca_res$x, 2, var)
      props <- vars / sum(vars)
      max(which(cumsum(props)<=covs_th))+2->last_axis
      if((last_axis-2)<2) {
        original_X_maps->mappa_pix;
        as(mappa_pix,"SpatialPixelsDataFrame")->mappa_pix;
        original_X_maps->mappa_ras
      }
      if((last_axis-2)>=2) {
        X_maps[complete.cases(X_maps),-c(1:2)]<-pca_res$x
        if(class(try(X_maps[,1:last_axis]->pca_covs))=="try-error") cat("\n","Please set a higher threshold for PCA axes selection") else pca_covs->pca_covs
        if(ncol(pca_covs)==3) {rasterFromXYZ(cbind(pca_covs[,1:2],pca_covs[,3]))->pca_ras} else if(ncol(pca_covs)>3) {
          apply(pca_covs[,-c(1:2)],2,function(x)rasterFromXYZ(cbind(pca_covs[,1:2],x)))->pca_ras
        }
        stack(pca_ras)->pca_stack
        names(pca_stack)<-paste("PC",1:length(names(pca_stack)),sep="")
        raster::crs(pca_stack)<-raster::crs(ra)
        pca_stack->mappa_ras
        as(pca_stack,"SpatialPixelsDataFrame")->mappa_pix
        raster::crs(mappa_pix)<-raster::crs(rast)
      }
      cat(paste("\n","Number of predictors after reduction:",length(names(mappa_ras))))
    }
  }
  if(!is.null(reduce_covs_by)) {
    if(reduce_covs_by=="corr"){
      cat(paste("\n","Number of predictors after reduction:",ncol(mappa_pix)))
    }
    if(reduce_covs_by=="variance") {
      cat(paste("\n","Number of predictors after reduction:",ncol(mappa_pix)))
    }
    if(reduce_covs_by=="pca") {
      cat(paste("\n","Number of predictors after reduction:",ncol(mappa_pix)))
    }
  if(!isTRUE(combine.covs)) {
    raster::unstack(mappa_ras)->res_ras
    res_ras[length(res_ras)+1]<-abiotic.covs
    raster::resample(res_ras[[length(res_ras)]],res_ras[[1]])->res_ras[[length(res_ras)]]
    the_min<-which.min(sapply(res_ras,function(z)length(na.omit(values(z)))))
    lapply(res_ras,function(z)raster::crop(z,res_ras[[the_min]]))->res_ras
    lapply(res_ras,function(z)raster::mask(z,res_ras[[the_min]]))->res_ras
    the_min<-which.min(sapply(res_ras,function(z)length(na.omit(values(z)))))
    lapply(res_ras,function(z)raster::mask(z,res_ras[[the_min]]))->res_ras
    raster::stack(raster::scale(raster::stack(res_ras)))->mappa_ras
    }
   }


  return(list(response=Y,predictors=mappa_ras,sig_species=valid_pairs))

  detachAllPackages <- function() {
    basic.packages <- c("package:stats","package:raster","package:sp","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
    package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
    package.list <- setdiff(package.list,basic.packages)
    if (length(package.list)>0)  for (package in package.list) detach(package, unload=TRUE,force=TRUE,character.only=TRUE)
  }
  detachAllPackages()
  detachAllPackages()
  R.utils::gcDLLs(gc=TRUE, quiet=TRUE)
}

