#' @title A ready to use function to compute species geographic range with default settings
#' @description This function performs \code{\link{minosse.data}} and \code{\link{minosse.target}} all at once for one or multiple target species.
#' @usage minosse(dat,species.name=NULL,domain=NULL,time.overlap=0.95,
#' coc_by="locality",min.occs=10,min.bkg=100,sampling.by.distance=TRUE,
#' n.sims=10,n.clusters=NULL,n.sims.clusters=NULL,prediction.ground=NULL,
#' abiotic.covs=NULL,combine.covs=FALSE,reduce_covs_by="pca",covs_th=0.95,crop.by.mcp=FALSE,constrain.predictors=FALSE,
#' temporal.tolerance=NULL,projection=NULL,lon_0=NULL,lat_0=NULL,seed=NULL)
#' @param dat A \eqn{n x m} dataframe where \eqn{n} are the single occurrences and \eqn{m} are the following columns: spec (the species name), x and y (longitude and latitude in decimal degrees, respectively) and loc_id (an id identifying the fossil locality).
#' @param species.name The character vector of the species names whose geographic ranges are to be estimated. If \code{NULL}, \code{minosse} runs for all the species in the fossil dataset.
#' @param domain Character or \code{NULL}. Only used if no prediction ground is provided. If set as \code{"land"}, then present day mainland portions are selected according to fossil data spatial distribution, if \code{"sea"}, marine domain portion is used as prediction ground. Default \code{NULL}.
#' @param time.overlap Numeric. The proportion of temporal intersection between the targer and the predictors' time span. Default is 0.95.
#' @param coc_by Character. This argument enables the cooccurrence analysis performed either at the locality level (then use \code{"locality"}) or at cell of the prediction ground level (then use \code{"cell"}). See details below.
#' @param min.occs Numeric. For both target and predictor species. The number occurrences below which \code{minosse.data} discards a species from being a valid target or elegible predictors. Default 10.
#' @param min.bkg Numeric. \code{minosse} function by default simulates as many pseudo absences as the presences, thereby this is the minimum number of pseudo absences to simulate if a species occurrence number is below this value.
#' @param sampling.by.distance Logical. If \code{TRUE} pseudo absences are simulated with an intensity proportional to the distance to the presence data. If \code{FALSE} a pure spatial random distribution is simulated.
#' @param n.sims Numeric or \code{NULL}. The number of pseudo absences simulations (see details).
#' @param n.clusters Numeric or \code{NULL}. The number of cores to use during spatial interpolations. If "automatic", the number of used cores is equal to the number of predictors. If predictors' number > the avaialble cores, all cores - 1 is then used. Default is \code{NULL}.
#' @param n.sims.clusters Numeric or \code{NULL}. The number of machine cores to use when setting multiple pseudo absences simulations. By setting "automatic", as many cores are used as the number of simulations. In the case number of simulations is higher than the available cores, then the number of cores = available cores -1. Default \code{NULL}.
#' @param prediction.ground Either a raster or a SpatialPolygons class object where to perform all the spatial interpolations. This will be the prediction ground used when running \code{minosse.target}.
#' @param abiotic.covs The raster or rasters' stack of additional environmental predictors.
#' @param combine.covs Logical. Should \code{minosse.data} collate species and abiotic predictors when performing variables' number reduction? Default \code{FALSE}. See details.
#' @param crop.by.mcp Logical. If \code{TRUE}, the interpoalation of the predictors species data are limited to the \code{prediction.grund} area delimited by the MCP of ALL the fossil occurrences. Default \code{FALSE}.
#' @param constrain.predictors Logical. Removing from the predictors' record all the localities not complying with spatial and temporal restrictions? Default is \code{FALSE}. See \code{minosse.data} function details.
#' @param temporal.tolerance Numeric. If \code{constrain.predictors} is \code{TRUE} this is the maximum difference (expressed in Million years unit, i.e. 0.1 = one kylo years) allowed between target and predictors species' age estimate of the localites. See \code{minosse.data} function details.
#' @param projection Character. This argument works only if \code{prediction.ground} is \code{NULL}. This is the euqual-area projection for spatial interpolations. A character string in the proj4 format or either \code{"moll"} (Mollweide) or \code{"laea"} (Lambert Azimuthal equal area) projections (see details in \code{minosse.data function}).
#' @param lon_0 Numeric. Only if \code{prediction.ground} is \code{NULL}. The longitude of the projection centre used when setting either \code{"moll"} or \code{"laea"} projections.  If \code{NULL} the mean longitude of the whole fossil record is used. Default \code{NULL}.
#' @param lat_0 Numeric. Only if \code{prediction.ground} is \code{NULL}. The latitude of the projection centre used when setting \code{"laea"} projection.  If \code{NULL} the mean latitude of whole fossil record is used. Default \code{NULL}.
#' @param seed Numeric. The \code{seed} number for experiment replication.
#' @details NULL
#' @export
#' @return A list of three objects where the first one is the polygon of the target species geographic range (a SpatialPolygons object), the second element is the output of \code{minosse.target} function (see \code{minosse.target} function for details) and the last one is the result (if available) of the cooccurrence analysis.
#' If \code{minosse} function is performed for multiple species all at once, then \code{minosse} output described above is replicated for each target species.
#' @author Francesco Carotenuto, francesco.carotenuto@unina.it
#' @examples
#'   \dontrun{
#'   library(raster)
#'   data(lgm)
#'   raster(system.file("exdata/prediction_ground.gri", package="EcoPast"))->prediction_ground
#'
#'   mam<-minosse(dat=lgm,species.name="Mammuthus_primigenius",domain=NULL,
#'   time.overlap=0.95,prediction.ground=prediction_ground,crop.by.mcp=FALSE,
#'   coc_by="locality",min.occs=3,min.bkg=100,sampling.by.distance=TRUE,
#'   constrain.predictors=FALSE, temporal.tolerance=NULL,n.sims=10,n.clusters=NULL,
#'   n.sims.clusters=NULL,projection="laea",lon_0 = NULL,lat_0 = NULL,seed=625)
#'
#'   }


minosse<-function(dat,
                  species.name=NULL,
                  domain=NULL,
                  time.overlap=time.overlap,
                  coc_by="locality",
                  min.occs=10,
                  min.bkg=100,
                  sampling.by.distance=TRUE,
                  n.sims=10,
                  n.clusters=NULL,
                  n.sims.clusters=NULL,
                  prediction.ground=NULL,
                  abiotic.covs=NULL,
                  combine.covs=FALSE,
                  reduce_covs_by="pca",
                  covs_th=0.95,
                  crop.by.mcp=FALSE,
                  constrain.predictors=FALSE,
                  temporal.tolerance=NULL,
                  projection=NULL,
                  lon_0=NULL,
                  lat_0=NULL,
                  seed=NULL){

  if(is.null(species.name)) unique(as.character(dat$spec))->spec_list else strsplit(species.name," ")->spec_list
  if(is.null(projection)) stop("MInOSSE requires equal area-projected coordinates reference system to properly perform spatial analyses")
  minosse_dat<-list()
  minosse_res_list<-list()
  for(i in 1:length(spec_list)){
    if(class(try(minosse_dat[[i]]<-minosse.data(
      obj=dat,
      species_name=spec_list[[i]],
      domain=domain,
      coc.by=coc_by,
      min.occs=min.occs,
      abiotic.covs=abiotic.covs,
      combine.covs=combine.covs,
      reduce_covs_by="pca",
      covs_th=0.95,
      c.size="mean",
      bkg.predictors="presence",
      min.bkg=min.bkg,
      sampling.by.distance=sampling.by.distance,
      prediction.ground=prediction.ground,
      projection=projection,
      lon_0=lon_0,
      lat_0=lat_0,
      n.clusters=n.clusters,
      seed=seed)
    ))=="try-error") {next} else {minosse_dat[[i]]<-minosse_dat[[i]];

    if(class(try(minosse_res_list[[i]]<-minosse.target(
      resp=minosse_dat[[i]][[1]],
      predictors=minosse_dat[[i]][[2]],
      bkg="presence",
      min.bkg = min.bkg,
      n.sims=n.sims,
      sampling.by.distance=sampling.by.distance,
      n.folds=1,
      n.sims.clusters=n.sims.clusters,
      seed=seed)))=="try-error") {next} else minosse_res_list[[i]]<-minosse_res_list[[i]]
    }
    print(i)
  }
  gc()

  minosse_res_list[which(!is.na(minosse_res_list))]->minosse_res_list
  spec_list[which(!is.na(minosse_res_list))]->spec_list
  minosse_dat[which(!is.na(minosse_res_list))]->minosse_dat
  lapply(minosse_res_list,function(x) minosse.poly(x,th_num=3))->minosse_poly_list
  minosse_res_list<-lapply(seq(1:length(spec_list)),function(x)list(minosse_poly=minosse_poly_list[[x]],minosse_res=minosse_res_list[[x]],sig_species=minosse_dat[[x]]$sig_species))
  names(minosse_res_list)<-spec_list
  if(length(minosse_res_list)==1) minosse_res_list[[1]]->minosse_res_list
  return(minosse_res_list)
}
