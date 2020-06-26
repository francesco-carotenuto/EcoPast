#' @title The species geographic range polygon
#' @description This function turns \code{\link{minosse}} geographic range raster map into a spatial polygon with the related surface area.
#' @usage minosse.poly(minosse_res, th_num=3)
#' @param minosse_res  the output of \code{\link{minosse.target}} function or the \code{minosse_res} element list in the minosse function output.
#' @param th_num Numeric. The Regression Kriging prediction map binarization threshold value index as reported in the \code{optimal.thresholds} function in the PresenceAbsence package. Default 3 = MaxSens+Spec,	maximizes (sensitivity+specificity)/2.
#' @importFrom raster rasterToPolygons
#' @importFrom rgeos gUnaryUnion gArea
#' @details NULL
#' @export
#' @return The target species geographic range as a SpatialPolygonsDataFrame object.
#' @author Francesco Carotenuto, francesco.carotenuto@unina.it
#' @examples
#'   \dontrun{
#'   library(raster)
#'   data(lgm)
#'   raster(system.file("exdata/prediction_ground.gri", package="EcoPast"))->prediction_ground
#'
#'   minosse_dat<-minosse.data(obj=lgm,species_name="Mammuthus_primigenius",
#'   domain="land",time.overlap=0.95,coc.by="locality",min.occs=3,abiotic.covs=NULL,
#'   combine.covs=FALSE,reduce_covs_by="pca",covs_th=0.95,c.size="mean",
#'   bkg.predictors="presence",min.bkg=100,sampling.by.distance=TRUE,
#'   crop.by.mcp=FALSE,constrain.predictors=FALSE,temporal.tolerance=NULL,
#'   prediction.ground=prediction_ground,projection="laea",lon_0=85,lat_0=45,
#'   n.clusters=NULL,seed=625)
#'
#'   minosse_res<-minosse.target(resp=minosse_dat[[1]],predictors=minosse_dat[[2]],
#'   bkg="presence", min.bkg = 100,n.sims=10,sampling.by.distance=TRUE,n.folds=1,
#'   n.sims.clusters=NULL,seed=625)
#'
#'   minosse.poly(minosse_res)
#'   }

minosse.poly<-function(minosse_res, th_num=3) {
  # library(raster)
  # library(rgeos)
  minosse_res[[1]][[th_num]]->selected
  rasterToPolygons(selected,function(x) x>0 ,dissolve=TRUE)->minosse_poly
  gUnaryUnion(minosse_poly)->minosse_poly
  as(minosse_poly,"SpatialPolygonsDataFrame")->minosse_poly
  names(minosse_poly@data)<-"area"
  gArea(minosse_poly)->minosse_poly@data$area
  return(minosse_poly)
}
