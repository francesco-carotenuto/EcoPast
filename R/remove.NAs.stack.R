#' @title Remove unpaired NA cells from a RasterStack
#'
#' @description The function removes unpaired NA cells from a RasterStack.
#' @param rast.stack a \code{RasterStack} object where to remove NAs.
#' @export
#' @importFrom raster calc stack
#' @return A \code{RasterStack} object without unpaired NA cells.
#' @author Mirko Di Febbraro
#' @examples \donttest{
#' library(raster)
#' raster(system.file("exdata/prediction_ground_thinned.gri",
#' package="EcoPast"))->prediction_ground_thinned
#' raster(system.file("exdata/prediction_ground.gri", package="EcoPast"))->prediction_ground
#' s1<-stack(prediction_ground, prediction_ground_thinned)
#' plot(s1)
#'
#' summary(rasterToPoints(s1))
#'
#' s2<-remove.NAs.stack(s1)
#'
#' summary(rasterToPoints(s2))
#' }

remove.NAs.stack<-function(rast.stack){
  nom<-names(rast.stack)
  test1<-calc(rast.stack, fun=sum)
  test1[!is.na(test1)]<-1
  test2<-rast.stack*test1
  test2<-stack(test2)
  names(test2)<-nom
  return(test2)
}
