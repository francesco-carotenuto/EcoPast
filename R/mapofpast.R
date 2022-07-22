#' @title Creates prediction ground and paleo dem
#' @description This function uses Scotese and Wright (2018) PALEOMAPS to generate the prediction ground and paleo dems of the requested temporal interal
#' @usage mapofpast(Ma,domain,ext=NULL,projection,pts=NULL)
#' @param Ma An integer value comprised in between 0 and 540
#' @param domain Character. Eithter \code{"land"}, for terrestrial species, or \code{"sea"}, for marine species.
#' @param ext Optional. A numeric vector with the coordinates (in wgs84 format) of the rectangle's vertices enclosing the domain. It must be c(min longitude, max longitude, min latitude, max latitude). See details
#' @param projection Character. This is the euqual-area projection of the final maps. A CRS object in the proj4 format or either \code{"moll"} (Mollweide) or \code{"laea"} (Lambert Azimuthal equal area) projections. See details.
#' @param pts Optional. A dataframe with the x and y coordinates of individual fossil localities. If specified the mean center of the fossil record spatial distribution is considered as projection's centre. See details.
#' @export
#' @importFrom utils unzip
#' @importFrom raster drawPoly
#' @importFrom grDevices dev.off
#' @importFrom raster plot
#' @importFrom sp plot
#' @importFrom graphics par
#' @details By setting a numeric vector for the argumento \code{ext}, the user can select a portion of the choosen domain.
#' If the argument is NULL, when running the function, the user is asket to draw the polygon encolosing the portion
#' of the domain to be used in all the predictions. In this case, the user will draw a polygon's vertices by clicking
#' on the map with the lef button of the mouse. Then, by pressing Esc, the function will crop the portion of the domain
#' enclosed by the drawn polygon. If the argument \code{projection} is either \code{"laea"} or \code{"moll"}, the projection's centre is computed
#' as the extent's centre either if set in the \code{ext} argument either if set by drawing a polygon. In the case the argument \code{pts}
#' is specified, then the mean centre of the fossil record distribution is used as projection centre.
#' @return A list of two rasters to be used in the \code{\link{minosse.data}} function: the first is the prediction domain to be used in the \code{prediction.ground} argument, and the second is a paleodem to be used in the \code{abiotic.covs} argument.
#' @author Francesco Carotenuto, francesco.carotenuto@unina.it
#' @examples
#'   \dontrun{
#'   mapofpast(Ma=540,domain="land",ext=c(-50,50,-30,30),projection="laea",pts=NULL)
#'   }



mapofpast<-function(Ma,domain,ext=NULL,projection,pts=NULL){

 
  if (!requireNamespace("chronosphere", quietly = TRUE)) {
    stop("Package \"chronosphere\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if(is.null(Ma)) stop("Please select a temporal interval between 0 and 540 Mya")
  if(is.null(projection)) stop("Please set a projected coordinate system. Available choices are: laea, moll or a CRS object")

 # library(raster)
 # library(ncdf4)
 # library(sp)
 # library(chronosphere)
 # library(utils)

  print("Warning! If you are asking for deep time world map you need to rotate the geographic coordinates of your fossil dataset. If you didn't rotate your coordinates, please use the chronosphere package")

  Ma->time.interval
  paste("_",time.interval,"Ma",sep="")->time.interval
  mps<-utils::unzip(system.file("exdata/paleomaps.zip", package="EcoPast"))
  original_dem<-raster::raster(mps[as.numeric(grep(time.interval,mps))])
  strsplit(mps[grep(time.interval,mps)],"Ma.nc")[[1]]->selected_time
  strsplit(selected_time,"_")[[1]][length(strsplit(selected_time,"_")[[1]])]->selected_time
  print(paste("The temporal interval you asked for is:",selected_time,"Mya"))
  print("Generating paleo world map. Citation for map source: Scotese, C. R., & Wright, N. (2018). PALEOMAP Paleodigital Elevation Models (PaleoDEMS) for the Phanerozoic.")
  original_dem->final_map

  raster::plot(final_map)
  if(!is.null(pts)){
    sp::coordinates(pts)<-~x+y
    raster::crs(pts)<-CRS("+proj=longlat")
    points(pts,pch=19,col=2,cex=0.5)
  }

  if(is.null(ext)) {
    print("Please select the a polygon's vertices by clicking the left-button mouse and then press Esc")
    raster::drawPoly()->poly_mask
    dev.off()
    }
  if(is.numeric(ext)) {
    raster::extent(ext)->poly_mask
    as(poly_mask,"SpatialPolygons")->poly_mask
    sp::plot(poly_mask,add=TRUE,border="red")
    }
  if(domain=="land") {raster::mask(final_map,poly_mask)->final_map}
  raster::crop(final_map,extent(poly_mask))->final_map


  if(projection%in%c("laea","moll")) {
    if(!is.null(pts)) {
      mean(pts$x)->lon_0
      mean(pts$y)->lat_0
    }
    if(is.null(pts)) {
      mean(extent(final_map)[1:2])->lon_0
      mean(extent(final_map)[3:4])->lat_0

    }
  if(projection=="laea") my.prj<-CRS(paste("+proj=laea +lat_0=",lat_0, " +lon_0=",lon_0," +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs",sep=""))
  if(projection=="moll") my.prj<-CRS(paste("+proj=moll +lon_0=",lon_0," +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs",sep=""))
  }
  if(class(projection)=="CRS")  my.prj<-projection

  raster::projectRaster(final_map,crs=my.prj)->final_map

  if(domain=="sea") {
      final_map<0->final_map
    final_map[final_map==0]<-NA
      }
  if(domain=="land") {
    final_map->final_map
    prj_mask<-raster::rasterToPolygons(final_map,function(x)x>0,dissolve=TRUE)
    raster::mask(final_map,prj_mask)->final_map
    raster::crop(final_map,prj_mask)->final_map
    !is.na(final_map)->final_map
    final_map[final_map==0]<-NA
    }

  raster::projectRaster(original_dem,crs=crs(final_map))->original_dem
  raster::crop(original_dem,final_map)->original_dem
  par(mfrow=c(1,2))
  raster::plot(final_map,col="lightgreen",legend=FALSE)
  if(!is.null(pts)) {
    points(spTransform(pts,CRSobj=raster::crs(final_map)),pch=19,col="red",cex=0.5)
  }
  chronosphere::mapplot(original_dem,col="earth",box = TRUE, legend = TRUE)

  list(ground=final_map,dem=original_dem)->mp
  return(mp)
}




