#' @title Shift points to the nearest non-NA cell of a base raster within a maximum distance
#' @usage fix.coastal.points(data,xy.cols,r,ncell,clip.by.poly=FALSE,poly,occ.desaggregation=FALSE)
#' @description The function moves points to the nearest non-NA cell of a base raster within a maximum distance.
#' @param data a matrix or data.frame containing the points coordinates.
#' @param xy.cols a numeric vector indicating the IDs of the x-y columns.
#' @param r a \code{RasterLayer} object to be used as base raster map.
#' @param ncell an integer indicating the maximum number of cells away from the points where they can be shifted. If the nearest non-NA cell fall farther than this distance from a give point, this will be dropped.
#' @param clip.by.poly a logical vector indicating if the base raster must be clipped along the boundaries of a polygon.
#' @param poly a \code{SpatialPolygons} object used to clip the base raster map.
#' @param occ.desaggregation a logical vector indicating if to apply the \code{\link{occ.desaggregation.RASTER}} function.
#' @export
#' @importFrom raster extent extend crop mask extract adjacent cellFromXY xyFromCell
#' @importFrom sp SpatialPoints spDists
#' @importFrom tcltk tkProgressBar setTkProgressBar
#' @return A matrix including the shifted points
#' @author Mirko Di Febbraro
#' @examples
#' \donttest{
#' library(raster)
#'raster(system.file("exdata/prediction_ground.gri", package="EcoPast"))->prediction_ground
#'data(lgm)
#'
#'data_clupus<-subset(lgm, spec=="Canis_lupus")
#'
#'coordinates(data_clupus)<-~x+y
#'proj4string(data_clupus)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#'data_clupus<-spTransform(data_clupus, proj4string(prediction_ground))
#'data_clupus<-as.data.frame(data_clupus)
#'
#'data_mammuth_reduced<-fix.coastal.points(data_clupus, xy.cols = 4:5,
#'r = prediction_ground, ncell = 1, occ.desaggregation = FALSE)
#' }

fix.coastal.points<-function(data,
                             xy.cols,
                             r,
                             ncell,
                             clip.by.poly=FALSE,
                             poly,
                             occ.desaggregation=FALSE)
{

  init<-nrow(data)

  if(clip.by.poly){
    ex<-extent(SpatialPoints(data[,xy.cols]))
    ex<-extend(ex, c(2,2))
    poly1<-crop(poly, ex)
    r1<-mask(crop(r, poly1), poly1, progress="text")
  }

  cc<-extract(r, SpatialPoints(data[,xy.cols]))
  vv<-which(is.na(cc))

  if(length(vv)!=0){
    if(clip.by.poly){
      cc<-extract(r1, SpatialPoints(data[,xy.cols]))
      vv1<-which(is.na(cc))
      vv<-sort(union(vv, vv1))
    }

    m<-matrix(rep(rep(1, ncell+1), 2*ncell+1), ncol=2*ncell+1, nrow=2*ncell+1)
    m[ncell+1,ncell+1]<-0

    shifted<-removed<-matrix(ncol=2)
    colnames(shifted)<-colnames(removed)<-colnames(data)[xy.cols]

    #pb_4<-tkProgressBar(min=0, max=length(vv), label="Progress (0 %)", width=600, title="MOVING OCCURRENCES")

    for(k in vv){
      #setTkProgressBar(pb_4,
      #                 value=which(vv==k),
      #                 label=paste("Progress (",
      #                             as.character(round(seq(1/length(vv)*100,100,length.out=length(vv)))[which(vv==k)]),
      #                             " %)",
      #                             sep=""))
      focal<-rbind(data[k,xy.cols])
      valori<-r[adjacent(r, cellFromXY(r, focal), directions=m, pairs=FALSE)]
      new.coords<-rbind(xyFromCell(r, adjacent(r, cellFromXY(r, focal), directions=m, pairs=FALSE))[!is.na(valori),])
      if(nrow(new.coords)!=0){
        new1<-spDists(new.coords, focal)
        minimo<-which.min(new1)
        data[k,xy.cols]<-new.coords[minimo[sample(1:length(minimo),1)],]
        shifted<-rbind(shifted, focal)
      } else {
        data[k,xy.cols]<-c(NA, NA)
        removed<-rbind(removed, focal)
      }
    }
    #close(pb_4)
    final1<-data[!is.na(data[,xy.cols[1]]),]
    if(occ.desaggregation&nrow(final1)>0){
      final2<-occ.desaggregation.RASTER(final1, xy.cols, r, F)
      } else {final2<-final1}
    if(occ.desaggregation==FALSE)final2<-final1
    cat("\n", "initial sample size: ", init, sep="")
    cat("\n","shifted points: ", nrow(shifted)-1, sep="")
    cat("\n","removed points: ", nrow(removed)-1, sep="")
    if(occ.desaggregation)cat("\n","removed duplicates: ", nrow(final1)-nrow(final2), sep="")
    cat("\n","final sample size: ", nrow(final2), sep="")
    cat("\n")
    return(final2)
  }
  if(length(vv)==0){
    if(occ.desaggregation){
      final2<-occ.desaggregation.RASTER(data, xy.cols, r, F)
      return(final2)
    }else return(data)
  }
}
