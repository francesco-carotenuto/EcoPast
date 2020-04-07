#' @title A function to plot the \code{\link{minosse}} results
#' @description This function provides two maps: one with a species occurrence probability and another with the binarized geographic range
#' @usage minosse.plot(minosse_res,th_num=3,pcex=1.5,cont_pol=NULL,title=NULL)
#' @param minosse_res  the output of \code{\link{minosse.target}} function or the \code{minosse_res} element list in the \code{minosse} function output.
#' @param th_num Numeric. The Regression Kriging prediction map binarization threshold value index as reported in the \code{optimal.thresholds} function in the PresenceAbsence package. Default 3 = MaxSens+Spec,	maximizes (sensitivity+specificity)/2.
#' @param pcex Numeric. The graphic parameter indicating the size of points showing target species fossil localities. Default 1.5.
#' @param cont_pol The spatial polygon object showing the land or sea contour on the maps. Default \code{NULL}.
#' @param title The title on the top of the plot. If \code{NULL}, then the name of the species is used. Default \code{NULL}.
#' @importFrom ggplot2 ggplot labs scale_fill_gradientn theme_minimal theme element_text scale_fill_gradient unit Stat
#' @importFrom gridExtra grid.arrange
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom ggspatial layer_spatial StatSpatialRasterDf GeomSpatialRaster StatSpatialRasterAnnotation StatSpatialRaster
#' @export
#' @details When running this function for the first time in your workspace it is necessary to firstly load \code{ggplo2}.
#' @return A plot showing the occurrence probability (on the left) and the threshold-based geographic range (on the right) maps of target species along with its fossil localities.
#' @author Francesco Carotenuto, francesco.carotenuto@unina.it
#' @examples
#'   \donttest{
#'
#'   library(raster)
#'   data(lgm)
#'   raster(system.file("exdata/prediction_ground.gri", package="EcoPast"))->prediction_ground
#'
#'   mam<-minosse(dat=lgm,species.name="Mammuthus_primigenius",domain=NULL,
#'   time.overlap=0.95,prediction.ground=prediction_ground,crop.by.mcp=FALSE,
#'   constrain.predictors=FALSE,temporal.tolerance=NULL,
#'   coc_by="locality",n.sims=10,projection="laea",lon_0 = 85,
#'   lat_0 = 45,seed=625)
#'
#'   library(ggplot2)
#'   minosse.plot(mam)
#'
#'   }



minosse.plot<-function(minosse_res,th_num=3,pcex=1.5,cont_pol=NULL,title=NULL){
  #   library(RColorBrewer)
  #   library(lattice)
  #   library(ggplot2)
  #   library(ggspatial)
  #   library(gridExtra)
  #   library(raster)


if(length(which((.packages())=="ggplot2"))==0) stop("Please load ggplot2 package the first time you run this function")
  if(any(names(minosse_res)%in%c("minosse_poly","minosse_res","sig_species"))){
    minosse_res$minosse_res->minosse_res
  } else minosse_res->minosse_res

  minosse_res$pts->occs
  prob<-minosse_res$prob_map
  binmap<-minosse_res$bin_maps[[th_num]]
  raster::crs(prob)->ras_crs
  if(class(try(raster::projectRaster(prob,crs=CRS("+proj=longlat +datum=WGS84")),silent=TRUE))=="try-error") {
    raster::crs(prob)<-NA
    raster::crs(binmap)<-NA
    raster::crs(occs)<-NA
  }
  minosse_col<-rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(1000))
  if(is.null(cont_pol)) spatial_mask<-cont_pol else spatial_mask<-ggspatial::layer_spatial(data=cont_pol,color="black",fill=NA)
  as.numeric(raster::values(binmap))->raster::values(binmap)
  prob_plot<-ggplot2::ggplot()+
    ggspatial::layer_spatial(data=prob,interpolate=FALSE,dpi = 600) +
    ggplot2::coord_sf(datum = as.character(ras_crs))+
    ggplot2::labs(title = "Probability map\n", x = "Longitude (meters)", y = "Latitude (meters)")+
    ggplot2::scale_fill_gradientn(name = "",colours=minosse_col,na.value="#FFFFFF00")+
    ggplot2::theme_minimal()+
    ggplot2::theme(plot.title = element_text(hjust = 0.5,size=8),legend.text=element_text(size=5),
          axis.title.x=element_text(size=7),
          axis.title.y=element_text(size=7),
          legend.key.width=unit(0.2,"cm"),
          axis.text=element_text(size=6))+
    spatial_mask

  bin_plot<-ggplot2::ggplot()+
    ggspatial::layer_spatial(data=binmap,interpolate=FALSE,dpi=600)  +
    #ggplot2::scale_fill_manual(values=c("gray88",minosse_col[1000]))+
    ggplot2::scale_fill_gradient(name = "",low="gray88",high=minosse_col[1000],na.value="#FFFFFF00")+
    ggplot2::labs(title = "Geographic range\n", x = "Longitude (meters)", y = "Latitude (meters)")  +
    ggspatial::layer_spatial(data=occs,show.legend=FALSE,size=pcex,color="lightskyblue1")+
    ggplot2::coord_sf(datum = as.character(ras_crs))+
    ggplot2::theme_minimal()+
    ggplot2::theme(plot.title = element_text(hjust = 0.5,size=8),legend.text=element_text(size=5),
                 axis.title.x=element_text(size=7),
                 axis.title.y=element_text(size=7),
                 legend.key.width=unit(0.2,"cm"),
                 axis.text=element_text(size=6))+
    spatial_mask
  if(is.null(title))  gridExtra::grid.arrange(prob_plot,bin_plot, layout_matrix=matrix(c(1,2), 1, 2, byrow=TRUE), top=minosse_res$nam) else gridExtra::grid.arrange(prob_plot,bin_plot, layout_matrix=matrix(c(1,2), 1, 2, byrow=TRUE), top=title)
  }
