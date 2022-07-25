#' @title A function to plot the \code{\link{minosse}} results
#' @description This function provides two maps: one with a species occurrence probability and another with the binarized geographic range
#' @usage minosse.plot(minosse_res,th_num=3,pcex=1.5,cont_pol=NULL,title=NULL)
#' @param minosse_res  the output of \code{\link{minosse.target}} function or the \code{minosse_res} element list in the \code{minosse} function output.
#' @param th_num Numeric. The Regression Kriging prediction map binarization threshold value index as reported in the \code{optimal.thresholds} function in the PresenceAbsence package. Default 3 = MaxSens+Spec,	maximizes (sensitivity+specificity)/2.
#' @param pcex Numeric. The graphic parameter indicating the size of points showing target species fossil localities. Default 1.5.
#' @param cont_pol The spatial polygon object showing the land or sea contour on the maps. Default \code{NULL}.
#' @param title The title on the top of the plot. If \code{NULL}, then the name of the species is used. Default \code{NULL}.
#' @importFrom sp sp.polygons coordinates
#' @importFrom gridExtra grid.arrange
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom grid textGrob gpar
#' @importFrom lattice levelplot
#' @importFrom raster cellStats
#' @importFrom latticeExtra layer_
#' @importFrom lattice panel.xyplot
#' @importFrom rasterVis rasterTheme
#' @importClassesFrom sp SpatialPoints SpatialPointsDataFrame
#' @importClassesFrom raster Raster RasterLayer
#' @export
#' @details When running this function for the first time in your workspace it is necessary to firstly load \code{ggplo2}.
#' @return A plot showing the occurrence probability (on the left) and the threshold-based geographic range (on the right) maps of target species along with its fossil localities.
#' @author Francesco Carotenuto, francesco.carotenuto@unina.it
#' @examples
#'   \dontrun{
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
#'   minosse.plot(mam)
#'
#'   }



minosse.plot<-function(minosse_res,th_num=3,pcex=1.5,cont_pol=NULL,title=NULL){
     #library(RColorBrewer)
     #library(lattice)
     #library(latticeExtra)
     #library(rasterVis)
     #library(gridExtra)
     #library(raster)
     #library(grid)

  if(any(names(minosse_res)%in%c("minosse_poly","minosse_res","sig_species"))){
    minosse_res$minosse_res->minosse_res
  } else minosse_res->minosse_res

  sp::coordinates(minosse_res$pts)->occs
  as.data.frame(occs)->occs
  rownames(occs)<-NULL
  prob<-minosse_res$prob_map
  binmap<-minosse_res$bin_maps[[th_num]]
  raster::crs(prob)->ras_crs
  
  minosse_col_prob<-rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(1000))
  my.theme.prob <-rasterVis::rasterTheme(region = minosse_col_prob)
  my.theme.bin <- rasterVis::rasterTheme(region = c("gray88",minosse_col_prob[1000]))
  
  
  if(!is.null(cont_pol)) spatial_mask<-latticeExtra::layer(sp.polygons(cont_pol,color="black",fill=NA),data=list(cont_pol=cont_pol)) else spatial_mask<-NULL
  my.at <- seq(from = raster::cellStats(prob, "min"),
               to = raster::cellStats(prob, "max"),
               length.out = 100 + 1)
  
  if(is.null(cont_pol)) {
    prob_plot<-rasterVis::levelplot(prob,par.settings = my.theme.prob,at=my.at,margin=FALSE,
                         main="Probability map",xlab="Longitude (meters)", ylab="Latitude (meters)")
    bin_plot<-rasterVis::levelplot(binmap,par.settings = my.theme.bin,at=my.at,margin=FALSE,
                        main="Geographic range",xlab = "Longitude (meters)", ylab= "Latitude (meters)")+
      latticeExtra::layer(panel.xyplot(x=occs[,1],y=occs[,2], pch=20, cex=pcex, col="lightskyblue1"),data=list(occs=occs,pcex=pcex))
  } else { if(!is.null(spatial_mask)) {
    prob_plot<-rasterVis::levelplot(prob,par.settings = my.theme.prob,at=my.at,margin=FALSE,
                                    main="Probability map",xlab="Longitude (meters)", ylab="Latitude (meters)")+
      spatial_mask
    bin_plot<-rasterVis::levelplot(binmap,par.settings = my.theme.bin,margin=FALSE,
                                   main="Geographic range",xlab = "Longitude (meters)", ylab= "Latitude (meters)")+
      latticeExtra::layer(panel.xyplot(x=occs[,1],y=occs[,2], pch=20, cex=pcex, col="lightskyblue1"),data=list(occs=occs,pcex=pcex))+
      spatial_mask
    } else {
    prob_plot<-rasterVis::levelplot(prob,par.settings = my.theme.prob,at=my.at,margin=FALSE,
                         main="Probability map",xlab="Longitude (meters)", ylab="Latitude (meters)")
    bin_plot<-rasterVis::levelplot(binmap,par.settings = my.theme.bin,margin=FALSE,
                        main="Geographic range",xlab = "Longitude (meters)", ylab= "Latitude (meters)")+
      latticeExtra::layer(panel.xyplot(x=occs[,1],y=occs[,2], pch=20, cex=pcex, col="lightskyblue1"),data=list(occs=occs,pcex=pcex))
    }
  }

  title1=grid::textGrob(minosse_res$nam, gp=grid::gpar(fontface="italic"))
  if(is.null(title))  gridExtra::grid.arrange(prob_plot,bin_plot, layout_matrix=matrix(c(1,2), 1, 2, byrow=TRUE), top=title1) else gridExtra::grid.arrange(prob_plot,bin_plot, layout_matrix=matrix(c(1,2), 1, 2, byrow=TRUE), top=title) 
  }
