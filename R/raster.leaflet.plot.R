#' @title Plot a \code{RasterLayer} or \code{RasterStack} object through the leaflet functionalities
#' @usage raster.leaflet.plot(RASTER.LAYER,PROJECT=FALSE,DOMAIN=NULL,PALETTE=NULL,
#' POSITION=NULL,REVERSE=TRUE,OPACITY=1,colNA=NA)
#' @description The function uses the \code{leaflet} package to plot raster maps, also allowing zoom and navigation functionalities.
#' @param RASTER.LAYER a \code{RasterLayer} or \code{RasterStack} object to be plotted.
#' @param PROJECT a logical vector indicating if the raster map must be projected to EPSG:3857.
#' @param DOMAIN a numerical vector indicating the range of raster values to be plotted.
#' @param PALETTE one of the palettes from PAL1 to PAL28.
#' @param POSITION a character vector indicating legend position ('topright' is the default).
#' @param REVERSE a logical vector indicating if to reverse the color palette.
#' @param OPACITY a numerical vector ranging from 0 to 1 indicating the opacity level.
#' @param colNA the color of NA cells.
#' @export
#' @import leaflet
#' @importFrom RColorBrewer brewer.pal
#' @importFrom raster projectRaster projectExtent values nlayers
#' @importFrom grDevices terrain.colors
#' @return A data.frame with the reduced points.
#' @author Mirko Di Febbraro
#' @examples
#' \donttest{
#'   library(raster)
#'   data(lgm)
#'   raster(system.file("exdata/prediction_ground.gri",
#'          package="EcoPast"))->prediction_ground
#'   raster.leaflet.plot(prediction_ground)

#' }

raster.leaflet.plot<-function(RASTER.LAYER,
                              PROJECT=FALSE,
                              DOMAIN=NULL,
                              PALETTE=NULL,
                              POSITION=NULL,
                              REVERSE=TRUE,
                              OPACITY=1,
                              colNA=NA){
  #require(leaflet)

  if(PROJECT){
  epsg3857<-"+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
  cat("Projecting Rasters...")
  RASTER.LAYER<-projectRaster(RASTER.LAYER,
                              projectExtent(RASTER.LAYER, crs = sp::CRS(epsg3857)),
                              method="ngb",
                              progress="text")
  cat("Done!")}


  if(is.null(PALETTE))PALETTE<-terrain.colors(100)
  if(is.null(POSITION))POSITION<-"topright"

  if(PROJECT){tile.opacity<-1}else{tile.opacity<-0}

  if(is.null(DOMAIN))DOMAIN<-values(RASTER.LAYER)

  if(class(RASTER.LAYER)=="RasterLayer"){
    colorNumeric(palette=PALETTE,
                 domain=DOMAIN,
                 na.color = colNA,
                 reverse=REVERSE)->COLORS

    m<-leaflet()%>%
      #addTiles(options=tileOptions(opacity=tile.opacity))%>%
        addRasterImage(RASTER.LAYER,
                       colors = COLORS,
                       project = FALSE,
                       opacity=OPACITY,
                       maxBytes = 16 * 1024 * 1024)%>%
        addLegend(position=POSITION,
                  pal=COLORS,
                  values=DOMAIN,
                  opacity=1,
                  title=names(RASTER.LAYER))
    m
  }else{
    base<-c("m<-leaflet()%>%")
        #addTiles(options=tileOptions(opacity=tile.opacity))%>%")


    X<-lapply(1:nlayers(RASTER.LAYER), function(j){
      paste(paste("addRasterImage(RASTER.LAYER[[",j,"]],
               colors =   colorNumeric(palette=PALETTE,
                                       domain=values(RASTER.LAYER[[",j,"]]),
                                       na.color = colNA,
                                       reverse=REVERSE),
               project = FALSE,
               opacity=OPACITY,
               maxBytes = 16 * 1024 * 1024,
              group=names(RASTER.LAYER[[",j,"]]))", sep=""),"%>% \n")})

    X<-do.call(paste, X)

    controls<-"addLayersControl(baseGroups=names(RASTER.LAYER),
                   options = layersControlOptions(collapsed=FALSE))"

    final<-paste(base, X, controls, sep="\n")
    eval(parse(text=final))
    m}
  return(m)
}




################################################################################
################################################################################

#' @title Import some user-defined palettes
#' @usage import.palettes()
#' @description The function imports some user-defined colour palettes.
#' @export
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @examples \dontrun{
#' import.palettes()
#' }

import.palettes<-function(){
  PAL1<-c("#F2F1A2", # Palette giallo, rosso, fucsia e viola
          "#FFFF00",
          "red",
          "#FC0345",
          "#B007ED",
          "#071DAD")


  PAL2<-c("#c2523c", # Palette marrone, giallo, verde e blu scuro
          "#eda113",
          "#ffff00",
          "#00db00",
          "#20998f",
          "#0b2c7a")

  PAL3<-c("#aff0e9", # Palette marrone, giallo, verde e blu scuro
          "#ffffb3",
          "#008040",
          "#fcba03",
          "#800000",
          "#69300d",
          "#ababab",
          "#fffcff")

  PAL4<-c("#fff700",
          "#f80702",
          "#85427b",
          "#4756c8",
          "#022979",
          "#26100b")


  sequential<-c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu",
                "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")

  X<-lapply(1:length(sequential), function(x)paste("PAL",
                                                   x+4,
                                                   "<<-RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info['",
                                                   sequential[x],
                                                   "',]$maxcolors,'",
                                                   sequential[x],"')", sep=""))
  for(x in 1:length(X))eval(parse(text=X[[x]]))

  diverging<-c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")

  X<-lapply(1:length(diverging), function(x)paste("PAL",
                                                  x+18,
                                                  "<<-RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info['",
                                                  diverging[x],
                                                  "',]$maxcolors,'",
                                                  diverging[x],"')", sep=""))
  for(x in 1:length(X))eval(parse(text=X[[x]]))

  PAL28<-c("#ff0000",
           "#ff3500",
           "#ff6b00",
           "#ffa100",
           "#ffd600",
           "#f2ff00",
           "#bcff00",
           "#87ff00",
           "#51ff00",
           "#1bff00",
           "#00ff1a",
           "#00ff50",
           "#00ff86",
           "#00ffbb",
           "#00fff1",
           "#00d7ff",
           "#00a2ff",
           "#006cff",
           "#0036ff",
           "#0000ff")

  return(c("diverging" ,"PAL1","PAL10","PAL11","PAL12","PAL13","PAL14","PAL15","PAL16" ,"PAL17",
    "PAL18","PAL19","PAL2","PAL20","PAL21","PAL22","PAL23","PAL24","PAL25","PAL26","PAL27",
    "PAL28","PAL3","PAL4","PAL5","PAL6","PAL7","PAL8","PAL9","sequential"))
}

################################################################################

#' @title Show the palettes from PAL1 to PAL28
#' @usage show.custom.palettes()
#' @description The function shows an example of the colour palettes imported through \code{import.palettes()}
#' @export
#' @importFrom graphics rect par plot text
#' @examples \dontrun{
#' show.custom.palettes()
#' }

show.custom.palettes<-function(){
  pal <- function(col, border = "transparent", ...){
    n <- length(col)
    plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
         axes = FALSE, xlab = "", ylab = "", ...)
    rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
  }

  ini.par<-par()
  par(mar=c(1,1,1,1), mfrow=c(4,7))
  lapply(paste("PAL", 1:28, sep=""), function(x){
    eval(parse(text=paste("pal(colorRampPalette(",x,")(100))",sep="")))
    text(0.5, 0.5, x)
  })
  suppressWarnings(par(ini.par))
}
