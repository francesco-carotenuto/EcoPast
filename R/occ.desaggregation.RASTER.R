#' @title Remove duplicated points falling into the same raster cell
#' @usage occ.desaggregation.RASTER(df, colxy, rast, plot=TRUE)
#' @description The function removes duplicated points falling into the same raster cell.
#' @param df a matrix or data.frame containing the points coordinates.
#' @param colxy a numeric vector indicating the IDs of the x-y columns.
#' @param rast a \code{RasterLayer} object to be used as base raster map.
#' @param plot a logical vector indicating if a plot must be displayed.
#' @export
#' @importFrom raster cellFromXY
#' @importFrom graphics points
#' @importFrom grDevices dev.new
#' @return A matrix including the reduced points
#' @author Mirko Di Febbraro
#' @examples \donttest{
#' library(raster)
#' raster(system.file("exdata/prediction_ground.gri", package="EcoPast"))->prediction_ground
#' data(lgm)
#'
#' data_mammuth<-subset(lgm, spec=="Mammuthus_primigenius")
#'
#' coordinates(data_mammuth)<-~x+y
#' proj4string(data_mammuth)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#' data_mammuth<-spTransform(data_mammuth, proj4string(prediction_ground))
#' data_mammuth<-as.data.frame(data_mammuth)
#'
#' data_mammuth_reduced<-occ.desaggregation.RASTER(df=data_mammuth, colxy=4:5,
#' rast=prediction_ground, plot=TRUE)
#'
#' }

occ.desaggregation.RASTER<-function(df, colxy, rast, plot=TRUE){
  #require(raster)
  df_ini<-df
  obj<-cellFromXY(rast, df[,colxy])
  if(any(is.na(obj))){stop("no NA admitted in species occurrences!")}
  obj1<-split(obj, obj)
  l<-sapply(obj1, length)
  if(max(l)>1){
    l1<-l[l>1]
    for(j in 1:length(l1)){
      w<-which(obj%in%as.numeric(names(l1[j])))
      found<-df[w,]
      s<-sample(w,1)
      df[w[w!=s],colxy]<-NA
    }
    df_final<-df[!is.na(df[,colxy[1]]),]
    if(plot==TRUE){
      #x11()
      dev.new()
      plot(df_ini[,colxy],main="distribution of occurences",sub=paste("# initial (black):",nrow(df_ini)," | # kept (red): ",nrow(df_final)),pch=19,col="black",cex=0.5)
      points(df_final[,colxy],pch=19,col="red",cex=0.2)
    }
    return(df_final)
  }
  if(max(l)==1)return(df)
}
