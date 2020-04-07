#' @importFrom automap autofitVariogram
#' @importClassesFrom raster Raster RasterStack
#' @importClassesFrom GSIF gstatModel
#' @importMethodsFrom plotKML summary

minosse.core<-Vectorize(function(dtset,
                                 mappa_ras,
                                 projection.ground,
                                 seed,
                                 predict_by=c("kriging"),
                                 configs) {

  # library(GSIF)
  # library(SDMTools)
  # library(sp)
  # library(gstat)
  # library(plotKML)
  # library(raster)
  # library(automap)

  if(is.null(projection.ground)) projection.ground<-mappa_ras else projection.ground<-projection.ground




  kriging_prediction<-function(x,mappa_ras,projection.ground,seed,configs){
    gimme.mod<-function(x,form,mappa_ras){
      form<-as.formula(paste("value~", paste(names(mappa_ras), collapse="+")));
      x[[1]]->x2
      ov<-sp::over(x[[1]],mappa_ras)
      x2@data<-cbind(x2@data,ov)
      x2[which(complete.cases(x2@data)),]->x2
      df_x2<-as.data.frame(x2)
      m<-minosse_fit.regModel(configs=configs,form, rmatrix=df_x2, mappa_ras,
                              method=configs$kriging$method,rvgm=NULL)
      vg_form<-as.formula(paste(parse(text=names(m@sp)[ncol(m@sp)]),"1",sep="~"))
      resid_vg<-automap::autofitVariogram(vg_form,m@sp[,ncol(m@sp)])
      if(class(try(mod<-minosse_fit.gstatModel(configs=configs,x[[1]], form,mappa_ras,method=configs$kriging$method,vgmFun=as.character(resid_vg$var_model[2,1])),silent = TRUE))=="try-error") {
        mod<-NA} else  mod<-mod
      return(mod)
    }

    set.seed(seed)

    if(is.na(mod<-gimme.mod(x,form,mappa_ras))) {
      pred<-NA;
      kr_MINOSSE<-NA;
      kr_auc<-NA;
      data_for_th<-NA;
      kr_th<-NA;
      RK_RMSE<-NA;
      RK_tvar<-NA
      xgb_RMSE<-NA;
      xgb_R2<-NA;
      MLA_mod<-NA} else {
        mod<-mod
        pred<-predict(mod,projection.ground,zmin=0,zmax=1,predict.method="RK")
        kr_MINOSSE<-raster(pred@predicted[,"value"])
        raster::crs(kr_MINOSSE)<-raster::crs(projection.ground)
        which(!is.na(extract(kr_MINOSSE,x[[2]])))->valid_points
        if(class(try(kr_auc<-SDMTools::auc(x[[2]][valid_points,]$value,extract(kr_MINOSSE,x[[2]][valid_points,]))))=="try-error") kr_auc<-NA else kr_auc<-kr_auc
        RK_RMSE<-summary(pred)$RMSE
        RK_tvar<-summary(pred)$tvar
        MLA_mod<-mod
        xgb_RMSE<-MLA_mod@regModel$results$RMSE
        xgb_R2<-MLA_mod@regModel$results$Rsquared
      }
    MINOSSE_list<-list(list(kr_MINOSSE=kr_MINOSSE,kr_auc=kr_auc,RK_RMSE=RK_RMSE,RK_tvar=RK_tvar,xgb_RMSE=xgb_RMSE,xgb_R2=xgb_R2,MLA_mod=MLA_mod))
    return(MINOSSE_list);
  }



  switch(as.character(predict_by),
         kriging={cat("\n","Performing RK prediction");
           kriging_prediction(dtset,as(mappa_ras,"SpatialPixelsDataFrame"),as(projection.ground,"SpatialPixelsDataFrame"),seed,configs=configs);
         });
}, "predict_by")
