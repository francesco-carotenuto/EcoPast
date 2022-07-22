#' @title The target species geographic range computation.
#' @description This function reconstructs the target species geographic range and probability distribution raster maps.
#' @usage minosse.target(resp,predictors,bkg="presence",min.bkg=NULL,n.sims=10,
#' n.folds=1,sampling.by.distance=TRUE,n.sims.clusters=NULL,seed=NULL)
#' @param resp The target species occurrence record. Usually the first element list in the \code{\link{minosse.data}} output.
#' @param predictors The predictor species raster stack. Usually the second element list in the \code{minosse.data} output.
#' @param bkg The number of pseudo absences to be simulated for each predictor species. If \code{"presence"}, \code{minosse.data} generates as many pseudo absences the presences in each species.
#' @param min.bkg Numeric. If \code{bkg.predictors} is set to \code{"presence"}, this is the minimum number of pseudo absences to simulate if a species occurrence number is below this value.
#' @param n.sims Numeric. The number of pseudo absences simulations (see details).
#' @param n.folds Numeric. The number of folds for AUC-based cross-validation of Regression Kriging predictions. Default 1.
#' @param sampling.by.distance Logical. If \code{TRUE} pseudo absences are simulated with an intensity proportional to the distance to the presence data. If \code{FALSE} a pure spatial random distribution is simulated.
#' @param n.sims.clusters Numeric or \code{NULL}. The number of machine cores to use when setting multiple pseudo absences simulations. Default \code{"automatic"}. Default is \code{NULL}.
#' @param seed Numeric. The \code{seed} to use for experiment replication.
#' @importFrom raster distanceFromPoints stack
#' @importFrom dismo randomPoints
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach
#' @importFrom foreach %dopar% foreach
#' @importFrom PresenceAbsence optimal.thresholds
#' @importClassesFrom raster Raster RasterStack
#' @importClassesFrom plotKML SpatialPredictions
#' @importMethodsFrom raster mean
#' 
#' @export
#' @details For Machine Learning Algorithms Barbet & Masin (https://doi.org/10.1111/j.2041-210X.2011.00172.x) suggest using as many pseudo absences as the number of the presences. For very poor occurrence datasets
#' (with less than 1000 presences), they suggest performing multiple (10 at least) pseudo absences simulations and then averaging the derived predictions. \code{minosse.target} function allows fully complying with these authors' suggestions.
#' @return A list of 7 objects is returned. The first item of the list is a list of 12 geographic range maps obtained by binarizing Regression Kriging output by means of an equivalent number of threshold values. The second element
#' of the list is the probability map yielded by Regression Kriging. The third element is a list containing eXtreme Gradient Boosting model and validation statistics for each of the replicates of the pseudo absences simulations (\code{n.sims}, see above).
#' The fourth element of the list shows the Regression Kriging binarization threshold values. The fifth element (validation) shows the Regression Kriging prediction (RMSE and total explained variance) and xgboost model validations (squared R and RMSE).
#' The two last elements show the target species name and the occurence record.
#' @author Francesco Carotenuto, francesco.carotenuto@unina.it
#' @examples
#'   \dontrun{
#'   library(raster)
#'   data(lgm)
#'   raster(system.file("exdata/prediction_ground.gri", package="EcoPast"))->prediction_ground
#'
#'   minosse_dat<-minosse.data(obj=lgm,species_name="Mammuthus_primigenius",
#'   domain=NULL,coc.by="locality",min.occs=3,abiotic.covs=NULL,
#'   combine.covs=FALSE,reduce_covs_by="pca",covs_th=0.95,c.size="mean",
#'   bkg.predictors="presence",min.bkg=100,sampling.by.distance=TRUE,
#'   crop.by.mcp=FALSE,prediction.ground=prediction_ground,projection=NULL,
#'   lon_0=NULL,lat_0=NULL,n.clusters="automatic",seed=625)
#'
#'   minosse_res<-minosse.target(resp=minosse_dat[[1]],predictors=minosse_dat[[2]],
#'   bkg="presence", min.bkg = 100,n.sims=10,n.folds=1,sampling.by.distance=TRUE,
#'   n.sims.clusters=NULL,seed=625)
#'
#'   }


minosse.target<-function(resp,
                         predictors,
                         bkg="presence",
                         min.bkg=NULL,
                         n.sims=10,
                         n.folds=1,
                         sampling.by.distance=TRUE,
                         n.sims.clusters=NULL,
                         seed=NULL){


  
  

    
  
  
     library(raster)
     library(dismo)
     library(parallel)
     library(doSNOW)
     library(doParallel)
     library(foreach)
     library(PresenceAbsence)



  length(getLoadedDLLs())
  Sys.getenv("R_MAX_NUM_DLLs", unset = 1000)
  configs<-list(kriging=list(method="xgboost",tuning.clusters=1,folds=5,repeats=5))
  mappa_ras<-predictors
  projection.ground<-NULL
  spatial.layer<-NULL
  predict_by<-"kriging"

  sim.bkg<-function(predict_by) {
    if(any(predict_by%in%c("kriging"))) {
      if(class(bkg)%in%c("Data.Frame","Matrix")) {
        backg<-bkg;
        dataset<-rbind(as.data.frame(resp[resp$value==1,])[,c("x","y","value")],backg)
      }

      if(is.numeric(bkg)){
        po_dis<-distanceFromPoints(mappa_ras[[1]], resp[resp$value==1,])
        po_dis[is.na(mappa_ras[[1]]),]<-NA
        dismo::randomPoints(po_dis,n=bkg,p=resp[resp$value==1,],prob=sampling.by.distance)->backg
        backg<-as.data.frame(backg)
        backg$value<-0
        dataset<-rbind(as.data.frame(resp[resp$value==1,])[,c("x","y","value")],backg)
        rownames(dataset)<-NULL
        sp::coordinates(dataset)<-~x+y
        raster::crs(dataset)<-raster::crs(resp)
      } else if(is.null(bkg)) {dataset<-resp[,"value"]}
      if(bkg=="presence") {
        po_dis<-distanceFromPoints(mappa_ras[[1]], resp[resp$value==1,])
        po_dis[is.na(mappa_ras[[1]]),]<-NA
        bkg_num<-length(resp[resp$value==1,])
        if(!is.null(min.bkg)) {if(bkg_num<min.bkg) bkg_num<-min.bkg} else bkg_num<-bkg_num
        backg<-dismo::randomPoints(po_dis,n=bkg_num,p=resp[resp$value==1,],prob=sampling.by.distance)
        backg<-as.data.frame(backg)
        backg$value<-0
        dataset<-rbind(as.data.frame(resp[resp$value==1,])[,c("x","y","value")],backg)
        rownames(dataset)<-NULL
        sp::coordinates(dataset)<-~x+y
        raster::crs(dataset)<-raster::crs(resp)
      }
    }
    if(exists("dataset")) return(dataset)
  }

  set.seed(seed)
  replicate(n.sims,sim.bkg(predict_by))->dataset_data

  if(any(predict_by%in%c("kriging"))) {
    dataset_data->dataset;
    dataset->original_dataset}

  data_subs<-function(dataset,n.folds){
    fold <- dismo::kfold(dataset,k = n.folds)
    as.data.frame(sapply(seq(1:n.folds), function(x)fold!=x))->DataSplitTable
    colnames(DataSplitTable) <- paste("RUN", 1:(n.folds),sep = "")
    DataSplitTable$Full<-rep(TRUE,nrow(DataSplitTable))
    return(DataSplitTable)
  }

  if(n.folds>1) {
    if(any(predict_by%in%c("kriging"))){
      set.seed(seed)
      lapply(dataset,function(x)data_subs(x,n.folds))->DataSplitTable
      dataset<-lapply(seq(1:n.sims),function(x){
        splitted_dataset<-apply(DataSplitTable[[x]],2,function(z){list(dataset[[x]][z==TRUE,],dataset[[x]][z==FALSE,])->splitted_dataset;
          return(splitted_dataset)})
        return(splitted_dataset)})
    } else dataset<-NULL

  } else if(n.folds==1) {
    if(any(predict_by%in%c("kriging"))){
      lapply(dataset,function(x)list(x,x))->dataset
    } else dataset<-NULL

  }

  ################################   Predictions of the response ###################################
  ################################   Predictions of the response ###################################
  ################################   Predictions of the response ###################################
  ################################   Predictions of the response ###################################

  if(n.folds==1) if(all(length(which(is.na(raster::extract(mappa_ras[[1]],dataset[[1]][[1]][which(dataset[[1]][[1]]$value==1),]))))>(length(which(dataset[[1]][[1]]$value==1))/2))) cat("\n","Warning, more than 50% of target species occurrences falls out of the prediction ground domain! Try to fix target species geographic coordinates by using fix.coastal.points function")

  if(n.folds>1) {
    if(length(dataset[[1]][[1]])==1) {
      if(all(length(which(is.na(raster::extract(mappa_ras[[1]],dataset[[1]][[1]][which(dataset[[1]][[1]]$value==1),]))))>(length(which(dataset[[1]][[1]]$value==1))/2))) cat("\n","Warning, more than 50% of target species occurrences falls out of the prediction ground domain! Try to fix target species geographic coordinates by using fix.coastal.points function")
    }
    if(length(dataset[[1]][[1]])>1) {
      if(all(length(which(is.na(raster::extract(mappa_ras[[1]],dataset[[1]][[1]][[1]][which(dataset[[1]][[1]][[1]]$value==1),]))))>(length(which(dataset[[1]][[1]][[1]]$value==1))/2))) cat("\n","Warning, more than 50% of target species occurrences falls out of the prediction ground domain! Try to fix target species geographic coordinates by using fix.coastal.points function")
    }
  }


  if(n.sims==1) n.sims.clusters<-NULL
  cat("\n","Performing prediction of the response")

  if(is.null(n.sims.clusters)) {
    MINOSSE_list<-list()
    for(i in 1:n.sims)  {
      if(n.folds==1) {
        if(class(try(MINOSSE_list[[i]]<-minosse.core(dtset=dataset[[i]],
                                                     mappa_ras=mappa_ras,
                                                     projection.ground = projection.ground,
                                                     seed=seed,
                                                     predict_by=predict_by,configs=configs),silent=TRUE))=="try-error") MINOSSE_list[[i]]<-NA else MINOSSE_list[[i]]<-MINOSSE_list[[i]]
      }

      if(n.folds>1){
        MINOSSE_list[[i]]<-lapply(seq(1:length(dataset[[i]])),function(z){
          if(class(try(MINOSSE_list<-minosse.core(dtset=dataset[[i]][[z]],
                                                  mappa_ras=scale(mappa_ras),
                                                  projection.ground = projection.ground,
                                                  seed=seed,
                                                  predict_by=predict_by,configs=configs)))=="try-error") MINOSSE_list<-NA else MINOSSE_list<-MINOSSE_list
                                                  return(MINOSSE_list)})
      }
    }
  };gc()

  if(!is.null(n.sims.clusters)) {

  if(n.sims.clusters=="automatic") {
    parallel::detectCores()-1->n.sims.clusters
    if(n.sims.clusters>n.sims) {n.sims.clusters<-n.sims}}

  if(n.sims.clusters>1) {
    cl <- parallel::makeCluster(n.sims.clusters, type="SOCK")

    i<-NULL
    z<-NULL
    parallel::clusterExport(cl,c("i","z","lapply","minosse_fit.gstatModel","minosse_fit.vgmModel","minosse_fit.regModel","predict.gstatModel","minosse.core","n.sims","n.folds","dataset",
                                 "mappa_ras","projection.ground","seed","predict_by","configs"),envir = environment())
    doSNOW::registerDoSNOW (cl)


    if(n.folds==1) {

       MINOSSE_list<- foreach::foreach(i=1:n.sims,.export =c("i","lapply","minosse_fit.gstatModel","minosse_fit.vgmModel","minosse_fit.regModel","predict.gstatModel","minosse.core","n.sims","n.folds","dataset", "mappa_ras","projection.ground","seed","predict_by","configs"),
                                       .packages = c("plotKML","gstat","automap","PresenceAbsence","snow", "doSNOW","raster","dismo","sp","caret","xgboost"),.verbose=FALSE) %dopar% {

                                         

        if(class(try(MINOSSE_list<-minosse.core(dtset=dataset[[i]],
                                         mappa_ras=mappa_ras,
                                         projection.ground = projection.ground,
                                         seed=seed,
                                         predict_by=predict_by,configs=configs)))=="try-error") MINOSSE_list<-NA else MINOSSE_list<-MINOSSE_list


    return(MINOSSE_list)}
    if(!is.null(cl))parallel::stopCluster(cl)
    };gc()


    if(n.folds>1) {
       MINOSSE_list<- foreach::foreach(i=1:n.sims,.export =c("i","z","lapply","minosse_fit.gstatModel","minosse_fit.vgmModel","minosse_fit.regModel","predict.gstatModel","minosse.core","n.sims","n.folds","dataset", "mappa_ras","projection.ground","seed","predict_by","configs"),
                                       .packages = c("plotKML","gstat","automap","PresenceAbsence","snow", "doSNOW","raster","dismo","sp","caret","xgboost"),.verbose=FALSE) %dopar% {

                                         
                                      MINOSSE_list<-list()
                                      for(z in 1:length(dataset[[i]])) {
       if(class(try(MINOSSE_list[[z]]<-minosse.core(dtset=dataset[[i]][[z]],
                                                    mappa_ras=mappa_ras,
                                                    projection.ground = projection.ground,
                                                    seed=seed,
                                                    predict_by=predict_by,configs=configs)))=="try-error") MINOSSE_list[[z]]<-NA else MINOSSE_list[[z]]<-MINOSSE_list[[z]]
                                      }
    return(MINOSSE_list)}
    if(!is.null(cl))parallel::stopCluster(cl)};gc()
  };gc()

    }

#


  MINOSSE_list<-MINOSSE_list[which(!is.na(MINOSSE_list))]
  MINOSSE_list<-MINOSSE_list[sapply(MINOSSE_list,function(x)any(!is.na(x)))]
  if(n.folds==1) {
    MINOSSE_list<-MINOSSE_list[which(sapply(MINOSSE_list,function(x)class(x$kriging$kr_MINOSSE)=="RasterLayer"))]
  }
  if(n.folds>1) {
    lapply(MINOSSE_list,function(x)sapply(x,function(z)all(is.na(unlist(z)))))->full_model_check
    MINOSSE_list[sapply(full_model_check,function(x)any(which(x==FALSE)==length(x)))]->MINOSSE_list
    lapply(MINOSSE_list, function(x)x[which(!is.na(x))])->MINOSSE_list
    }


  ################################   Predictions of the response END ###################################
  ################################   Predictions of the response END ###################################
  ################################   Predictions of the response END ###################################
  ################################   Predictions of the response END ###################################

  if(n.sims==1){
    if(n.folds==1) {
      lapply(MINOSSE_list,function(x){
        if(class(try(MINOSSE_obj<-list(MINOSSE_stack=x[[1]][[1]],auc_values=x[[1]][[2]])))=="try-error") MINOSSE_obj<-NA else MINOSSE_obj<-MINOSSE_obj;
        return(MINOSSE_obj)
      })->MINOSSE_obj
    }

    if(n.folds>1) {
      lapply(MINOSSE_list,function(x){
        MINOSSE_obj<-lapply(seq(1:length(x)),function(z){
          if(class(try(MINOSSE_obj<-list(MINOSSE_stack=x[[z]][[1]][[1]],auc_values=x[[z]][[1]][[2]])))=="try-error") MINOSSE_obj<-NA else MINOSSE_obj<-MINOSSE_obj;
          return(MINOSSE_obj)})
      })->MINOSSE_obj
    }

  }
  if(n.sims>1) {
    if(n.folds==1){
      lapply(MINOSSE_list,function(x){
        if(class(try(MINOSSE_obj<-list(MINOSSE_stack=x[[1]][[1]],auc_values=x[[1]][[2]]),silent=TRUE))=="try-error") MINOSSE_obj<-NA else MINOSSE_obj<-MINOSSE_obj;
        return(MINOSSE_obj)
      })->MINOSSE_obj
  }
    if(n.folds>1) {
      lapply(MINOSSE_list,function(x){
        MINOSSE_obj<-lapply(seq(1:length(x)),function(z){
          if(class(try(MINOSSE_obj<-list(MINOSSE_stack=x[[z]][[1]][[1]],auc_values=x[[z]][[1]][[2]]),silent=TRUE))=="try-error") MINOSSE_obj<-NA else MINOSSE_obj<-MINOSSE_obj;
          return(MINOSSE_obj)})
      })->MINOSSE_obj
    }
  }


  if(length(MINOSSE_obj)==0) stop("MINOSSE stopped because of no valid model")

  if(n.folds==1) {
    lapply(MINOSSE_obj,function(x)x[[1]])->MINOSSE_stack;
    lapply(MINOSSE_obj,function(x)x[[2]])->auc_values;
    total_th<-lapply(seq(1:length(MINOSSE_obj)), function(z){data.frame(plotID=raster::extract(MINOSSE_stack[[z]],original_dataset[[z]],cellnumbers=TRUE)[,1],Observed=original_dataset[[z]]$value,Predicted1=raster::extract(MINOSSE_stack[[z]],original_dataset[[z]]))->total_data_for_th;
      suppressWarnings(PresenceAbsence::optimal.thresholds(total_data_for_th[complete.cases(total_data_for_th),])->total_th);
      rownames(total_th)<-NULL;
      return(total_th)})
    total_th<-do.call(rbind,total_th)
    unique(as.character(total_th$Method))->met_total_names
    total_th<-aggregate(Predicted1~Method,data=total_th,function(x) mean(x)  ,simplify = FALSE)
    total_th<-total_th[match(met_total_names,total_th$Method),]
    mean(stack(MINOSSE_stack))->MINOSSE
  }

  if(n.folds>1) {
    lapply(MINOSSE_obj,function(x) x[[length(x)]][[1]])->MINOSSE_stack;
    mean(unlist(lapply(MINOSSE_obj,function(x) {
      if(class(try(sapply(x[1:(length(x)-1)],function(z)z[[2]])->y,silent=TRUE))=="try-error") y<-NA else y<-y
      return(y)} )),na.rm=TRUE)->auc_values;
    total_th<-lapply(seq(1:length(MINOSSE_obj)), function(z){data.frame(plotID=raster::extract(MINOSSE_stack[[z]],original_dataset[[z]],cellnumbers=TRUE)[,1],Observed=original_dataset[[z]]$value,Predicted1=raster::extract(MINOSSE_stack[[z]],original_dataset[[z]]))->total_data_for_th;
      suppressWarnings(PresenceAbsence::optimal.thresholds(total_data_for_th[complete.cases(total_data_for_th),])->total_th);
      rownames(total_th)<-NULL;
      return(total_th)})
    total_th<-do.call(rbind,total_th)
    unique(as.character(total_th$Method))->met_total_names
    total_th<-aggregate(Predicted1~Method,data=total_th,function(x) mean(x)  ,simplify = FALSE)
    total_th<-total_th[match(met_total_names,total_th$Method),]
    mean(stack(MINOSSE_stack))->MINOSSE
  }

  if(n.sims==1) {
    if(n.folds==1){
      RK_RMSE<-MINOSSE_list[[1]]$kriging$RK_RMSE;
      RK_tvar<-MINOSSE_list[[1]]$kriging$RK_tvar;
      xgb_R2<-mean(MINOSSE_list[[1]]$kriging$xgb_R2);
      xgb_RMSE<-mean(MINOSSE_list[[1]]$kriging$xgb_RMSE);
      full_model_stats<-list(RK_AUC=auc_values,RK_RMSE=RK_RMSE,RK_tvar=RK_tvar,xgb_R2=xgb_R2,xgb_RMSE=xgb_RMSE)
    }
    if(n.folds>1){
      RK_RMSE<-MINOSSE_list[[1]][[n.folds+1]]$kriging$RK_RMSE;
      RK_tvar<-MINOSSE_list[[1]][[n.folds+1]]$kriging$RK_tvar;
      xgb_R2<-mean(MINOSSE_list[[1]][[n.folds+1]]$kriging$xgb_R2);
      xgb_RMSE<-mean(MINOSSE_list[[1]][[n.folds+1]]$kriging$xgb_RMSE);
      full_model_stats<-list(RK_AUC=auc_values,RK_RMSE=RK_RMSE,RK_tvar=RK_tvar,xgb_R2=xgb_R2,xgb_RMSE=xgb_RMSE)
    }
  }

  if(n.sims>1){
    if(n.folds==1){
      lapply(MINOSSE_list,function(x){
        RK_RMSE<-x$kriging$RK_RMSE;
        RK_tvar<-x$kriging$RK_tvar;
        xgb_R2<-mean(x$kriging$xgb_R2,na.rm=TRUE);
        xgb_RMSE<-mean(x$kriging$xgb_RMSE,na.rm=TRUE);
        return(full_model_stats<-data.frame(RK_AUC=mean(unlist(auc_values),na.rm=TRUE),RK_RMSE=RK_RMSE,RK_tvar=RK_tvar,xgb_R2=xgb_R2,xgb_RMSE=xgb_RMSE))
      })->full_model_stats
      do.call(rbind,full_model_stats)->full_model_stats
      full_model_stats<-as.list(apply(full_model_stats,2,function(x)mean(x,na.rm=TRUE)))
    }
    if(n.folds>1){
      lapply(MINOSSE_list,function(x){
        RK_RMSE<-x[[length(x)]]$kriging$RK_RMSE;
        RK_tvar<-x[[length(x)]]$kriging$RK_tvar;
        xgb_R2<-mean(x[[length(x)]]$kriging$xgb_R2,na.rm=TRUE);
        xgb_RMSE<-mean(x[[length(x)]]$kriging$xgb_RMSE,na.rm=TRUE);
        return(full_model_stats<-data.frame(RK_AUC=auc_values,RK_RMSE=RK_RMSE,RK_tvar=RK_tvar,xgb_R2=xgb_R2,xgb_RMSE=xgb_RMSE))
      })->full_model_stats
      do.call(rbind,full_model_stats)->full_model_stats
      full_model_stats<-as.list(apply(full_model_stats,2,function(x)mean(x,na.rm=TRUE)))
    }
  }
  if(n.folds==1) lapply(MINOSSE_list, function(x)x$kriging$MLA_mod)->MINOSSE_list
  if(n.folds>1) lapply(MINOSSE_list, function(x)x[[length(x)]]$kriging$MLA_mod)->MINOSSE_list # only full model is kept

  bin_maps<-lapply(total_th[,2],function(x)MINOSSE>=x)
  names(bin_maps)<-total_th$Method
  print("Prediction done")
  return(list(bin_maps=bin_maps,prob_map=MINOSSE,MINOSSE_STATS=MINOSSE_list,thresholds=total_th,predictors=mappa_ras,validation=list(RK_AUC=full_model_stats$RK_AUC,RK_RMSE=full_model_stats$RK_RMSE,RK_tvar=full_model_stats$RK_tvar,xgb_R2=full_model_stats$xgb_R2,xgb_RMSE=full_model_stats$xgb_RMSE),nam=unique(resp$spec),pts=resp))
}

