#' @importFrom sp over
#' @importFrom doParallel registerDoParallel
#' @importFrom GSIF fit.vgmModel
#' @importFrom methods new
#' @importFrom stats aggregate as.formula gaussian glm na.omit na.pass predict resid shapiro.test step weights


minosse_fit.gstatModel<-function(configs=configs,observations, formulaString, covariates, method = list("GLM", "rpart", "randomForest", "quantregForest", "xgboost", "ranger"), dimensions = list("2D", "3D", "2D+T", "3D+T"), fit.family = gaussian(), stepwise = TRUE, vgmFun = "Exp", subsample = 5000, subsample.reg = 10000, ...){

  method <- method[[1]]
  dimensions <- dimensions[[1]]
  ## TH: the function only works with 2D maps at the moment:
  if(length(attr(sp::coordinates(observations), "dimnames")[[2]])>2){
    warning("This method uses only 2D coordinates of the points. For 3D data consider using the 'geosamples-class'.")
  }
  ## TH: the model has to include at least one covariate:
  if(length(all.vars(formulaString)[-1])==0){
    stop("No covariates have been specified in the 'formulaString'")
  }
  ## check argument 'fit.family'
  if(!missing(fit.family) & !method=="GLM"){ warning("'fit.family' argument will be ignored. Using 'method=\"GLM\"' instead.")  }

  ## overlay:
  ov <- over(observations, covariates)
  ## all variables of interest:
  tv <- all.vars(formulaString)[1]
  seln <- names(covariates) %in% all.vars(formulaString)[-1]
  ## check if all covariates are available:
  if(sum(!is.na(seln))==0){
    stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
  }
  ov <- cbind(data.frame(observations[,tv]), ov)

  ## copy coordinate column names for consistency:
  xyn <- which(names(ov) %in% attr(observations@bbox, "dimnames")[[1]])
  if(is.null(attr(covariates@bbox, "dimnames")[[1]])) {
    dim.names = attr(observations@bbox, "dimnames")[[1]]
  } else {
    dim.names = attr(covariates@bbox, "dimnames")[[1]]
  }
  if(length(xyn)==2) {
    names(ov)[xyn] <- dim.names[1:2]
  } else {
    names(ov)[xyn] <- dim.names
  }

  ## check the size of output:
  if(nrow(ov)==0|is.null(ov[,tv])) {
    stop("The 'over' operations resulted in an empty set.")
  }

  ## fit/filter the regression model:
  m <- minosse_fit.regModel(configs=configs,formulaString = formulaString, rmatrix = ov, predictionDomain = covariates[seln], method = method, dimensions = dimensions, fit.family = fit.family, stepwise = stepwise, vgmFun = vgmFun, subsample = subsample, subsample.reg = subsample.reg, ...)
  return(m)

}
