#' @importFrom intamap createIntamapObject getIntamapParams checkSetup preProcess postProcess estimateParameters methodParameters spatialPredict
#' @importFrom methods is
#' @importFrom utils capture.output modifyList
#' @author Francesco Carotenuto, francesco.carotenuto@unina.it

minosse_interpolate<-function (observations, predictionLocations, outputWhat = list(mean = TRUE,
                                                                                    variance = TRUE), obsChar = NA, methodName = "automatic",
                               optList = list(), cv = FALSE)
{
  
  
  
  obsChar = NA
  startTime = Sys.time()
  npred = ifelse(is.numeric(predictionLocations), predictionLocations,
                 nrow(coordinates(predictionLocations)))
  msg = paste("R", startTime, "interpolating", nrow(coordinates(observations)),
              "observations,", npred, "prediction locations\n")
  cat(msg)
  if ("formulaString" %in% names(optList)) {
    formulaString = as.formula(optList$formulaString)
  }  else {
    formulaString = as.formula("value~1")
  }
  localParams = list(confProj = TRUE, debug.level = 0)
  localParams = modifyList(localParams, optList)
  debug.level = localParams$debug.level
  if (debug.level >= 1) {
    print("sending debug info to processDescription")
    debugOutput = 1
    rm(debugOutput)
    debugOut = textConnection("debugOutput", "w")
    sink(debugOut, split = TRUE)
  }
  # if (methodName == "automatic" & !(sum(names(optList) %in%
  #                                       c("methodParameters", "variogramModel", "copulaParams"))) >=
  #     1) {
  #   methodName = chooseMethod(observations, predictionLocations,
  #                             formulaString, obsChar, maximumTime, outputWhat = outputWhat)
  # }
  krigingObject = createIntamapObject(observations = observations,
                                      predictionLocations = predictionLocations, formulaString = formulaString,
                                      outputWhat = outputWhat, obsChar = obsChar, params = getIntamapParams(localParams),
                                      class = methodName)
  debug.level = krigingObject$params$debug.level
  checkSetup(krigingObject)
  krigingObject = intamap::preProcess(krigingObject)
  if (is.null(krigingObject$variogramModel) && is.null(krigingObject$copulaParams) &&
      is.null(krigingObject$inverseDistancePower)) {
    krigingObject = estimateParameters(krigingObject)
  }
  krigingObjectMp = try(methodParameters(krigingObject))
  if (!is(krigingObjectMp, "try-error"))
    krigingObject = krigingObjectMp
  nsim = ifelse("nsim" %in% names(outputWhat), outputWhat$nsim,
                0)
  if (cv) {
    kObj = krigingObject
    predictions = krigingObject$observations
    depVar = as.character(krigingObject$formulaString[[2]])
    predictions@data = data.frame(var1.pred = NA, var1.var = NA,
                                  observed = observations@data[, depVar], residual = NA,
                                  zscore = NA)
    for (i in 1:dim(krigingObject$observations)[1]) {
      kObj$predictionLocations = krigingObject$observations[i,
                                                            ]
      kObj$observations = krigingObject$observations[-i,
                                                     ]
      if (debug.level == 0) {
        tmp = capture.output(kObj <- spatialPredict(kObj))
      }      else kObj <- spatialPredict(kObj)
      if ("var1.pred" %in% names(kObj$predictions) & "var1.var" %in%
          names(kObj$predictions)) {
        predictions@data[i, 1:2] = kObj$predictions@data[,
                                                         c("var1.pred", "var1.var")]
      }      else predictions@data[i, 1:2] = kObj$predictions@data[,
                                                                   c("mean", "variance")]
    }
    predictions$residual = predictions$observed - predictions$var1.pred
    predictions$zscore = predictions$residual/sqrt(predictions$var1.var)
    krigingObject$predictions = predictions
  }  else krigingObject = spatialPredict(krigingObject, nsim = nsim)
  krigingObject = postProcess(krigingObject)
  if (!is.null(krigingObject$returnPlot) && krigingObject$returnPlot)
    krigingObject$processPlot = ""
  else krigingObject$processPlot = ""
  krigingObject$processDescription = paste("Spatial prediction using the method ",
                                           class(krigingObject))
  if (debug.level >= 1) {
    sink()
    close(debugOut)
    krigingObject$processDescription = c(krigingObject$processDescription,
                                         debugOutput)
  }
  return(krigingObject)
}
