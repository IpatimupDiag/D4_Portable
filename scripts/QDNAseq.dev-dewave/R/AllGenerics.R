setGeneric("createNormalCalibration", function(object, force=FALSE)
    standardGeneric("createNormalCalibration"))
setGeneric("dewaveBins", function(object, user.calibration=FALSE, force=FALSE)
    standardGeneric("dewaveBins"))

#NoWaves
setGeneric("CorrectTumors", function(object, user.calibration, ...)
  standardGeneric("CorrectTumors"))
setGeneric("CorTumNorm", function(object, user.calibration)
  standardGeneric("CorTumNorm"))
setGeneric("SmoothNormals", function(object, ...)
  standardGeneric("SmoothNormals"))
