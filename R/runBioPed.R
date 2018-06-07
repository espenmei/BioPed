#' runBioPed
#'
#' Fits a biometric model based on parent identifiers.
#'
#' @param vrs string vector with names of all \eqn{p} outcome variables.
#' @param AModel list with OpenMx components for additive genetic model.
#' Must contain mxMatrix/mxAlgebra with dim \eqn{pxp} named VA.
#' @param DModel list with OpenMx components for dominance genetic model.
#' Must contain mxMatrix/mxAlgebra with dim \eqn{pxp} named VD.
#' @param EmModel list of OpenMx components for mother environmental model.
#' Must contain mxMatrix/mxAlgebra with dim \eqn{pxp} named VEm.
#' @param EfModel list of OpenMx components for father environmental model.
#' Must contain mxMatrix/mxAlgebra with dim \eqn{pxp} named VEf.
#' @param EnModel list of OpenMx components for nuclear family environmental model.
#' Must contain mxMatrix/mxAlgebra with dim \eqn{pxp} named VEn.
#' @param EpModel list of OpenMx components for twin/pregnancy environmental model.
#' Must contain mxMatrix/mxAlgebra with dim \eqn{pxp} named VEp.
#' @param EcModel list of OpenMx components for unique environmental model.
#' Must contain mxMatrix/mxAlgebra with dim \eqn{pxp} named VEc.
#' @param MModel list of OpenMx components for mean model.
#' Must contain mxMatrix/mxAlgebra with dim \eqn{px1} named My.
#' @param data_M dataframe with mother identifiers.
#' @param data_F dataframe with father identifiers.
#' @param data_PP dataframe with parent-pair identifiers.
#' @param data_P dataframe with pregnancy identifiers and mz dummy.
#' @param data_C dataframe with pregnancy identifiers, mz dummy, covariates and variables to be modelled (vrs).
#' @param ... Arguments to be passed to \code{MxTryhard()}.
#' @export
#' @examples
#'library(OpenMx)
#'library(dplyr)
#'library(labOpenMx)
#'
#'data(D_5var) # Load example dataframe
#'
#'nv = 5
#'vrs = paste0("Y", 1:nv)
#'
#'D_M = data.frame(Mid = unique(D_5var$Mid))
#'D_F = data.frame(Fid = unique(D_5var$Fid))
#'D_PP = data.frame(PPid = unique(D_5var$PPid))
#'D_P = D_5var %>%
#'  group_by(Pid) %>%
#'  summarise(PPid = unique(PPid),
#'            Mid = unique(Mid),
#'            Fid = unique(Fid),
#'            mz = unique(mz)) %>%
#'  ungroup() %>%
#'  data.frame()
#'
#'D_C = D_5var %>%
#'  dplyr::select(Y1:Y5, female, mz, Pid) %>%
#'  data.frame()
#'
#'# Fit an "independent pathway" model
#'# Additive genetic model
#'Amod = list(mxMatrix("Full", nv, 1, T, 0.5, labF("lCA", nv, 1), name = "lCA"),
#'            mxMatrix("Diag", nv, nv, T, 0.5, labD("lSA", nv), name = "lSA"),
#'            mxAlgebra(lCA%*%t(lCA) + lSA%*%t(lSA), name = "VA"))
#'
#'# pregnancy/twin specific environmental model
#'Tmod = list(mxMatrix("Full", nv, 1, T, 0.5, labF("lCT", nv, 1), name = "lCT"),
#'            mxMatrix("Diag", nv, nv, T, 0.5, labD("lST", nv), name = "lST"),
#'            mxAlgebra(lCT%*%t(lCT) + lST%*%t(lST), name = "VEp"))
#'
#'# Unique environmental model
#'Emod = list(mxMatrix("Full", nv, 1, T, 0.5, labF("lCE", nv, 1), name = "lCE"),
#'            mxMatrix("Diag", nv, nv, T, 0.5, labD("lSE", nv), name = "lSE"),
#'            mxAlgebra(lCE%*%t(lCE) + lSE%*%t(lSE), name = "VEc"))
#'
#'# Mean model
#'Mmod = list(mxMatrix("Full", nv, 2, T, 0, labF("b", nv, 2), name = "B"),
#'            mxMatrix("Full", 2, 1, F, 1, c("x1", "data.female"), name = "X"),
#'            mxAlgebra(B%*%X, name = "My"))
#'
#'fit = runBioPed(vrs, AModel = Amod, EpModel = Tmod, EcModel = Emod, MModel = Mmod,
#'                data_M = D_M, data_F = D_F, data_PP = D_PP, data_P = D_P, data_C = D_C)
#'summary(fit)
#'
#' @details
#' The arguments \code{vrs, data_M, data_F, data_PP, data_P} and \code{data_C} must all be supplied.
#' In addition, at least one of the arguments \code{AModel, DModel, EpModel, EnModel, EcModel} or \code{MModel} must be supplied.
runBioPed = function(vrs,
                     AModel, DModel, EmModel, EfModel, EnModel, EpModel, EcModel, MModel,
                     data_M, data_F, data_PP, data_P, data_C, ...) {

  # Check arguments
  if(missing(vrs)) {
    stop("Variable names must be supplied")
  }

  if(any(c(missing(data_M), missing(data_F), missing(data_PP), missing(data_P), missing(data_C)))) {
    stop("A dataframe must be supplied for each level")
  }

  if(all(c(missing(AModel), missing(DModel), missing(EnModel), missing(EpModel), missing(EcModel)))) {
    stop("Some kind of covariance-model must be specified")
  }

  # Define variables
  nv = length(vrs)

  if(missing(AModel)) {
    AModel = mxMatrix("Zero", nv, nv, name = "VA")
  }

  if(missing(DModel)) {
    DModel = mxMatrix("Zero", nv, nv, name = "VD")
  }

  if(missing(EmModel)) {
    EmModel = mxMatrix("Zero", nv, nv, name = "VEm")
  }

  if(missing(EfModel)) {
    EfModel = mxMatrix("Zero", nv, nv, name = "VEf")
  }

  if(missing(EnModel)) {
    EnModel = mxMatrix("Zero", nv, nv, name = "VEn")
  }

  if(missing(EpModel)) {
    EpModel = mxMatrix("Zero", nv, nv, name = "VEp")
  }

  if(missing(EcModel)) {
    EcModel = mxMatrix("Zero", nv, nv, name = "VEc")
  }

  if(missing(MModel)) {
    MModel = mxMatrix("Zero", nv, 1, name = "My")
  }

  # Create model
  mod = modBioPed(nv, vrs,
                  AModel, DModel, EmModel, EfModel, EnModel, EpModel, EcModel, MModel,
                  data_M, data_F, data_PP, data_P, data_C)
  # Fit the model
  fit = mxTryHard(mod, ...)
  # Return fitted model
  return(fit)

}
