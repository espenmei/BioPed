#' The actual OpenMx model
#' @param nv number of variables to be analysed
#' @param vrs string vector with names of all \eqn{p} outcome variables.
#' @param AModel list with OpenMx components for additive genetic model.
#' Must contain mxMatrix/mxAlgebra with dim \eqn{pxp} named VAf.
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
#' @return fitted OpenMx model
modBioPed = function(nv, vrs,
                     AModel, DModel, EmModel, EfModel, EnModel, EpModel, EcModel, MModel,
                     data_M, data_F, data_PP, data_P, data_C) {

  # Variable labels
  #lat_PA = paste0("A", 1:nv, "_PA")
  lat_M = c(paste0("A", 1:nv, "_PA"),
            paste0("E", 1:nv, "_M"))
  lat_F = c(paste0("A", 1:nv, "_PA"),
            paste0("E", 1:nv, "_F"))
  lat_PP = paste0("D", 1:nv, "_PP")
  lat_P = c(paste0("A", 1:nv, "_P"),
            paste0("D", 1:nv, "_P"),
            paste0("E", 1:nv, "_P"),
            paste0("Edummy", 1:nv, "_P"))
  man_C = vrs
  lat_C = c(paste0("A", 1:nv, "_C"),
            paste0("D", 1:nv, "_C"),
            paste0("E", 1:nv, "_C"))
  all_C = c(man_C, lat_C)

  # Dimensions
  #nlat_PA = length(lat_PA)
  nlat_M = length(lat_M)
  nlat_F = length(lat_F)
  nlat_PP = length(lat_PP)
  nlat_P = length(lat_P)
  nlat_C = length(lat_C)
  nall_C = length(all_C)

  # Constants
  nvS00 = mxMatrix("Zero", nv, nv, name = "nvS00")
  nvD10 = mxMatrix("Iden", nv, nv, name = "nvD10")
  nvD05 = mxMatrix("Diag", nv, nv, F, 0.5, name = "nvD05")
  const = list(nvS00, nvD10, nvD05)

  # Parent models
  # A
  #mod_PA = mxModel("mod_PA", type = "RAM", latentVars = lat_PA,
  #                 AModel,
  #                 mxMatrix("Zero", nlat_PA, nlat_PA, dimnames = list(lat_PA, lat_PA), name = "A"),
  #                 mxMatrix("Zero", 0, nlat_PA, dimnames = list(NULL, lat_PA), name = "F"),
  #                 mxMatrix("Zero", 1, nlat_PA, dimnames = list(NA, lat_PA), name = "M"),
  #                 mxExpectationRAM("A", "VA", "F", "M"))
  #mod_M = mxModel(mod_PA, mxData(D_M, "raw", primaryKey = "Mid"), name = "mod_M")
  #mod_F = mxModel(mod_PA, mxData(D_F, "raw", primaryKey = "Fid"), name = "mod_F")

  mod_M = mxModel("mod_M", type = "RAM", latentVars = lat_M,
                  mxData(D_M, "raw", primaryKey = "Mid"),
                  AModel, EmModel, const,
                  mxAlgebra(rbind(cbind(VA, nvS00),
                                  cbind(nvS00, VEm)), dimnames = list(lat_M, lat_M), name = "S_M"),
                  mxMatrix("Zero", nlat_M, nlat_M, dimnames = list(lat_M, lat_M), name = "A"),
                  mxMatrix("Zero", 0, nlat_M, dimnames = list(NULL, lat_M), name = "F"),
                  mxMatrix("Zero", 1, nlat_M, dimnames = list(NA, lat_M), name = "M"),
                  mxExpectationRAM("A", "S_M", "F", "M"))

  mod_F = mxModel("mod_F", type = "RAM", latentVars = lat_F,
                  mxData(D_F, "raw", primaryKey = "Fid"),
                  AModel, EfModel, const,
                  mxAlgebra(rbind(cbind(VA, nvS00),
                                  cbind(nvS00, VEf)), dimnames = list(lat_F, lat_F), name = "S_F"),
                  mxMatrix("Zero", nlat_F, nlat_F, dimnames = list(lat_F, lat_F), name = "A"),
                  mxMatrix("Zero", 0, nlat_F, dimnames = list(NULL, lat_F), name = "F"),
                  mxMatrix("Zero", 1, nlat_F, dimnames = list(NA, lat_F), name = "M"),
                  mxExpectationRAM("A", "S_F", "F", "M"))

  # Parent-pair/nuclear model
  # D
  mod_PP = mxModel("mod_PP", type = "RAM", latentVars = lat_PP,
                   mxData(D_PP, "raw", primaryKey = "PPid", sort = F),
                   DModel,
                   mxAlgebra(VD, dimnames = list(lat_PP, lat_PP), name = "S_PP"),
                   mxMatrix("Zero", nlat_PP, nlat_PP, dimnames = list(lat_PP, lat_PP), name = "A"),
                   mxMatrix("Zero", 0, nlat_PP, dimnames = list(NULL, lat_PP), name = "F"),
                   mxMatrix("Zero", 1, nlat_PP, dimnames = list(NA, lat_PP), name = "M"),
                   mxExpectationRAM("A", "S_PP", "F", "M"))

  # Pregnancy model
  # A, D, Ep, Edummy
  mod_P = mxModel("mod_P", type = "RAM", latentVars = lat_P, mod_PP, mod_M, mod_F,
                  mxData(D_P, "raw", primaryKey = "Pid", sort = F), const,
                  AModel, DModel, EpModel,
                  mxAlgebra(rbind(cbind(data.mz*0.50*VA, nvS00, nvS00, nvS00),
                                  cbind(nvS00, data.mz*0.75*VD, nvS00, nvS00),
                                  cbind(nvS00, nvS00, VEp, nvS00),
                                  cbind(nvS00, nvS00, nvS00, nvS00)),
                            dimnames = list(lat_P, lat_P),
                            name = "S_P"),
                  mxMatrix("Zero", nlat_P, nlat_P, dimnames = list(lat_P, lat_P), name = "A"),
                  mxMatrix("Zero", 0, nlat_P, dimnames = list(NULL, lat_P), name = "F"),
                  mxMatrix("Zero", 1, nlat_P, dimnames = list(NA, lat_P), name = "M"),
                  mxAlgebra(rbind(nvS00, nvD05, nvS00, nvS00),
                            dimnames = list(lat_P, lat_PP),
                            name = "T_PonPP",
                            joinKey = "PPid",
                            joinModel = "mod_PP"),
                  mxAlgebra(rbind(cbind(nvD05, nvS00),
                                  cbind(nvS00, nvS00),
                                  cbind(nvS00, nvS00),
                                  cbind(nvS00, nvD10)),
                            dimnames = list(lat_P, lat_M),
                            name = "T_PonM",
                            joinKey = "Mid",
                            joinModel = "mod_M"),
                  mxAlgebra(rbind(cbind(nvD05, nvS00),
                                  cbind(nvS00, nvS00),
                                  cbind(nvS00, nvS00),
                                  cbind(nvS00, nvD10)),
                            dimnames = list(lat_P, lat_F),
                            name = "T_PonF",
                            joinKey = "Fid",
                            joinModel = "mod_F"),
                  mxExpectationRAM("A", "S_P", "F", "M", between = c("T_PonPP", "T_PonM", "T_PonF")))

  # Child model
  # Y, A, D, Ec
  mod_C = mxModel("mod_C", type = "RAM", latentVars = lat_C, manifestVars = man_C, mod_P,
                  mxData(D_C, "raw",  sort = F), const, MModel,
                  AModel, DModel, EcModel,
                  mxAlgebra(rbind(cbind(nvS00, nvS00, nvS00, nvS00),
                                  cbind(nvS00, (1 - data.mz)*0.50*VA, nvS00, nvS00),
                                  cbind(nvS00, nvS00, (1 - data.mz)*0.75*VD, nvS00),
                                  cbind(nvS00, nvS00, nvS00, VEc)),
                            dimnames = list(all_C, all_C),
                            name = "S_C"),
                  mxAlgebra(rbind(cbind(nvS00, nvD10, nvD10, nvD10),
                                  cbind(nvS00, nvS00, nvS00, nvS00),
                                  cbind(nvS00, nvS00, nvS00, nvS00),
                                  cbind(nvS00, nvS00, nvS00, nvS00)),
                            dimnames = list(all_C, all_C),
                            name = "A_C"),
                  mxMatrix("Full", nv, nall_C, F,
                           cbind(diag(1, nv, nv), matrix(0, nv, nlat_C)),
                           dimnames = list(man_C, all_C),
                           name = "F"),
                  mxAlgebra(cbind(t(My), 0*t(My), 0*t(My), 0*t(My)), dimnames = list(NULL, all_C), name = "M"),

                  mxAlgebra(rbind(cbind(nvS00, nvS00, nvD10, nvD10),
                                  cbind(nvD10, nvS00, nvS00, nvS00),
                                  cbind(nvS00, nvD10, nvS00, nvS00),
                                  cbind(nvS00, nvS00, nvS00, nvS00)),
                            dimnames = list(all_C, lat_P),
                            name = "T_ConP",
                            joinKey = "Pid",
                            joinModel = "mod_P"),

                  mxExpectationRAM("A_C", "S_C", "F", "M", between = "T_ConP"))

  return(mod_C)

}
