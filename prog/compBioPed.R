rm(list = ls())

library(dplyr)
library(tidyr)
library(labOpenMx)

data(D_5var) # Load example dataframe

nv = 5
vrs = paste0("Y", 1:nv)

# Data
D_M = data.frame(Mid = unique(D_5var$Mid))
D_F = data.frame(Fid = unique(D_5var$Fid))
D_PP = data.frame(PPid = unique(D_5var$PPid))
D_P = D_5var %>%
  group_by(Pid) %>%
  summarise(PPid = unique(PPid),
            Mid = unique(Mid),
            Fid = unique(Fid),
            mz = unique(mz)) %>%
  ungroup() %>%
  data.frame()
D_C = D_5var %>%
  dplyr::select(Y1:Y5, female, mz, Pid) %>%
  data.frame()

D_C = D_C %>%
  group_by(Pid) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  data.frame()
D_C$Pid = droplevels(D_C$Pid)


D_P = data.frame(D_P[D_P$Pid %in% D_C$Pid, ])
D_P$Pid = droplevels(D_P$Pid)
D_P$PPid = droplevels(D_P$PPid)
D_P$Mid = droplevels(D_P$Mid)
D_P$Fid = droplevels(D_P$Fid)

D_PP = data.frame(D_PP[D_PP$PPid %in% D_P$PPid, ])
colnames(D_PP) = "PPid"
D_PP$PPid = droplevels(D_PP$PPid)

D_F = data.frame(D_F[D_F$Fid %in% D_P$Fid, ])
colnames(D_F) = "Fid"
D_F$Fid = droplevels(D_F$Fid)

D_M = data.frame(D_M[D_M$Mid %in% D_P$Mid, ])
colnames(D_M) = "Mid"
D_M$Mid = droplevels(D_M$Mid)

# Additive genetic model
lA = mxMatrix("Lower", nv, nv, T, 0.5, labL("la", nv), name = "lA")
VA = mxAlgebra(lA%*%t(lA), name = "VA")

# Dominance genetic model
lD = mxMatrix("Lower", nv, nv, F, 0.0, labL("ld", nv), name = "lD")
VD = mxAlgebra(lD%*%t(lD), name = "VD")

# Mother environmental model
lEm = mxMatrix("Lower", nv, nv, T, 0.7, labL("lem", nv), name = "lEm")
VEm = mxAlgebra(lEm%*%t(lEm), name = "VEm")

# Nuclear environmental model
lEn = mxMatrix("Lower", nv, nv, F, 0.0, labL("len", nv), name = "lEn")
VEn = mxAlgebra(lEn%*%t(lEn), name = "VEn")

# Twin specific environmental model
lEp = mxMatrix("Lower", nv, nv, F, 0.0, labL("lep", nv), name = "lEp")
VEp = mxAlgebra(lEp%*%t(lEp), name = "VEp")

# Unique environmental model
lEc = mxMatrix("Lower", nv, nv, T, 0.5, labL("lec", nv), name = "lEc")
VEc = mxAlgebra(lEc%*%t(lEc), name = "VEc")

# Mean model
B = mxMatrix("Full", nv, 1, T, 0, labF("b", nv, 2), name = "B")
X = mxMatrix("Full", 2, 1, F, 1, c("x1", "data.female"), name = "X")
My = mxAlgebra(B, name = "My")

AMod = list(lA, VA)
DMod = list(lD, VD)
EmMod = list(lEm, VEm)
EnMod = list(lEn, VEn)
EpMod = list(lEp, VEp)
EcMod = list(lEc, VEc)
MMod = list(B, X, My)
fit = runBioPed(vrs,
                AModel = AMod,
                DModel = DMod,
                EmModel = EmMod,
                EnModel = EnMod,
                EpModel = EpMod,
                EcModel = EcMod,
                MModel = MMod,
                data_M = D_M,
                data_F = D_F,
                data_PP = D_PP,
                data_P = D_P,
                data_C = D_C,
                extraTries = 2)
summary(fit)
mxEval(VEm, fit)
# Compare
# ---------------------------------------------------------------------
Dtw = D_C %>%
  group_by(Pid) %>%
#  filter(n() == 2) %>%
  mutate(tw = 1:n()) %>%
  ungroup() %>%
  gather(var, val, -mz, -tw, -Pid) %>%
  unite(twvar, tw, var, sep = "_") %>%
  spread(twvar, val, sep = "_") %>%
  mutate(rA = ifelse(mz == 1, 1, 0.50),
         rD = ifelse(mz == 1, 1, 0.25)) %>%
  data.frame()

mT = mxModel("Trad",
             AMod, DMod, EnMod, EpMod, EcMod,
             mxAlgebra(rbind(cbind(VA + VD + VEn + VEp + VEc, data.rA%x%VA + data.rD%x%VD + VEn + VEp),
                             cbind(data.rA%x%VA + data.rD%x%VD + VEn + VEp, VA + VD + VEn + VEp + VEc)), name = "V"),
             mxMatrix("Full", 1, nv*2, T, 0, c(paste0("m", 1:nv), paste0("m", 1:nv)), name = "M"),
             mxExpectationNormal("V", "M", colnames(Dtw)[grep("Y", colnames(Dtw))]),
             mxFitFunctionML(),
             mxData(Dtw, "raw"))
fmT = mxRun(mT)
summary(fmT)
summary(fit)
logLik(fit)
logLik(fmT)

mxEval(B, fit)
mxEval(M, fmT)
mxEval(lEc, fit)
mxEval(lEc, fmT)
mxEval(lA, fit)
mxEval(lA, fmT)
