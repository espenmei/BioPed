rm(list = ls())

library(MASS)
library(tidyr)
library(devtools)

# Parameters
a = diag(c(0.8, 0.7, 0.6, 0.5, 0.4))
rA = matrix(c(1.0, 0.5, 0.3, 0.0, 0.0,
              0.5, 1.0, 0.5, 0.3, 0.0,
              0.3, 0.5, 1.0, 0.5, 0.3,
              0.0, 0.3, 0.5, 1.0, 0.5,
              0.0, 0.0, 0.3, 0.5, 1.0), 5, 5, byrow = T)

d = diag(c(0.8, 0.7, 0.6, 0.5, 0.4))
rD = matrix(c(1.0, 0.5, 0.3, 0.0, 0.0,
              0.5, 1.0, 0.5, 0.3, 0.0,
              0.3, 0.5, 1.0, 0.5, 0.3,
              0.0, 0.3, 0.5, 1.0, 0.5,
              0.0, 0.0, 0.3, 0.5, 1.0), 5, 5, byrow = T)

en = diag(c(0.0, 0.0, 0.0, 0.0, 0.0))
rEn = matrix(c(1.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 1.0), 5, 5, byrow = T)

ep = diag(c(0.0, 0.0, 0.0, 0.0, 0.0))
rEp = matrix(c(1.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 1.0), 5, 5, byrow = T)

ec = diag(c(0.3, 0.3, 0.3, 0.3, 0.3))
rEc = matrix(c(1.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 1.0), 5, 5, byrow = T)


nv = 5
m = c(2, 2, 2, 2, 2)
b_fem = c(1.5, 1.5, 1.5, 1.5, 1.5)
G = diag(0.5, nv, nv)
N_M = 10000
N_F = 10000
N_P = 2000
set.seed(1234)

# Mothers
Mpop = 1:N_M
A_Mpop = mvrnorm(N_M, rep(0, nv), rA)

# Fathers
Fpop = 1:N_F
A_Fpop = mvrnorm(N_F, rep(0, nv), rA)

# Pregnancies
P_M = sample(Mpop, N_P, replace = T)
P_F = sample(Fpop, N_P, replace = T)
P_PP = as.numeric(factor(paste0(P_M, "_", P_F)))

# Parent pairs
PP = unique(P_PP)
N_PP = length(PP)
D_PP = mvrnorm(N_PP, rep(0, nv), rD)
E_PP = mvrnorm(N_PP, rep(0, nv), rEn)

NO_P = sample(1:2, N_P, replace = T, prob = c(0.7, 0.3))
mz_P = ifelse(NO_P == 2, sample(0:1, sum(NO_P == 2), replace = T), 0)

# Children
Y_C = matrix(NA, sum(NO_P), nv)
colnames(Y_C) = paste0("Y", 1:nv)
female_C = vector("numeric", sum(NO_P))
ij = 1
for(i in 1:N_P) {
  # A
  rA_P = mvrnorm(1, rep(0, nv), G%*%rA)
  A_P = G%*%(A_Mpop[P_M[i], ] + A_Fpop[P_F[i], ]) + mz_P[i]*rA_P

  # D
  rD_P = mvrnorm(1, rep(0, nv), 0.75*rD)
  D_P = G%*%D_PP[P_PP[i], ] + mz_P[i]*rD_P

  # Ep
  E_P = mvrnorm(1, rep(0, nv), rEp)
  female_p = sample(0:1, 1)
  for(j in 1:NO_P[i]) {
    female_C[ij] = ifelse(mz_P[i] == 1, female_p, sample(0:1, 1))

    rA_C = mvrnorm(1, rep(0, nv), G%*%rA)
    A_C = A_P + (1 - mz_P[i])*rA_C

    rD_C = mvrnorm(1, rep(0, nv), 0.75*rD)
    D_C = D_P + (1 - mz_P[i])*rD_C

    # E_c
    E_C = mvrnorm(1, rep(0, nv), rEc)
    Y_C[ij, ] = m + b_fem*female_C[ij] + a%*%A_C + d%*%D_C +
      en%*%E_PP[P_PP[i], ] + ep%*%E_P + ec%*%E_C
    ij = ij + 1
  }
}

# Data
D_M = data.frame(Mid = factor(unique(P_M)))
D_F = data.frame(Fid = factor(unique(P_F)))
D_PP = data.frame(PPid = factor(PP))
D_P = data.frame(Pid = factor(1:N_P),
                 PPid = factor(P_PP),
                 Mid = factor(P_M),
                 Fid = factor(P_F),
                 mz = mz_P)
D_C = data.frame(Y_C,
                 Pid = factor(rep(1:N_P, times = NO_P)),
                 female = female_C)

D_5var = left_join(D_C, D_P, "Pid")
devtools::use_data(D_5var, pkg = "C:/programmer/r/packages/BioPed", overwrite = T)
