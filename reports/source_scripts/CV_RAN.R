# Additional information for the following article:

# Multi-environment genomic selection in rice elite breeding lines
# 
# Nguyen V. H.1,2,3,4, Morantte R. I. Z.3, Lopena V.3, Verdeprado H.3, Murori R.3, Ndayiragije A.3, Katiyar S. 3, Islam M. R. 3, Juma R.3, Galvez H.4, Glaszmann J-C.1,2, Cobb J. N. 3,5, Bartholomé J.2,6,7*
#   
# 1 CIRAD, UMR AGAP Institut, F-34398, Montpellier, France
# 2 UMR AGAP Institut, Univ Montpellier, CIRAD, INRAE, Institut Agro, Montpellier, France
# 3 Rice Breeding Platform, International Rice Research Institute, DAPO Box7777, Metro
# Manila, Philippines 
# 4 University of the Philippines Los Baños, College, Los Baños, Laguna 4031, Philippines
# 5 RiceTec. Inc, PO Box 1305, Alvin, TX 77512
# 6 CIRAD, UMR AGAP Institut, Cali, Colombia
# 7 Alliance Bioversity-CIAT, Cali, Colombia
# *Corresponding author: jerome.bartholome@cirad.fr


# Script for one repetition

rm(list = ls())
# load packages ----------------------------------------------------------------
library(magrittr)
library(BGLR)
library(dplyr)

print(Sys.time())

# load data ---------------------------------------------------------------
load("data_CV_RAN.RData")

# function for covar matrices ---------------------------------------------

covarMat <- function(pheno, geno, EC) {
  
  # Co-var structure of E 
  Z <- as.matrix(model.matrix(~ tibble::deframe(pheno[, "ENV"]) - 1))
  G.E <- tcrossprod(Z)
  
  # Co-var structure of W 
  covarW <- function(W) {
    W <- as.matrix(W)
    G.W <- tcrossprod(scale(W)) / ncol(W)
    return(G.W)
  }
  G.W <- covarW(EC[, -1])
  
  # Co-var structure G
  IDs <- tibble::deframe(pheno[, "GID"])
  X <- as.matrix(geno)
  
  colNames <- colnames(X)[-1]
  IDx <- X[, 1]
  X <- X[, -1]
  
  rownames(X) <- IDx
  colnames(X) <- colNames
 
  G.G <- tcrossprod(scale(X)) / ncol(X)
  
  IDs <- factor(IDs, levels = rownames(G.G))
  Z <- as.matrix(model.matrix( ~ IDs - 1))
  G.G <- tcrossprod(tcrossprod(Z, G.G), Z)  #This is ZGZ’
  
  # Co-var structure G*E 
  G.GxE <- G.G * G.E
  
  # Co-var structure G*W 
  G.GxW <- G.G * G.W
  
  return(list(G.G, G.E, G.W, G.GxE, G.GxW))
}

# set parameters ----------------------------------------------------------
# param <- commandArgs(trailingOnly = TRUE)
param <- c(3, 4)

rep <- as.numeric(param[1])
nEnv <- as.numeric(param[2])

# CV - random -------------------------------------------------------------

# number of iterations
nIter <- 25000
# burnin period
burnIn <- 5000
# thin
thin <- 10

traits <- c("DTF", "HT", "YLD")

model <-
  c("G", "GE", "GW", "GEW", "GE.GxE", "GW.GxW", "GEW.GxE.GxW")

res <- data.frame(
  rep = NA,
  nEnv = NA,
  trait = NA,
  env = NA,
  model = NA,
  PA = NA,
  PA_r = NA
)
a = 0

for (trait in traits) {
  cat("\n")
  print("==========================================")
  print(Sys.time())
  print(trait)
  for (Env in unique(pheno$ENV)) {
    cat("\n")
    print(Sys.time())
    print(Env)
    
    # randomly select env for training set
    set.seed(rep)
    vect_env <- unique(pheno$ENV)[-match(Env, unique(pheno$ENV))]
    rEnv <- sample(vect_env, nEnv)
    pheno_int <- pheno[pheno$ENV %in% c(rEnv, Env),]
    
    # subset phenotypic data
    y <- tibble::deframe(pheno_int[, trait])
    yNA <- y
    valSet <- pheno_int$ENV == Env
    yNA[valSet] <- NA
    
    # create relationship matrices 
    geno_int <- geno[geno$GID %in% unique(pheno_int$GID), ]
    
    EC_int <- pheno_int %>% select(ENV) %>%
      left_join(EC_24)

    G <- covarMat(pheno = pheno_int,
             geno = geno_int,
             EC = EC_int[,-1])
    
    G.G <- G[[1]]
    G.E <- G[[2]]
    G.W <- G[[3]]
    G.GxE <- G[[4]]
    G.GxW <- G[[5]]
    
    # model choice
    for (m in model) {
      print(m)
      ETA <- list()
      if (m == "G") {
        ETA[[1]] <- list(K = G.G, model = 'RKHS')
      } else if (m == "GE") {
        ETA[[1]] <- list(K = G.G, model = 'RKHS')
        ETA[[2]] <- list(K = G.E, model = 'RKHS')
      } else if (m == "GW") {
        ETA[[1]] <- list(K = G.G, model = 'RKHS')
        ETA[[2]] <- list(K = G.W, model = 'RKHS')
      } else if (m == "GEW") {
        ETA[[1]] <- list(K = G.G, model = 'RKHS')
        ETA[[2]] <- list(K = G.E, model = 'RKHS')
        ETA[[3]] <- list(K = G.W, model = 'RKHS')
      } else if (m == "GE.GxE") {
        ETA[[1]] <- list(K = G.G, model = 'RKHS')
        ETA[[2]] <- list(K = G.E, model = 'RKHS')
        ETA[[3]] <- list(K = G.GxE, model = 'RKHS')
      } else if (m == "GW.GxW") {
        ETA[[1]] <- list(K = G.G, model = 'RKHS')
        ETA[[2]] <- list(K = G.W, model = 'RKHS')
        ETA[[3]] <- list(K = G.GxW, model = 'RKHS')
      } else if (m == "GEW.GxE.GxW") {
        ETA[[1]] <- list(K = G.G, model = 'RKHS')
        ETA[[2]] <- list(K = G.E, model = 'RKHS')
        ETA[[3]] <- list(K = G.W, model = 'RKHS')
        ETA[[4]] <- list(K = G.GxE, model = 'RKHS')
        ETA[[5]] <- list(K = G.GxW, model = 'RKHS')
      }
      
      fm <- BGLR(
        y = yNA,
        ETA = ETA,
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        verbose = F,
        saveAt = rep
      )
      a = a + 1
      res[a, "rep"] <- rep
      res[a, "nEnv"] <- nEnv
      res[a, "trait"] <- trait
      res[a, "env"] <- Env
      res[a, "model"] <- m
      res[a, "PA"] <-
        round(cor(fm$yHat[valSet], y[valSet], use = "na.or.complete"), 2)
      res[a, "PA_r"] <-
        round(cor(fm$yHat[valSet],
                  y[valSet],
                  method = "spearman",
                  use = "na.or.complete"),
              2)
    }
  }
}

write.csv(res, paste(rep, "Acc_CV_Ran.csv", sep = "_"))
