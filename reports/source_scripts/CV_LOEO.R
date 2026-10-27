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



rm(list = ls())
# load packages ----------------------------------------------------------------
library(magrittr)
library(BGLR)

# load data ---------------------------------------------------------------

load("data_CV_LOEO.RData")

# CV - LOEO ---------------------------------------------------------------

G.W <- G.W24
G.GxW <- G.GxW24

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
  trait = NA,
  env = NA,
  model = NA,
  PA = NA,
  PA_r = NA
)
a = 0

for (trait in traits) {
  print(trait)
  for (Env in unique(pheno$ENV)) {
    print(Env)
    y <- tibble::deframe(pheno[, trait])
    yNA <- y
    valSet <- pheno$ENV == Env
    yNA[valSet] <- NA
    
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
        verbose = F
      )
      a = a + 1
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



