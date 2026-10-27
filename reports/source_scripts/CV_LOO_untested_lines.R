rm(list = ls())
# load packages ----------------------------------------------------------------
library(magrittr)
library(BGLR)

# load data ---------------------------------------------------------------

load("data_CV_LOEO.RData")

# CV - LOEO ---------------------------------------------------------------

# Common lines
cl <- names(table(pheno$GID) >= 15)
cl <- as.numeric(cl[which(table(pheno$GID) >= 15)])

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
  line = NA,
  model = NA,
  PA = NA,
  PA_r = NA
)

pred <- list()

a = 0

for (trait in traits) {
  print("-----------------------")
  print(trait)
  for (L in cl) {
    print(L)
    y <- tibble::deframe(pheno[, trait])
    yNA <- y
    valSet <- pheno$GID == L
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
      res[a, "line"] <- L
      res[a, "model"] <- m
      res[a, "PA"] <-
        round(cor(fm$yHat[valSet], y[valSet], use = "na.or.complete"), 2)
      res[a, "PA_r"] <-
        round(cor(fm$yHat[valSet],
                  y[valSet],
                  method = "spearman",
                  use = "na.or.complete"),
              2)
      pred[[a]] <- pheno[valSet,] %>%
        select(GID,
               EnvID,
               ENV,
               trait) %>%
        rename(y = trait) %>%
        mutate(
          ebv = fm$yHat[valSet],
          trait = trait,
          line = L,
          model = m
        )
    }
  }
}
