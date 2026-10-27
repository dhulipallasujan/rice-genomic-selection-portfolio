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



# load packages ----------------------------------------------------------------
library(magrittr)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# load data ---------------------------------------------------------------

# The conversion of the additional file 3 into numeric format (allele count)
# was done with Tassel (Data --> Numerical Genotype)
# The first column of the output file from Tassel was renamed GID

geno <- read_tsv("Additional file 3.num.txt")
pheno <- read_csv("Additional file 4.csv")
EC_24 <- read_csv("Additional file 5.csv")

# format data -------------------------------------------------------------

#Genotypic information
geno <- geno %>%
  filter(GID %in% unique(pheno$GID)) %>%
  arrange(GID)

gid <- intersect(unique(pheno$GID), geno$GID)

#Phenotypic information
pheno <- pheno %>%
  filter(!is.na(BLUPs)) %>%
  pivot_wider(
    names_from = Trait,
    values_from = BLUPs,
    values_fill = NA,
    values_fn = mean
  ) %>%
  mutate(EnvID = paste(Country, Location, Year, Season, sep = "_"))

# Check
length(unique(pheno$GID))
unique(pheno$GID) %in% geno$GID
table(pheno$EnvID)


# Env data
Env <- data.frame(
  EnvID =
    c(
      "Bangladesh_Gazipur_2019_Wet",
      "Bangladesh_Nizmawna_2019_Wet",
      "India_Hyderabad_2018_Wet",
      "India_Cuttack_2019_Wet",
      "India_Hyderabad_2019_Dry",
      "India_Hyderabad_2019_Wet",
      "India_Maruteru_2019_Wet",
      "India_Raipur_2019_Wet",
      "Kenya_Ahero_2019_Dry",
      "Kenya_Mwea_2020_Wet",
      "Mozambique_Chokwe_2020_Wet",
      "Mozambique_Maputo_2020_Wet",
      "Philippines_Los Banos_2019_Dry",
      "Philippines_Los Banos_2019_Wet",
      "Tanzania_Dakawa_2020_Wet"
    ),
  ENV =
    c(
      "BD_GZ_19W",
      "BD_NM_19W",
      "IN_HY_18W",
      "IN_CU_19W",
      "IN_HY_19D",
      "IN_HY_19W",
      "IN_MA_19W",
      "IN_RP_19W",
      "KE_AH_19D",
      "KE_MW_20W",
      "MZ_CK_20W",
      "MZ_MP_20W",
      "PH_LB_19D",
      "PH_LB_19W",
      "TZ_DK_20W"
    )
)

pheno <- pheno %>% left_join(Env)

EC_24l <- pheno %>% select(ENV) %>%
  left_join(EC_24)




# Co-var structure of E  --------------------------------------------------

Z <- as.matrix(model.matrix(~ tibble::deframe(pheno[, "ENV"]) - 1))
d <- colSums(Z)
V <- Z
for (i in 1:ncol(Z)) {
  V[, i] <- V[, i] / sqrt(d[i])
}
EVD.E <- list(vectors = V, values = d)
G.E <- tcrossprod(Z)

# Co-var structure of W ---------------------------------------------------

## Imputing, centering, standarizing.
# Parameters for centering and standarizing
covarW <- function(W) {
  W <- as.matrix(W)
  G.W <- tcrossprod(scale(W)) / ncol(W)
  EVD.W <- eigen(G.W)
  
  return(list(G.W, EVD.W))
}

covar24 <- covarW(EC_24l[, -1])

# Co-var structure G ------------------------------------------------------


IDs <- tibble::deframe(pheno[, "GID"])
X <- as.matrix(geno)

colNames <- colnames(X)[-1]
IDx <- X[, 1]
X <- X[, -1]

rownames(X) <- IDx
colnames(X) <- colNames

## Imputing, centering, standarizing.
# Parameters for centering and standarizing

G.G <- tcrossprod(scale(X)) / ncol(X)

IDs <- factor(IDs, levels = rownames(G.G))
Z <- as.matrix(model.matrix( ~ IDs - 1))
G.G <- tcrossprod(tcrossprod(Z, G.G), Z)  #This is ZGZ’

EVD.G <- eigen(G.G)
rownames(EVD.G$vectors) <- rownames(G.G)


# Co-var structure G*E ----------------------------------------------------

G.GxE <- G.G * G.E
EVD.GxE <- eigen(G.GxE)


# Co-var structure G*W ----------------------------------------------------
G.W24 <- covar24[[1]]
G.GxW24 <- G.G * G.W24

EVD.GxW24 <- eigen(G.GxW24)


# save data ---------------------------------------------------------------

save(geno, pheno, G.E, G.W24, G.G, G.GxE, G.GxW24, file = "data_CV_LOEO.RData")
save(geno, pheno, EC_24, file = "data_CV_RAN.RData")

