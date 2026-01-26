# # 1 loading packages ---------------------------------------
packages <- c("RandomFields"
              , "fields"
              , "Rcpp"
              , "data.table"
              , "dplyr"
              , "Matrix"
              , "tidyr"
              , "plyr"
              , "MASS"
              , "Hmisc"
              , "parallel"
              , "knnWLLR"
              , "tensorflow"
              , "lubridate"
              # , "latex2exp"
              , "np"
              # , "tensor"
              # , "CompRandFld"
              # , "PLRModels"
              # , "profvis"
              # , "RODBC"
              , "mvtnorm"
              # , "gstat"
              # , "RColorBrewer"
              , "reticulate"
              # , "SpecsVerification"
              # , "scoringRules"
              # , "verification"
            ) 
# 2  library
for(i in 1:length(packages))
{
  if(!lapply(packages[i], require,
             character.only = TRUE)[[1]])
  {
    install.packages(packages[i])
    # library(packages[i])
    lapply(packages[i], require,
           character.only = TRUE)
  }else{lapply(packages[i], require,
               character.only = TRUE)}
}
# x=lapply(packages, require, character.only = TRUE)
# rm(list=ls())
rm(i, packages)


