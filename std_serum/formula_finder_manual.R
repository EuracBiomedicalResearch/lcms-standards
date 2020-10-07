k <- 1 # MIX nº k
z <- 1 # POS = 1 ; NEG = 2
cmp <- "dimethylglycine" # abbreviation

### Start auto #########################################################

MZML_PATH <- "C:/Users/mgarciaaloy/Documents/mzML/"
#MZML_PATH <- "C:/Users/vveri/Documents/stds_ID_files_trial/"
polarity.all <- c("POS", "NEG") 

library(xcms)
library(Rdisop)
library(CompoundDb)
library(doParallel)
library(magrittr)
library(SummarizedExperiment)
library(MsCoreUtils)
library(MetaboCoreUtils)

C_rule <- function(x, error = 0.1){
  C <- round(c((x/1.1) - ((x/1.1)*error), (x/1.1) + ((x/1.1)*error)))
  return(C)
}
H_rule <- function(C, N){
  H <- 2*C + N + 2
  return(H)
}
RPU_rule <- function(C=0, H=0, N=0, P=0, Si=0, Cl=0, Fl=0, I=0){
  RPU <- (C+Si) - (H+Cl+Fl+I)/2 + (N+P)/2 + 1
  return(RPU)
}

std_info <- read.table(
  "../data/standards_dilution.txt",
  sep = "\t", header = TRUE, as.is = TRUE)

std_info$mz_POS <- NA
std_info$mz_NEG <- NA
for(i in 1:nrow(std_info)){
  std_info$mz_POS[i] <- unlist(
    CompoundDb::mass2mz(getMolecule(as.character(std_info$formula[i]))$exactmass, 
            adduct = as.character(std_info$POS[i])))
  std_info$mz_NEG[i] <- unlist(
    CompoundDb::mass2mz(getMolecule(as.character(std_info$formula[i]))$exactmass, 
            adduct = as.character(std_info$NEG[i])))
}
std_info$form_pred_POS <- NA
std_info$form_pred_NEG <- NA

mix.levels <- levels(factor(std_info$mix))

injections <- read.table("../data/std_serum_files.txt", 
                         sep = "\t", header = TRUE, as.is = TRUE)
injections <- injections[injections$mode == "FS", ]
injections <- injections[grep("Water_Mix", injections$class), ]
injections <- injections[grep("_High", injections$class), ]
injections <- injections[grep("_2_", injections$mzML), ]

cwp <- CentWaveParam(
  ppm = 50,
  peakwidth = c(2, 20),
  snthresh = 10,
  mzdiff = 0.001,
  prefilter = c(3, 500),
  noise = 100,
  integrate = 2)

# Select STDs used in mix "k"
std_info.i <- std_info[std_info$mix == mix.levels[k], ]
if(mix.levels[k] == 5){
  std_info.i <- std_info.i[std_info.i$name != 
                             "Glycerophospho-inositol", ]
} else if(mix.levels[k] == 9){
  std_info.i <- std_info.i[std_info.i$name != 
                             "Eplerone", ]
} else if(mix.levels[k] == 17){
  std_info.i <- std_info.i[std_info.i$name != 
                             "3-Hydroxy-DL-kynurenine", ]
}

if(polarity.all[z] == "NEG"){ 
  pol <- -1
  ion <- -1.007276 
} else if(polarity.all[z] == "POS"){ 
  pol <- 1
  ion <- 1.007276 
}
idx <- injections$type == paste0(
  "Mix",
  sprintf(paste0("%0", ceiling(log10(10 + 1L)), "d"), 1:20)[k])
filename <- paste0(injections$folder[idx], "/", injections$mzML[idx])
#filename <- injections$mzML[idx]
filename <- filename[grep(polarity.all[z], filename)]
data_raw <- readMSData(paste0(MZML_PATH, filename),
                       mode = "onDisk")
std_info.i$mz <- std_info.i[, which(colnames(std_info.i) == 
                                      paste0("mz_",polarity.all[z]))]
std_info.i <- std_info.i[!is.na(std_info.i$mz),]

i <- which(std_info.i$abbreviation == cmp)
eic <- chromatogram(data_raw, aggregationFun = "max", 
                    mz = std_info.i$mz[i] + 0.01 * c(-1, 1))
pks <- data.frame(chromPeaks(findChromPeaks(eic, param = cwp)))

### Finish auto #######################################################



# Is there any detected peak?
any(which((pks$rtmin - 10) < std_info.i$RT[i] & (pks$rtmax + 10) > std_info.i$RT[i]) ==  which.max(pks$maxo))

# Is the peak with a close RT the one with the highest height?
which((pks$rtmin - 10) < std_info.i$RT[i] & (pks$rtmax + 10) > std_info.i$RT[i]) ==  which.max(pks$maxo)



### Start auto #########################################################

pks <- pks[which((pks$rtmin - 30) < std_info.i$RT[i] & (pks$rtmax + 30) > std_info.i$RT[i]), ]
if(nrow(pks) > 1){pks <- pks[which.max(pks$maxo), ]}
sps <- data_raw[[closest(pks$rt, rtime(data_raw),
                         duplicates = "closest")]]
sps2 <- as.data.frame(sps)
sps2$i100 <- (sps2$i / max(sps2$i))*100
idx <- c(which.min((abs(std_info.i$mz[i] - sps2$mz) / std_info.i$mz[i])*1e6),
         which.min((abs((std_info.i$mz[i] + 1.003355) - sps2$mz) / (std_info.i$mz[i] + 1.003355))*1e6),
         which.min((abs((std_info.i$mz[i] + 1.003355*2) - sps2$mz) / (std_info.i$mz[i] + 1.003355*2))*1e6))
sps2 <-sps2[idx,]
masses <- sps2$mz
intensities <- sps2$i

if(grepl("Na", std_info.i[i, which(colnames(std_info.i) == 
                                   polarity.all[z])])){
  molecules <- decomposeIsotopes(masses, intensities, ppm = 20,
                                 elements = initializeElements(
                                   c("C", "H", "N", "O", 
                                     "P", "S", "Na")))
}else{
    molecules <- decomposeIsotopes(masses, intensities, ppm = 20)
    }
formulas <- data.frame(formula = getFormula(molecules))
tmp <- lapply(as.character(formulas$formula), countElements)
if(class(tmp) == "numeric"){
  formulas$C <- tmp["C"]
  formulas$H <- tmp["H"]
  formulas$O <- tmp["O"]
  formulas$N <- tmp["N"]
  formulas$S <- tmp["S"]
  formulas$P <- tmp["P"]
  formulas$Na <- tmp["Na"]
  formulas[is.na(formulas)] <- 0
} else{
  formulas$C <- NA
  formulas$H <- NA
  formulas$O <- NA
  formulas$N <- NA
  formulas$S <- NA
  formulas$P <- NA
  formulas$Na <- NA
  for(j in 1:nrow(formulas)){
    if(length(tmp[[j]][names(tmp[[j]]) == "C"]) == 1){
      formulas$C[j] <- tmp[[j]][names(tmp[[j]]) == "C"]  
    } else {formulas$C[j] <- 0}
    if(length(tmp[[j]][names(tmp[[j]]) == "H"]) == 1){
      formulas$H[j] <- tmp[[j]][names(tmp[[j]]) == "H"]
    } else {formulas$H[j] <- 0}
    if(length(tmp[[j]][names(tmp[[j]]) == "O"]) == 1){
      formulas$O[j] <- tmp[[j]][names(tmp[[j]]) == "O"]
    } else {formulas$O[j] <- 0}
    if(length(tmp[[j]][names(tmp[[j]]) == "N"]) == 1){
      formulas$N[j] <- tmp[[j]][names(tmp[[j]]) == "N"]
    } else {formulas$N[j] <- 0}
    if(length(tmp[[j]][names(tmp[[j]]) == "S"]) == 1){
      formulas$S[j] <- tmp[[j]][names(tmp[[j]]) == "S"]
    } else {formulas$S[j] <- 0}
    if(length(tmp[[j]][names(tmp[[j]]) == "P"]) == 1){
      formulas$P[j] <- tmp[[j]][names(tmp[[j]]) == "P"]
    } else {formulas$P[j] <- 0}
    if(length(tmp[[j]][names(tmp[[j]]) == "Na"]) == 1){
      formulas$Na[j] <- tmp[[j]][names(tmp[[j]]) == "Na"]
    } else {formulas$Na[j] <- 0}
  } # close formula "j"
}


formulas$C_cntrb <- (formulas$C * 1.07) / ((formulas$C * 1.07) + (formulas$H * 0.012) + 
  (formulas$O * 0.038) + (formulas$N * 0.37) + (formulas$S * 0.76))
formulas$C_exp <- (((intensities[2] / intensities[1])*100) * formulas$C_cntrb) / 1.1
C_error <- 0.25
formulas$C_min <- formulas$C_exp - (formulas$C_exp * C_error)
formulas$C_max <- formulas$C_exp + (formulas$C_exp * C_error)
formulas$C_rule <- (formulas$C > formulas$C_min) & (formulas$C < formulas$C_max)

formulas$A1 <- ((formulas$C * 1.07) + (formulas$H * 0.012) + 
                        (formulas$O * 0.038) + (formulas$N * 0.37) + (formulas$S * 0.76))
formulas$A1_min <- formulas$A1 - (formulas$A1 * C_error)
formulas$A1_max <- formulas$A1 + (formulas$A1 * C_error)
A1_x100 <- ((intensities[2] / intensities[1])*100)
formulas$A1_rule <- (A1_x100 > formulas$A1_min) & (A1_x100 < formulas$A1_max)

formulas$S_cntrb <- (formulas$S * 4.295776) / (
  (formulas$C * 0.011449) + (formulas$H * 1.44e-06) + 
    (formulas$O * 1.444e-05) + (formulas$N * 0.001369) + (formulas$S * 4.295776))
formulas$S_exp <- (((intensities[3] / intensities[1])*100) * formulas$S_cntrb) / 4.295776
S_error <- 0.5
formulas$S_min <- formulas$S_exp - (formulas$S_exp * S_error)
formulas$S_max <- formulas$S_exp + (formulas$S_exp * S_error)
formulas$S_rule <- (formulas$S >= formulas$S_min) & (formulas$S <= formulas$S_max)

formulas$H_rule <- FALSE
formulas$H_rule[formulas$H <= H_rule(formulas$C, formulas$N)] <- TRUE
formulas$N_rule <- formulas$N %% 2 == round(masses[1] - ion) %% 2
formulas$RPU_rule <- FALSE
formulas$RPU_rule[RPU_rule(C = formulas$C, H = formulas$H - pol, 
                           N = formulas$N, P = formulas$P) >= 0] <- TRUE
if(grepl("Na", std_info.i[i, which(colnames(std_info.i) == 
                                   polarity.all[z])])){
  formulas$RPU_rule[RPU_rule(C = formulas$C, H = formulas$H, 
                             N = formulas$N, P = formulas$P) %% 1 > 0] <- FALSE
} else{
  formulas$RPU_rule[RPU_rule(C = formulas$C, H = formulas$H - pol, 
                             N = formulas$N, P = formulas$P) %% 1 > 0] <- FALSE 
}
if(grepl("Na", std_info.i[i, which(colnames(std_info.i) == 
                                   polarity.all[z])])){
  formulas.ok <- formulas[formulas$C_rule & formulas$H_rule & 
                            formulas$N_rule & formulas$RPU_rule & 
                            formulas$A1_rule & formulas$S_rule & 
                            formulas$Na == 1, ]
}else{
  formulas.ok <- formulas[formulas$C_rule & formulas$H_rule & 
                            formulas$N_rule & formulas$RPU_rule & 
                            formulas$A1_rule & formulas$S_rule, ]
}


### Finish auto #######################################################

# How many formulas "OK" there are?
nrow(formulas.ok)


formulas.ok$isodev <- NA
for(j in 1:nrow(formulas.ok)){
  formulas.ok$isodev[j] <- mean(abs(
    ((getIsotope(getMolecule(as.character(formulas.ok$formula[j])), 
                seq(1,3))[2,] / 
       max(getIsotope(getMolecule(as.character(formulas.ok$formula[j])), 
                      seq(1,3))[2,]))*100)[2] - 
      (intensities / max(intensities)*100)[2]))
} # close formula "j"
formulas.ok <- formulas.ok[order(formulas.ok$isodev), ]

if(nrow(formulas.ok) < 6){
  std_info[
    which(std_info$name == std_info.i$name[i]),
    which(colnames(std_info) == 
            paste0("form_pred_", polarity.all[z]))] <- 
    paste(formulas.ok$formula, collapse = "; ") 
} else{
  std_info[
    which(std_info$name == std_info.i$name[i]),
    which(colnames(std_info) == 
            paste0("form_pred_", polarity.all[z]))] <- 
    paste0(paste(formulas.ok$formula[1:5], collapse = "; "), "; etc.")
}



