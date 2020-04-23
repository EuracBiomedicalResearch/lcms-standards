library(MetaboCoreUtils)
library(Rdisop)
load("isotopic_patterns_tables.RData")
rm(startpoint)
samples <- ls()
startpoint <- Sys.time()

std_info <- read.table(
  "../data/standards_dilution.txt",
  sep = "\t", header = TRUE, as.is = TRUE)
std_info <- std_info[, c("mix", "name", "abbreviation", "formula",
                         "POS", "NEG", "quality_POS", "quality_NEG")]
idx <- which(std_info$name == "3-Hydroxy-DL-kynurenine" & 
               std_info$mix == 17)
std_info <- std_info[-idx, ]
std_info <- std_info[-which(std_info$name == "Eplerone"), ]
idx <- which(std_info$name == "Glycerophospho-inositol")
std_info <- std_info[-idx, ]
std_info <- std_info[!is.na(std_info$quality_POS), ]
std_info <- std_info[std_info$quality_POS == "0_perfect" | 
                       std_info$quality_POS == "1_good", ]
std_info <- std_info[std_info$quality_NEG == "0_perfect" | 
                       std_info$quality_NEG == "1_good", ]
std_info <- std_info[std_info$POS == "[M+H]+" & 
                       std_info$NEG == "[M-H]-", ]
std_info$mzneut <- NA
for(i in 1:nrow(std_info)){
  std_info$mzneut[i] <- getMolecule(as.character(std_info$formula[i]))$exactmass
}
tmp <- lapply(as.character(std_info$formula), countElements)
tmp <- do.call(dplyr::bind_rows, tmp)
tmp[is.na(tmp)] <- 0
std_info <- cbind(std_info, tmp)


# Problematic samples
#### Paracetamol in sample PLS3: we've seen that the real A+1 ion is 
#### 153.1060, which mz deviation is 0.028
PLS3$i2[which(PLS3$abrv == "paracetamol")] <- NA
PLS3$mz2[which(PLS3$abrv == "paracetamol")] <- NA

# Calculate the experimental relative intensity of ion A+1
for(i in 1:length(samples)){
  tmp <- get(samples[i])
  tmp$A1 <- (tmp$i2 / tmp$i1)*100
  tmp$A1[tmp$i2 < 500] <- NA
  tmp$A2 <- (tmp$i3 / tmp$i1)*100
  tmp$A2[tmp$i3 < 100] <- NA
  tmp$A1[tmp$i1 > 5e5] <- NA
  if(substring(samples[i], 1, 1) == "P"){ion <- 1.007276
  } else if(substring(samples[i], 1, 1) == "N"){ion <- -1.007276}
  tmp$mz1dev <- (abs(tmp$mz1 - (std_info$mzneut + ion)) / 
                   (std_info$mzneut + ion))*1e6
  assign(samples[i], tmp)
}

# Put all data in a list
data <- c()
for(i in 1:length(samples)){
  data[[i]] <- get(samples[i])
  names(data)[[i]] <- samples[i]
}

# Calculate the theorical relative intensity of ion A+1
A1 <- do.call(cbind, lapply(data, function(x){x$A1}))
std_info$A1 <- ((std_info$C * 1.07) + (std_info$H * 0.012) + 
                  (std_info$O * 0.038) + (std_info$N * 0.37) + 
                  (std_info$S * 0.76))

std_info$A1_min <- apply(A1, 1, function(x){min(x, na.rm = TRUE)})
std_info$A1_max <- apply(A1, 1, function(x){max(x, na.rm = TRUE)})
for(i in 1:nrow(std_info)){
  std_info$A1_error_max[i] <- 
    max(abs((std_info$A1_min[i]-std_info$A1[i])*100) / std_info$A1[i], 
        abs((std_info$A1_max[i]-std_info$A1[i])*100)/std_info$A1[i])
}


# Calculate the theorical C contribution in A+1
std_info$C_cntrb <- (std_info$C * 1.07) / std_info$A1
C1 <- (A1*std_info$C_cntrb)/1.1
std_info$C1_min <- apply(C1, 1, function(x){min(x, na.rm = TRUE)})
std_info$C1_max <- apply(C1, 1, function(x){max(x, na.rm = TRUE)})
for(i in 1:nrow(std_info)){
  std_info$C1_error_max[i] <- 
    max((abs(std_info$C1_min[i] - std_info$C[i]) / std_info$C[i])*100, 
        (abs(std_info$C1_max[i] - std_info$C[i]) / std_info$C[i])*100)
}

# Calculate the ppm deviations
ppm <- do.call(cbind, lapply(data, function(x){x$mz1dev}))
std_info$ppm1_max <- apply(ppm, 1, function(x){max(x, na.rm = TRUE)})
plot(std_info$mzneut, std_info$ppm1_max, xlab = "neutral mz", 
     ylab = "max(ppm deviation)")
abline(lm(std_info$ppm1_max~std_info$mzneut), col = "grey", lty=2)

aggregate(std_info$ppm1_max, list(std_info$mzneut > 600), max)
aggregate(std_info$ppm1_max, list(std_info$mzneut > 250), max)


# Calculate the theorical relative intensity of ion A+2
A2 <- do.call(cbind, lapply(data, function(x){x$A2}))
std_info$A2 <- (
    (std_info$C * ((1.070*1.070)/100)) + 
    (std_info$H * ((0.012*0.012)/100)) + 
    (std_info$O * ((0.038*0.038)/100 + 0.2)) + 
    (std_info$N * ((0.370*0.370)/100)) + 
    (std_info$S * ((0.760*0.760)/100 + 4.29)) +
    (std_info$Cl * 24.22)
    )
std_info$A2_min <- apply(A2, 1, function(x){min(x, na.rm = TRUE)})
std_info$A2_max <- apply(A2, 1, function(x){max(x, na.rm = TRUE)})
for(i in 1:nrow(std_info)){
  std_info$A2_error_max[i] <- 
    max(abs((std_info$A2_min[i]-std_info$A2[i])*100) / std_info$A2[i], 
        abs((std_info$A2_max[i]-std_info$A2[i])*100)/std_info$A2[i])
}

# Calculate the theorical S contribution in A+2




# Formula finder
z <- 13
tmp <- PLW1
if(tmp$mz1[z] <= 250){myppm <- 60
} else if(tmp$mz1[z] > 250 & tmp$mz1[z] < 600){myppm <- 25
} else if(tmp$mz1[z] >= 600){myppm<- 15}
if(tmp$i1[z] < 5e5){myerror <- 0.5} else {myerror <- 1.5}

tmp$mz1[z]
tmp.frm1 <- data.frame(formula = decomposeMass(tmp$mz1[z], ppm = myppm)$formula) 
tmp.frm1 <- cbind(
  tmp.frm1, 
  do.call(dplyr::bind_rows, 
          lapply(as.character(tmp.frm1$formula), countElements)))
tmp.frm1[is.na(tmp.frm1)] <- 0

### Theoretical rules
tmp.frm1$H_rule <- tmp.frm1$H <= tmp.frm1$C*2 + tmp.frm1$N + 2
tmp.frm1$N_rule <- tmp.frm1$N %% 2 == round(tmp$mz1[z] - 1.007276) %% 2
tmp.frm1$RPU_rule <- # RPU <- (C+Si) - (H+Cl+Fl+I)/2 + (N+P)/2 + 1
  tmp.frm1$C - 
  (tmp.frm1$H - 1)/2 + 
  (tmp.frm1$N + tmp.frm1$P)/2  >= 0 

### Experimental rules
tmp.frm1$A1 <- ((tmp.frm1$C * 1.07) + (tmp.frm1$H * 0.012) + 
                  (tmp.frm1$O * 0.038) + (tmp.frm1$N * 0.37) + 
                  (tmp.frm1$S * 0.76))
tmp$i1[z] 
formatC(tmp$i1[z], format = "e", digits = 1) # < 5e5
tmp$i2[z] # > 500
tmp.A1 <- (tmp$i2[z] / tmp$i1[z]) * 100
(A1_range <- tmp.A1 + (tmp.A1 * myerror * c(-1, 1)))
tmp.frm1$A1_rule <- tmp.frm1$A1 >= A1_range[1] & tmp.frm1$A1 <= A1_range[2]

tmp.frm1$C_cntrb <- (tmp.frm1$C * 1.07) / tmp.frm1$A1
tmp.frm1$C_thr <- (tmp.A1*tmp.frm1$C_cntrb)/1.1
tmp.frm1$C_thr_min <- tmp.frm1$C_thr - (tmp.frm1$C_thr * myerror)
tmp.frm1$C_thr_max <- tmp.frm1$C_thr + (tmp.frm1$C_thr * myerror)
tmp.frm1$C_rule <- 
  tmp.frm1$C >= tmp.frm1$C_thr_min & 
  tmp.frm1$C <= tmp.frm1$C_thr_max


nrow(tmp.frm1)
sum((tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule))

table(tmp.frm1$A1_rule, (tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule))

table(tmp.frm1$C_rule, (tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule))

table(tmp.frm1$A1_rule, tmp.frm1$C_rule, (tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule))

tmp.frm1$formula[tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule &
                   tmp.frm1$A1_rule & tmp.frm1$C_rule]
std_info$name[z]
std_info$formula[z]
#z=1 (glutamic) C5H10NO4 vs C4H10N3OS 
#z=2 (xanthine) C5H5N4O2 vs C2H9N4O2S
#z=13 (glutathione oxidized)