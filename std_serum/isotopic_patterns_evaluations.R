library(MetaboCoreUtils)
library(Rdisop)
library(RColorBrewer)
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

#### Acetylornithine in sample NLS1: the code is taking the peak at 
#### 186 seconds, but the correct peak is the one at 170 seconds
NLS1$i1[which(NLS1$abrv == "acetylornithine")] <- NA
NLS1$mz1[which(NLS1$abrv == "acetylornithine")] <- NA

#### Tyrosine in sample NLS1: in reallity there it is not ionizating, 
#### the ion the code is taking is the 13C of a co-eluting compound with mz 179
NLS1$i1[which(NLS1$abrv == "tyrosine")] <- NA
NLS1$mz1[which(NLS1$abrv == "tyrosine")] <- NA


##### 3-Hydroxykynurenine
PLS1$i2[which(PLS1$abrv == "hydroxykynurenine_3")] <- NA
PLS1$mz2[which(PLS1$abrv == "hydroxykynurenine_3")] <- NA
PLS2$i2[which(PLS2$abrv == "hydroxykynurenine_3")] <- NA
PLS2$mz2[which(PLS2$abrv == "hydroxykynurenine_3")] <- NA
PLS3$i2[which(PLS3$abrv == "hydroxykynurenine_3")] <- NA
PLS3$mz2[which(PLS3$abrv == "hydroxykynurenine_3")] <- NA

for(i in 1:length(samples)){
  tmp <- get(samples[i])
  tmp$A1 <- (tmp$i2 / tmp$i1)*100 # Calculate the experimental relative intensity of ion A+1
  tmp$A1[tmp$i2 < 500] <- NA      # Exclude info of ions A+1 with an intensity < 500
  tmp$A2 <- (tmp$i3 / tmp$i1)*100 # Calculate the experimental relative intensity of ion A+2
  tmp$A2[tmp$i3 < 500] <- NA      # Exclude info of ions A+2 with an intensity < 500
  #tmp$A1[tmp$i1 > 5e5] <- NA      # Exclude saturated compounds
  #tmp$A2[tmp$i1 > 5e5] <- NA      # Exclude saturated compounds
  if(substring(samples[i], 1, 1) == "P"){ion <- 1.007276
  } else if(substring(samples[i], 1, 1) == "N"){ion <- -1.007276}
  tmp$mz1dev <- (abs(tmp$mz1 - (std_info$mzneut + ion)) / 
                   (std_info$mzneut + ion))*1e6
  tmp$mz1dalt <- abs(tmp$mz1 - (std_info$mzneut + ion))
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
std_info$A1 <- NA
for(i in 1:nrow(std_info)){
  tmp <- getMolecule(std_info$formula[i])$isotopes[[1]]
  std_info$A1[i] <- (tmp[2,2]/tmp[2,1])*100
}
#std_info$A1 <- ((std_info$C * 1.07) + (std_info$H * 0.012) + 
#                  (std_info$O * 0.038) + (std_info$N * 0.37) + 
#                  (std_info$S * 0.76))
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
dalton <- do.call(cbind, lapply(data, function(x){x$mz1dalt}))
i1 <- do.call(cbind, lapply(data, function(x){x$i1}))
mzvals <- do.call(cbind, lapply(data, function(x){x$mz1}))

ppm <- ppm[,substring(samples, 1, 1) != "N" &
             substring(samples, 3, 3) != "S"]

std_info$ppm1_max <- apply(ppm, 1, function(x){max(x, na.rm = TRUE)})
std_info$dalton_max <- apply(dalton, 1, function(x){max(x, na.rm = TRUE)})
plot(std_info$mzneut, std_info$ppm1_max, xlab = "neutral mass", 
     ylab = "max(ppm deviation)")
abline(lm(std_info$ppm1_max~std_info$mzneut), col = "grey", lty=2)

col_tmp <- brewer.pal(6, name = "Set1")
col_mtx <- matrix(ncol = ncol(ppm), nrow = nrow(ppm))
pch_mtx <- matrix(ncol = ncol(ppm), nrow = nrow(ppm))

par(mfrow=c(2,2))

col_mtx[,grep("N", substring(samples, 1, 1))] <- col_tmp[1]
col_mtx[,grep("P", substring(samples, 1, 1))] <- col_tmp[2]

col_mtx[,grep("H", substring(samples, 2, 2))] <- col_tmp[3]
col_mtx[,grep("L", substring(samples, 2, 2))] <- col_tmp[4]

col_mtx[,grep("S", substring(samples, 3, 3))] <- col_tmp[5]
col_mtx[,grep("W", substring(samples, 3, 3))] <- col_tmp[6]

pch_mtx[,grep("S", substring(samples, 3, 3))] <- 16
pch_mtx[,grep("W", substring(samples, 3, 3))] <- 8


plot(mzvals, ppm, col = paste0(col_mtx, 90), pch=16)
plot(mzvals, dalton, col = paste0(col_mtx, 90), pch=16)
plot(i1, ppm, col = paste0(col_mtx, 90), pch=16)
plot(i1, dalton, col = paste0(col_mtx, 90), pch=16)



legend("topright", col = col_tmp[1:2], legend = c("NEG", "POS"), pch=16)
legend("topright", col = col_tmp[3:4], legend = c("High", "Low"), pch=16)
legend("topright", col = col_tmp[5:6], legend = c("Serum", "Water"), pch=16)

plot(mzvals, ppm, col = paste0(col_mtx, 90), pch=pch_mtx)
plot(mzvals, dalton, col = paste0(col_mtx, 90), pch=pch_mtx)
plot(i1, ppm, col = paste0(col_mtx, 90), pch=pch_mtx)
legend("topright", col = col_tmp[1:2], legend = c("NEG", "POS"), pch=16)
plot(i1, dalton, col = paste0(col_mtx, 90), pch=pch_mtx)
legend("topright", col = "black", legend = c("Serum", "Water"), pch=c(16,8))


col_tmp <- brewer.pal(4, name = "Paired")
col_mtx[,] <- "#808080"
col_mtx[std_info$N >= 1 & std_info$N <= 4,] <- col_tmp[1]
col_mtx[std_info$N >= 5,] <- col_tmp[2]
par(mfrow=c(1,2))
plot(A1, C1, col = paste0(col_mtx, 90), pch=16)
legend("bottomright", col = c("#808080", col_tmp[1:2]), 
       legend = c("N = 0", "N = 1-4", "N >= 5"), pch=16)
col_mtx[,] <- "#808080"
col_mtx[std_info$S == 1,] <- col_tmp[3]
col_mtx[std_info$S == 2,] <- col_tmp[4]
plot(A1, C1, col = paste0(col_mtx, 90), pch=16)
legend("bottomright", col = c("#808080", col_tmp[3:4]), 
       legend = c("S = 0", "S = 1", "S = 2"), pch=16)

mzvals <- mzvals[,substring(samples, 3, 3) == "W"]
ppm <- ppm[,substring(samples, 3, 3) == "W"]
dalton <- dalton[,substring(samples, 3, 3) == "W"]
i1 <- i1[,substring(samples, 3, 3) == "W"]
col_mtx <- col_mtx[,substring(samples, 3, 3) == "W"]


aggregate(std_info$ppm1_max, list(std_info$mzneut > 600), max)
aggregate(std_info$ppm1_max, list(std_info$mzneut > 250), max)
max(std_info$ppm1_max)

# Calculate the theorical relative intensity of ion A+2
A2 <- do.call(cbind, lapply(data, function(x){x$A2}))
std_info$A2 <- NA
for(i in 1:nrow(std_info)){
  tmp <- getMolecule(std_info$formula[i])$isotopes[[1]]
  std_info$A2[i] <- (tmp[2,3]/tmp[2,1])*100
}
#std_info$A2 <- (
#    (std_info$C * ((1.070*1.070)/100)) + 
#    (std_info$H * ((0.012*0.012)/100)) + 
#    (std_info$O * ((0.038*0.038)/100 + 0.2)) + 
#    (std_info$N * ((0.370*0.370)/100)) + 
#    (std_info$S * ((0.760*0.760)/100 + 4.29)) +
#    (std_info$Cl * 24.22)
#    )
std_info$A2_min <- apply(A2, 1, function(x){min(x, na.rm = TRUE)})
std_info$A2_max <- apply(A2, 1, function(x){max(x, na.rm = TRUE)})
for(i in 1:nrow(std_info)){
  std_info$A2_error_max[i] <- 
    max(abs((std_info$A2_min[i]-std_info$A2[i])*100) / std_info$A2[i], 
        abs((std_info$A2_max[i]-std_info$A2[i])*100)/std_info$A2[i])
}

# Calculate the theorical S contribution in A+2




# Formula finder -----------------------------------------------------
z <- 1
tmp <- PLW2

formula_finder <- matrix(nrow=nrow(std_info), ncol = length(samples))
rownames(formula_finder) <- std_info$abbreviation
colnames(formula_finder) <- samples

formula_finder <- matrix(nrow=nrow(std_info), ncol = 7)
rownames(formula_finder) <- std_info$abbreviation
colnames(formula_finder) <- c("start", "thr_filters", "A1_filter", "C1_filter", "A2_filter", "rank", "rankp")

for(z in 1:10#:nrow(std_info)
    ){
  if(tmp$mz1[z] <= 250){myppm <- 40
  } else if(tmp$mz1[z] > 250 & tmp$mz1[z] < 600){myppm <- 25
  } else if(tmp$mz1[z] >= 600){myppm<- 15}
  if(tmp$i1[z] < 5e5){
    myerror_A1 <- 0.4
    myerror_A2 <- 2
  } else {
    myerror_A1 <- 1.5
    myerror_A2 <- 5
  }
  
  #tmp$mz1[z]
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
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1$A1 <- NA
  tmp.frm1$A2 <- NA
  for(i in 1:nrow(tmp.frm1)){
    tmpx <- getMolecule(tmp.frm1$formula[i])$isotopes[[1]]
    tmp.frm1$A1[i] <- (tmpx[2,2]/tmpx[2,1])*100
    tmp.frm1$A2[i] <- (tmpx[2,3]/tmpx[2,1])*100
  }
  #tmp.frm1$A1 <- ((tmp.frm1$C * 1.07) + (tmp.frm1$H * 0.012) + 
  #                  (tmp.frm1$O * 0.038) + (tmp.frm1$N * 0.37) + 
  #                  (tmp.frm1$S * 0.76))
  
  #tmp$i1[z] 
  #formatC(tmp$i1[z], format = "e", digits = 1) # < 5e5
  #tmp$i2[z] # > 500
  tmp.A1 <- (tmp$i2[z] / tmp$i1[z]) * 100
  (A1_range <- tmp.A1 + (tmp.A1 * myerror_A1 * c(-1, 1)))
  tmp.frm1$A1_rule <- tmp.frm1$A1 >= A1_range[1] & tmp.frm1$A1 <= A1_range[2]
  
  tmp.frm1$C_cntrb <- (tmp.frm1$C * 1.07) / tmp.frm1$A1
  tmp.frm1$C_thr <- (tmp.A1*tmp.frm1$C_cntrb)/1.07
  tmp.frm1$C_thr_min <- tmp.frm1$C_thr - (tmp.frm1$C_thr * myerror_A1)
  tmp.frm1$C_thr_max <- tmp.frm1$C_thr + (tmp.frm1$C_thr * myerror_A1)
  tmp.frm1$C_rule <- 
    tmp.frm1$C >= tmp.frm1$C_thr_min & 
    tmp.frm1$C <= tmp.frm1$C_thr_max
  
  tmp.A2 <- (tmp$i3[z] / tmp$i1[z]) * 100
  (A2_range <- tmp.A2 + (tmp.A2 * myerror_A2 * c(-1, 1)))
  if(A2_range[1] < 0){A2_range[1] <- 0}
  tmp.frm1$A2_rule <- tmp.frm1$A2 >= A2_range[1] & tmp.frm1$A2 <= A2_range[2]
  
  
  formula_finder[z, "start"] <- nrow(tmp.frm1)
  formula_finder[z, "thr_filters"] <- sum(tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule)
  formula_finder[z, "A1_filter"] <- sum(tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule & tmp.frm1$A1_rule)
  formula_finder[z, "C1_filter"] <- sum(tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule & tmp.frm1$A1_rule & tmp.frm1$C_rule)
  formula_finder[z, "A2_filter"] <- 
    sum(tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule & 
          tmp.frm1$A1_rule & tmp.frm1$C_rule & tmp.frm1$A2_rule)
  
  tmpx <- tmp.frm1[tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule & 
                    tmp.frm1$A1_rule & tmp.frm1$C_rule & tmp.frm1$A2_rule,]
  tmpx$A1x <- abs(tmpx$A1 - tmp.A1)
  tmpx$A2x <- abs(tmpx$A2 - tmp.A2)
  tmpx$A <- NA
  #tmpx$Aw <- NA
  for(i in 1:nrow(tmpx)){
    tmpx$A[i] <- mean(c(tmpx$A1x[i], tmpx$A2x[i]))
    #tmpx$Aw[i] <- weighted.mean(c(tmpx$A1x[i], tmpx$A2x[i]), c(2,1))
  }
  tmpx <- tmpx[order(tmpx$A),]
  tmpx$H <- tmpx$H - 1
  tmpx$fml <- NA
  tmpx$A1p <- NA
  tmpx$A2p <- NA
  tmpx$Ap <- NA
  for(i in 1:nrow(tmpx)){
    tmpx$fml[i] <- pasteElements(tmpx[i,2:7])
    tmpx$A1p[i] <- abs(100 - ((tmp.A1*100)/tmpx$A1[i]))
    tmpx$A2p[i] <- abs(100 - ((tmp.A2*100)/tmpx$A2[i]))
    tmpx$Ap[i] <- mean(c(tmpx$A1p[i], tmpx$A2p[i]))
  }
  tmpx <- tmpx[order(tmpx$A), ]
  formula_finder[z, "rank"] <- which(tmpx$fml == std_info$formula[z])
  tmpx <- tmpx[order(tmpx$Ap),]
  formula_finder[z, "rankp"] <- which(tmpx$fml == std_info$formula[z])
}


#table(tmp.frm1$A1_rule, (tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule))

#table(tmp.frm1$C_rule, (tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule))

#table(tmp.frm1$A1_rule, tmp.frm1$C_rule, (tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule))

#tmp.frm1$formula[tmp.frm1$H_rule & tmp.frm1$N_rule & tmp.frm1$RPU_rule &
#                   tmp.frm1$A1_rule & tmp.frm1$C_rule]
#std_info$name[z]
#std_info$formula[z]
#z=1 (glutamic) C5H10NO4 vs C4H10N3OS 
#z=2 (xanthine) C5H5N4O2 vs C2H9N4O2S
#z=13 (glutathione oxidized)


