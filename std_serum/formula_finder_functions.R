library(Rdisop)
library(MetaboCoreUtils)

getmzneut <- function(mz = numeric(0),
                      adduct = c("[M+H]+", "[M+Na]+", 
                                 "[M-H]-", "[M-H+HCOOH]-", "[M+Cl]-", "[M-H+HCOONa]-",
                                 "[2M-H]-", "[2M+H]+",
                                 "[M+H-H2O]+")){
  addtb <- rbind(
    c("[M+H]+", 1.007276),
    c("[M+Na]+", 22.98980),
    c("[M-H]-", -1.007276),
    c("[M-H+HCOOH]-", 44.99820),
    c("[M+Cl]-", 34.96885),
    c("[M-H+HCOONa]-", 66.98017),
    c("[2M+H]+", 1.007276),
    c("[2M-H]-", -1.007276),
    c("[M+H-H2O]+", -17.00328)
  )
  if(grepl("2M", adduct)){
    (mz - as.numeric(addtb[addtb[,1] == adduct,2]))/2
  } else {
    mz - as.numeric(addtb[addtb[,1] == adduct,2])
  }
}


getform <- function(mzval = numeric(0),
                    elements = "CHNOPS",
                    ppm = 10){
  Rdisop::decomposeMass(mzval, 
                        ppm = ppm,
                        elements = initializeElements(
                          unlist(strsplit(elements, 
                                          "(?<=.)(?=[[:upper:]])", 
                                          perl=TRUE)))
  )$formula
}



update_form <- function(formulas = character(0),
                        adduct = c("[M+H]+", "[M+Na]+", 
                                   "[M-H]-", "[M-H+HCOOH]-", "[M+Cl]-", "[M-H+HCOONa]-",
                                   "[2M-H]-", "[2M+H]+",
                                   "[M+H-H2O]+"),
                        action = c("add", "remove")){
  addtb <- data.frame(rbind(
    c("[M+H]+", "H"),
    c("[M+Na]+", "Na"),
    c("[M-H]-", "H"),
    c("[M-H+HCOOH]-", "CH2O2"),
    c("[M+Cl]-", "Cl"),
    c("[M-H+HCOONa]-", "CHO2Na"),
    c("[2M+H]+", "H"),
    c("[2M-H]-", "H"),
    c("[M+H-H2O]+", "H2O")
  ))
  colnames(addtb) <- c("adduct", "formulas")
  addtb <- cbind(
    addtb, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(addtb$formula), 
                   MetaboCoreUtils::countElements)))
  addtb$H[grep("M-H", addtb$adduct)] <- 
    addtb$H[grep("M-H", addtb$adduct)] - 1
  addtb$H[addtb$adduct == "[M-H]-"] <- -1
  addtb$H[addtb$adduct == "[2M-H]-"] <- -1
  addtb$H[addtb$adduct == "[M+H-H2O]+"] <- -1
  addtb$O[addtb$adduct == "[M+H-H2O]+"] <- -1
  
  tmp.frm1 <- data.frame(formulas)
  tmp.frm1$formulas <- as.character(tmp.frm1$formulas)
  #tmp.frm1 <- rbind(addtb[addtb[, 1] == adduct, 2], tmp.frm1)
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  
  addtb[setdiff(names(tmp.frm1), names(addtb))] <- NA
  tmp.frm1[setdiff(names(addtb), names(tmp.frm1))] <- NA
  tmp.frm1 <- rbind(addtb[addtb$adduct == adduct,], tmp.frm1)
  tmp.frm1 <- tmp.frm1[,-1]
  tmp.frm1[is.na(tmp.frm1)] <- 0
  idx <- ncol(tmp.frm1)
  tmp.frm1$formok <- NA
  if(action == "add"){
    if(adduct %in% c("[2M+H]+", "[2M-H]-")){
      for(i in 2:nrow(tmp.frm1)){
        tmp.frm1[i,2:idx] <- tmp.frm1[i,2:idx]*2
      }
    }
    for(i in 2:nrow(tmp.frm1)){
      tmp.frm1[i, 2:idx] <- tmp.frm1[i, 2:idx] + tmp.frm1[1, 2:idx]
      tmp.frm1$formok[i] <- pasteElements(tmp.frm1[i, 2:idx])
    }
    } else if(action == "remove"){
      for(i in 2:nrow(tmp.frm1)){
        tmp.frm1[i, 2:idx] <- tmp.frm1[i, 2:idx] - tmp.frm1[1, 2:idx]
      }
      if(adduct %in% c("[2M+H]+", "[2M-H]-")){
        for(i in 2:nrow(tmp.frm1)){
          tmp.frm1[i,2:idx] <- tmp.frm1[i,2:idx]/2
        }
      }
      for(i in 2:nrow(tmp.frm1)){
        tmp.frm1$formok[i] <- pasteElements(tmp.frm1[i, 2:idx])
      }
    }
  tmp.frm1$formok[-1]
}


adduct_filter <- function(formulas = character(0),
                          adduct = c("[M+H]+", "[M+Na]+", 
                                     "[M-H]-", "[M-H+HCOOH]-", "[M+Cl]-", "[M-H+HCOONa]-",
                                     "[2M-H]-", "[2M+H]+",
                                     "[M+H-H2O]+")
                          ){
  addtb <- rbind(
    c("[M+H]+", "H"),
    c("[M+Na]+", "Na"),
    c("[M-H]-", "H"),
    c("[M-H+HCOOH]-", "CH2O2"),
    c("[M+Cl]-", "Cl"),
    c("[M-H+HCOONa]-", "CHO2Na"),
    c("[2M+H]+", "H"),
    c("[2M-H]-", "H"),
    c("[M+H-H2O]+", "H")
  )
  tmp.frm1 <- data.frame(formulas)
  tmp.frm1$formulas <- as.character(tmp.frm1$formulas)
  tmp.frm1 <- rbind(addtb[addtb[, 1] == adduct, 2], tmp.frm1)
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  tmp.frm1[is.na(tmp.frm1)] <- 0
  if(adduct == "[2M+H]+"){
    tmp.frm1 <- tmp.frm1[-1,]
    tmp.frm1$H <- tmp.frm1$H - 1
    for(i in 2:ncol(tmp.frm1)){
      tmp.frm1 <- tmp.frm1[tmp.frm1[,i] %% 2 == 0, ]
    }
  } else if(adduct == "[2M-H]-"){
    tmp.frm1 <- tmp.frm1[-1,]
    tmp.frm1$H <- tmp.frm1$H + 1
    for(i in 2:ncol(tmp.frm1)){
      tmp.frm1 <- tmp.frm1[tmp.frm1[,i] %% 2 == 0, ]
    }
  } else {
    for(i in 2:ncol(tmp.frm1)){
      tmp.frm1 <- tmp.frm1[tmp.frm1[,i] >= tmp.frm1[1,i] , ]
      tmp.frm1 <- tmp.frm1[-1,]
    }
  }
  tmp.frm1$formulas
}


H_rule <- function(formulas = character(0)){
  
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1 <- rbind(tmp.frm1, "CHNOPS")
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  tmp.frm1[is.na(tmp.frm1)] <- 0
  tmp.frm1 <- tmp.frm1[-nrow(tmp.frm1),]
  
  tmp.frm1$rule <- tmp.frm1$H <= tmp.frm1$C*2 + tmp.frm1$N + 2
  tmp.frm1$rule[tmp.frm1$P > 0] <- TRUE
  tmp.frm1$formula[tmp.frm1$rule]
}



H_rule2 <- function(mzval = numeric(0),
                   formulas = character(0)){
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1 <- rbind(tmp.frm1, "CHNOPS")
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  tmp.frm1[is.na(tmp.frm1)] <- 0
  tmp.frm1 <- tmp.frm1[-nrow(tmp.frm1),]
  
  tmp.frm1$rule <- (trunc(mzval) + tmp.frm1$H) %% 4 == 0
  tmp.frm1$rule[tmp.frm1$P > 0] <- TRUE
  tmp.frm1$formula[tmp.frm1$rule]
}



N_rule <- function(mzval = numeric(0),
                   formulas = character(0)){
  if(mzval < 500){
    tmp.frm1 <- data.frame(formula = formulas)
    tmp.frm1$formula <- as.character(tmp.frm1$formula)
    tmp.frm1 <- rbind(tmp.frm1, "N")
    tmp.frm1 <- cbind(
      tmp.frm1, 
      do.call(dplyr::bind_rows, 
              lapply(as.character(tmp.frm1$formula), 
                     MetaboCoreUtils::countElements)))
    tmp.frm1[is.na(tmp.frm1)] <- 0
    tmp.frm1 <- tmp.frm1[-nrow(tmp.frm1),]
    
    tmp.frm1$formula[tmp.frm1$N %% 2 == round(mzval) %% 2]
  } else{
    #return(formulas)
    print("Warning! Rule not applied because mzval > 500 Da! :P")
  }
}


RPU_rule <- function(formulas = character(0)){
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1 <- rbind(tmp.frm1, "CHNOPSSiClFlI")
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  tmp.frm1[is.na(tmp.frm1)] <- 0
  tmp.frm1 <- tmp.frm1[-nrow(tmp.frm1),]
  
  # RPU <- (C+Si) - (H+Cl+Fl+I)/2 + (N+P)/2 + 1
  tmp.frm1$formula[(
    (tmp.frm1$C + tmp.frm1$Si) - 
      ((tmp.frm1$H + tmp.frm1$Cl + tmp.frm1$Fl + tmp.frm1$I)/2) + 
      ((tmp.frm1$N + tmp.frm1$P)/2) + 
      1
  )  >= 0]
}



isotope_rule <- function(formulas = character(0),
                         intensity = numeric(0),
                         range = 0.4,
                         isotope = 1){
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1$A1 <- NA
  for(i in 1:nrow(tmp.frm1)){
    tmpx <- Rdisop::getMolecule(tmp.frm1$formula[i])$isotopes[[1]]
    tmp.frm1$A1[i] <- (tmpx[2,isotope+1]/tmpx[2,1])*100
  }
  A1_range <- intensity + (intensity * range * c(-1, 1))
  tmp.frm1$formula[tmp.frm1$A1 >= A1_range[1] & 
                     tmp.frm1$A1 <= A1_range[2]]
}



rank_form <- function(formulas = character(0),
                      i1 = numeric(0),
                      i2 = numeric(0)){
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1$A1 <- NA
  tmp.frm1$A2 <- NA
  for(i in 1:nrow(tmp.frm1)){
    tmpx <- Rdisop::getMolecule(tmp.frm1$formula[i])$isotopes[[1]]
    tmp.frm1$A1[i] <- (tmpx[2, 2] / tmpx[2, 1])*100
    tmp.frm1$A2[i] <- (tmpx[2, 3] / tmpx[2, 1])*100
  }
  tmp.frm1$A1d <- abs(tmp.frm1$A1 - i1)
  tmp.frm1$A2d <- abs(tmp.frm1$A2 - i2)
  tmp.frm1$Ad <- NA
  for(i in 1:nrow(tmp.frm1)){
    tmp.frm1$Ad[i] <- mean(c(tmp.frm1$A1d[i], tmp.frm1$A2d[i]))
  }
  tmp.frm1 <- tmp.frm1[order(tmp.frm1$Ad), ]
  tmp.frm1$rank <- seq(nrow(tmp.frm1))
  return(tmp.frm1[,c("formula", "rank")])
}




