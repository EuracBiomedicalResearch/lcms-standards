library(Rdisop)
library(MetaboCoreUtils)


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



H_rule <- function(formulas = character(0)){
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  tmp.frm1[is.na(tmp.frm1)] <- 0
  
  tmp.frm1$formula[tmp.frm1$H <= tmp.frm1$C*2 + tmp.frm1$N + 2]
}



N_rule <- function(mzval = numeric(0),
                   formulas = character(0)){
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  tmp.frm1[is.na(tmp.frm1)] <- 0
  
  tmp.frm1$formula[tmp.frm1$N %% 2 == round(mzval) %% 2]
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
