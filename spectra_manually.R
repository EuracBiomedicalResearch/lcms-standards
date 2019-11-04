polarity <- "POS" # specify "POS" or "NEG"

study <- "standards_dilution" # specify "internal_standards" OR 
                              # "standards_dilution"
mixnum <- 17

source("R/which_within.R")
library(xcms)
library(magrittr)
library(CompoundDb)
library(Rdisop)

MZML_PATH <- "C:/Users/mgarciaaloy/Documents/mzML/"

# Get the information regarding the injection sequence:
injections <- read.table(paste0("data/", study, "_files.txt"), #import the file
                         sep = "\t", header = TRUE, as.is = TRUE)
myfiles <- injections$mzML # get file names
myfiles <- myfiles[grep(polarity, myfiles)]
if(study == "internal_standards"){
  (myfiles <- myfiles[grep(paste(10:15, collapse="|"), myfiles)])
} else if(study == "standards_dilution"){
  (myfiles <- c(myfiles[grep(paste0("MIX ", mixnum, "A"), myfiles)], 
                myfiles[grep(paste0("MIX ", mixnum, "B"), myfiles)]))
}


da = 0.01
RTd = 10

exclusion_list <- read.table("data/exclusion_mz.txt",
                             sep = "\t", header = TRUE, as.is = TRUE)
exclusion_list <- exclusion_list[, grep(polarity, colnames(exclusion_list))]


# Get the information regarding the standards that are in the samples:
std_info <- read.table(paste0("data/", study, ".txt"),
                       sep = "\t", header = TRUE, as.is = TRUE)
std_info$mzneut = NA
for(i in seq(nrow(std_info))){
  if(grepl("C", std_info$formula[i])){std_info$mzneut[i] = 
    getMolecule(as.character(std_info$formula[i]))$exactmass}else{
      std_info$mzneut[i] = as.numeric(std_info$formula[i])}
}
std_info$name <- c(substring(std_info$name, 1, 33))
if(study == "standards_dilution"){
  std_info <- subset(std_info, mix == mixnum)
}
std_info$mzneut = NA
for(i in seq(nrow(std_info))){
  if(grepl("C", std_info$formula[i])){std_info$mzneut[i] = 
    getMolecule(as.character(std_info$formula[i]))$exactmass}else{
      std_info$mzneut[i] = as.numeric(std_info$formula[i])}
}

mycompound <- "Epinephrine"
mycompound <- std_info[grep(mycompound, std_info$name),]
mzvalue <- unlist(mass2mz(mycompound$mzneut, 
                          adduct=
                            as.character(
                              mycompound[, 
                                         grep(polarity, colnames(mycompound))])
))
rtvalue <- mycompound$RT

for(j in seq(length(myfiles))){
  # Import data:
  raw_data <- readMSData(paste0(MZML_PATH, myfiles[j]), mode = "onDisk")
  
  # Get EIC:
  chr <- chromatogram(raw_data, aggregationFun = "max",
                      mz = c(mzvalue - da, mzvalue + da))
  #plot(chr)
  #plot(chr, xlim = c(chr@.Data[[1]]@rtime[which.max(
  #  chr@.Data[[1]]@intensity)] - RTd, chr@.Data[[1]]@rtime[which.max(
  #    chr@.Data[[1]]@intensity)] + RTd))
  
  # Get the spectrum:
  if(abs(chr@.Data[[1]]@rtime[which.max(chr@.Data[[1]]@intensity)] - rtvalue) 
     < 30){
    rtvalue <- chr@.Data[[1]]@rtime[which.max(chr@.Data[[1]]@intensity)]
  }
  sps <- raw_data %>%
    filterRt(rt = c(rtvalue - 0.5,
                    rtvalue + 0.5)) %>%
    spectra
  sps.df <- as.data.frame(sps[[2]])
  
  #plot(sps.df$mz, sps.df$i, type="h", xlab="mz", ylab="intensity", 
  #     main=paste0(gsub(".mzML", "", substring(myfiles[j], 26)), " - RT: ", round(sps[[2]]@rt,3)))
  #text(sps.df$mz[order(sps.df$i, decreasing = TRUE)[1:5]], 
  #   sps.df$i[order(sps.df$i, decreasing = TRUE)[1:5]], 
  #   round(sps.df$mz[order(sps.df$i, decreasing = TRUE)[1:5]], 4), 
  #   cex=0.7)
  #points(sps.df$mz[unlist(which_within(mzvalue, sps.df$mz))], 
  #       sps.df$i[unlist(which_within(mzvalue, sps.df$mz))],
  #       pch = 8, col = "red")
  
  if(study == "internal_standards"){
    mytitle <- paste0(gsub(".mzML", "", substring(myfiles[j], 26)), " - RT: ", 
                      round(sps[[2]]@rt,3))
  }else if(study == "standards_dilution"){
    mytitle <- paste0(gsub(".mzML", "", substring(myfiles[j], 14)), " - RT: ", 
                      round(sps[[2]]@rt,3))
  }
  
  sps.df.clean <- sps.df[-unlist(which_within(exclusion_list, sps.df$mz, mzd = 0.02)),]
  plot(sps.df.clean$mz, sps.df.clean$i, type="h", xlab="mz", ylab="intensity", 
       main=paste0(mytitle, " - RT: ", round(sps[[2]]@rt,3)))
  text(sps.df.clean$mz[order(sps.df.clean$i, decreasing = TRUE)[1:5]], 
       sps.df.clean$i[order(sps.df.clean$i, decreasing = TRUE)[1:5]], 
       round(sps.df.clean$mz[order(sps.df.clean$i, decreasing = TRUE)[1:5]], 4), 
       cex=0.7)
  points(sps.df.clean$mz[unlist(which_within(mzvalue, sps.df.clean$mz))], 
         sps.df.clean$i[unlist(which_within(mzvalue, sps.df.clean$mz))],
         pch = 8, col = "red")
  if(polarity=="POS" & mycompound$POS!="[M+H]+"){
    mzmol <- unlist(mass2mz(mycompound$mzneut, adduct = "[M+H]+"))
    points(sps.df.clean$mz[unlist(which_within(mzmol, sps.df.clean$mz))], 
           sps.df.clean$i[unlist(which_within(mzmol, sps.df.clean$mz))],
           pch = 8, col = "blue")
  }else if(polarity=="NEG" & mycompound$NEG!="[M-H]-"){
    mzmol <- unlist(mass2mz(mycompound$mzneut, adduct = "[M-H]-"))
    points(sps.df.clean$mz[unlist(which_within(mzmol, sps.df.clean$mz))], 
           sps.df.clean$i[unlist(which_within(mzmol, sps.df.clean$mz))],
           pch = 8, col = "blue")
  }
  points(sps.df.clean$mz[unlist(which_within(0, sps.df.clean$mz))], 
         sps.df.clean$i[unlist(which_within(0, sps.df.clean$mz))],
         pch = 8, col = "green")
}
