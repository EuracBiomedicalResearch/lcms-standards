mycompound <- "Isoleucine"
polarity <- "NEG" # specify "POS" or "NEG"

study <- "standards_dilution" # specify "internal_standards" OR 
                              # "standards_dilution"
mixnum <- 3

source("R/which_within.R")
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



# Get the information regarding the standards that are in the samples:
std_info <- read.table(paste0("data/", study, ".txt"),
                       sep = "\t", header = TRUE, as.is = TRUE)
if(study == "internal_standards"){
  colnames(std_info) <- c("name", "molecular_weight", "POS", "NEG", "RT")
} else if(study == "standards_dilution"){
  colnames(std_info) <- c("mix", "name", "molecular_weight", "POS", "NEG", "RT")
}
std_info$name <- c(substring(std_info$name, 1, 33))
if(study == "standards_dilution"){
  std_info <- subset(std_info, mix == mixnum)
}


mycompound <- std_info[grep(mycompound, std_info$name),]
mzvalue <- mycompound[, grep(polarity, colnames(mycompound))]
da = 0.01
RTd = 10

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
  sps <- raw_data %>%
    filterRt(rt = c(chr@.Data[[1]]@rtime[which.max(
      chr@.Data[[1]]@intensity)] - 0.5,
      chr@.Data[[1]]@rtime[which.max(
        chr@.Data[[1]]@intensity)] + 0.5)) %>%
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
  
  sps.df.clean <- sps.df[-unlist(which_within(exclusion_list, sps.df$mz)),]
  plot(sps.df.clean$mz, sps.df.clean$i, type="h", xlab="mz", ylab="intensity", 
       main=paste0(mytitle, " - RT: ", round(sps[[2]]@rt,3)))
  text(sps.df.clean$mz[order(sps.df.clean$i, decreasing = TRUE)[1:5]], 
       sps.df.clean$i[order(sps.df.clean$i, decreasing = TRUE)[1:5]], 
       round(sps.df.clean$mz[order(sps.df.clean$i, decreasing = TRUE)[1:5]], 4), 
       cex=0.7)
  points(sps.df.clean$mz[unlist(which_within(mzvalue, sps.df.clean$mz))], 
         sps.df.clean$i[unlist(which_within(mzvalue, sps.df.clean$mz))],
         pch = 8, col = "red")
}
