MZML_PATH <- "C:/Users/mgarciaaloy/Documents/mzML/"

study <- "standards_dilution" # specify "internal_standards" OR 
                              # "standards_dilution"
mixnum <- 17 # specify which MIX

library(xcms)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(Rdisop)
library(CompoundDb)

# the information regarding the injection sequence:
injections <- read.table(paste0("data/", study, "_files.txt"), #import the file
                         sep = "\t", header = TRUE, as.is = TRUE)
myfiles <- injections$mzML # get file names
if(study == "standards_dilution"){
  (myfiles <- myfiles[grep(paste0("MIX ", mixnum, "[A-Z]"), myfiles)])
}
myfiles_pos <- myfiles[grep("POS", myfiles)]
myfiles_neg <- myfiles[grep("NEG", myfiles)]


# Get the information regarding the standards that are in the samples:
std_info <- read.table(paste0("data/", study, ".txt"),
                       sep = "\t", header = TRUE, as.is = TRUE)
std_info$name <- c(substring(std_info$name, 1, 33))
if(study == "standards_dilution"){
  std_info <- subset(std_info, mix == mixnum)
}
std_info$mzneut = NA
for(i in seq(nrow(std_info))){
  if(grepl("C", std_info$formula[i])){std_info$mzneut[i] = 
    getMolecule(as.character(std_info$formula[i]))$exactmass}else{
      std_info$mzneut[i] = paste(std_info$formula[i])}
}
write.csv(std_info, "x.csv", row.names = FALSE)
std_info = read.csv("x.csv")


# EIC:
da = 0.01
if(study == "internal_standards"){
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
} else if(study == "standards_dilution"){
  getPalette = colorRampPalette(brewer.pal(9, "YlOrRd"))
}
mycols = getPalette(length(myfiles_pos))  
data_pos <- readMSData(paste0(MZML_PATH, myfiles_pos), mode = "onDisk")
data_neg <- readMSData(paste0(MZML_PATH, myfiles_neg), mode = "onDisk")


mycompound <- "Xanthosine"
mycompound <- std_info[grep(mycompound, std_info$name),]
mycompound

mzpos <- unlist(mass2mz(mycompound$mzneut, 
                        adduct = as.character(mycompound$POS)))
mzneg <- unlist(mass2mz(mycompound$mzneut, 
                        adduct = as.character(mycompound$NEG)))

chr_pos = chromatogram(data_pos, mz = c(mzpos - da, mzpos + da),
                       aggregationFun = "max")
chr_neg = chromatogram(data_neg, 
                       mz = c(mzneg - da, mzneg + da),
                       aggregationFun = "max")
mycompound$RT
RTd = 200
newRT = 151
par(mfrow=c(1,2))
plot(chr_pos, col=mycols, ylab="",
     xlim=c(newRT - RTd, newRT + RTd))
abline(v=newRT, lty=2)
plot(chr_neg, col=mycols, ylab="",
     xlim=c(newRT - RTd, newRT + RTd))
abline(v=newRT, lty=2)
abline(v=newRT, lty=2)
