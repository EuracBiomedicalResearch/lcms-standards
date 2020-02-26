# Preliminaries -------------------------------------------------------
MZML_PATH <- "/data/massspec/mzML/" 
#MZML_PATH <- "C:/Users/mgarciaaloy/Documents/mzML/"
polarity <- c("POS", "NEG")

library(Rdisop)
library(CompoundDb)
library(RColorBrewer)
library(xcms)
cwp <- CentWaveParam(
  peakwidth = c(2, 20), 
  ppm = 50, 
  snthresh = 5,
  mzdiff = 0.001,
  prefilter = c(3, 1000),
  noise = 100,
  integrate = 2)

ncores <- Sys.getenv("SLURM_JOB_CPUS_PER_NODE", 3)
register(bpstart(MulticoreParam(ncores)))

# Injections
injections <- read.table(
  "../data/std_serum_files.txt", 
  sep = "\t", header = TRUE, as.is = TRUE)
injections <- injections[injections$mode == "FS", ]
inj_blankQC <- injections[injections$type == "Blank_QC", ]
inj_blankQC <- inj_blankQC[grep("QC_5", inj_blankQC$mzML), ]
inj_blankQC$class <- "QC"
injections <- injections[grep("Mix", injections$mzML), ]

# Import table of standards:
std_info <- read.table(
  "../data/standards_dilution.txt",
  sep = "\t", header = TRUE, as.is = TRUE)

# Coloring factors
col_class <- brewer.pal(6, "Paired")
col_class <- col_class[-4]
names(col_class) <- c("Water_Low", "Water_High", "QC", "QC_Low", "QC_High")

# Isomers by formula --------------------------------------------------
# Get the repeated formulas
fmls <- unique(std_info$formula[duplicated(std_info$formula)])

for(i in seq(length(fmls))){
  
  # standards with formula "i"
  tmp <- std_info[std_info$formula == fmls[i], ]
  
  # calculate the mz neutral
  mzneut <- getMolecule(as.character(fmls[i]))$exactmass
  
  # start loop by ionization mode:
  for(j in 1:2){
    polarity.j <- polarity[j]
    injections.j <- injections[injections$polarity == polarity.j, ]
    
    filename <- paste0("images/isomers/", fmls[i], "_", polarity.j, ".png")
    
    png(file = filename, width = 1000, height = 500*nrow(tmp))
    par(mfrow = c(nrow(tmp), 2))
    
    # start loop by compound
    for(k in seq(nrow(tmp))){
      
      mzval <- unlist(mass2mz(
        mzneut,
        adduct = as.character(tmp[i, 
                                  which(colnames(tmp) == polarity.j)])))
      
      # Read the mzML files
      injections.k <- injections.j[
        as.numeric(gsub("Mix", "", injections.j$type)) == tmp$mix[k], ]
      injections.k <- rbind(inj_blankQC, injections.k)
      injections.k$class <- gsub(paste0(injections.k$type[10], "_"), "", 
                                 injections.k$class)
      myfiles <- paste(injections.k$folder, injections.k$mzML, sep = "/") 
      data_raw <- readMSData(paste0(MZML_PATH, myfiles), 
                             pdata = new("NAnnotatedDataFrame", 
                                         injections.k),
                             mode = "onDisk")
      # EIC
      chr <- chromatogram(data_raw, 
                          mz = mzval + 0.01 * c(-1, 1),
                          aggregationFun = "max")
      
      chr2 <- findChromPeaks(chr, param = cwp)
      delta_rt <- abs(chromPeaks(chr2)[, "rt"] - tmp$RT[k])
      pks <- data.frame(chromPeaks(chr2)[delta_rt < 60, , drop=FALSE])
      if (nrow(pks)) {
        pks2 <- pks[0, ]
        for(j in levels(factor(pks$column))){
          pks3 <- pks[pks$column == j, ]
          pks2 <- rbind(pks2, pks3[which.max(pks3$maxo),])
        }
        rtmin <- min(pks2$rtmin)
        rtmax <- max(pks2$rtmax)
        ylim <- c(0, max(pks$maxo, na.rm = TRUE))
      } else {
        rtmin <- tmp$RT[k] - 10
        rtmax <- tmp$RT[k] + 10
        ylim <- c(0, max(sapply(chr, intensity), na.rm = TRUE))
      }
      plot(chr2, peakType = "none", 
           col = col_class[chr2$class],
           main = paste0(tmp$name[k], ": ", 
                                tmp[k, colnames(tmp) == polarity.j], 
                                " (mix ", tmp$mix[k], ")"))
      plot(chr2, peakType = "none", 
           col = col_class[chr2$class],
           xlim = c(rtmin, rtmax))
      abline(v = tmp$RT[k], lty = 2)
      
      dev.off()
      
    } # send compound "k"
  } # close mode "j"
} # close formula "i"

# Isomers by ion ------------------------------------------------------