#' Extract ion chromatograms and return the corresponding peaks
#' 
#' @param file A character with file name to be read
#' @param mz A numeric vector of mass-to-charge (mz) values for the 
#' chromatogram(s). If it mz value is named, the output will also associate 
#' each detected peak with the corresponding compound name. 
#' See examples below for details.
#' @param mzd A numeric vector with the desidered mz range (in Da) for calculate
#' the mz range (default 0.01) that will be used to extract each ion.
#' @param param List of parameters used in CentWaveParam for peak detection
#' @return The list of detected peaks corresponding to mz values
#' @examples
#' listpeaks(file = system.file("sciex/20171016_POOL_POS_1_105-134.mzML", package = "msdata"),
#' mz = c(C = 117.0856, V = 124.1006), 
#' param = CentWaveParam(peakwidth = c(2, 10), ppm = 30))
listpeaks <- function(file = character(0), mz = numeric(0), mzd = 0.01, 
                param = CentWaveParam()) {
  data <- readMSData(file, mode = "onDisk")
  mzmat <- cbind(mz - mzd, mz + mzd)
  chr <- chromatogram(data, mz = mzmat, aggregationFun = "max")
  chr <- findChromPeaks(chr, param = param)
  chr <- data.frame(chromPeaks(chr))
  if (length(names(mz))){
    chr$name <- names(mz)[chr$row]
  }
  chr
}
