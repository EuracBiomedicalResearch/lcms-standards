#' Title Write MS data to a .ms file
#'
#' @param outputms path to the .ms file that will be created
#' @param compound string specifying the name of the measured compound (or any 
#' placeholder). This field is mandatory for SIRIUS.
#' @param ms1 data frame of MS peaks (m/z and intensities)
#' @param ms2 data frame of MS/MS peaks (m/z and intensities)
#' @param formula string specifying the molecular formula of the compound. 
#' This information is helpful if you already know the correct molecular formula 
#' and just want to compute a fragmentation tree or recalibrate the spectrum
#' @param ionization the ionization mode.
#' @param charge is redundant if you already provided the ion mode. Otherwise, 
#' it gives the charge of the ion (1 or -1).
#' @param collision collision energy if available
#'
#'
#' @examples
#' 
#' ms1 <- data.frame(mz = 1:4,intensity = 1:4)
#' ms2 <- data.frame(mz = 2:4,intensity = 2:4)
#' writetoms(outputms = "example.ms", compound = "Bicuculline",
#'           formula = "C20H17NO6", ion = "[M+H]+", charge = "+1" , 
#'           ms1 = ms1, ms2 = ms2, collision = 35)
#' writetoms(outputms = "example.ms", compound = "Bicuculline", ion = "[M+H]+",
#'           ms1 = ms1)

writetoms <- function(outputms, ms1, ms2 = NULL, compound = "_", formula = NULL, 
                      ionization = NULL, charge = NULL, collision = NULL) {
  lines <- paste0(">compound ", compound)
  if(length(d <- dim(ms1)) != 2 || d[1] == 0 || d[2] != 2) 
    stop("ms1 must be 2-dimensional with 2 columns and number of rows > 0")
  parentmass <- ms1[1, 1]
  if(!is.null(formula))
    lines <- c(lines, paste0(">formula ", formula))
  lines <- c(lines, paste0(">parentmass ", parentmass))
  if(!is.null(ionization))
    lines <- c(lines, paste0(">ionization ", ionization))
  if(!is.null(ms2)) {
    lines <- c(lines, ifelse(!is.null(collision), 
                             paste0(">collision ", collision), ">ms2"))
    for (row in 1:dim(ms2)[1])
      lines <- c(lines, paste0(ms2[row, 1], " ", ms2[row, 2]))
  }
  lines <- c(lines, ">ms1")
  for (row in 1:dim(ms1)[1])
    lines <- c(lines, paste0(ms1[row, 1], " ", ms1[row, 2]))
  # write to file
  fileConn <- file(outputms)
  writeLines(lines, outputms)
  close(fileConn)
}
