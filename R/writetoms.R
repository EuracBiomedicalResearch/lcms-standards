

#' Title Write MS data to a .ms file
#'
#' @param outputms path to the .ms file that will be created
#' @param compound string specifying the name of the measured compound (or any placeholder). This field is mandatory for SIRIUS.
#' @param formula string specifying the molecular formula of the compound. This information is helpful if you already know the correct molecular formula and just want to compute a fragmentation tree or recalibrate the spectrum
#' @param parentmass  the mass of the parent peak
#' @param ionization the ionization mode.
#' @param charge is redundant if you already provided the ion mode. Otherwise, it gives the charge of the ion (1 or -1).
#' @param ms1 data frame of MS peaks (m/z and intensities)
#' @param ms2 data frame of MS/MS peaks (m/z and intensities)
#' @param collision collision energy if available
#'
#'
#' @examples
#' ms1=data.frame(mz=1:4,intensity=1:4)
#' ms2=data.frame(mz=2:4,intensity=2:4)
#' writetoms(outputms = "~/example.ms", compound = "Bicuculline", formula = "C20H17NO6", parentmass = "368.113616943359", ion = "[M+H]+", charge ="+1" , ms1 = ms1, ms2 = ms2, collision = 35)
#' writetoms(outputms = "example.ms", compound = "Bicuculline", parentmass = "368.113616943359", ion = "[M+H]+", ms1 = ms1)
#' 
writetoms<-function(outputms, compound=NULL, formula=NULL, parentmass=NULL, ionization=NULL, charge=NULL, ms1=NULL, ms2=NULL, collision=NULL)
{
  if(!is.null(compound))
    lines=paste0(">compound ",compound)
  if(!is.null(formula))
    lines=c(lines,paste0(">formula ",formula))
  if(!is.null(parentmass))
    lines=c(lines,paste0(">parentmass ",parentmass))
  if(!is.null(ionization))
    lines=c(lines,paste0(">ionization ",ionization))
  if(!is.null(ms2))
  {
    lines=c(lines,ifelse(!is.null(collision),paste0(">collision ",collision),">ms2"))
    for (row in 1:dim(ms2)[1])
      lines=c(lines,paste0(ms2[row,1]," ",ms2[row,2]))
  }
  if(!is.null(ms1))
  {
    lines=c(lines,">ms1")
    for (row in 1:dim(ms1)[1])
      lines=c(lines,paste0(ms1[row,1]," ",ms1[row,2]))
  }
  # write to file
  fileConn<-file(outputms)
  writeLines(lines, outputms)
  close(fileConn)
}

