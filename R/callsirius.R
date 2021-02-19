#' Title Call SIRIUS and get the resulting ranking
#'
#' @param input path to the input .ms file
#' @param diroutput directory where the output of SIRIUS will be saved
#'
#' @return data frame with the ranking of the chemical formulas
#'
#' @examples
#' callsirius(input = "~/Documents/demo-data/ms/Kaempferol.ms",diroutput = "~/direxample")

callsirius<-function(input, diroutput)
{
  system2('/Applications/sirius.app/Contents/MacOS/sirius', args = c('--input', input, '--output', diroutput,'formula'))
  ranking=read.table(paste0(diroutput,"/formula_identifications.tsv"),sep="\t",header=T)
  return(ranking)
}