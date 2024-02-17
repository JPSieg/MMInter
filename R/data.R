#' @title df.NISTSRD46
#' @description A tidy formatted version of the NIST-SRD-46 critically selected stability constants database
#' @format A data frame with 89824 rows and 10 variables:
#' \describe{
#'   \item{\code{M}}{character Temperature in degCelcius}
#'   \item{\code{L}}{character The absorbance measured using a spectrometer}
#'   \item{\code{Type}}{character The experiment name}
#'   \item{\code{constant}}{character The RNA sequence}
#'   \item{\code{temperature}}{character The buffer}
#'   \item{\code{ionicstrength}}{character The file path}
#'   \item{\code{ligandenNr}}{interger The cell in a run}
#'   \item{\code{metalNr}}{interger The blank the spectrometer uses}
#'   \item{\code{Lit}}{character The pathlength of the cuvette in cm}
#'}
#' @source \url{https://doi.org/10.18434/M32154}
"df.NISTSRD46"

#' @title df.conc
#' @description Total biological concentration of 15 metabolites, Mg, Mn, Zn, and Ca in E. coli
#' @format A data frame with 19 rows and 2 variables:
#' \describe{
#'   \item{\code{Metabolites}}{character The chemical compound in the E. coli}
#'   \item{\code{Concentration}}{numeric The total molar concetration of each chemical compound in E. coli}
#'}
"df.conc"
