#' Ocean temperatures from Argo profiling floats
#'
#' A dataset containing ocean temperature measurements from
#' three pressure levels (depths), measured by profiling
#' floats from the Argo program. Data collected in Jan, Feb,
#' and March of 2016.
#' @format A data frame with 32436 rows and 6 columns
#' \describe{
#'   \item{lon}{longitude in degrees between 0 and 360}
#'   \item{lat}{latitude in degrees between -90 and 90}
#'   \item{day}{time in days}
#'   \item{temp100}{Temperature at 100 dbars (roughly 100 meters)}
#'   \item{temp150}{Temperature at 150 dbars (roughly 150 meters)}
#'   \item{temp200}{Temperature at 200 dbars (roughly 200 meters)}
#' }
#' @source Mikael Kuusela. Argo program: \url{https://argo.ucsd.edu/}
"argo2016"
