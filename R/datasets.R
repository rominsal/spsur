#' A classical Spatial Phillips-Curve
#'
#' A data set from Anselin (1988, p. 203-2011) used to estimate a Spatial
#' Phillips-Curve for 25 counties in South-West Ohio for two time
#' periods (1981 and 1983).
#'
#' @format A data frame with 25 rows and 10 variables:
#' \describe{
#'   \item{COUNTY}{County coded as a name.}
#'   \item{WAGE83}{Changes in wage rates for 1983.}
#'   \item{UN83}{Inverse unemployment rate in 1983.}
#'   \item{NMR83}{Net migration rate 1983.}
#'   \item{SMSA}{Dummy variable to identify counties defined as
#'     Standard Metropolitan Statistical Areas (SMSA = 1).}
#'   \item{WAGE82}{Changes in wage rates for 1982.}
#'   \item{WAGE81}{Changes in wage rates for 1981.}
#'   \item{UN80}{Inverse unemployment rate in 1980.}
#'   \item{NMR80}{Net migration rate 1983.}
#'   \item{WAGE80}{changes in wage rates.}
#' }
#'
#' @source Anselin (1988, p. 203-211)
#'
#' @references
#'   \itemize{
#'     \item Anselin, L. (1988). \emph{Spatial Econometrics:
#'       Methods and Models}. Springer Science & Business Media.
#'   }
"spc"


#' Spatial weight matrix for South-West Ohio Counties to estimate
#' Spatial Phillips-Curve
#'
#' A spatial weight matrix row-standardized based on first order
#' contiguity criterium.
#'
#' @format A row-standardized squared matrix with 25 rows and columns.
#' The rows and columns follow the same order than provinces included in
#' \emph{spc} data frame.
#'
#' @source Anselin (1988, p. 207)
#'
#' @references
#'   \itemize{
#'     \item Anselin, L. (1988). \emph{Spatial Econometrics:
#'       Methods and Models}. Springer Science & Business Media.
#'   }
"Wspc"

#' Homicides in U.S. counties
#'
#' Homicides and selected socio-economic characteristics for continental
#' U.S. counties. Data for four decennial census years: 1960, 1970, 1980
#' and 1990.
#'
#' @format An spatial feature (sf) object with 3085 rows and 41 variables:
#' \describe{
#'   \item{NAME}{County coded as a name (factor)}
#'   \item{STATE_NAME}{state fips code (factor)}
#'   \item{FIPS}{state fips code (factor)}
#'   \item{SOUTH}{dummy variable for Southern counties (South = 1)}
#'   \item{HR60, HR70, HR80, HR90}{homicide rate per 100,000 (1960, 1970,
#'                                1980, 1990)}
#'   \item{HC60, HC70, HC80, HC90}{homicide count, three year average centered
#'                                 on 1960, 1970, 1980, 1990}
#'   \item{PO60, PO70, PO80, PO90}{county population, 1960, 1970, 1980, 1990}
#'   \item{RD60, RD70, RD80, RD90}{resource deprivation 1960, 1970, 1980,
#'                                 1990 (principal component, see Codebook for details)}
#'   \item{PS60, PS70, PS80, PS90}{population structure 1960, 1970, 1980,
#'                                1990 (principal component, see Codebook for details)}
#'   \item{UE60, UE70, UE80, UE90}{unemployment rate 1960, 1970, 1980, 1990}
#'   \item{DV60, DV70, DV80, DV90}{divorce rate 1960, 1970, 1980, 1990
#'                                 (\% males over 14 divorced)}
#'   \item{MA60, MA70, MA80, MA90}{median age 1960, 1970, 1980, 1990}
#'   \item{FP59, FP69, FP79, FP89}{\% families below poverty 1960, 1970, 1980,
#'                                 1990 (see Codebook for details)}
#'   \item{geometry}{Multipolygon geometry of the spatial feature object}   
#' }
#'
#' @source S. Messner, L. Anselin, D. Hawkins, G. Deane, S. Tolnay, R. Baller
#'   (2000). An Atlas of the Spatial Patterning of County-Level
#'   Homicide, 1960-1990. Pittsburgh, PA, National Consortium on
#'   Violence Research (NCOVR)
#'   \url{https://geodacenter.github.io/data-and-lab/ncovr/}
#'
#' @references
#'   \itemize{
#'     \item Baller, R., L. Anselin, S. Messner, G. Deane and D. Hawkins
#'       (2001). Structural covariates of US county homicide rates:
#'       incorporating spatial effects. \emph{Criminology} 39, 561-590.
#'   }
#'
"NCOVR.sf"
