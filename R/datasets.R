

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


#' Regional unemployment rates Italian provinces
#'
#' A panel dataset containing unemployment rates and other economic
#' variables for Italian NUTS-3 provinces during the years 2010-2014.
#'
#' @format A data frame with 515 rows and 17 variables:
#' \describe{
#'   \item{prov}{province (NUTS-3) coded as an integer}
#'   \item{name}{province (NUTS-3) coded as a name}
#'   \item{reg}{region (NUTS-2) coded as a name}
#'   \item{area}{area of the province (km~2~)}
#'   \item{year}{year}
#'   \item{unrate}{unemployment rate (percentage)}
#'   \item{agri}{share of employment in agriculture (percentage)}
#'   \item{ind}{share of employment in industry (percentage)}
#'   \item{cons}{share of employment in construction (percentage)}
#'   \item{serv}{share of employment in services (percentage)}
#'   \item{popdens}{population density}
#'   \item{ln_popdens}{population density in logarithms}
#'   \item{empgrowth}{employment growth rate (percentage)}
#'   \item{partrate}{labor force participation rate, i.e. the
#'                   ratio between the total labor force and the
#'                   working population}
#'   \item{long}{longitude of the centroid of the province}
#'   \item{lat}{latitude of the centroid of the province}
#'   \item{South}{dummy variable with unit value for southern provinces}
#' }
#'
#' @source Italian National Institute of Statistics (ISTAT)
#'         \url{https://www.istat.it/}
"unemp_it"


#' Spatial weight matrix for Italian provinces
#'
#' A spatial weight matrix row-standardized for Italian NUTS-3 provinces
#'
#' @format A row-standardized squared matrix with 103 rows and columns.
#' The rows and columns follow the same order than provinces included in
#' \code{\link{unemp_it}} data frame.
#'
#' @source Italian National Institute of Statistics (ISTAT)
#'         \url{https://www.istat.it/}
"W_italy"


#' Homicides in U.S. counties
#'
#' Homicides and selected socio-economic characteristics for continental
#' U.S. counties. Data for four decennial census years: 1960, 1970, 1980 and 1990.
#'
#' @format A data frame with 3085 rows and 69 variables:
#' \describe{
#'   \item{NAME}{County coded as a name (factor)}
#'   \item{STATE_NAME}{state fips code (factor)}
#'   \item{STATE_FIPS}{state fips code (factor)}
#'   \item{CNTY_FIPS}{county fips code (character)}
#'   \item{FIPS}{combined state and county fips code (character)}
#'   \item{STFIPS}{State fips code (numeric)}
#'   \item{COFIPS}{county fips code (numeric)}
#'   \item{FIPSNO}{fips code as numeric variable}
#'   \item{SOUTH}{dummy variable for Southern counties (South = 1)}
#'   \item{HR60}{homicide rate per 100,000 (1960, 1970, 1980, 1990)}
#'   \item{HR70}{homicide rate per 100,000 (1960, 1970, 1980, 1990)}
#'   \item{HR80}{homicide rate per 100,000 (1960, 1970, 1980, 1990)}
#'   \item{HR90}{homicide rate per 100,000 (1960, 1970, 1980, 1990)}
#'   \item{HC60}{homicide count, three year average centered on 1960,
#'     1970, 1980, 1990}
#'   \item{HC70}{homicide count, three year average centered on 1960,
#'     1970, 1980, 1990}
#'   \item{HC80}{homicide count, three year average centered on 1960,
#'     1970, 1980, 1990}
#'   \item{HC90}{homicide count, three year average centered on 1960,
#'     1970, 1980, 1990}
#'   \item{PO60}{county population, 1960, 1970, 1980, 1990}
#'   \item{PO70}{county population, 1960, 1970, 1980, 1990}
#'   \item{PO80}{county population, 1960, 1970, 1980, 1990}
#'   \item{PO90}{county population, 1960, 1970, 1980, 1990}
#'   \item{RD60}{resource deprivation 1960, 1970, 1980, 1990 (principal
#'     component, see Codebook for details)}
#'   \item{RD70}{resource deprivation 1960, 1970, 1980, 1990 (principal
#'     component, see Codebook for details)}
#'   \item{RD80}{resource deprivation 1960, 1970, 1980, 1990 (principal
#'     component, see Codebook for details)}
#'   \item{RD90}{resource deprivation 1960, 1970, 1980, 1990 (principal
#'     component, see Codebook for details)}
#'   \item{PS60}{population structure 1960, 1970, 1980, 1990 (principal
#'     component, see Codebook for details)}
#'   \item{PS70}{population structure 1960, 1970, 1980, 1990 (principal
#'     component, see Codebook for details)}
#'   \item{PS80}{population structure 1960, 1970, 1980, 1990 (principal
#'     component, see Codebook for details)}
#'   \item{PS90}{population structure 1960, 1970, 1980, 1990 (principal
#'     component, see Codebook for details)}
#'   \item{UE60}{unemployment rate 1960, 1970, 1980, 1990}
#'   \item{UE70}{unemployment rate 1960, 1970, 1980, 1990}
#'   \item{UE80}{unemployment rate 1960, 1970, 1980, 1990}
#'   \item{UE90}{unemployment rate 1960, 1970, 1980, 1990}
#'   \item{DV60}{divorce rate 1960, 1970, 1980, 1990 (\% males
#'     over 14 divorced)}
#'   \item{DV70}{divorce rate 1960, 1970, 1980, 1990 (\% males
#'     over 14 divorced)}
#'   \item{DV80}{divorce rate 1960, 1970, 1980, 1990 (\% males
#'     over 14 divorced)}
#'   \item{DV90}{divorce rate 1960, 1970, 1980, 1990 (\% males
#'     over 14 divorced)}
#'   \item{MA60}{median age 1960, 1970, 1980, 1990}
#'   \item{MA70}{median age 1960, 1970, 1980, 1990}
#'   \item{MA80}{median age 1960, 1970, 1980, 1990}
#'   \item{MA90}{median age 1960, 1970, 1980, 1990}
#'   \item{POL60}{log of population 1960, 1970, 1980, 1990}
#'   \item{POL70}{log of population 1960, 1970, 1980, 1990}
#'   \item{POL80}{log of population 1960, 1970, 1980, 1990}
#'   \item{POL90}{log of population 1960, 1970, 1980, 1990}
#'   \item{DNL60}{log of population density 1960, 1970, 1980, 1990}
#'   \item{DNL70}{log of population density 1960, 1970, 1980, 1990}
#'   \item{DNL80}{log of population density 1960, 1970, 1980, 1990}
#'   \item{DNL90}{log of population density 1960, 1970, 1980, 1990}
#'   \item{MFIL59}{log of median family income 1960, 1970, 1980, 1990}
#'   \item{MFIL69}{log of median family income 1960, 1970, 1980, 1990}
#'   \item{MFIL79}{log of median family income 1960, 1970, 1980, 1990}
#'   \item{MFIL89}{log of median family income 1960, 1970, 1980, 1990}
#'   \item{FP59}{\% families below poverty 1960, 1970, 1980, 1990 (see
#'     Codebook for details)}
#'   \item{FP69}{\% families below poverty 1960, 1970, 1980, 1990 (see
#'     Codebook for details)}
#'   \item{FP79}{\% families below poverty 1960, 1970, 1980, 1990 (see
#'     Codebook for details)}
#'   \item{FP89}{\% families below poverty 1960, 1970, 1980, 1990 (see
#'     Codebook for details)}
#'   \item{BLK60}{\% black 1960, 1970, 1980, 1990}
#'   \item{BLK70}{\% black 1960, 1970, 1980, 1990}
#'   \item{BLK80}{\% black 1960, 1970, 1980, 1990}
#'   \item{BLK90}{\% black 1960, 1970, 1980, 1990}
#'   \item{GI59}{Gini index of family income inequality 1960, 1970, 1980, 1990}
#'   \item{GI69}{Gini index of family income inequality 1960, 1970, 1980, 1990}
#'   \item{GI79}{Gini index of family income inequality 1960, 1970, 1980, 1990}
#'   \item{GI89}{Gini index of family income inequality 1960, 1970, 1980, 1990}
#'   \item{FH60}{\% female headed households 1960, 1970, 1980, 1990}
#'   \item{FH70}{\% female headed households 1960, 1970, 1980, 1990}
#'   \item{FH80}{\% female headed households 1960, 1970, 1980, 1990}
#'   \item{FH90}{\% female headed households 1960, 1970, 1980, 1990}
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
"NAT"


#' Homicides in U.S. counties in panel format
#'
#' Homicides and selected socio-economic characteristics for continental
#' U.S. counties in panel format for years: 1960, 1970, 1980 and 1990.
#'
#' @format A data frame with 12340 rows and 25 variables:
#' \describe{
#'   \item{NAME}{County coded as a name (factor)}
#'   \item{STATE_NAME}{state fips code (factor)}
#'   \item{STATE_FIPS}{state fips code (factor)}
#'   \item{CNTY_FIPS}{county fips code (character)}
#'   \item{FIPS}{combined state and county fips code (character)}
#'   \item{STFIPS}{State fips code (numeric)}
#'   \item{COFIPS}{county fips code (numeric)}
#'   \item{FIPSNO}{fips code as numeric variable}
#'   \item{SOUTH}{dummy variable for Southern counties (South = 1)}
#'   \item{time}{year of the census observation: 1960, 1970, 1980, 1990}
#'   \item{HR}{homicide rate per 100,000}
#'   \item{HC}{homicide count, three year average centered}
#'   \item{PO}{county population}
#'   \item{RD}{resource deprivation (principal component, see Codebook
#'     for details)}
#'   \item{PS}{population structure (principal component, see Codebook
#'     for details)}
#'   \item{UE}{unemployment rate }
#'   \item{DV}{divorce rate (\% males over 14 divorced)}
#'   \item{MA}{median age}
#'   \item{POL}{log of population}
#'   \item{DNL}{log of population density}
#'   \item{MFIL}{log of median family income}
#'   \item{FP}{\% families below poverty (see Codebook for details)}
#'   \item{BLK}{\% black}
#'   \item{GI}{Gini index of family income inequality}
#'   \item{FH}{\% female headed households}
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
"NATpanel"


#' Spatial weight matrix for U.S. Counties
#'
#' A spatial weight matrix row-standardized based on first order
#'   contiguity criterium.
#'
#' @format A row-standardized squared matrix with 6085 rows and columns.
#' The rows and columns follow the same order than Counties included in
#' \code{\link{NAT}} data frame.
#'
#' @source \url{https://geodacenter.github.io/data-and-lab/ncovr/}
#'
#' @references
#'   \itemize{
#'     \item Baller, R., L. Anselin, S. Messner, G. Deane and D. Hawkins
#'       (2001). Structural covariates of US county homicide rates:
#'       incorporating spatial effects. \emph{Criminology} 39, 561-590.
#' }
"W"


#' Simulated data from a spsur-SAR DGP
#'
#' Dataset simulated from a spsur-SAR data generating process with:
#' Number of periods: nT = 4
#' Number of equations: nG = 3
#' Number of spatial units: nR = 49
#' Number of independent variables for equation = 5
#' Spatial parameters: \eqn{\lambda_1 = \lambda_2 = \lambda_3 = 0.5}
#' Variances: \eqn{\sigma_1^2 = \sigma_2^2 = \sigma_3^2 = 1}
#' Covariances: \eqn{\sigma_{12} = \sigma_{13} = \sigma_{23} = 0.5}
#'
#' @format
#' \describe{
#'   \item{Ysar}{Vector of dependent variables (588x1)}
#'   \item{XXsar}{Matrix of independent variables (588x588)}
#' }
#'
#' @seealso
#' \itemize{
#'   \code{\link{dgp_spsur}}: Function to simulate DGP from spatial SUR
#'     processes.
#' }
"XXsar"


#' Simulated data from a spsur-SAR DGP
#'
#' Dataset simulated from a spsur-SAR data generating process with:
#' Number of periods: nT = 4
#' Number of equations: nG = 3
#' Number of spatial units: nR = 49
#' Number of independent variables for equation = 5
#' Spatial parameters: \eqn{\lambda_1 = \lambda_2 = \lambda_3 = 0.5}
#' Variances: \eqn{\sigma_1^2 = \sigma_2^2 = \sigma_3^2 = 1}
#' Covariances: \eqn{\sigma_{12} = \sigma_{13} = \sigma_{23} = 0.5}
#'
#' @format
#' \describe{
#'   \item{Ysar}{Vector of dependent variables (588x1)}
#'   \item{XXsar}{Matrix of independent variables (588x588)}
#' }
#'
#' @seealso
#' \itemize{
#'   \code{\link{dgp_spsur}}: Function to simulate DGP from spatial SUR
#'     processes.
#' }
"Ysar"


#' Spatial weight matrix for simulated spsur-SAR DGP
#'
#' A spatial weight matrix row-standardized based on first order
#'   contiguity criterium.
#'
#' @format A row-standardized squared matrix with 49 rows and columns.
#'
"Ws"



