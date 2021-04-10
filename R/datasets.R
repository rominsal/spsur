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
#'   \item{geometry}{geometry of sf object.}
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

#' Within/Exit mobility index and incidence COVID-19 at Spain provinces
#'
#' Weekly within/exit mobility indices and COVID-19 incidence in Spain at provincial level. A total of 17 weeks from February 21 to 
#' May 21. Every week starts on a Friday and ends on the following Thursday. All travels 
#' are expressed with respect to the pre-COVID week of February 14th-20th, 2020 (week0). A value lower than 1
#' indicate a reduction of the mobility. Upper values than 1 indicate an increase in the mobility with respect
#' to the reference week (week0). 
#'
#' @format A data frame with 850 rows and 11 variables:
#' \describe{
#'   \item{province}{Province name coded as a factor}
#'   \item{indiv}{Province coded as a number}
#'   \item{time}{Number of week afther pre-COVID (week0: February 14th-20th, 2020). 
#'   time=1 (week1, February 21th-27th); time=2 (week2,February 28th-March 5th);....}
#'   \item{Within}{Mobility index within of the province. See details for a formal definition} 
#'   \item{Exits}{Mobility index of exits of the province. See details for a formal definition}
#'   \item{Emergence}{Dummy variable. 1 if the Emergence State ("Estado de Alarma") is active in the week "time" but economic
#'    activity of essential services are allowed. 0 in anoter case.
#'    The Emergence State was active in Spain from February 14th, 2020 to May 21th, 2020}
#'   \item{EmergenceTotal}{Dummy variable. 1 if the Emergence State ("Estado de Alarma") is active in the week "time" and economic 
#'   activity of essential services are not allowed. 0 in anoter case}
#'   \item{Old65}{Percentage of population aged 65 and older in the province}
#'   \item{Density}{Inhabitants (in thousands) per km^2 in the province}
#'   \item{Essential}{Percentage of firms in the province "indiv" with essential activities (food, health) over the total of firms in the province "indiv" 
#'   in the province (e.g., food, heath and some economic subsectors of industry and construction)}
#'   \item{Incidence}{Weekly incidence in the week "time-1" in logs}
#' }
#' @section Details:
#' \itemize{
#' \item Mobility
#'   \cr
#'   \cr
#' The mobility indices Within and Exits has been obtain as ratio of the total number of weekly travels
#' in reference to the total number of travels in the week of reference (week0). \cr
#' In particular, 
#' \cr 
#' \cr 
#'  \eqn{
#'  Within = \frac{\textnormal{Number of travels within the province 
#'  'indiv' in the week 'time'}} 
#'   {\textnormal{Number of travels within the province 'indiv' in the 
#'   reference week (week0)}} }
#'   \cr 
#'   \cr 
#'   \eqn{
#'   Exits = \frac{ 
#'     \begin{tabular}{c}
#'       \textnormal{Number of travels with origin in the 
#'                   province 'indiv' and arrival to another} \\
#'       \textnormal{province in the week 'time'} 
#'     \end{tabular} } 
#'   { \begin{tabular}{c}
#'       \textnormal{Number of travels with origin in the province 
#'                   'indiv' and arrival to another} \\
#'        \textnormal{province in the reference week (week0)} 
#'     \end{tabular} } 
#'    }
#'   \cr 
#'   \cr
#'   A ‘travel’ is a displacement from an origin to a destination of at least 500m. A travel can have several stages.
#'   These stages of the same travel are calculated based on the duration of the intermediate stop. For example
#'   , if I move from origin A to destination B with a stop at a point C of long duration, it
#'   is considered two travels, but if the stop is short, it is considered a single travel. For example,
#'   I can go by train from Madrid to Alicante, and there takes a bus to Benidorm. The travel will be one
#'   (Madrid-Benidorm). If I do the same but the stop in Alicante is long (or an overnight stay, for example),
#'   two travels will be considered (Madrid-Alicante, and Alicante-Benidorm). Similarly, if I go by car from
#'   Madrid to Alicante and stop for 15 minutes to take a coffee, it is also considered only one travel and not two.
#'   The travels considered in this study are always from 500m, due to the limitation of source data that is based
#'   on mobile telephony and its antennas. But one travel can be 600 meters or 600km.
#'   \cr
#'   \cr
#'   \item Incidence
#'   \cr
#'   \cr
#'   \eqn{
#'   Incidence = \log \left( 
#'     \frac{ 
#'       \begin{tabular}{c}
#'          \textnormal{total diagnostic cases of 
#'                      COVID-19, PCR test in the week 'time'} \\
#'          \textnormal{at the province 'indiv'} 
#'       \end{tabular} }  
#'   {\textnormal{total population in the province 'indiv'}}
#'   \right) }
#'   \cr
#'   \cr
#'   \item Essential activities
#'   \cr
#'   \cr
#'   Essential activities. Economic activities whose activities are essential for the population. By example, essential activities are
#'   healthcare, food supply, State security, media and communication, refuse collection, management and public transport, etc.
#'   Non essential activities are by example, restaurants, hotels, hairdressers, etc. A full list in the Spanish official bulletin
#'   \cr
#'   \cr
#'   \url{https://www.boe.es/boe/dias/2020/03/29/pdfs/BOE-A-2020-4166.pdf}
#' }
#' @source Ministerio de Transportes, Movilidad y Agenda Urbana. Spain Government.
#'   \cr
#'   \url{https://www.mitma.gob.es/ministerio/covid-19/evolucion-movilidad-big-data}
#'   \cr
#'   The National Statistics Institute 
#'   \cr
#'   \url{https://www.ine.es/en/index.htm}
#'   \cr
#'   Instituto de Salud Carlos III
#'   \cr
#'   \url{https://cnecovid.isciii.es/covid19/}
#'
#'
"spain.covid"

#' Spain geometry
#'
#' A sf object with the Spanish geometry
#'
#' @format A sf object with the Spanish geometry
#' 
#' \describe{
#'   \item{PROVINCIA}{province name coded as a factor}
#'   \item{ID_INE}{province coded as a number}
#'   }
"spain.covid.sf"