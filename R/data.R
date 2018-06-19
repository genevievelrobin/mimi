#' Excerpt of the 2016 Public Use American Census Survey (Alabama only)
#'
#' A dataset containing answers of 24614 Alabama households to 20 questions
#'
#' @format survey A data frame with 24614 rows and 20 columns:
#' \describe{
#'   \item{NP}{Number of persons in household}
#'   \item{ACCESS}{Access to the internet. 1 yes 0 no.}
#'   \item{AGS}{Sales of agriculture products ($, yearly)}
#'   \item{BATH}{Bathtub or shower. 0 yes 1 no.}
#'   \item{BDSP}{Number of bedrooms in household.}
#'   \item{BROADBND}{Cellular data plan for a smartphone or other mobile device}{1 yes 2 no}
#'   \item{COMPOTHX}{Other computer equipment. 1 yes 2 no}
#'   \item{CONP}{Condo fee ($, monthly)}
#'   \item{ELEP}{Electricity ($, monthly)}
#'   \item{FS}{Food Stamps. 0 no 1 yes}
#'   \item{FULP}{Fuel cost ($, yearly)}
#'   \item{GASP}{Gas ($, monthly)}
#'   \item{MHP}{Mobile home costs}{$, yearly}
#'   \item{REFR}{Refrigerator, 1 yes, 2 no.}
#'   \item{RMSP}{Number of rooms in household}
#'   \item{RWAT}{Hot and cold running water. 1 yes 2 no}
#'   \item{SATELLITE}{Satellite internet service. 1 yes 2 no.}
#'   \item{WATP}{Water ($, yearly)}
#'   \item{FFINCP}{Family income allocation flag (past 12 months) 0 No 1 yes.}
#' }
#' @format fes A factor indicating the Family and Employment Status with the following levels:
#' \describe{
#' \item{1}{Married-couple family: Husband and wife in Labor Force}
#' \item{2}{Married-couple family: Husband in labor force, wife.not in LF}
#' \item{3}{Married-couple family: Husband not in LF,wife in LF}
#' \item{4}{Married-couple family: Neither husband nor wife in LF}
#' \item{5}{Other family: Male householder, no wife present, in LF}
#' \item{6}{Other family: Male householder, no wife present, not in LF}
#' \item{7}{Other family: Female householder, no husband present, in LF}
#' \item{8}{Other family: Female householder, no husband present, not in LF}
#' \item{9}{Other}
#' }
#' @source \url{https://factfinder.census.gov/faces/nav/jsf/pages/searchresults.xhtml?refresh=t}
"acs2016"
