#' Age adjusted death rates for six causes, United States 1900 to 1998
#'
#' Age adjusted death rates per 100,000 for accidents, cancer, heart
#' disease, influenza and pneumonia, stroke and tuberculosis in the Death
#' Registration States (1900 to 1932) and the United States (1933 to 1998),
#' with the resident population.
#'
#' @format A data frame with 594 rows and 4 columns:
#' \describe{
#'   \item{year}{Calendar year.}
#'   \item{cause}{Cause of death (factor).}
#'   \item{rate}{Age adjusted death rate per 100,000.}
#'   \item{population}{United States resident population.}
#' }
#' @source CDC National Center for Health Statistics, National Vital
#'   Statistics System, Age Adjusted Death Rates for Selected Causes, Death
#'   Registration States 1900 to 1932 and United States 1933 to 1998
#'   (HIST293), as compiled by van Doren T (2024), Estimating time points of
#'   significant change in cause-specific mortality: joinpoint regression in
#'   R, The Denominator, \doi{10.6069/JK1Y-7Z16}.
#' @examples
#' head(us_mortality)
#' table(us_mortality$cause)
"us_mortality"

#' Standard populations for age adjustment
#'
#' Standard population distributions by age group for direct age
#' adjustment with [age_adjust()] and [jp_cluster()].
#'
#' @format A data frame with columns `standard`, `age_group` and
#'   `population`. Standards: `"us2000"` (2000 United States standard
#'   population, Census P25-1130, 19 age groups), `"segi1960"` (Segi world
#'   standard, 18 groups, per 100,000) and `"esp2013"` (2013 European
#'   standard population, 19 groups, per 100,000).
#' @source United States Census Bureau, Current Population Reports
#'   P25-1130 (2000 standard); Segi M (1960), Cancer mortality for selected
#'   sites in 24 countries; Eurostat (2013), Revision of the European
#'   Standard Population.
#' @examples
#' head(std_populations)
#' tapply(std_populations$population, std_populations$standard, sum)
"std_populations"

#' Prepregnancy body mass index in Ontario, 2012 to 2021
#'
#' Prepregnancy body mass index (BMI) of 1,170,997 singleton births in
#' Ontario, Canada, by fiscal year (April to March) from 2012-2013 to
#' 2021-2022: the number of births in each World Health Organization BMI
#' category, and the yearly mean and standard deviation of BMI.
#'
#' @format A data frame with 60 rows and 6 columns:
#' \describe{
#'   \item{year}{Start of the fiscal year (2012 is April 2012 to March
#'     2013).}
#'   \item{category}{BMI category (ordered from underweight to obesity
#'     class III).}
#'   \item{count}{Births in the category.}
#'   \item{births}{All births in the fiscal year.}
#'   \item{mean_bmi,sd_bmi}{Yearly mean and standard deviation of BMI in
#'     kg/m^2, rounded to one decimal and repeated over categories.}
#' }
#' @source Khademioore S, Darling EK, D'Souza R, Muraca GM, Guo Y,
#'   Anderson LN (2026). Prepregnancy body mass index in Ontario from 2012
#'   to 2022: a population-based descriptive study. Canadian Journal of
#'   Public Health. \doi{10.17269/s41997-026-01214-w}. Table 3, licensed
#'   under CC BY 4.0.
#' @examples
#' head(ontario_bmi)
"ontario_bmi"
