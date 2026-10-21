# Prepregnancy body mass index of singleton births in Ontario, Canada, fiscal
# years 2012-2013 to 2021-2022, from Table 3 (Appendix A) of
# Khademioore S, Darling EK, D'Souza R, Muraca GM, Guo Y, Anderson LN (2026).
# Prepregnancy body mass index in Ontario from 2012 to 2022: a
# population-based descriptive study. Canadian Journal of Public Health.
# doi:10.17269/s41997-026-01214-w. Licensed under CC BY 4.0; reshaped to
# long format.
year <- 2012:2021
births <- c(107457, 113350, 115375, 118943, 123328, 121714, 120177, 119566,
            115034, 116053)
mean_bmi <- c(25.3, 25.3, 25.3, 25.3, 25.3, 25.6, 25.6, 25.7, 25.8, 26.0)
sd_bmi <- c(6.2, 6.1, 6.1, 6.0, 6.0, 6.2, 6.2, 6.2, 6.2, 6.3)
counts <- rbind(
  Underweight = c(9262, 9619, 10078, 10011, 10574, 9853, 9386, 9023, 8256,
                  7741),
  `Normal weight` = c(53904, 56842, 57480, 59051, 60770, 58600, 57341, 56455,
                      53770, 53332),
  Overweight = c(21955, 22967, 23525, 24423, 25394, 25612, 25713, 25745,
                 25282, 25972),
  `Obesity class I` = c(14541, 15576, 15678, 16871, 17469, 17813, 17943,
                        18529, 18038, 18827),
  `Obesity class II` = c(3957, 4415, 4547, 4581, 5010, 5165, 5218, 5122, 5129,
                         5317),
  `Obesity class III` = c(3838, 3931, 4067, 4006, 4111, 4671, 4576, 4692, 4559,
                          4864)
)
stopifnot(all(colSums(counts) == births))
ontario_bmi <- data.frame(
  year = rep(year, each = nrow(counts)),
  category = factor(rep(rownames(counts), length(year)),
                    levels = rownames(counts)),
  count = as.vector(counts),
  births = rep(births, each = nrow(counts)),
  mean_bmi = rep(mean_bmi, each = nrow(counts)),
  sd_bmi = rep(sd_bmi, each = nrow(counts))
)
usethis::use_data(ontario_bmi, overwrite = TRUE)
