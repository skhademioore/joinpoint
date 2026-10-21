# Age adjusted death rates for six leading causes of death, Death
# Registration States 1900 to 1932 and United States 1933 to 1998
# (CDC National Vital Statistics System, historical table HIST293),
# with the United States resident population, as compiled by van Doren
# (2024, doi:10.6069/JK1Y-7Z16); `ref` is the folder of the compiled files.
ref <- "documentation/refs/joinpoint-regression-tutorial-master/"
md <- utils::read.csv(file.path(ref, "mortality data.csv"))
pop <- utils::read.csv(file.path(ref, "usa pop.csv"))
us_mortality <- merge(md, pop, by = "year")
names(us_mortality) <- c("year", "cause", "rate", "population")
us_mortality$cause <- factor(us_mortality$cause)
us_mortality <- us_mortality[order(us_mortality$cause, us_mortality$year), ]
rownames(us_mortality) <- NULL
usethis::use_data(us_mortality, overwrite = TRUE)

# Fitted values of linear joinpoint models with 0 to 3 joinpoints published
# with the same tutorial, kept as a validation fixture.
file.copy(file.path(ref, "model fits.csv"),
          "tests/testthat/fixtures/tutorial_model_fits.csv", overwrite = TRUE)
