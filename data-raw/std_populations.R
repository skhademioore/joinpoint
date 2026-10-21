# Standard populations for direct age adjustment.
us2000 <- c(3794901, 15191619, 19919840, 20056779, 19819518, 18257225,
            17722067, 19511370, 22179956, 22479229, 19805793, 17224359,
            13307234, 10654272, 9409940, 8725574, 7414559, 4900234, 4259173)
us_lab <- c("00", "01-04", "05-09", "10-14", "15-19", "20-24", "25-29",
            "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64",
            "65-69", "70-74", "75-79", "80-84", "85+")
five <- function(lo) sprintf("%02d-%02d", lo, lo + 4)
segi <- c(12000, 10000, 9000, 9000, 8000, 8000, 6000, 6000, 6000, 6000,
          5000, 4000, 4000, 3000, 2000, 1000, 500, 500)
segi_lab <- c(five(seq(0, 80, 5)), "85+")
esp <- c(5000, 5500, 5500, 5500, 6000, 6000, 6500, 7000, 7000, 7000, 7000,
         6500, 6000, 5500, 5000, 4000, 2500, 1500, 1000)
esp_lab <- c(five(seq(0, 85, 5)), "90+")
mk <- function(std, lab, pop) {
  data.frame(standard = std, age_group = lab, population = pop)
}
std_populations <- rbind(
  mk("us2000", us_lab, us2000),
  mk("segi1960", segi_lab, segi),
  mk("esp2013", esp_lab, esp)
)
stopifnot(sum(segi) == 1e5, sum(esp) == 1e5, sum(us2000) == 274633642)
usethis::use_data(std_populations, overwrite = TRUE)
