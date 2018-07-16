library(fitdistrplus)
library(ggplot2)

# print(summary(catalogue_2))
# print(summary(stars_catalogue))
# print(summary(galaxies_catalogue))

fit <- lm(catalogue_2$PM ~ catalogue_2$B_IMAGE + catalogue_2$A_IMAGE)
print(summary(fit))
print(confint(fit))
