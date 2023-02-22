# Loading packages 
require(ggplot2)
require(MASS)
require(ggpubr)
require(usdm)
require(psych)
require(lmerTest)
require(sjPlot)
require(car)
require(ggeffects)

# Loading data
worms<- read.csv("Earthworms.csv", header= TRUE)
str(worms)


# Checking for outliers
## No outliers were found
boxplot(data= worms, Worm.abundance ~ Habitat)


# Plotting the spread of earthworms across soil parameters to infer the relationships for future comparison with the GLM results
## Earthworms and soil hydraulic conductivity (Ksat)
plot(Worm.abundance ~ Ksat, data= worms)
abline(lm(Worm.abundance ~ Ksat, data = worms))

## Earthworms and soil bulk density 
plot(Worm.abundance ~ Bulk.density, data= worms)
abline(lm(Worm.abundance ~ Bulk.density, data = worms))

## Earthworms and soil porosity
plot(Worm.abundance ~ Porosity, data= worms)
abline(lm(Worm.abundance ~ Porosity, data = worms))

# Checking the centrality and spread of the data
hist(worms$Worm.abundance)
var(worms$Worm.abundance, na.rm = TRUE)
sd(worms$Worm.abundance,na.rm = TRUE)

## Calculating the median and range since earthworm abundance is a quantitative discrete variable (count) and therefore using means and standard deviation would not be appropriate
### CC
median(worms$Worm.abundance [worms$Habitat == "CC"])
range(worms$Worm.abundance [worms$Habitat == "CC"])
### TSA
median(worms$Worm.abundance [worms$Habitat == "TSA"])
range(worms$Worm.abundance [worms$Habitat == "TSA"])
### RF
median(worms$Worm.abundance [worms$Habitat == "RF"])
range(worms$Worm.abundance [worms$Habitat == "RF"])
### SSA
median(worms$Worm.abundance [worms$Habitat == "SSA"])
range(worms$Worm.abundance [worms$Habitat == "SSA"])
### TM
median(worms$Worm.abundance [worms$Habitat == "TM"])
range(worms$Worm.abundance [worms$Habitat == "TM"])


# Transforming the categorical explanatory variable "Habitat" into a factor
worms$Habitat <- as.factor(worms$Habitat)


# Fitting the model
m1 <- lm(Worm.abundance ~ Habitat + scale(Ksat) + scale(Bulk.density) + scale(Porosity), data= worms)
summary(m1)

## Since the diagnostic plots are not perfect and suggest that the residuals violate the assumption of normal distribution, I chose to do a GLM
## Additionally, the residuals are slightly right-skewed and earthworm abundance is count data. Therefore, Poisson distribution could be appropriate for the model
par(mfrow=c(2,2))
plot(m1)
hist(resid(m1), main="", xlab="Earthworm abundance", col="grey")


# Fitting the maximized GLM with poisson distribution
## The response variable is earthworm abundance and the main explanatory variable is habitat type. All other variables were added as covariates.
## The covariates were z-standardized to account for different scales.
## No interactions between the explanatory variables were added to simplify model interpretation and reduce the risk of overfitting the data since the sample size is limited
m2 <- glm(Worm.abundance ~ Habitat + scale(Ksat) + scale(Bulk.density) + scale(Porosity), family = "poisson", data= worms)
summary(m2)

## Plotting the residuals to check homogeneity of variance
### The diagnostic plots suggest that there are unequal variances
par(mfrow=c(2,2))
plot(m2)
sum(cooks.distance(m2)>1) # no outliers were removed because they represent real observations and are biologically relevant 


# Since unequal variances can be due to collinearity between the explanatory variables, the Variance Inflation Factor was used to test for collinearity among the covariates
## The results indicate that there is no collinearity between these variables since all VIFs are below the threshold of 3
## Therefore, all variables were included in the model because they are biologically relevant
car::vif(m2)
pairs.panels(worms)

## Pseudo R-squared of m2 [1-(366.41/506.47)= 0.27]

## Checking the distribution of m2
### Dispersion Parameter= 366.41/39= 9.40
### The model is overdispersed, which suggests that the model is not capturing all of the variance in the data


# Therefore, to account for overdispersion, the family was changed to the quasipoisson distribution
m3 <- glm(Worm.abundance ~ Habitat + scale(Ksat) + scale(Porosity) + scale(Bulk.density), data = worms, family = "quasipoisson")
summary(m3)
plot(m3)

## Dispersion Parameter= 8.95. Overdispersion is still high
## Pseudo R-squared= 0.27


# To account for overdispersion, the family was changed to the negative binomial error distribution
## The log link function was used to account for a non-linear and more complex relationship between the response, explanatory variable and covariates
m4 <- glm.nb(Worm.abundance ~ Habitat + scale(Ksat) + scale(Porosity) + scale(Bulk.density), data = worms, link = log)
summary(m4)
plot(m4) # It looks fine
sum(cooks.distance(m4)>1)# less outliers than m2

## The overdispersion was corrected: Dispersion Parameter= 1.80 (as it is approximately 1, I accepted it)
## Even though the Pseudo R-squared is lower [1-(55.12/69.59)= 0.21] comparatively to "m3", "m4" is not overdispersed
## The overdispersion in "m3" was affecting the interpretation of the model due to the biased estimation of parameter estimates and significance, which would lead to invalid conclusions.
## Therefore, "m4" was selected as the final model 


## Writing the model equation to better visualize it and understand the effect size
### Main equation: Earthworm abundance = 2.48 + 0.72 * HabitatRF - 0.21 * HabitatSSA + 0.82 * HabitatTM - 0.26 * HabitatTSA + 0.07 * Ksat + 0.001 * Porosity + 0.08 * Bulk density

### Since only Habitat TM was significant, I can't assume that the coefficients for all of the other variables are significantly different from zero
### Therefore, I only focused on the variable TM:    
#### ln(Earthworm abundance) = 2.48 + 0.82 * HabitatTM

### The effect size = 15 earthworms
#### e^(2.48+0.82) - e^(2.48) 
exp(2.48+0.82) - exp(2.48)

### Predicted earthworm abundance on TM= e^(2.48) * e^(0.82)
exp(2.48) * exp(0.82)

### This matches the prediction plot
test<-ggpredict(m4, term="Habitat") 
plot(test)

ggplot(test, aes(x = x, y = predicted)) +
  geom_point(size = 2) + 
  scale_x_discrete(limits = c("CC", "TSA", "RF", "SSA", "TM")) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(worms, mapping = aes(x = Habitat, y = Worm.abundance), color = "darkgrey") +
  labs(x= "Habitat type", y = "Earthworm abundance") +
  theme_classic()


# Estimating the goodness-of-fit
## The variable "Habitat" is significant because the model fits the data well
anova(m4, test="Chisq")


# Comparing the two models
## Creating a null model for comparison
null_m <- glm.nb(Worm.abundance ~ 1, data = worms, link = log)
summary(null_m)

### The AIC of the null model (356.19) is only slightly lower than the AIC of the final model (357.68)
### As the AIC is a measure of the model's relative quality, it provides a trade-off between the its complexity and the goodness of its fit
### Since the aim was to select a biologically relevant model that is based on observed data, regardless of its quality of fit, the final model ("m4") was selected
