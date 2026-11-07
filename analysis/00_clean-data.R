library(tidyverse)

heart_data <- read.csv("rmr_data.csv")

### Explortory Analysis ###

### Resting heart rate ###

mean(heart_data$Height..cm.)
sd(heart_data$Height..cm.)

## Resting metabolic rate ###

mean(heart_data$RMR..kcal.day.)
sd(heart_data$RMR..kcal.day.)


## Box plot ###


### BASED on scatter plot there is no correlation between resting Hr and rmr 
ggplot(heart_data, aes(x = Height..cm., y = RMR..kcal.day.)) +
  geom_point()


hr_fit <- lm(RMR..kcal.day. ~ Resting_HR..bpm., data = heart_data)
summary(hr_fit)

plot(hr_fit)

heihgt_fit <- lm(RMR..kcal.day. ~ Height..cm., data = heart_data)

summary(heihgt_fit)

plot(heihgt_fit)
