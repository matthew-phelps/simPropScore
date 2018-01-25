# Propensity score matching simulation
rm(list = ls())
library(lava)
library(survival)
library(prodlim)

# CREATE LVM OBJECT -------------------------------------------------------
# age is std norm dist. Will set dist for sex later
m <- lvm(A ~ age + sex)

# Set dist for A, sex as binom
distribution(m, ~A) <- binomial.lvm()
distribution(m, ~sex) <- binomial.lvm(p = 0.4)

# Add time to event and time to censor variable with Cox-Weibull dist
distribution(m, ~event_time) <- coxWeibull.lvm()
distribution(m, ~cen_time) <- coxWeibull.lvm()


# Can't observe both event and censor time for a person, so chose min. of two and
# # Add "status" variable
m <- eventTime(m, time ~ min(event_time, cen_time), "status")


regression(m, A ~ age + sex) <- c(0.3, 0.7)
regression(m, event_time ~ A + age + sex) <- c(-0.4, -0.4, 0)


# SIMULATE ----------------------------------------------------------------
set.seed(13)
x <- sim(m, 1000)

coxph(Surv(time, status) ~ A, data = x)

plot(prodlim(Surv(time, status) ~ A, data = x))
