# Propensity score matching simulation
rm(list = ls())
library(lava)
library(survival)
library(prodlim)
library(riskRegression)
source("load-functions.R")
# CREATE LVM OBJECT -------------------------------------------------------
# age is std norm dist. Will set dist for sex later
m <- lvm(A ~ age + sex)

# Set dist for A, sex as binom
distribution(m, ~A) <- binomial.lvm()
distribution(m, ~sex) <- binomial.lvm(p = 0.4)

# Add time to event and time to censor variable with Cox-Weibull dist
distribution(m, ~event_time) <- coxWeibull.lvm()
distribution(m, ~cen_time) <- coxWeibull.lvm()


# Can't observe both event and censor time for a person, so chose min. of two
# and # Add "status" variable
m <- eventTime(m, time ~ min(event_time=1, cen_time=0), "status")

# Specify "true" beta for log reg
regression(m, A ~ age + sex) <- c(0.3, 0.7)

# specify "true" beta for cox model 
regression(m, event_time ~ A + age + sex) <- c(0,-0.4, -0.4)


# SIMULATE ----------------------------------------------------------------
set.seed(13)
x <- sim(m, 1000)
head(x)
glm(A~age+sex, data = x,family=binomial)
coxph(Surv(time, status) ~ A , data = x)

coxph(Surv(time, status) ~ A + age + sex , data = x)

plot(prodlim(Surv(time, status) ~ A, data = x))



# FUNCTION ----------------------------------------------------------------

run <- function(..., n = 1000) {
  # takes lvm object "m"\
  # browser()
  d <- simulate(m, n = n)
  setDT(d)
  ps.fit <- glm(A~age+sex,data=d, family = binomial)
  ps.score <- predictRisk(ps.fit,newdata=d)
  matched.x <- match.nn(A=d$A,dt=d,ps=ps.score, M=1)
  f <- coxph(Surv(time,status)~A,data=matched.x)
  summary(f)
  structure(exp(coef(f)["A"]),
            names = c("A"))
}
x <- run()
x


x <- sim(run, 100)
x

print(x)
summary(x,
        estimate = c('A'),
        true = c(exp(0))
density(x)

