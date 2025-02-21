# CALCULATIONS FOR ÖRESUND WORKSHOP PRESENTATION
#
# Author: Sebastian Moretto Krog (sebastian.moretto.krog@regionh.dk)
# Presented at the Öresund Workshop 2025
#
# Please don't hesitate to contact me if you are interested 
# in discussing the contents
#
# Date: 14.02.2025


# Libraries
library(tidyverse)
library(pwr)
library(reshape2)
library(survival)
library(survminer)
library(ggsurvfit)
library(flexsurv)
library(simsurv)


# ----------------------------------------
# TCP model
# ----------------------------------------

# Define the dose-fractionations for the comparisons
interventions <- tibble(
  dose = c(9, 12, 7, 10, 14, 8, 12.5),
  fractions = c(3, 2, 5, 3, 2, 5, 3),
  total_dose = dose * fractions,
  Schedule = paste(dose, fractions, sep = " Gy x ")
)
#View(interventions)

# Calculate BED at some alpha beta ratio
BED <- function(fractions, dose_per_fraction, alpha_beta) {
  return(fractions * dose_per_fraction * (1 + (dose_per_fraction / alpha_beta)))
}

# Calculate the N-fraction equivalent dose
NfxED <- function(target_fractions, fractions, dose_per_fraction, alpha_beta) {
  total_dose <- fractions * dose_per_fraction
  BED_original <- BED(fractions, dose_per_fraction, alpha_beta)
  NfxED_output <- 
    0.5 * (
      sqrt(4*BED_original * target_fractions * alpha_beta + target_fractions^2*alpha_beta^2)
      - target_fractions * alpha_beta)
  
  ## Check calculations
  if (identical(
    round(BED_original, 2), 
    round(BED(target_fractions, NfxED_output/target_fractions, alpha_beta), 2))) {
    return(NfxED_output)
  }
  
  warning("Dose calculation does not match BED calculation")
  return(NA)
}

# Calculate the 3fxED
interventions <- interventions |>
  mutate(`3fxED_6` = round(NfxED(3, fractions, dose, 6), 2))

# TCP model from the HyTEC group, Soltys et al., 2021
# https://doi.org/10.1016/j.ijrobp.2020.11.021
#TCP <- function(D, D50, g50) {
#  exponent <- 4 * g50 * (D / D50 - 1)
#  tcp <- exp(exponent) / (1 + exp(exponent))
#  return(tcp)
#}
# g50 <- 0.814
# D50 <- 19.44

# We use a shortened version for the 3fxED
TCP_3fxED <- function(D) {
  return(1 / (25.9455 * exp(-0.16759 * D) + 1))
}

# Calculate the tcp for the interventions
interventions <- mutate(interventions, tcp = TCP_3fxED(`3fxED_6`))

# Datapoints for the TCP Plot
plot_data <- tibble(
  dose = seq(15, 40, by = 0.1),
  tcp = TCP_3fxED(dose) 
)

# TCP model plot without lines
ggplot(plot_data, aes(x = dose, y = tcp)) +
  geom_line(linewidth = 1) +
  labs(title = "HyTEC TCP model",
       x = "3fxED, Gy",
       y = "2-Year Tumor Control Probability") +
  scale_x_continuous(
    limits = c(0, 50), 
    breaks = seq(0, 50, by = 5),
    labels = c(seq(0, 40, by = 5), "", "")
  ) +
  scale_y_continuous(
    limits = c(0, 1), 
    breaks = 0:10/10,
    labels = c("0", as.character(1:9/10), "1")
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(-5, "pt"))

# TCP model plot with lines
ggplot(plot_data, aes(x = dose, y = tcp)) +
  geom_line(linewidth = 1) +
  geom_point(data = interventions, 
             aes(x = `3fxED_6`, y = tcp), shape = 4, size = 3) +
  geom_segment(data = interventions, 
               aes(x = `3fxED_6`, xend = `3fxED_6`, y = 0, yend = tcp), 
               linetype = "dashed", color = "black") +
  geom_segment(data = interventions, 
               aes(x = 0, xend = `3fxED_6`, y = tcp, yend = tcp), 
               linetype = "dashed", color = "black") +
#  geom_text(data = interventions, # Add labels
#            aes(x = `3fxED_6`, y = tcp, label = Schedule), 
#            hjust = -0.1, vjust = 1, color = "black", size = 4) +
  labs(title = "HyTEC TCP model",
       x = "3fxED, Gy",
       y = "2-Year Tumor Control Probability") +
  scale_x_continuous(
    limits = c(0, 50), 
    breaks = seq(0, 50, by = 5),
    labels = c(seq(0, 40, by = 5), "", "")
  ) +
  scale_y_continuous(
    limits = c(0, 1), 
    breaks = 0:10/10,
    labels = c("0", as.character(1:9/10), "1")
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(-5, "pt"))


# ----------------------------------------
# Sample size calculations
# ----------------------------------------

alpha <- 0.05 # Assuming a two-sided endpoint
sig_level <- alpha # It's the same thing
power <- 0.8 # For a phase 3 trial a power of 0.9 should be considered
beta <- 1 - power

# Prepare the intervention list
interventions <- interventions |> arrange(`3fxED_6`)

# Calculate n per group (sample size) depending on prop_a and prop_b
# Uses a z-test (which, for a 2 by 2 table is not so different from a chi-squared)
# See https://stats.stackexchange.com/questions/2391/what-is-the-relationship-between-a-chi-squared-test-and-test-of-equal-proportion
n_power_prop_a_b <- function(a, b) {
  h <- ES.h(a, b) # This calculates the H-statistic.
  if (h > 0.01) { # It can fail when the doses are too close or the same.
    calc <- pwr.2p.test(h, sig.level = alpha, power = power)
    return(calc$n)
  }
  return(NA)
}

# Create the pairwise comparison matrix for TCP values
tcp_values <- interventions$tcp

pairwise_matrix <- sapply(tcp_values, function(a) {
  sapply(tcp_values, function (b) { return(n_power_prop_a_b(a, b)) })
})

# Set row and column names for better readability
rownames(pairwise_matrix) <- gsub(" ", "", interventions$Schedule, fixed = TRUE)
colnames(pairwise_matrix) <- gsub(" ", "", interventions$Schedule, fixed = TRUE)

# Show the pairwise comparison matrix (raw values)
print(pairwise_matrix)

melted_matrix <- melt(pairwise_matrix)

# Plot the matrix
ggplot(data = melted_matrix, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient(low = "green", high = "red", trans="log",
                      space = "Lab", 
                      name="n per Group") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Var2, Var1, label = ceiling(value)), color = "black", size = 4) +
  ggtitle("Sample size for comparison based on TCP") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())


# ----------------------------------------
# Survival statistics
# ----------------------------------------

# Source: Shariq Mohammed 
# (https://shariq-mohammed.github.io/files/cbsa2019/2-power-and-sample-size.html)

# We use our previous calculations for expected 2 year control
tcp_a <- interventions %>% filter(Schedule == "12 Gy x 2") %>% pull(tcp)
tcp_b <- interventions %>% filter(Schedule == "12.5 Gy x 3") %>% pull(tcp)

# We can derive the HR from the proportions/risk at a given time point
HR_calc <- function(prop_a, prop_b) {
  (1-prop_a)/(1-prop_b)
}

# The estimated HR from the TCP model for our chosen dose-fractions
HR_calc(tcp_a, tcp_b) # 4.1
# Very optimistic. Firstly, the model might not be accurate for higher doses
# Secondly, because of how doses are prescribed differently, we
#   might not actually be using so widely different doses

# We select a smaller hazard ratio for the calculations
HR <- 2.5

# And adjust the expected proportion for the upper dose:
tcp_c <- 1-(1-tcp_a)/HR

# Such that it gives a HR of 2.5
HR_calc(tcp_a, tcp_c) # = 2.5

# We can calculate the number of events as follows
surv_events <- function(HR, alpha = 0.05, beta = 0.2) {
  Z_a <- qnorm(alpha/2, lower.tail = FALSE)
  Z_b <- qnorm(beta, lower.tail = FALSE)
  return (4*(Z_a + Z_b)^2/(log(HR))^2)
}

# For a HR of 2.5
surv_events(HR)

# If we assume a exponential survival function then we can estimate the
# median survival from a given proportion at a time point as follows:
median_survival <- function(proportion, time) {
  time * log(1/2)/log(proportion)
}

# And thus the median survival as predicted from the TCP model (24 months)
t_m <- median_survival(c(tcp_a, tcp_c), 24)

# To calculate the probability of events given median survival,
#   the time of accrual, and the followup time afterwards,
#   we can use the following function:
prob_event = function(median_survival, accrual, followup) {
  lambda <- log(2)/median_survival
  term1 <- exp(-lambda*followup) - exp(-lambda*(accrual + followup))
  probEvent <- 1 - (term1/(accrual*lambda))
  probEvent
}

# We can calculate the probability of events for each arm
# 36 month accrual, 12 month follow-up
p <- prob_event(t_m, 36, 12)

# Different ways to calculate means...
#p_mean <- 1/(0.5/p[[1]] + 0.5/p[[2]])
p_mean <- mean(p)

# Events needed given HR = 2.5
sample_size_total <- ceiling(surv_events(HR)/p_mean)
sample_size_total


# ----------------------------------------
# Extracting IPD and estimating survival
# ----------------------------------------

# Use IPDfromKM to estimate IPD https://biostatistics.mdanderson.org/shinyapps/IPDfromKM/
# This has been done for both curves from the Zeng et al. 2023 study
# https://doi.org/10.1016/j.ijrobp.2022.09.076
# The extracted individual patient data been saved in "data/zeng_ipd.csv"
zeng_ipd <- read_csv("data/zeng_ipd.csv") %>% 
  transmute(
    rectime = `Survival time`,
    status = Status,
    group = as.character(`Treatment group`),
    n3fxED = case_when(
      group == "24" ~ NfxED(3, 2, 12, 6),
      group == "28" ~ NfxED(3, 2, 14, 6)
    )
  )

# We can recreate the plotted survival curves
survfit2(Surv(rectime, status) ~ group, data = zeng_ipd) |>
  ggsurvfit(linewidth = 1, type = "risk") +
  coord_cartesian(xlim = c(0, 48), ylim = c(0, 1))

# Fit the flexible parametric survival model
# Here we use a 3 knot spline.
mod_flex <- flexsurvspline(Surv(rectime, status) ~ group, data = zeng_ipd, k = 3)

# Overlay it with the Zeng data
ggflexsurvplot(mod_flex, fun="cumhaz", ylim=c(0,1), conf.int=F, censor=F,
               title="3-knot spline fit on extracted IPD from Zeng. et al.")

# And we can make some predictions using our model
pred <- predict(mod_flex, type="cumhaz",
                times=c(6,12,24), 
                newdata = zeng_ipd %>% distinct(group))

# Gather the data for a nice table
pred_table <- pred %>% mutate(group=c(24,28)) %>% 
  unnest(.pred) %>% 
  transmute(group, .eval_time, cumhaz=.pred_cumhaz) %>% 
  pivot_wider(names_from=.eval_time, values_from=cumhaz, names_prefix = "time=") %>% 
  mutate(`TCP@time=24` = 1-TCP_3fxED(NfxED(3,2,group/2,6)))

# Hazard ratio for the Zeng fitted model
max(pred_table$`time=24`)/min(pred_table$`time=24`) # 1.81


# ----------------------------------------
# Simulation
# ----------------------------------------

# Redefine the model using the 3fxED as a continous variable
mod_flex_var <- flexsurv::flexsurvspline(Surv(rectime, status) ~ n3fxED, 
                                         data = zeng_ipd, k = 3)

# The HR is likely overestimated... we limit it to 2.5
simulation_doses <- tibble(n3fxED = NfxED(3,c(2,3),c(12,11.806), 6))

# And check our HR calculation
pred_trial <- predict(mod_flex_var, type="cumhaz",
                      times=c(24), 
                      newdata = simulation_doses)
HR_trial <- max(pred_trial$.pred_cumhaz)/min(pred_trial$.pred_cumhaz)
HR_trial # 2.5

# To simulate the trial, we use the sim surv functions.
# Source: simsurv package
mod_logcumhaz <- function(t, x, betas, knots) {
  # Obtain the basis terms for the spline-based log
  # cumulative hazard (evaluated at time t)
  basis <- flexsurv::basis(knots, log(t))
  
  res <- 
    betas[["gamma0"]] * basis[[1]] + 
    betas[["gamma1"]] * basis[[2]] +
    betas[["gamma2"]] * basis[[3]] +
    betas[["gamma3"]] * basis[[4]] +
    betas[["gamma4"]] * basis[[5]] +
    betas[["n3fxED"]] * x[["n3fxED"]]
  res
}

simulate_ipd <- function(N=250, h0=F) {
  # Create a data frame with the subject IDs and treatment covariates
  cov <- tibble(id = 1:N,
                group = rbinom(N, 1, 0.5),
                n3fxED = case_when(
                  group == 1 ~ simulation_doses$n3fxED[1],
                  group == 0 ~ simulation_doses$n3fxED[2]))
  
  # To test the null hypothesis we assume no difference in effect
  if (h0) cov$n3fxED <- mean(cov$n3fxED)
  
  # Simulate the event times for local failure
  dat <- simsurv(betas = mod_flex_var$coefficients, # "true" parameter values
                 x = cov,                   # covariate data for n individuals
                 knots = mod_flex_var$knots,    # knot locations for splines
                 logcumhazard = mod_logcumhaz,  # definition of log cum hazard
                 maxt = 48,               # up to 48 months of total follow-up
                 interval = c(1E-8,100000)) # interval for root finding
  
  # Simulate accrual
  acc <- dat |> mutate(
    acc_maxt = 1:N/N*36+12,
    acc_status = case_when(acc_maxt < eventtime ~ 0, T ~ status),
    acc_eventtime = case_when(acc_maxt < eventtime ~ acc_maxt, T ~ eventtime),
  )
  
  ipd <- left_join(acc, cov, by="id")
  ipd
}

# We will calculate the Cox PH, and report the result for the logrank test
logrank_test <- function(ipd) {
  logrank <- survdiff(formula = Surv(acc_eventtime, acc_status) ~ group, data = ipd)
  c(pval=logrank$pvalue, events=sum(logrank$obs), chisq=logrank$chisq)
}

cox_ph_test <- function(ipd) {
  c_lf <- coxph(Surv(acc_eventtime, acc_status) ~ group, data=ipd)
  c(hr = exp(c_lf$coefficients[["group"]]), p=summary(c_lf)$logtest["pvalue"])
}

simulate_trial <- function(N=234, h0=F, plot=F) {
  ipd <- simulate_ipd(N, h0)
  
  if (plot) {
    survfit2(Surv(acc_eventtime, acc_status) ~ group, data = ipd) |>
      ggsurvfit(linewidth = 1, type = "risk") +
      coord_cartesian(xlim = c(0, 48), ylim = c(0, 1)) +
      add_risktable() +
      add_censor_mark() +
      ggtitle(paste("Simulation", plot)) +
      theme(plot.title = element_text(size=18))
  } else {
    c(logrank=logrank_test(ipd), cph=cox_ph_test(ipd))
  }
}

run_1000_sims <- function() {
  # Set a seed for the simulations
  set.seed(543543)
  
  if (!file.exists("data/sim_true.Rds")) {
    print(system.time({ # Takes around 85+85 seconds for 1000 sims
      sim_true <- replicate(1000, simulate_trial())
      
      # For kicks, we can also simulate the h0 trials
      # This is not strictly necessary as we are testing all the other simtrials
      #   using p = 0.05
      sim_h0 <- replicate(1000, simulate_trial(h0=T))
    }))
    
    saveRDS(sim_true, "data/sim_true.Rds")
    saveRDS(sim_h0, "data/sim_h0.Rds")
    
  }
}

# Care.. it takes a long time!
run_1000_sims()
sim_true <- readRDS("data/sim_true.Rds")
sim_h0 <- readRDS("data/sim_h0.Rds")

sim <- tibble(
    p_logrank=sim_true[1,], 
    events=sim_true[2,], 
    chisq=sim_true[3,],
    hr=sim_true[4,],
    p_cph=sim_true[5,],
    h0=F
  ) |> add_row(tibble(
    p_logrank=sim_h0[1,], 
    events=sim_h0[2,], 
    chisq=sim_h0[3,], 
    hr=sim_h0[4,],
    p_cph=sim_h0[5,],
    h0=T
  ))

power <- filter(sim, h0==F) |> 
  summarize(power = sum((p_logrank < 0.05))/n()) |> pull()

alpha <- filter(sim, h0==T) |>
  summarize(alpha = sum((p_logrank< 0.05))/n()) |> pull()

binom.test(sum(filter(sim, h0==F)$p_logrank <= 0.05), length(filter(sim, h0==F)$p_logrank))

sim_plot <- sim |> mutate(p05 = p_cph >= 0.05)

# Volcanoplot
ggplot(data=filter(sim), aes(x=log10(hr), y=-log10(p_cph), color=h0)) + geom_point()
