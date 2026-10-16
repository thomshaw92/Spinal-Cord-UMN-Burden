# Bayesian longitudinal and brain-cord coupling models (Supplementary Table S2, Figures S1 and S2).
suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
})
options(brms.backend = "cmdstanr", mc.cores = parallel::detectCores())

here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
data <- Sys.getenv("DATA", file.path(here, "..", "data"))
out <- file.path(here, "..", "output")
dir.create(out, showWarnings = FALSE)

s <- read_csv(file.path(data, "sessions.csv"), show_col_types = FALSE)
demo <- read_csv(file.path(data, "demographics.csv"), show_col_types = FALSE)
lev <- read_csv(file.path(data, "levels.csv"), show_col_types = FALSE)

d <- s %>% left_join(demo, by = c("participant", "session")) %>%
  mutate(age_c = age_at_scan - mean(age_at_scan[session == 1]),
         time_mo = months_since_enrolment,
         group = factor(phenotype, levels = c("classic-ALS", "slow-ALS", "UMN-dominant", "PLS")),
         spectrum = factor(ifelse(phenotype %in% c("PLS", "UMN-dominant"), "UMN", "ALS"), levels = c("ALS", "UMN")),
         subses = paste(participant, session))

priors <- c(prior(normal(0, 1), class = "b"), prior(normal(0, 1), class = "Intercept"),
            prior(normal(0, 1), class = "sd"), prior(normal(0, 1), class = "sigma"))

fit <- function(form, dl, adapt) {
  brm(form, data = dl, prior = priors, chains = 4, iter = 4000, warmup = 2000,
      control = list(adapt_delta = adapt, max_treedepth = 12), seed = 1729, refresh = 0)
}

summarise_fit <- function(f, ...) {
  post <- posterior::as_draws_df(f)
  keep <- posterior::variables(f)[grepl("^b_|^sd_|^sigma$", posterior::variables(f))]
  posterior::summarise_draws(posterior::subset_draws(post, variable = keep), mean, sd,
                             ~quantile(.x, c(0.025, 0.5, 0.975)), rhat = posterior::rhat,
                             ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail) %>%
    rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`) %>%
    mutate(pd = map_dbl(variable, ~ mean(sign(post[[.x]]) == sign(median(post[[.x]])))), ...)
}

marginal <- function(f, term, model, metric) {
  nd <- expand_grid(time_mo = seq(0, 12, length.out = 40), !!term := levels(f$data[[term]]),
                    age_c = 0, sex = levels(factor(f$data$sex))[1])
  if ("vertlevel" %in% names(f$data)) nd$vertlevel <- levels(f$data$vertlevel)[1]
  bind_cols(nd, as_tibble(fitted(f, newdata = nd, re_formula = NA, probs = c(0.025, 0.975)))) %>%
    rename(level = all_of(term), fit = Estimate, lo = Q2.5, hi = Q97.5) %>%
    mutate(model = model, metric = metric)
}

LEVELS <- list(mtr = 2:5, fa = 2:5, gma = 2:5, csa = 2:3, fdc = NULL)
post <- list()
me <- list()
for (term in c("group", "spectrum")) {
  for (m in names(LEVELS)) {
    if (is.null(LEVELS[[m]])) {
      dl <- d %>% select(participant, time_mo, all_of(term), age_c, sex, val = all_of(m)) %>% drop_na()
      form <- as.formula(paste("val_z ~ time_mo *", term, "+ age_c + sex + (1 | participant)"))
      adapt <- if (term == "group") 0.95 else 0.97
    } else {
      cols <- paste0(m, "_c", LEVELS[[m]])
      dl <- d %>% select(participant, session, subses, time_mo, all_of(term), age_c, sex) %>%
        left_join(lev %>% select(participant, session, all_of(cols)), by = c("participant", "session")) %>%
        pivot_longer(all_of(cols), names_to = "vertlevel", values_to = "val") %>%
        mutate(vertlevel = factor(toupper(sub("^.*_", "", vertlevel)))) %>% drop_na(val)
      form <- as.formula(paste("val_z ~ time_mo *", term,
                               "+ age_c + sex + (1 | participant) + (1 | subses) + (1 | vertlevel)"))
      adapt <- 0.97
    }
    dl$val_z <- as.numeric(scale(dl$val))
    f <- fit(form, dl, adapt)
    post[[paste(term, m)]] <- summarise_fit(f, model = term, metric = m, n_obs = nrow(dl))
    if (term == "spectrum" || m %in% c("mtr", "fa")) me[[paste(term, m)]] <- marginal(f, term, term, m)
  }
}
write_csv(bind_rows(post), file.path(out, "longitudinal_posteriors.csv"))
write_csv(bind_rows(me), file.path(out, "longitudinal_marginal_effects.csv"))

PAIRS <- tribble(~brain, ~cord, "fdc", "mtr", "fdc", "fa", "fd", "fa", "logfc", "csa")
post <- list()
me <- list()
for (i in seq_len(nrow(PAIRS))) {
  b <- PAIRS$brain[i]
  c <- PAIRS$cord[i]
  dl <- d %>% select(participant, group, age_c, sex, brain = all_of(b), cord = all_of(c)) %>% drop_na() %>%
    mutate(brain_z = as.numeric(scale(brain)), cord_z = as.numeric(scale(cord)))
  f <- fit(cord_z ~ brain_z * group + age_c + sex + (1 | participant), dl, 0.95)
  post[[i]] <- summarise_fit(f, brain = b, cord = c, n_obs = nrow(dl))
  r <- range(dl$brain_z)
  nd <- expand_grid(brain_z = seq(r[1], r[2], length.out = 40), group = levels(dl$group),
                    age_c = 0, sex = levels(factor(dl$sex))[1])
  me[[i]] <- bind_cols(nd, as_tibble(fitted(f, newdata = nd, re_formula = NA, probs = c(0.025, 0.975)))) %>%
    rename(fit = Estimate, lo = Q2.5, hi = Q97.5) %>% mutate(brain = b, cord = c)
}
write_csv(bind_rows(post), file.path(out, "coupling_posteriors.csv"))
write_csv(bind_rows(me), file.path(out, "coupling_marginal_effects.csv"))
