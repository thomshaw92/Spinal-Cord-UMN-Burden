# Within-ALS sensitivity analysis: Penn UMN Scale against cord FA and brain FDC (Supplementary Table S2D).

suppressPackageStartupMessages(library(brms))

arg <- grep("--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(arg)) dirname(normalizePath(sub("--file=", "", arg))) else getwd()
data <- Sys.getenv("DATA", file.path(here, "..", "data"))
out <- file.path(here, "..", "output")

s <- read.csv(file.path(data, "sessions.csv"))
dm <- read.csv(file.path(data, "demographics.csv"))[c("participant", "session", "age_at_scan", "sex",
                                                       "months_since_onset")]
stopifnot(!anyDuplicated(s[c("participant", "session")]), !anyDuplicated(dm[c("participant", "session")]))
d <- merge(s, dm, by = c("participant", "session"), all.x = TRUE)
stopifnot(nrow(d) == nrow(s))
enrol_age <- mean(d$age_at_scan[d$session == 1])
d <- d[d$phenotype %in% c("slow-ALS", "classic-ALS"), ]

prep <- function(m) {
  x <- d[complete.cases(d[c(m, "pumns_total", "age_at_scan", "sex", "months_since_onset")]), ]
  stopifnot(all(x$months_since_onset > 0))
  x$y <- (x[[m]] - mean(x[[m]])) / sd(x[[m]])
  pm <- ave(x$pumns_total, x$participant)
  x$penn_between <- (pm - mean(x$pumns_total)) / sd(x$pumns_total)
  x$penn_within <- (x$pumns_total - pm) / sd(x$pumns_total)
  x$age_c <- x$age_at_scan - enrol_age
  x$sex <- factor(x$sex)
  ld <- log(x$months_since_onset)
  x$logdur_z <- (ld - mean(ld)) / sd(ld)
  x
}

priors <- c(prior(normal(0, 1), class = b), prior(normal(0, 1), class = Intercept),
            prior(normal(0, 1), class = sd), prior(normal(0, 1), class = sigma))
forms <- list(
  "1 age, sex" = y ~ penn_between + penn_within + age_c + sex + (1 | participant),
  "2 age, sex, log duration" = y ~ penn_between + penn_within + age_c + sex + logdur_z + (1 | participant))

rows <- list()
for (m in c("fa", "fdc")) {
  x <- prep(m)
  if (Sys.getenv("DRY_RUN") == "1") next
  for (v in names(forms)) {
    fit <- brm(forms[[v]], data = x, family = gaussian(), prior = priors,
               chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 1729,
               control = list(adapt_delta = 0.95, max_treedepth = 12),
               backend = "cmdstanr", refresh = 0, silent = 2)
    dr <- as_draws_df(fit)
    np <- nuts_params(fit)
    for (tm in c("b_penn_between", "b_penn_within")) {
      z <- dr[[tm]]
      rows[[length(rows) + 1]] <- data.frame(
        measure = toupper(m), model = v, term = sub("^b_", "", tm), mean = mean(z),
        lower95 = unname(quantile(z, 0.025)), upper95 = unname(quantile(z, 0.975)),
        pd = max(mean(z > 0), mean(z < 0)), n_sessions = nrow(x),
        n_participants = length(unique(x$participant)),
        divergent = sum(np$Value[np$Parameter == "divergent__"]), max_rhat = max(rhat(fit), na.rm = TRUE))
    }
  }
}
if (length(rows)) {
  dir.create(out, showWarnings = FALSE)
  write.csv(do.call(rbind, rows), file.path(out, "als_sensitivity_posteriors.csv"), row.names = FALSE)
}
