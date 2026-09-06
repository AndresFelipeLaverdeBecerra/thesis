options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {stop("Usage: 02_run_cts5_clinical_benchmark.R INPUT_CSV OUTPUT_DIR")}
SRC <- normalizePath(args[[1]], mustWork = TRUE)
BASE <- args[[2]]
RES <- file.path(BASE, "results")
FIG <- file.path(BASE, "figures")
if (dir.exists(BASE) && length(list.files(BASE, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Output directory must be empty or absent")}

dir.create(RES, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("survival", "cmprsk", "riskRegression", "prodlim", "survAUC")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_packages) > 0) {stop(paste0("Missing R package(s): ",
                                               paste(missing_packages, collapse = ", "),
                                               "\nDo not install packages inside every batch job. Install them once in the ",
                                               "persistent R library/environment and rerun."))}

suppressPackageStartupMessages({
  library(survival)
  library(cmprsk)
  library(riskRegression)
  library(prodlim)
  library(survAUC)})

d0 <- read.csv(SRC, check.names = FALSE)

required_cols <- c("POST5_TIME_MONTHS", "POST5_EVENT_TYPE", "CTS5_COMPLETE_CASE", "CTS5_SCORE", 
                   "CTS5_RISK_GROUP_LITERATURE", "CTS5_STRICT_SCOPE_PROXY", "NPI",
                   "AGE_AT_DIAGNOSIS", "HORMONE_THERAPY", "CHEMOTHERAPY", "RADIO_THERAPY")

missing_cols <- setdiff(required_cols, names(d0))
if (length(missing_cols) > 0) {stop("Missing required columns: ", paste(missing_cols, collapse = ", "))}

num <- function(x) suppressWarnings(as.numeric(x))

d0$POST5_TIME_MONTHS <- num(d0$POST5_TIME_MONTHS)
d0$POST5_EVENT_TYPE  <- num(d0$POST5_EVENT_TYPE)
d0$CTS5_COMPLETE_CASE <- num(d0$CTS5_COMPLETE_CASE)
d0$CTS5_SCORE <- num(d0$CTS5_SCORE)
d0$CTS5_STRICT_SCOPE_PROXY <- num(d0$CTS5_STRICT_SCOPE_PROXY)
d0$NPI <- num(d0$NPI)
d0$AGE_AT_DIAGNOSIS <- num(d0$AGE_AT_DIAGNOSIS)

d <- d0[d0$CTS5_COMPLETE_CASE == 1 & is.finite(d0$POST5_TIME_MONTHS) & d0$POST5_TIME_MONTHS > 0 &
          d0$POST5_EVENT_TYPE %in% c(0, 1, 2) & is.finite(d0$CTS5_SCORE) & is.finite(d0$NPI) &
          is.finite(d0$AGE_AT_DIAGNOSIS),]

d$CTS5_GROUP <- factor(d$CTS5_RISK_GROUP_LITERATURE, levels = c("LOW", "INTERMEDIATE", "HIGH"))

if (any(is.na(d$CTS5_GROUP))) {stop("CTS5 risk-group coding contains unexpected/missing values.")}

normalize_binary_factor <- function(x) {
  z <- trimws(toupper(as.character(x)))
  out <- ifelse(z %in% c("YES", "Y", "TRUE", "1", "POSITIVE", "RECEIVED"), "YES",
                ifelse(z %in% c("NO", "N", "FALSE", "0", "NEGATIVE", "NOT RECEIVED"), "NO", z))
  factor(out)}

d$HORMONE_THERAPY_F <- normalize_binary_factor(d$HORMONE_THERAPY)
d$CHEMOTHERAPY_F <- normalize_binary_factor(d$CHEMOTHERAPY)
d$RADIO_THERAPY_F <- normalize_binary_factor(d$RADIO_THERAPY)

d$RFS_EVENT_CS <- as.integer(d$POST5_EVENT_TYPE == 1)
d$CTS5_Z <- as.numeric(scale(d$CTS5_SCORE))


# Descriptive distribution and outcome counts
group_counts <- do.call(rbind, lapply(levels(d$CTS5_GROUP), 
                                      function(g) {x <- d[d$CTS5_GROUP == g, ]
                                      data.frame(CTS5_GROUP = g, N = nrow(x), 
                                                 COMPOSITE_RFS_EVENTS = sum(x$POST5_EVENT_TYPE == 1),
                                                 COMPETING_DEATHS = sum(x$POST5_EVENT_TYPE == 2), 
                                                 OTHER_CENSOR = sum(x$POST5_EVENT_TYPE == 0),
                                                 EVENT_PCT = 100 * mean(x$POST5_EVENT_TYPE == 1), 
                                                 CTS5_MEAN = mean(x$CTS5_SCORE), CTS5_SD = sd(x$CTS5_SCORE), 
                                                 CTS5_MEDIAN = median(x$CTS5_SCORE), NPI_MEAN = mean(x$NPI), 
                                                 stringsAsFactors = FALSE)}))

write.csv(group_counts, file.path(RES, "02_cts5_group_event_counts.csv"), row.names = FALSE)

png(file.path(FIG, "01_cts5_score_distribution.png"), width = 1200, height = 800, res = 130)
hist(d$CTS5_SCORE, breaks = 35, main = "CTS5 score distribution — METABRIC post-5y complete cases", xlab = "CTS5 score")
abline(v = c(3.13, 3.86), lty = 2, lwd = 2)
dev.off()


# Kaplan-Meier: cause-specific RFS endpoint (competing deaths censored)
km_fit <- survfit(Surv(POST5_TIME_MONTHS, RFS_EVENT_CS) ~ CTS5_GROUP, data = d)

png(file.path(FIG, "02_km_cause_specific_rfs_by_cts5.png"), width = 1200, height = 800, res = 130)
plot(km_fit, xlab = "Months after 5-year landmark", ylab = "Cause-specific RFS survival probability",
     xlim = c(0, 120), lwd = 2, mark.time = FALSE, main = "Post-5y cause-specific RFS by CTS5 group")
legend("bottomleft", legend = levels(d$CTS5_GROUP), lty = 1, lwd = 2, bty = "n")
dev.off()


# Cumulative incidence: event 1 with other-cause death as competing event
ci_all <- cmprsk::cuminc(ftime = d$POST5_TIME_MONTHS, fstatus = d$POST5_EVENT_TYPE, group = d$CTS5_GROUP,
                         cencode = 0)

cause1_names <- grep(" 1$", names(ci_all), value = TRUE)
ci_rfs <- ci_all[cause1_names]
class(ci_rfs) <- class(ci_all)

png(file.path(FIG, "03_cif_composite_rfs_by_cts5.png"), width = 1200, height = 800, res = 130)
plot(ci_rfs, xlab = "Months after 5-year landmark", ylab = "Cumulative incidence of composite RFS event",
     xlim = c(0, 120), lwd = 2, main = "Composite RFS cumulative incidence by CTS5 group")
dev.off()

cuminc_at <- function(data, t0, failcode = 1) {
  z <- cmprsk::cuminc(ftime = data$POST5_TIME_MONTHS, fstatus = data$POST5_EVENT_TYPE, cencode = 0)
  nm <- grep(paste0(" ", failcode, "$"), names(z), value = TRUE)[1]
  if (is.na(nm) || length(nm) == 0) return(NA_real_)
  obj <- z[[nm]]
  idx <- which(obj$time <= t0)
  if (length(idx) == 0) return(0)
  obj$est[max(idx)]}

horizons <- c(24, 60, 120)
obs_cif <- do.call(rbind,lapply(levels(d$CTS5_GROUP), 
                                function(g) {x <- d[d$CTS5_GROUP == g, ]
                                do.call(rbind, lapply(horizons, 
                                                      function(t0) {data.frame(CTS5_GROUP = g, HORIZON_MONTHS_POST5 = t0, 
                                                                               YEARS_FROM_DIAGNOSIS = (60 + t0) / 12, 
                                                                               N_GROUP = nrow(x),
                                                                               OBSERVED_CIF_RFS = cuminc_at(x, t0, 1),
                                                                               OBSERVED_CIF_COMPETING_DEATH = cuminc_at(x, t0, 2))}))}))

write.csv(obs_cif, file.path(RES, "03_observed_cumulative_incidence_by_cts5.csv"), row.names = FALSE)


# Cause-specific Cox models
clinical_formula_cs <- as.formula("Surv(POST5_TIME_MONTHS, RFS_EVENT_CS) ~ 
                                  AGE_AT_DIAGNOSIS + NPI + HORMONE_THERAPY_F + CHEMOTHERAPY_F + RADIO_THERAPY_F")

cox_cts5 <- coxph(Surv(POST5_TIME_MONTHS, RFS_EVENT_CS) ~ CTS5_SCORE, data = d, x = TRUE, y = TRUE)
cox_cts5_z <- coxph(Surv(POST5_TIME_MONTHS, RFS_EVENT_CS) ~ CTS5_Z, data = d, x = TRUE, y = TRUE)
cox_cts5_group <- coxph(Surv(POST5_TIME_MONTHS, RFS_EVENT_CS) ~ CTS5_GROUP, data = d, x = TRUE, y = TRUE)
cox_npi <- coxph(Surv(POST5_TIME_MONTHS, RFS_EVENT_CS) ~ NPI, data = d, x = TRUE, y = TRUE)
cox_clinical <- coxph(clinical_formula_cs, data = d, x = TRUE, y = TRUE)
cox_table <- function(fit, model_name) {
  s <- summary(fit)
  cf <- as.data.frame(s$coefficients)
  ci <- as.data.frame(s$conf.int)
  data.frame(MODEL = model_name, TERM = rownames(cf), BETA = cf[, "coef"], HR = cf[, "exp(coef)"],
             HR_LOWER_95 = ci[, "lower .95"], HR_UPPER_95 = ci[, "upper .95"], P_VALUE = cf[, "Pr(>|z|)"],
             stringsAsFactors = FALSE, row.names = NULL)}

cox_results <- rbind(cox_table(cox_cts5, "CTS5_CONTINUOUS_PER_UNIT"), cox_table(cox_cts5_z, "CTS5_CONTINUOUS_PER_SD"),
                     cox_table(cox_cts5_group, "CTS5_CATEGORICAL"), cox_table(cox_npi, "NPI_CONTINUOUS"),
                     cox_table(cox_clinical, "CLINICAL_BASELINE"))

write.csv(cox_results, file.path(RES, "04_cause_specific_cox_results.csv"), row.names = FALSE)

ph_table <- function(fit, model_name) {
  z <- cox.zph(fit)
  tab <- as.data.frame(z$table)
  tab$TERM <- rownames(tab)
  tab$MODEL <- model_name
  rownames(tab) <- NULL
  tab[, c("MODEL", "TERM", setdiff(names(tab), c("MODEL", "TERM")))]}

ph_results <- rbind(ph_table(cox_cts5, "CTS5_CONTINUOUS"), ph_table(cox_cts5_group, "CTS5_CATEGORICAL"),
                    ph_table(cox_npi, "NPI_CONTINUOUS"), ph_table(cox_clinical, "CLINICAL_BASELINE"))

write.csv(ph_results, file.path(RES, "05_proportional_hazards_tests.csv"), row.names = FALSE)


# Fine-Gray competing-risk models
fg_table <- function(fit, model_name, terms) {
  beta <- as.numeric(fit$coef)
  se <- sqrt(diag(fit$var))
  z <- beta / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  data.frame(MODEL = model_name, TERM = terms, BETA = beta, SHR = exp(beta), 
             SHR_LOWER_95 = exp(beta - 1.96 * se), SHR_UPPER_95 = exp(beta + 1.96 * se), P_VALUE = p,
             stringsAsFactors = FALSE)}

fg_cts5 <- cmprsk::crr(ftime = d$POST5_TIME_MONTHS, fstatus = d$POST5_EVENT_TYPE,
                       cov1 = as.matrix(d["CTS5_SCORE"]), failcode = 1, cencode = 0)
mm_group <- model.matrix(~ CTS5_GROUP, data = d)[, -1, drop = FALSE]
fg_group <- cmprsk::crr(ftime = d$POST5_TIME_MONTHS, fstatus = d$POST5_EVENT_TYPE, cov1 = mm_group,
                        failcode = 1, cencode = 0)
fg_npi <- cmprsk::crr(ftime = d$POST5_TIME_MONTHS, fstatus = d$POST5_EVENT_TYPE, cov1 = as.matrix(d["NPI"]),
                      failcode = 1, cencode = 0)

mm_clin <- model.matrix(
  ~ AGE_AT_DIAGNOSIS + NPI + HORMONE_THERAPY_F + CHEMOTHERAPY_F + RADIO_THERAPY_F, data = d)[, -1, drop = FALSE]

fg_clinical <- cmprsk::crr(ftime = d$POST5_TIME_MONTHS, fstatus = d$POST5_EVENT_TYPE, cov1 = mm_clin, failcode = 1,
                           cencode = 0)

fg_results <- rbind(fg_table(fg_cts5, "CTS5_CONTINUOUS", "CTS5_SCORE"), 
                    fg_table(fg_group, "CTS5_CATEGORICAL", colnames(mm_group)), fg_table(fg_npi, "NPI_CONTINUOUS", "NPI"),
                    fg_table(fg_clinical, "CLINICAL_BASELINE", colnames(mm_clin)))

write.csv(fg_results, file.path(RES, "06_fine_gray_results.csv"), row.names = FALSE)


# Repeated 5x5-fold OOF clinical-baseline linear predictor
set.seed(20260831)
n_repeats <- 5
n_folds <- 5
oof_sum <- rep(0, nrow(d))
oof_count <- rep(0, nrow(d))

strata <- interaction(d$POST5_EVENT_TYPE, d$CTS5_GROUP, drop = TRUE)

make_stratified_folds <- function(strata, k) {fold <- integer(length(strata))
for (s in levels(strata)) {idx <- which(strata == s)
idx <- sample(idx, length(idx), replace = FALSE)
fold[idx] <- rep(seq_len(k), length.out = length(idx))}
fold}

for (r in seq_len(n_repeats)) {
  fold_id <- make_stratified_folds(strata, n_folds)
  for (f in seq_len(n_folds)) {
    tr <- which(fold_id != f)
    te <- which(fold_id == f)
    
    fit_oof <- coxph(clinical_formula_cs, data = d[tr, , drop = FALSE], x = TRUE, y = TRUE)
    pred <- predict(fit_oof, newdata = d[te, , drop = FALSE], type = "lp", reference = "zero")
    oof_sum[te] <- oof_sum[te] + as.numeric(pred)
    oof_count[te] <- oof_count[te] + 1}}

if (any(oof_count != n_repeats)) {stop("OOF prediction count failed: each patient must have one prediction per repeat.")}

d$CLINICAL_BASELINE_OOF_LP <- oof_sum / oof_count


# Harrell C and Uno C
harrell_from_score <- function(score, label) {
  cc <- survival::concordance(Surv(d$POST5_TIME_MONTHS, d$RFS_EVENT_CS) ~ score, reverse = TRUE)
  data.frame(MODEL = label, HARRELL_C = unname(cc$concordance), SE = unname(sqrt(cc$var)), stringsAsFactors = FALSE)}

# CTS5 and NPI are predefined scores. The clinical baseline is fitted, therefore its discrimination is evaluated with repeated OOF predictions.
score_cts5 <- d$CTS5_SCORE
score_npi <- d$NPI
score_clin_oof <- d$CLINICAL_BASELINE_OOF_LP

harrell_results <- rbind(harrell_from_score(score_cts5, "CTS5"), harrell_from_score(score_npi, "NPI"),
                         harrell_from_score(score_clin_oof, "CLINICAL_BASELINE_OOF_5X5"))

write.csv(harrell_results, file.path(RES, "07_harrell_c.csv"), row.names = FALSE)
surv_rsp <- Surv(d$POST5_TIME_MONTHS, d$RFS_EVENT_CS)
scores_for_uno <- list(CTS5 = score_cts5, NPI = score_npi, CLINICAL_BASELINE_OOF_5X5 = score_clin_oof)

uno_results <- do.call(rbind,
                       lapply(horizons, function(t0) {do.call(rbind, lapply(names(scores_for_uno), 
                                                                            function(nm) {val <- survAUC::UnoC(Surv.rsp = surv_rsp, Surv.rsp.new = surv_rsp,
                                                                                                               lpnew = scores_for_uno[[nm]], time = t0)
                                                                            data.frame(MODEL = nm, HORIZON_MONTHS_POST5 = t0, 
                                                                                       YEARS_FROM_DIAGNOSIS = (60 + t0) / 12, 
                                                                                       UNO_C = as.numeric(val), stringsAsFactors = FALSE)}))}))

write.csv(uno_results, file.path(RES, "08_uno_c.csv"), row.names = FALSE)


# Apparent competing-risk prediction metrics. These metrics are calculated from models fitted on the full eligible cohort.
# The fitted clinical baseline uses OOF predictions for concordance evaluation.
hist_formula <- Hist(POST5_TIME_MONTHS, POST5_EVENT_TYPE) ~ 1

csc_cts5 <- riskRegression::CSC(Hist(POST5_TIME_MONTHS, POST5_EVENT_TYPE) ~ CTS5_SCORE, data = d)
csc_npi <- riskRegression::CSC(Hist(POST5_TIME_MONTHS, POST5_EVENT_TYPE) ~ NPI, data = d)
csc_clinical <- riskRegression::CSC(Hist(POST5_TIME_MONTHS, POST5_EVENT_TYPE) ~ AGE_AT_DIAGNOSIS + NPI + 
                                      HORMONE_THERAPY_F + CHEMOTHERAPY_F + RADIO_THERAPY_F, data = d)

eval_times <- sort(unique(c(horizons, seq(6, 120, by = 6))))

score_obj <- riskRegression::Score(object = list(CTS5 = csc_cts5, NPI = csc_npi, CLINICAL_BASELINE = csc_clinical),
                                   formula = hist_formula, data = d, cause = 1, metrics = c("auc", "brier"),
                                   times = eval_times, conf.int = TRUE, null.model = FALSE, cens.method = "ipcw")

auc_df <- as.data.frame(score_obj$AUC$score)
brier_df <- as.data.frame(score_obj$Brier$score)

write.csv(auc_df, file.path(RES, "09_time_dependent_auc_apparent.csv"), row.names = FALSE)
write.csv(brier_df, file.path(RES, "10_time_dependent_brier_apparent.csv"), row.names = FALSE)

detect_time_col <- function(x) {
  hit <- intersect(c("times", "time"), names(x))
  if (length(hit) == 0) stop("Could not identify time column in Score output.")
  hit[1]}

detect_brier_col <- function(x) {
  hit <- intersect(c("Brier", "brier"), names(x))
  if (length(hit) == 0) stop("Could not identify Brier column in Score output.")
  hit[1]}

time_col <- detect_time_col(brier_df)
brier_col <- detect_brier_col(brier_df)

trapz_mean <- function(t, y, tau) {keep <- is.finite(t) & is.finite(y) & t <= tau
t <- t[keep]
y <- y[keep]
o <- order(t)
t <- t[o]
y <- y[o]
if (length(t) == 0) return(NA_real_)
t <- c(0, t)
y <- c(0, y)
sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2) / tau}

ibs_rows <- list()
ix <- 1
for (model_name in unique(brier_df$model)) {
  z <- brier_df[brier_df$model == model_name, ]
  for (tau in c(60, 120)) {
    ibs_rows[[ix]] <- data.frame(MODEL = model_name, TAU_MONTHS_POST5 = tau, YEARS_FROM_DIAGNOSIS = (60 + tau) / 12,
                                 IBS = trapz_mean(z[[time_col]], z[[brier_col]], tau), stringsAsFactors = FALSE)
    ix <- ix + 1}}

ibs_df <- do.call(rbind, ibs_rows)
write.csv(ibs_df, file.path(RES, "11_integrated_brier_score_apparent.csv"), row.names = FALSE)


# Calibration at primary 60-month post-landmark horizon
pred_risk <- list(CTS5 = as.numeric(riskRegression::predictRisk(csc_cts5, newdata = d, times = 60, cause = 1)),
                  NPI = as.numeric(riskRegression::predictRisk(csc_npi, newdata = d, times = 60, cause = 1)),
                  CLINICAL_BASELINE = as.numeric(riskRegression::predictRisk(csc_clinical, newdata = d, times = 60, cause = 1)))

make_decile <- function(x, nbin = 10) {r <- rank(x, ties.method = "first")
pmin(nbin, pmax(1, ceiling(nbin * r / length(r))))}

calibration_rows <- list()
ix <- 1

for (model_name in names(pred_risk)) {
  p <- pred_risk[[model_name]]
  bin <- make_decile(p, 10)
  for (b in sort(unique(bin))) {z <- d[bin == b, ]
  calibration_rows[[ix]] <- data.frame(MODEL = model_name, DECILE = b, N = nrow(z),
                                       MEAN_PREDICTED_RISK_60M = mean(p[bin == b]),
                                       OBSERVED_CIF_RFS_60M = cuminc_at(z, 60, 1),
                                       stringsAsFactors = FALSE)
  ix <- ix + 1}}

calibration_df <- do.call(rbind, calibration_rows)
write.csv(calibration_df, file.path(RES, "12_calibration_60m_apparent.csv"), row.names = FALSE)
png(file.path(FIG, "04_calibration_60m_apparent.png"), width = 1200, height = 900, res = 130)
plot(c(0, 1), c(0, 1), type = "n", xlab = "Mean predicted risk at 5 years post-landmark",
     ylab = "Observed cumulative incidence at 5 years post-landmark", main = "Calibration at primary 5-year post-landmark horizon")
abline(0, 1, lty = 2)
pch_values <- c(16, 17, 15)
i <- 1
for (model_name in unique(calibration_df$MODEL)) {z <- calibration_df[calibration_df$MODEL == model_name, ]
points(z$MEAN_PREDICTED_RISK_60M, z$OBSERVED_CIF_RFS_60M, type = "b",
       pch = pch_values[i])
i <- i + 1}
legend("topleft", legend = unique(calibration_df$MODEL), pch = pch_values, lty = 1, bty = "n")
dev.off()


# Strict-scope proxy sensitivity
ds <- d[d$CTS5_STRICT_SCOPE_PROXY == 1, ]

ds$RFS_EVENT_CS <- as.integer(ds$POST5_EVENT_TYPE == 1)
ds$CTS5_GROUP <- factor(ds$CTS5_RISK_GROUP_LITERATURE, levels = c("LOW", "INTERMEDIATE", "HIGH"))
cox_strict <- coxph(Surv(POST5_TIME_MONTHS, RFS_EVENT_CS) ~ CTS5_SCORE, data = ds)
fg_strict <- cmprsk::crr(ftime = ds$POST5_TIME_MONTHS, fstatus = ds$POST5_EVENT_TYPE, cov1 = as.matrix(ds["CTS5_SCORE"]),
                         failcode = 1, cencode = 0)

strict_summary <- rbind(data.frame(ANALYSIS = "CAUSE_SPECIFIC_COX", N = nrow(ds), EFFECT = "HR_PER_CTS5_UNIT",
                                   ESTIMATE = exp(coef(cox_strict)[1]), LOWER_95 = exp(confint(cox_strict)[1, 1]),
                                   UPPER_95 = exp(confint(cox_strict)[1, 2]), P_VALUE = summary(cox_strict)$coefficients[1, "Pr(>|z|)"]),
                        data.frame(ANALYSIS = "FINE_GRAY", N = nrow(ds), EFFECT = "SHR_PER_CTS5_UNIT", ESTIMATE = exp(fg_strict$coef[1]),
                                   LOWER_95 = exp(fg_strict$coef[1] - 1.96 * sqrt(fg_strict$var[1, 1])), 
                                   UPPER_95 = exp(fg_strict$coef[1] + 1.96 * sqrt(fg_strict$var[1, 1])), P_VALUE = 2 * pnorm(
                                     abs(fg_strict$coef[1] / sqrt(fg_strict$var[1, 1])), lower.tail = FALSE)))

write.csv(strict_summary, file.path(RES, "13_strict_scope_sensitivity.csv"), row.names = FALSE)


# Compact model-comparison summary
find_auc_col <- function(x) {hit <- intersect(c("AUC", "auc"), names(x))
if (length(hit) == 0) stop("Could not identify AUC column.")
hit[1]}

auc_time_col <- detect_time_col(auc_df)
auc_value_col <- find_auc_col(auc_df)

selected_auc <- auc_df[auc_df[[auc_time_col]] %in% horizons, c("model", auc_time_col, auc_value_col), drop = FALSE]
names(selected_auc) <- c("MODEL", "HORIZON_MONTHS_POST5", "AUC")

harrell_compact <- harrell_results[, c("MODEL", "HARRELL_C")]
uno60 <- uno_results[uno_results$HORIZON_MONTHS_POST5 == 60, c("MODEL", "UNO_C")]
names(uno60)[2] <- "UNO_C_60M"
ibs60 <- ibs_df[ibs_df$TAU_MONTHS_POST5 == 60, c("MODEL", "IBS")]
names(ibs60)[2] <- "IBS_60M_APPARENT"
auc60 <- selected_auc[selected_auc$HORIZON_MONTHS_POST5 == 60, c("MODEL", "AUC")]
names(auc60)[2] <- "AUC_60M_APPARENT"

comparison_summary <- Reduce(function(x, y) merge(x, y, by = "MODEL", all = TRUE), list(harrell_compact, uno60, auc60, ibs60))
write.csv(comparison_summary, file.path(RES, "14_model_comparison_summary.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(RES, "sessionInfo.txt"))
message("Analysis completed. Results written to: ", RES)
message("Figures written to: ", FIG)