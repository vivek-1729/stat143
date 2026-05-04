# analysis.R
# ----------
# Replicates Nevo & Ritov (2012) "Around the goal".
# Uses Premier League data for 2007-08, 2009-10, 2010-11 (3 seasons).
# Requires: survival_data.csv produced by process_data.py
#
# Models fitted (matching paper notation):
#   III  -- null: ProbWin + Season
#   II   -- full with frailty: + Goal + TimeOfFirstGoal + FirstGoalTeam +
#            log(TimeFromFirstGoal) + frailty(game_id)
#   IV   -- ProbWin + Season + Goal:TimeOfFirstGoal
#   V    -- ProbWin + Season + Goal + Goal:TimeOfFirstGoal + Goal:TimeFromFirstGoal
#   VI   -- ProbWin + Season + Goal + Goal:TimeOfFirstGoal + Goal:log(TimeFromFirstGoal)
#   VII  -- stratified baseline (no interaction), for Figure 1

library(survival)

# ── Load data ────────────────────────────────────────────────────────────────

df <- read.csv("survival_data.csv", stringsAsFactors = FALSE)
df$date   <- as.Date(df$date)
# Season as factor (0=2007-08, 1=2009-10, 2=2010-11); paper used binary for 2 seasons.
# We include as factor so each season gets its own indicator vs. the baseline.
df$SeasonF <- factor(df$Season)

cat("Observations:", nrow(df), "\n")
cat("obs_type=1:", sum(df$obs_type == 1), " | obs_type=2:", sum(df$obs_type == 2), "\n")
cat("Events (obs1):", sum(df$obs_type == 1 & df$event == 1),
    " | Events (obs2):", sum(df$obs_type == 2 & df$event == 1), "\n")
cat("Mean ProbWin:", round(mean(df$ProbWin, na.rm = TRUE), 3), "\n")
cat("Mean TimeOfFirstGoal (obs2):",
    round(mean(df$TimeOfFirstGoal[df$obs_type == 2]), 1), "\n")
cat("Seasons:", paste(levels(df$SeasonF), collapse = ", "), "\n\n")

# Drop rows with missing ProbWin
df <- df[!is.na(df$ProbWin), ]

# Survival object (counting-process / start-stop format)
surv_obj <- with(df, Surv(tstart, tstop, event))

# ── Helper: print a formatted coefficient table ───────────────────────────────

print_table <- function(model, title, digits = 3) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat(" ", title, "\n")
  cat(strrep("=", 60), "\n", sep = "")
  s <- summary(model)$coefficients
  # Remove frailty rows for cleaner display
  s <- s[!grepl("frailty", rownames(s)), , drop = FALSE]
  coefs <- s[, "coef"]
  ses   <- s[, "se(coef)"]
  # Frailty models return Chisq instead of z; regular models return z (= coef/se)
  if ("z" %in% colnames(s)) {
    zvals <- s[, "z"]
    pvals <- s[, "Pr(>|z|)"]
  } else {
    zvals <- sign(coefs) * sqrt(s[, "Chisq"])
    pvals <- s[, "p"]
  }
  out <- data.frame(
    coef       = round(coefs,        digits),
    exp_coef   = round(exp(coefs),   digits),
    se_coef    = round(ses,          digits),
    z          = round(zvals,        digits),
    p          = signif(pvals,       digits)
  )
  colnames(out) <- c("coef", "exp(coef)", "se(coef)", "z", "p")
  print(out)
  invisible(s)
}

# ── tt() function for log(t - TimeOfFirstGoal) ───────────────────────────────
# For obs_type=1 rows TimeOfFirstGoal=0 → return 0 (no contribution).
# For obs_type=2 rows TimeOfFirstGoal>0 → return log(t - T1), floored at 0.5
# to avoid log(0) issues when t = T1.
log_time_from_first <- function(x, t, ...) {
  ifelse(x == 0, 0, log(pmax(t - x, 0.5)))
}

# ── Model III: null (ProbWin + Season) ──────────────────────────────────────

m3 <- coxph(surv_obj ~ ProbWin + SeasonF,
            data = df, ties = "efron")
print_table(m3, "Model III (null): ProbWin + Season")

# ── Model II: full model with frailty ───────────────────────────────────────
# Covariates for I{j=2} are zero-padded for obs_type=1 rows,
# so they contribute 0 to the linear predictor for those rows.

m2 <- coxph(surv_obj ~ ProbWin + SeasonF +
              Goal + TimeOfFirstGoal + FirstGoalTeam +
              tt(TimeOfFirstGoal) +
              frailty(game_id, distribution = "gamma"),
            data = df, ties = "efron",
            tt = log_time_from_first)
print_table(m2, "Model II (full + frailty): paper Table 2 equivalent")

# Frailty variance test (H0: theta = 0)
frailty_key <- names(m2$history)[1]
theta_val   <- m2$history[[frailty_key]]$theta
cat("\nFrailty variance (theta):", formatC(theta_val, digits = 4, format = "g"),
    "  (near 0 → frailty not needed, consistent with paper)\n")

# ── Model IV: ProbWin + Season + TimeOfFirstGoal for j=2 ────────────────────

m4 <- coxph(surv_obj ~ ProbWin + SeasonF + TimeOfFirstGoal,
            data = df, ties = "efron")
print_table(m4, "Model IV: paper Table 3 equivalent")

# ── Model V: + Goal + TimeFromFirstGoal (linear, time-dependent) ─────────────

m5 <- coxph(surv_obj ~ ProbWin + SeasonF + Goal + TimeOfFirstGoal +
              tt(TimeOfFirstGoal),
            data = df, ties = "efron",
            tt = function(x, t, ...) ifelse(x == 0, 0, t - x))
print_table(m5, "Model V: paper Table 4 equivalent")

# ── Model VI: chosen — log(TimeFromFirstGoal) time-dependent ─────────────────

m6 <- coxph(surv_obj ~ ProbWin + SeasonF + Goal + TimeOfFirstGoal +
              tt(TimeOfFirstGoal),
            data = df, ties = "efron",
            tt = log_time_from_first)
print_table(m6, "Model VI (chosen): paper Table 5 equivalent")

# ── Likelihood ratio tests ───────────────────────────────────────────────────

cat("\n", strrep("=", 60), "\n", sep = "")
cat(" Likelihood Ratio Tests\n")
cat(strrep("=", 60), "\n", sep = "")

lrt <- function(m_full, m_null, df_diff, label) {
  stat <- 2 * (m_full$loglik[2] - m_null$loglik[2])
  p    <- pchisq(stat, df = df_diff, lower.tail = FALSE)
  cat(sprintf("  %-35s chi2=%.3f  df=%d  p=%.3f\n", label, stat, df_diff, p))
}

lrt(m4, m3, 1,  "Model IV vs III (TimeOfFirstGoal):")
lrt(m6, m3, 3,  "Model VI vs III (Goal+T1+logTfT1):")
lrt(m6, m4, 2,  "Model VI vs IV  (Goal+logTfT1):")
lrt(m5, m4, 2,  "Model V  vs IV  (Goal+TfT1):")
cat(sprintf("  %-35s AIC(IV)=%.1f  AIC(V)=%.1f  AIC(VI)=%.1f\n",
            "AIC comparison (IV / V / VI):",
            AIC(m4), AIC(m5), AIC(m6)))

# ── Model VII: stratified baseline ──────────────────────────────────────────

m7 <- coxph(surv_obj ~ ProbWin + SeasonF + strata(obs_type),
            data = df, ties = "efron")

# ── Figures ──────────────────────────────────────────────────────────────────
# Note: basehaz() does not work on models with tt() terms. We use:
#   - m7 (stratified, no tt) for Figures 1 & 4
#   - m3 (null) baseline + m6 coefficients for manual survival curve computation

cat("\nGenerating figures...\n")

mean_pw <- mean(df$ProbWin, na.rm = TRUE)
t_grid  <- seq(1, 89, by = 1)

# Helper: interpolate cumulative baseline hazard onto t_grid
interp_H0 <- function(bh, t) {
  H <- approx(bh$time, bh$hazard, xout = t, rule = 2)$y
  H[is.na(H)] <- 0
  H
}

# --- Figure 1: Cumulative baseline hazards for Model VII (stratified) --------

bh7 <- basehaz(m7, centered = FALSE)
bh7_1 <- bh7[bh7$strata == "obs_type=1", ]
bh7_2 <- bh7[bh7$strata == "obs_type=2", ]

pdf("fig1_cumulative_baseline_hazards.pdf", width = 7, height = 5)
par(mar = c(4, 4.5, 2, 4.5))
ylim_max <- max(bh7$hazard) * 1.1
plot(bh7_1$time, bh7_1$hazard, type = "s", lwd = 2, col = "black",
     xlim = c(0, 90), ylim = c(0, ylim_max),
     xlab = "Time (minutes)", ylab = "Cumulative Baseline Hazard",
     main = "Figure 1: Cumulative Baseline Hazard (Model VII)")
lines(bh7_2$time, bh7_2$hazard, type = "s", lwd = 2, col = "steelblue", lty = 2)
# Hazard ratio overlaid on right axis
t_com <- seq(5, 88, by = 1)
h1_a  <- interp_H0(bh7_1, t_com)
h2_a  <- interp_H0(bh7_2, t_com)
ratio <- ifelse(h1_a > 0, h2_a / h1_a, NA)
par(new = TRUE)
plot(t_com, ratio, type = "l", lwd = 1.5, col = "firebrick", lty = 3,
     axes = FALSE, xlab = "", ylab = "", ylim = c(0, max(ratio, na.rm = TRUE) * 1.1))
axis(4, col = "firebrick", col.axis = "firebrick", las = 1)
mtext("Hazard Ratio (2nd/1st)", side = 4, line = 3, col = "firebrick", cex = 0.8)
legend("topleft",
       legend = c("First Goal", "Second Goal", "Hazard Ratio"),
       col    = c("black", "steelblue", "firebrick"),
       lwd    = c(2, 2, 1.5), lty = c(1, 2, 3), bty = "n")
dev.off()

# --- Figure 2: Survival curves for different ProbWin (FirstGoalTime) ---------
# Use m7 stratum 1 baseline (no tt, so survfit works); apply m6 ProbWin coef.

bh_base <- bh7_1   # baseline for obs_type=1 (FirstGoal)
H0_base <- interp_H0(bh_base, t_grid)

cols    <- c("navy", "steelblue", "darkorange", "firebrick")
pw_vals <- c(0.15, 0.35, 0.55, 0.75)

pdf("fig2_survival_by_probwin.pdf", width = 7, height = 5)
par(mar = c(4, 4, 2, 1))
plot(NA, xlim = c(0, 90), ylim = c(0, 1),
     xlab = "Time (minutes)", ylab = "Survival Function",
     main = "Figure 2: Model VI Survival Curves by ProbWin")
for (i in seq_along(pw_vals)) {
  # For obs_type=1: Goal=0, T1=0, SeasonF=0 (baseline season) → those terms = 0
  lp <- coef(m6)["ProbWin"] * pw_vals[i]
  S  <- exp(-H0_base * exp(lp))
  lines(t_grid, S, col = cols[i], lwd = 2)
}
legend("topright",
       legend = paste0("ProbWin=", pw_vals),
       col = cols, lwd = 2, bty = "n")
dev.off()

# --- Figure 3: SecondGoalTime survival curves for different T1 ---------------
# Use m7 stratum 2 baseline + m6 coefficients; integrate numerically because
# the log(t-T1) term makes the LP time-varying.

bh2_step <- bh7_2   # baseline cumulative hazard for SecondGoalTime

T1_vals <- c(25, 45, 65)
cols3   <- c("steelblue", "darkorange", "firebrick")

pdf("fig3_survival_by_T1.pdf", width = 7, height = 5)
par(mar = c(4, 4, 2, 1))
plot(NA, xlim = c(0, 90), ylim = c(0, 1),
     xlab = "Time (minutes)", ylab = "Survival Function",
     main = "Figure 3: Model VI Second Goal Survival by TimeOfFirstGoal")

# "No first goal" reference curve: obs_type=1 at mean ProbWin
lp_ref <- coef(m6)["ProbWin"] * mean_pw
S_ref  <- exp(-H0_base * exp(lp_ref))
lines(c(0, t_grid), c(1, S_ref), col = "black", lwd = 2)

for (i in seq_along(T1_vals)) {
  T1      <- T1_vals[i]
  t_after <- t_grid[t_grid > T1]

  # H0 at T1 (truncation point) and at all later times
  H0_T1    <- interp_H0(bh2_step, T1)
  H0_after <- interp_H0(bh2_step, t_after)

  # Time-varying linear predictor at each t in t_after
  lp_t <- coef(m6)["ProbWin"] * mean_pw +
    coef(m6)["Goal"] * 1 +
    coef(m6)["TimeOfFirstGoal"] * T1 +
    coef(m6)["tt(TimeOfFirstGoal)"] * log(pmax(t_after - T1, 0.5))

  # Conditional survival: exp(-∫_{T1}^{t} h0(s) exp(lp(s)) ds)
  # Approximate by summing hazard increments × exp(lp at right endpoint)
  dH0 <- diff(c(H0_T1, H0_after))
  cum_haz <- cumsum(dH0 * exp(lp_t))
  S_cond  <- c(1, exp(-cum_haz))

  lines(c(T1, t_after), S_cond, col = cols3[i], lwd = 2)
}
legend("topright",
       legend = c("No First Goal", paste0("T1=", T1_vals)),
       col    = c("black", cols3), lwd = 2, bty = "n")
dev.off()

# --- Figure 4: Overall cumulative baseline hazard (use m7 pooled / m3) -------
# Paper uses Model VI baseline; we use m7's obs_type=1 as the closest proxy.

pdf("fig4_baseline_hazard.pdf", width = 7, height = 5)
par(mar = c(4, 4, 2, 1))
bh3 <- basehaz(m3, centered = FALSE)
plot(bh3$time, bh3$hazard, type = "s", lwd = 2,
     xlab = "Time (minutes)", ylab = "Cumulative Baseline Hazard",
     main = "Figure 4: Cumulative Baseline Hazard (Model III)")
dev.off()

cat("Figures saved: fig1–fig4 (PDF).\n")

# ── Export coefficient tables to CSV ─────────────────────────────────────────

extract_coefs <- function(model, model_name) {
  s     <- summary(model)$coefficients
  s     <- s[!grepl("frailty", rownames(s)), , drop = FALSE]
  coefs <- s[, "coef"]
  ses   <- s[, "se(coef)"]
  if ("z" %in% colnames(s)) {
    zvals <- s[, "z"];          pvals <- s[, "Pr(>|z|)"]
  } else {
    zvals <- sign(coefs) * sqrt(s[, "Chisq"]);  pvals <- s[, "p"]
  }
  data.frame(
    model    = model_name,
    variable = rownames(s),
    coef     = round(coefs,       3),
    exp_coef = round(exp(coefs),  3),
    se_coef  = round(ses,         3),
    z        = round(zvals,       3),
    p        = signif(pvals,      3),
    row.names = NULL
  )
}

tables <- rbind(
  extract_coefs(m3, "Model_III"),
  extract_coefs(m2, "Model_II"),
  extract_coefs(m4, "Model_IV"),
  extract_coefs(m5, "Model_V"),
  extract_coefs(m6, "Model_VI")
)
write.csv(tables, "model_tables.csv", row.names = FALSE)
cat("Coefficient tables saved → model_tables.csv\n")

# ── RQ3: Team-specific equalizer speed (home conceded first) ─────────────────
#
# Subset: obs_type=2 rows where the away team scored first (FirstGoalTeam=0).
# In these 422 observations the home team is trailing and we track whether
# (and how quickly) they score next — their "equalizer."
#
# Cox model: Surv(tstart, tstop, event) ~ ProbWin + SeasonF + TimeOfFirstGoal
#            + home_team
#
# team HRs are relative to Wigan (reference: most obs, all 3 seasons → stable).
# Reading is excluded: 0 equalization events → coefficient undefined.
# Teams present in only 1 season are flagged; their estimates conflate team
# and season effects (wider CIs make this visible).

library(ggplot2)

eq_df <- df[df$obs_type == 2 & df$FirstGoalTeam == 0, ]
eq_df <- eq_df[eq_df$home_team != "Reading", ]   # 0 events → inestimable

# Per-team summary for labels and reliability flag
team_info <- aggregate(
  cbind(n_obs = event, n_events = event) ~ home_team + Season,
  data   = eq_df,
  FUN    = length
)
team_summary <- aggregate(
  cbind(n_obs = event, n_events = event) ~ home_team,
  data = eq_df,
  FUN  = function(x) c(n = length(x), ev = sum(x))
)
# Rebuild cleanly
team_summary <- data.frame(
  home_team = names(table(eq_df$home_team)),
  n_obs     = as.integer(table(eq_df$home_team)),
  n_events  = as.integer(tapply(eq_df$event, eq_df$home_team, sum)),
  n_seasons = as.integer(tapply(eq_df$Season, eq_df$home_team,
                                function(x) length(unique(x))))
)

# Reference team: Wigan (most observations, all 3 seasons)
eq_df$team_f <- relevel(factor(eq_df$home_team), ref = "Wigan")

eq_cox <- coxph(
  Surv(tstart, tstop, event) ~ ProbWin + SeasonF + TimeOfFirstGoal + team_f,
  data  = eq_df,
  ties  = "efron"
)

cat("\n", strrep("=", 60), "\n", sep = "")
cat(" RQ3: Team equalizer speed (home conceded first)\n")
cat(strrep("=", 60), "\n", sep = "")
cat("Observations:", nrow(eq_df), " | Events:", sum(eq_df$event),
    " | Teams:", nlevels(eq_df$team_f), "\n\n")

# Global test for team effects
lr_null <- coxph(Surv(tstart, tstop, event) ~ ProbWin + SeasonF + TimeOfFirstGoal,
                 data = eq_df, ties = "efron")
lr_stat <- 2 * (eq_cox$loglik[2] - lr_null$loglik[2])
lr_df   <- nlevels(eq_df$team_f) - 1
lr_p    <- pchisq(lr_stat, df = lr_df, lower.tail = FALSE)
cat(sprintf("Global LRT for team effects: chi2 = %.2f, df = %d, p = %.3f\n",
            lr_stat, lr_df, lr_p))

# Extract team coefficients + CIs
s_eq   <- summary(eq_cox)$coefficients
team_r <- grepl("^team_f", rownames(s_eq))

team_coefs <- data.frame(
  home_team = c("Wigan", sub("^team_f", "", rownames(s_eq)[team_r])),
  coef      = c(0,    s_eq[team_r, "coef"]),
  se        = c(0,    s_eq[team_r, "se(coef)"]),
  HR        = c(1,    exp(s_eq[team_r, "coef"])),
  CI_low    = c(1,    exp(s_eq[team_r, "coef"] - 1.96 * s_eq[team_r, "se(coef)"])),
  CI_high   = c(1,    exp(s_eq[team_r, "coef"] + 1.96 * s_eq[team_r, "se(coef)"])),
  p_val     = c(NA,   s_eq[team_r, "Pr(>|z|)"])
)

# Merge in team counts and reliability flag
team_coefs <- merge(team_coefs, team_summary, by = "home_team")
team_coefs$label <- paste0(team_coefs$home_team,
                           " (", team_coefs$n_events, "/", team_coefs$n_obs, ")")
team_coefs$reliable <- team_coefs$n_seasons == 3   # in all 3 seasons

# Sort by HR for plot
team_coefs <- team_coefs[order(team_coefs$HR), ]
team_coefs$y_pos <- seq_len(nrow(team_coefs))

# Print table
print_cols <- c("home_team", "HR", "CI_low", "CI_high", "p_val", "n_events", "n_obs", "n_seasons")
out_tbl <- team_coefs[order(team_coefs$HR, decreasing = TRUE), print_cols]
out_tbl[, c("HR","CI_low","CI_high")] <- lapply(out_tbl[, c("HR","CI_low","CI_high")], round, 3)
out_tbl$p_val <- ifelse(is.na(out_tbl$p_val), "ref", signif(out_tbl$p_val, 3))
print(out_tbl, row.names = FALSE)

write.csv(out_tbl, "figures/team_equalizer_table.csv", row.names = FALSE)

# ── Forest plot with 95% CIs ─────────────────────────────────────────────────
p_eq <- ggplot(team_coefs, aes(x = HR, y = y_pos)) +
  # CI bars — thinner for 1-season teams to signal lower reliability
  geom_segment(aes(x = CI_low, xend = CI_high, y = y_pos, yend = y_pos,
                   alpha = ifelse(reliable, 1, 0.45)), linewidth = 0.7) +
  # Point estimate — filled for reliable (3-season), hollow for 1–2 season
  geom_point(aes(shape = reliable, fill = HR > 1),
             size = 3, stroke = 0.9) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 19),
                     labels = c("FALSE" = "1–2 seasons", "TRUE" = "All 3 seasons"),
                     name   = "Data coverage") +
  scale_fill_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#2980b9"),
                    guide  = "none") +
  scale_alpha_identity() +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(breaks = team_coefs$y_pos, labels = team_coefs$label) +
  labs(
    title    = "Team response after conceding first: equalizer hazard multipliers",
    subtitle = paste0(
      "Cox model (n=", nrow(eq_df), " obs, ", sum(eq_df$event), " events) | ",
      "reference = Wigan | controls: ProbWin, Season, T1\n",
      "Label format: team (equalizations / games conceding first) | ",
      "hollow = 1–2 seasons only (wider CI expected)\n",
      sprintf("Global test for team effects: chi2=%.1f, df=%d, p=%.3f",
              lr_stat, lr_df, lr_p)
    ),
    x = "Hazard multiplier (>1 = faster equalizer than Wigan, controlling for strength)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    axis.text.y   = element_text(size = 9),
    legend.position = "bottom"
  )

ggsave("figures/team_equalizer_ranking.png", p_eq,
       width = 10, height = 12, dpi = 150)
cat("\nSaved → figures/team_equalizer_ranking.png\n")
