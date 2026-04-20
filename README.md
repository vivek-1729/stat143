# Around the Goal — Replication Study

Replication of Nevo & Ritov (2012), *"Around the goal: Examining the effect of the first goal on the second goal in soccer using survival analysis methods"*, using Premier League data.

## Overview

The paper models the time to the first two home-team goals in a soccer match as a survival problem. Each game produces two observations:

- **FirstGoalTime** (`obs_type=1`): time until the home team scores its first goal, censored by the first away goal or minute 90
- **SecondGoalTime** (`obs_type=2`): left-truncated at the time of the first goal in the game (regardless of team), tracking when the home team scores next

The key finding: a first goal makes the next goal *less* likely immediately, but this effect reverses as the game progresses — captured by a log(time since first goal) term in a Cox proportional hazards model.

## Data

| File | Description |
|------|-------------|
| `eng-premier-league.csv` | Minute-level goal data for all PL seasons (1990–2025) |
| `E0-2.csv` | Betting odds, 2009-10 season |
| `E0-3.csv` | Betting odds, 2007-08 season |
| `E0-4.csv` | Betting odds, 2010-11 season |
| `E0-25.csv` | Betting odds, 2024-25 season (goal data incomplete for this season) |

**Seasons used:** 2007-08, 2009-10, 2010-11 — 1,140 games, 2,180 observations.

`ProbWin` (home team win probability) is computed from BetBrain average odds (`BbAvH/D/A`), the same source as the original paper.

> Note: `eng-premier-league.csv` only has goal data through November 2024 for the 2024-25 season (96/380 games), so that season is excluded.

## Pipeline

### Step 1 — Build the survival dataset

```bash
python3 process_data.py
```

Outputs `survival_data.csv` with columns:

| Column | Description |
|--------|-------------|
| `game_id` | Unique game identifier |
| `obs_type` | 1 = FirstGoalTime, 2 = SecondGoalTime |
| `tstart` / `tstop` | Start-stop interval (counting process format) |
| `event` | 1 if home team scored, 0 if censored |
| `ProbWin` | Implied home win probability from avg bookmaker odds |
| `Season` | Season indicator (0=2007-08, 1=2009-10, 2=2010-11) |
| `Goal` | 0 for obs_type=1, 1 for obs_type=2 |
| `TimeOfFirstGoal` | Minute of first goal in game (0 for obs_type=1) |
| `FirstGoalTeam` | 1=home scored first, 0=away (only meaningful for obs_type=2) |

### Step 2 — Fit models and generate output

```bash
Rscript analysis.R
```

Requires the `survival` package. Outputs:

- `model_tables.csv` — coefficient tables for all models
- `fig1_cumulative_baseline_hazards.pdf` — stratified baseline hazards (paper Figure 1)
- `fig2_survival_by_probwin.pdf` — FirstGoalTime survival by ProbWin (paper Figure 2)
- `fig3_survival_by_T1.pdf` — SecondGoalTime survival by TimeOfFirstGoal (paper Figure 3)
- `fig4_baseline_hazard.pdf` — overall cumulative baseline hazard (paper Figure 4)

## Models

All models use a Cox PH hazard (counting-process format). The indicator `I{j=2}` terms are implemented by zero-padding covariates for `obs_type=1` rows.

| Model | Formula | Notes |
|-------|---------|-------|
| III | `ProbWin + Season` | Null model |
| II | `+ Goal + TimeOfFirstGoal + FirstGoalTeam + log(t−T1) + frailty` | Full model (paper Table 2) |
| IV | `+ TimeOfFirstGoal` | Parsimonious (paper Table 3) |
| V | `+ Goal + TimeOfFirstGoal + (t−T1)` | Linear time-from-first-goal (paper Table 4) |
| VI | `+ Goal + TimeOfFirstGoal + log(t−T1)` | **Chosen model** (paper Table 5) |
| VII | Stratified baseline by obs_type | For baseline hazard figures |

`log(t−T1)` is a time-dependent covariate handled via `tt()` in R's `coxph`.

## Key Results (Model VI)

Compared to the paper (our results → paper results):

| Variable | coef | exp(coef) | p |
|----------|------|-----------|---|
| ProbWin | 1.644 → 1.915 | 5.17 → 6.75 | <0.001 |
| Goal | −0.536 → −0.594 | 0.585 → 0.552 | 0.009 |
| TimeOfFirstGoal | 0.004 → 0.011 | 1.004 → 1.011 | 0.241 |
| **log(TimeFromFirstGoal)** | **0.158 → 0.160** | **1.17 → 1.17** | **0.006** |

The log(TimeFromFirstGoal) coefficient is nearly identical. Frailty variance θ ≈ 5×10⁻⁷ (effectively 0), replicating the paper's conclusion that game-level random effects are not needed.

## Reference

Nevo, D. & Ritov, Y. (2012). Around the goal: Examining the effect of the first goal on the second goal in soccer using survival analysis methods. *arXiv:1207.6796*.
