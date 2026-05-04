#!/usr/bin/env python3
"""
process_data.py
---------------
Converts minute-level goal data (eng-premier-league.csv) and betting odds
(E0-*.csv) into a survival analysis dataset matching the methodology of
Nevo & Ritov (2012) "Around the goal".

Seasons used: 2007-08, 2009-10, 2010-11
(E0-25.csv covers 2024-25 but eng-premier-league.csv goal data for that
season is only partial — 96/380 games — so it is excluded here.)

Each game produces up to 2 rows:
  obs_type=1  FirstGoalTime  : time until home team scores first goal,
                               censored by first away goal or minute 90
  obs_type=2  SecondGoalTime : left-truncated at TimeOfFirstGoal (first goal
                               in game regardless of team); tracks next home
                               goal, censored by next away goal or minute 90

Output columns:
  game_id, obs_type, tstart, tstop, event
  ProbWin          -- implied win prob from avg bookmaker odds (home team)
  Season           -- integer season label (for multi-season models)
  Goal             -- 0 for obs_type=1, 1 for obs_type=2  (= I{j=2})
  TimeOfFirstGoal  -- 0 for obs_type=1, actual minute for obs_type=2
  FirstGoalTeam    -- 0 for obs_type=1; 1=home/0=away for obs_type=2
  home_team, away_team, date
"""

import re
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path

BASE = Path(__file__).parent

# ---------------------------------------------------------------------------
# Seasons: (goal-data season string, betting CSV, Season indicator for model)
# ---------------------------------------------------------------------------
SEASONS = [
    ("2007-2008", "E0-3.csv",         0),
    ("2009-2010", "E0-2.csv",         1),
    ("2010-2011", "E0-4.csv",         2),
    ("2011-2012", "E0-1112-E0.csv",   3),
    ("2012-2013", "E0-1213-E0.csv",   4),
    ("2013-2014", "E0-1314-E0.csv",   5),
    ("2014-2015", "E0-1415-E0.csv",   6),
    ("2015-2016", "E0-1516-E0.csv",   7),
]

# ---------------------------------------------------------------------------
# Team name mapping: goal-data names → betting-data names
# We build this dynamically since different E0 files may have slight variations;
# the static map covers all known mismatches between the two sources.
# ---------------------------------------------------------------------------
TEAM_MAP = {
    # 2007-08 / 2009-10 / 2010-11 era
    "Arsenal FC":             "Arsenal",
    "Aston Villa":            "Aston Villa",
    "Birmingham City":        "Birmingham",
    "Blackburn Rovers":       "Blackburn",
    "Blackpool FC":           "Blackpool",
    "Bolton Wanderers":       "Bolton",
    "Burnley FC":             "Burnley",
    "Chelsea FC":             "Chelsea",
    "Derby County":           "Derby",
    "Everton FC":             "Everton",
    "Fulham FC":              "Fulham",
    "Hull City":              "Hull",
    "Liverpool FC":           "Liverpool",
    "Manchester City":        "Man City",
    "Manchester United":      "Man United",
    "Middlesbrough FC":       "Middlesbrough",
    "Newcastle United":       "Newcastle",
    "Portsmouth FC":          "Portsmouth",
    "Reading FC":             "Reading",
    "Sunderland AFC":         "Sunderland",
    "Stoke City":             "Stoke",
    "Tottenham Hotspur":      "Tottenham",
    "West Bromwich Albion":   "West Brom",
    "West Ham United":        "West Ham",
    "Wigan Athletic":         "Wigan",
    "Wolverhampton Wanderers":"Wolves",
    # 2011-16 era additions
    "Cardiff City":           "Cardiff",
    "Norwich City":           "Norwich",
    "Queens Park Rangers":    "QPR",
    "Swansea City":           "Swansea",
    "Watford FC":             "Watford",
    # 2024-25 era (included for completeness)
    "AFC Bournemouth":        "Bournemouth",
    "Brighton & Hove Albion": "Brighton",
    "Brentford FC":           "Brentford",
    "Crystal Palace":         "Crystal Palace",
    "Ipswich Town":           "Ipswich",
    "Leicester City":         "Leicester",
    "Nottingham Forest":      "Nott'm Forest",
    "Southampton FC":         "Southampton",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_minute(t: str) -> int:
    """
    Convert goal time string to integer minute, capping stoppage time.
    '67' → 67,  '45+2' → 45,  '90+6' → 90
    Paper convention: stoppage-time goals recorded as last minute of half.
    """
    m = re.match(r"^(\d+)", str(t).strip())
    if m:
        return min(int(m.group(1)), 90)
    raise ValueError(f"Cannot parse goal time: {t!r}")


def parse_game_string(game: str):
    """'Arsenal FC vs. Chelsea FC 2:1' → ('Arsenal FC', 'Chelsea FC')"""
    without_score = game.rsplit(" ", 1)[0]
    home, away = without_score.split(" vs. ", 1)
    return home.strip(), away.strip()


def parse_bet_date(s: str) -> datetime:
    for fmt in ("%d/%m/%Y", "%d/%m/%y"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            pass
    raise ValueError(f"Cannot parse date: {s!r}")


def make_row(game_id, obs_type, tstart, tstop, event,
             prob_win, season_label, goal, t1, fgt,
             home_std, away_std, date):
    return {
        "game_id":         game_id,
        "obs_type":        obs_type,
        "tstart":          tstart,
        "tstop":           tstop,
        "event":           event,
        "ProbWin":         prob_win,
        "Season":          season_label,
        "Goal":            goal,
        "TimeOfFirstGoal": t1,
        "FirstGoalTeam":   fgt,
        "home_team":       home_std,
        "away_team":       away_std,
        "date":            date,
    }

# ---------------------------------------------------------------------------
# Load goal data (all seasons at once, filter per season in loop)
# ---------------------------------------------------------------------------

all_goals = pd.read_csv(BASE / "eng-premier-league.csv")
all_goals["minute"] = all_goals["time"].apply(parse_minute)
all_goals["date"] = pd.to_datetime(all_goals["date"])

home_raw, away_raw = zip(*all_goals["game"].apply(parse_game_string))
all_goals["home_raw"] = home_raw
all_goals["away_raw"] = away_raw

# ---------------------------------------------------------------------------
# Process each season
# ---------------------------------------------------------------------------

all_rows = []
game_id = 0
season_stats = []

for season_str, bet_file, season_label in SEASONS:
    print(f"\nProcessing {season_str} ({bet_file}) ...")

    goals = all_goals[all_goals["season"] == season_str].copy()

    # Load betting data
    # Older files (pre-2014) use BbAvH/D/A (BetBrain avg — same source as paper).
    # Newer files use AvgH/D/A. Fall back gracefully.
    bet = pd.read_csv(BASE / bet_file)
    bet = bet.dropna(subset=["Date"])
    bet["date"] = pd.to_datetime(bet["Date"].apply(parse_bet_date))
    h_col = "BbAvH" if "BbAvH" in bet.columns else "AvgH"
    d_col = "BbAvD" if "BbAvD" in bet.columns else "AvgD"
    a_col = "BbAvA" if "BbAvA" in bet.columns else "AvgA"
    bet["p_H"] = 1.0 / pd.to_numeric(bet[h_col], errors="coerce")
    bet["p_D"] = 1.0 / pd.to_numeric(bet[d_col], errors="coerce")
    bet["p_A"] = 1.0 / pd.to_numeric(bet[a_col], errors="coerce")
    bet["ProbWin"] = bet["p_H"] / (bet["p_H"] + bet["p_D"] + bet["p_A"])

    bet_lookup = bet.set_index(["date", "HomeTeam"])["ProbWin"].to_dict()
    all_bet_games = bet[["date", "HomeTeam", "AwayTeam", "ProbWin"]].copy()
    all_bet_games = all_bet_games.rename(columns={"HomeTeam": "home_std",
                                                   "AwayTeam": "away_std"})

    # Games with goal data
    goal_game_keys = set(
        zip(goals["date"],
            goals["home_raw"].map(lambda x: TEAM_MAP.get(x, x)))
    )
    bet_game_keys = set(zip(all_bet_games["date"], all_bet_games["home_std"]))
    no_goal_keys  = bet_game_keys - goal_game_keys   # 0-0 draws

    rows = []
    unmatched = 0

    # ── Games with at least one goal ─────────────────────────────────────
    for (game_str, date), grp in goals.groupby(["game", "date"]):
        game_id += 1
        home_raw_ = grp["home_raw"].iloc[0]
        away_raw_ = grp["away_raw"].iloc[0]
        home_std  = TEAM_MAP.get(home_raw_, home_raw_)
        away_std  = TEAM_MAP.get(away_raw_, away_raw_)

        grp = grp.sort_values("minute")
        home_goals = sorted(grp[grp["scoring_team"] == home_raw_]["minute"].tolist())
        away_goals = sorted(grp[grp["scoring_team"] == away_raw_]["minute"].tolist())

        prob_win = bet_lookup.get((date, home_std), np.nan)
        if np.isnan(prob_win):
            unmatched += 1

        fh = home_goals[0] if home_goals else None
        fa = away_goals[0] if away_goals else None

        # Obs 1: FirstGoalTime
        if fh is not None and (fa is None or fh <= fa):
            stop1, ev1 = fh, 1
        elif fa is not None:
            stop1, ev1 = fa, 0
        else:
            stop1, ev1 = 90, 0

        rows.append(make_row(game_id, 1, 0, stop1, ev1,
                             prob_win, season_label, 0, 0, 0,
                             home_std, away_std, date))

        # Obs 2: SecondGoalTime
        all_sorted = sorted(home_goals + away_goals)
        t1 = all_sorted[0]
        if t1 >= 90:
            continue

        home_first = (fh is not None) and (fa is None or fh <= fa)

        if home_first:
            next_h = home_goals[1] if len(home_goals) > 1 else None
            next_a = away_goals[0] if away_goals else None
        else:
            next_h = home_goals[0] if home_goals else None
            next_a = away_goals[1] if len(away_goals) > 1 else None

        INF   = 9999
        nh    = next_h if next_h is not None else INF
        na    = next_a if next_a is not None else INF
        stop2 = min(nh, na, 90)

        if stop2 <= t1:   # degenerate interval (two goals same minute)
            continue

        ev2 = 1 if (next_h is not None and next_h <= 90 and next_h < na) else 0

        rows.append(make_row(game_id, 2, t1, stop2, ev2,
                             prob_win, season_label, 1, t1, int(home_first),
                             home_std, away_std, date))

    # ── 0-0 games (no goals → obs_type=1 censored at 90) ────────────────
    for (date_, home_std_) in no_goal_keys:
        game_id += 1
        row_bet  = all_bet_games[
            (all_bet_games["date"] == date_) &
            (all_bet_games["home_std"] == home_std_)
        ].iloc[0]
        rows.append(make_row(game_id, 1, 0, 90, 0,
                             row_bet["ProbWin"], season_label, 0, 0, 0,
                             home_std_, row_bet["away_std"], date_))

    n1  = sum(r["obs_type"] == 1 for r in rows)
    n2  = sum(r["obs_type"] == 2 for r in rows)
    ev1 = sum(r["event"] == 1 and r["obs_type"] == 1 for r in rows)
    ev2 = sum(r["event"] == 1 and r["obs_type"] == 2 for r in rows)
    print(f"  obs_type=1: {n1} ({n1 - ev1} censored, {ev1} events)")
    print(f"  obs_type=2: {n2} ({n2 - ev2} censored, {ev2} events)")
    print(f"  Unmatched ProbWin: {unmatched}")
    season_stats.append((season_str, n1, n2, ev1, ev2))
    all_rows.extend(rows)

# ---------------------------------------------------------------------------
# Assemble and save
# ---------------------------------------------------------------------------

df = pd.DataFrame(all_rows).sort_values(["game_id", "obs_type"]).reset_index(drop=True)

n_games    = df["game_id"].nunique()
n_obs      = len(df)
n_obs1     = (df["obs_type"] == 1).sum()
n_obs2     = (df["obs_type"] == 2).sum()
ev1_all    = df[(df["obs_type"] == 1) & (df["event"] == 1)].shape[0]
ev2_all    = df[(df["obs_type"] == 2) & (df["event"] == 1)].shape[0]
missing_pw = df["ProbWin"].isna().sum()
mean_pw    = df["ProbWin"].mean()
mean_t1    = df[df["obs_type"] == 2]["TimeOfFirstGoal"].mean()
pct_home_first = df[(df["obs_type"] == 2)]["FirstGoalTeam"].mean()

print("\n" + "=" * 55)
print("  Combined survival dataset summary")
print("=" * 55)
print(f"  Seasons:                 {', '.join(s[0] for s in SEASONS)}")
print(f"  Games (unique):          {n_games}  (paper: 760)")
print(f"  Total observations:      {n_obs}   (paper: 1433)")
print(f"  obs_type=1 (FirstGoal):  {n_obs1}  ({n_obs1 - ev1_all} censored, {ev1_all} events)")
print(f"  obs_type=2 (SecondGoal): {n_obs2}  ({n_obs2 - ev2_all} censored, {ev2_all} events)")
print(f"  Missing ProbWin:         {missing_pw}")
print(f"  Mean ProbWin:            {mean_pw:.3f}  (paper: 0.448)")
print(f"  Mean TimeOfFirstGoal:    {mean_t1:.1f}  (paper: 30.8)")
print(f"  % first goals by home:   {pct_home_first*100:.1f}%  (paper: 60%)")
print("=" * 55)

df.to_csv(BASE / "survival_data.csv", index=False)
print("  Saved → survival_data.csv")
