"""Intervals between dosing and the two MRI sessions (Methods, Participants).

    python src/scan_intervals.py                 # spreadsheet from config
    python src/scan_intervals.py path/to.xlsx    # or give it explicitly

Reports, overall and by arm, the days from MRI 1 to dosing, dosing to MRI 2,
and MRI 1 to MRI 2, and flags implausible orderings.

Outputs
  results/tables/scan_intervals_summary.csv   aggregate only, safe to commit
  <spreadsheet dir>/scan_intervals_per_participant.csv
                                              per-participant dates, NEVER
                                              written inside the repository
"""
import sys
from pathlib import Path

import pandas as pd

try:
    from scipy.stats import mannwhitneyu
except ImportError:            # descriptives still work without scipy
    mannwhitneyu = None

REPO = Path(__file__).resolve().parent.parent
DATE_COLS = ["Dosingdate", "Incldate", "MRI_date_1", "MRI_date_2"]
INTERVALS = [("mri1_to_dose", "MRI 1 -> dosing (days before dose)"),
             ("dose_to_mri2", "Dosing -> MRI 2 (days after dose)"),
             ("mri1_to_mri2", "MRI 1 -> MRI 2 (days between scans)")]


def parse_mixed(s):
    """Parse a column containing both YYYY-MM-DD and DD/MM/YYYY.

    Slash-formatted dates in this dataset are day-first (values >12 in the
    first field confirm this). Parsing each format explicitly avoids silent
    month/day transposition, which would be invisible for days 1-12.
    """
    s = s.astype(str).str.strip()
    iso = pd.to_datetime(s, format="%Y-%m-%d", errors="coerce")
    eur = pd.to_datetime(s, format="%d/%m/%Y", errors="coerce")
    out = iso.fillna(eur)
    if out.isna().any():
        raise ValueError(f"Unparsed dates in {s.name}: {s[out.isna()].tolist()}")
    return out


def describe(x):
    x = x.dropna()
    return dict(n=len(x), median=x.median(), q1=x.quantile(.25),
                q3=x.quantile(.75), min=x.min(), max=x.max())


def fmt(d, label):
    return (f"{label}: median {d['median']:.0f}, IQR {d['q1']:.0f}-{d['q3']:.0f}, "
            f"range {d['min']:.0f}-{d['max']:.0f}, n={d['n']}")


def audit(df):
    df = df.copy()
    for c in DATE_COLS:
        df[c] = parse_mixed(df[c])

    df["mri1_to_dose"] = (df["Dosingdate"] - df["MRI_date_1"]).dt.days  # + = pre-dose
    df["dose_to_mri2"] = (df["MRI_date_2"] - df["Dosingdate"]).dt.days  # + = post-dose
    df["mri1_to_mri2"] = (df["MRI_date_2"] - df["MRI_date_1"]).dt.days

    flags = []
    if (df["mri1_to_dose"] < 0).any():
        flags.append(f"MRI 1 AFTER dosing: rows {df.index[df['mri1_to_dose'] < 0].tolist()}")
    if (df["dose_to_mri2"] <= 0).any():
        flags.append(f"MRI 2 NOT after dosing: rows {df.index[df['dose_to_mri2'] <= 0].tolist()}")
    if (df["Incldate"] > df["Dosingdate"]).any():
        flags.append(f"Inclusion AFTER dosing: rows {df.index[df['Incldate'] > df['Dosingdate']].tolist()}")

    rows = []
    print("=== Overall ===")
    for col, lab in INTERVALS:
        d = describe(df[col]); print(" ", fmt(d, lab))
        rows.append(dict(interval=col, arm="All", **d, U=None, p=None))

    if "Treatment" in df.columns:
        print("\n=== By arm ===")
        for col, lab in INTERVALS:
            print(f"\n {col}")
            groups = []
            for arm, g in df.groupby("Treatment"):
                d = describe(g[col]); print("  ", fmt(d, str(arm)))
                rows.append(dict(interval=col, arm=str(arm), **d, U=None, p=None))
                groups.append(g[col].dropna().values)
            if len(groups) == 2 and mannwhitneyu is not None:
                u, p = mannwhitneyu(groups[0], groups[1], alternative="two-sided")
                print(f"   Mann-Whitney U = {u:.0f}, p = {p:.3f}")
                rows.append(dict(interval=col, arm="between-arm test",
                                 n=None, median=None, q1=None, q3=None,
                                 min=None, max=None, U=u, p=p))

    if flags:
        print("\n=== FLAGS ===")
        for f in flags:
            print("  !", f)

    return df, pd.DataFrame(rows)


def main():
    sys.path.insert(0, str(REPO))
    import config as cfg

    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(cfg.SCAN_DATES_XLSX)
    if not src.exists():
        sys.exit(f"spreadsheet not found: {src}\nset SCAN_DATES_XLSX in config_local.py")

    per, summary = audit(pd.read_excel(src))

    out_summary = Path(cfg.SAVETABLE) / "scan_intervals_summary.csv"
    summary.round(3).to_csv(out_summary, index=False)

    out_private = src.parent / "scan_intervals_per_participant.csv"
    if REPO in out_private.resolve().parents:
        sys.exit("refusing to write per-participant dates inside the repository")
    per.to_csv(out_private, index=False)

    print(f"\nWrote {out_summary}  (aggregate)")
    print(f"Wrote {out_private}  (per participant, private)")


if __name__ == "__main__":
    main()
