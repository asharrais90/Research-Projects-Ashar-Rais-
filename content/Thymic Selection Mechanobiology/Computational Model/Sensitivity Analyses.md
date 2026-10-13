"""

sensitivity_analysis.py

Requires: Thymic_selection_mechanobiology.py in the same directory.

Run:  python sensitivity_analysis.py

      python sensitivity_analysis.py --sa 1 2 3 4 5 6 7 8 9

"""

  

from __future__ import annotations

import argparse, os, sys, time

from typing import Any, Dict, List, Tuple

import numpy as np

import matplotlib.pyplot as plt

  

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Thymic_selection_mechanobiology_model import (

    Params, calibrate_thresholds, ThymicEnvironment, Thymocyte

)

  

# ── Output directory ──────────────────────────────────────────────────

FIG_DIR = r"C:/Academics/Coding projects/Research projects/Thymic mechanobiology/Sensitivity Analysis figures"

os.makedirs(FIG_DIR, exist_ok=True)

  

# ── Plot style ────────────────────────────────────────────────────────

C_POS = "#378ADD"; C_NEG = "#E24B4A"; C_NGL = "#1D9E75"; C_SEC = "#888780"

plt.rcParams.update({

    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,

    "legend.fontsize": 8, "figure.dpi": 150, "savefig.bbox": "tight",

})

  
  

# ══════════════════════════════════════════════════════════════════════

# SHARED HELPERS

# ══════════════════════════════════════════════════════════════════════

  

# SA run defaults — smaller than full model, large enough for

# ±3 percentage-point resolution on fate rates.

# NOTE: grid 100×100, n_initial=500 gives ~5% density, similar to

# the full 100×100 / n_initial=1000 density at 10%.

# Increase n_initial to 1000 if you need tighter resolution.

_SA_DEFAULTS: Dict[str, Any] = dict(

    n_initial  = 1000,

    grid_w     = 100,

    grid_h     = 100,

    T_run_hr   = 288.0,  

)

  
  

DEFAULT_SEEDS: Tuple[int, ...] = (42, 123, 456, 789, 1024)

  
  

def _make_params(

    scenario: str = "uniform",

    seed: int = 42,

    recalibrate: bool = True,

    frozen_thresholds: Tuple[float, float] | None = None,

    **overrides,

) -> Params:

    """Build a Params object for a sensitivity run.

  

    Parameters are applied in this order:

      1. model defaults,

      2. SA runtime defaults,

      3. caller overrides.

  

    recalibrate=True reproduces the original SA behaviour.

    recalibrate=False requires explicit frozen_thresholds so that a

    parameter perturbation does not silently move the selection boundaries.

    """

    p = Params(scenario=scenario, seed=seed)

    for k, v in _SA_DEFAULTS.items():

        setattr(p, k, v)

    for k, v in overrides.items():

        if not hasattr(p, k):

            raise AttributeError(f"Params has no attribute '{k}'")

        setattr(p, k, v)

  

    if frozen_thresholds is not None:

        p.theta_pos, p.theta_neg = frozen_thresholds

    elif recalibrate:

        p.theta_pos, p.theta_neg = calibrate_thresholds(p)

    else:

        raise ValueError(

            "recalibrate=False requires explicit frozen_thresholds."

        )

  

    return p

  
  

def _reference_thresholds(scenario: str = "uniform",) -> Tuple[float, float]:

    """Compute the baseline thresholds once for an OAT sensitivity group."""

    p_ref = _make_params(scenario=scenario, recalibrate=True)

    return float(p_ref.theta_pos), float(p_ref.theta_neg)

  
  

def _resolve_seeds(

    seed: int | None = None,

    seeds: Tuple[int, ...] = DEFAULT_SEEDS,

) -> Tuple[int, ...]:

    """Keep backward compatibility with the old single-seed API."""

    return (int(seed),) if seed is not None else tuple(int(s) for s in seeds)

  
  

def _run_silent(p: Params) -> ThymicEnvironment:

    """Run the full simulation with no console output."""

    n_steps = int(p.T_run_hr * 60.0 / p.dt)

    model   = ThymicEnvironment(p)

    for _ in range(n_steps):

        model.step()

        # ── Early stop: no point running empty timesteps ──────────────

        if not any(

            isinstance(a, Thymocyte) and a.alive

            for a in model.agents

        ):

            break

    return model

  
  

def _extract(model: ThymicEnvironment) -> Dict[str, Any]:

    """Pull outcome metrics from a completed run into a flat dict."""

    fc    = model.fate_counts

    total = max(sum(fc.values()), 1)

    pos   = fc["positive_selection"]

    neg   = fc["negative_selection"]

    ngl   = fc["death_by_neglect"]

  

    r: Dict[str, Any] = dict(

        pos_pct        = 100.0 * pos / total,

        neg_pct        = 100.0 * neg / total,

        neglect_pct    = 100.0 * ngl / total,

        neg_pos_ratio  = neg / max(pos, 1),

        pos_n=pos, neg_n=neg, neglect_n=ngl,

    )

  

    cs = getattr(model, "contact_stats", {})

    for fate_key, sh in [

        ("positive_selection", "pos"),

        ("negative_selection", "neg"),

        ("death_by_neglect",   "ngl"),

    ]:

        recs = cs.get(fate_key, [])

        if recs:

            r[f"{sh}_cfrac"]  = float(np.mean([x["contact_frac"]  for x in recs]))

            r[f"{sh}_csteps"] = float(np.mean([x["contact_steps"] for x in recs]))

            r[f"{sh}_meanS"]  = float(np.mean([x["mean_S"]        for x in recs]))

        else:

            r[f"{sh}_cfrac"] = r[f"{sh}_csteps"] = r[f"{sh}_meanS"] = np.nan

    return r

  
  

def _ms(rows: List[Dict], key: str) -> Tuple[float, float]:

    vals = [r[key] for r in rows if not np.isnan(r.get(key, np.nan))]

    if len(vals) > 1:

        return float(np.mean(vals)), float(np.std(vals, ddof=1))

    if len(vals) == 1:

        return float(vals[0]), 0.0

    return np.nan, np.nan

  
  

def _range(rows: List[Dict], key: str) -> float:

    vals = [r[key] for r in rows if not np.isnan(r.get(key, np.nan))]

    return float(max(vals) - min(vals)) if vals else np.nan

  
  

def _sensitivity_flag(rows: List[Dict], key: str) -> Tuple[float, str]:

    """

    Compare the range across parameter values to the mean within-condition SD.

    Returns (range, annotation_string).

    Falls back gracefully when SD data is absent (single-seed runs).

    """

    param_range = _range(rows, key)

  

    sd_vals = [r.get(f"{key}_sd", np.nan) for r in rows]

    sd_vals  = [v for v in sd_vals if not np.isnan(v)]

  

    if not sd_vals:

        return param_range, "[no SD]"

  

    mean_sd = float(np.mean(sd_vals))

    if mean_sd == 0.0:

        return param_range, "[SD=0]"

  

    snr  = param_range / mean_sd

    flag = f"SNR={snr:.1f}  ◄ {'REAL' if snr > 2.0 else 'NOISE'}"

    return param_range, flag

  
  

# ══════════════════════════════════════════════════════════════════════

# SA-1: STOCHASTIC REPRODUCIBILITY

# ══════════════════════════════════════════════════════════════════════

  

def sa1_seed_reproducibility(

    seeds: Tuple[int, ...] = DEFAULT_SEEDS,

    scenario: str = "uniform",

) -> Dict:

    """

    Five independent seeds, all other parameters fixed.

    Pass criterion: SD < 2 percentage points per fate class.

    """

    print(f"\n{'─'*60}\nSA-1  Seed reproducibility  |  scenario={scenario}\n{'─'*60}")

    rows = []

    for s in seeds:

        p  = _make_params(scenario=scenario, seed=s)

        t0 = time.time()

        m  = _run_silent(p)

        r  = _extract(m); r["seed"] = s

        rows.append(r)

        print(f"  seed={s:5d}  pos={r['pos_pct']:5.1f}%  neg={r['neg_pct']:5.1f}%"

              f"  neglect={r['neglect_pct']:5.1f}%  ({time.time()-t0:.0f}s)")

  

    all_pass = True

    print("\n  Summary:")

    for key, label in [("pos_pct","Positive"),("neg_pct","Negative"),("neglect_pct","Neglect")]:

        mu, sd = _ms(rows, key)

        flag   = "✓ PASS" if sd < 2.0 else "✗ WARN — increase n_initial"

        all_pass = all_pass and (sd < 2.0)

        print(f"    {label:12s}  mean={mu:5.1f}%  SD={sd:4.1f}%  {flag}")

    if not all_pass:

        print("  → Increase n_initial before trusting SA-2–6.")

  

    return dict(seeds=list(seeds), rows=rows, scenario=scenario)

  
  

def plot_sa1(res: Dict, save: bool = True, show: bool = False) -> None:

    seeds = res["seeds"]; rows = res["rows"]

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))

    fig.suptitle("SA-1 — Stochastic Reproducibility (seed sweep)", fontweight="bold")

  

    for ax, key, ylabel, color, fate in zip(

        axes,

        ["pos_pct", "neg_pct", "neglect_pct"],

        ["Positive selection (%)", "Negative selection (%)", "Death by neglect (%)"],

        [C_POS, C_NEG, C_NGL],

        ["pos", "neg", None],

    ):

        vals = [r[key] for r in rows]

        mu   = np.mean(vals); sd = np.std(vals, ddof=1) if len(vals) > 1 else 0.0

        ax.plot(range(len(seeds)), vals, "o-", color=color, lw=1.5, ms=7, zorder=3)

        ax.axhline(mu, ls="--", color=color, alpha=0.7, lw=1.2, label=f"mean={mu:.1f}%")

        ax.fill_between(range(len(seeds)), mu-sd, mu+sd,

                        alpha=0.18, color=color, label=f"±1 SD ({sd:.1f}%)")

        ax.set_xticks(range(len(seeds)))

        ax.set_xticklabels(seeds, rotation=45, ha="right", fontsize=8)

        ax.set_xlabel("Random seed"); ax.set_ylabel(ylabel)

        ax.set_title(ylabel.split(" (")[0]); ax.legend(fontsize=7)

  

    plt.tight_layout()

    if save:

        path = os.path.join(FIG_DIR, "SA1_seed_reproducibility.png")

        fig.savefig(path)

        print(f"  → saved: {path}")

    if show:

        plt.show()

    plt.close(fig)

  
  

# ══════════════════════════════════════════════════════════════════════

# SHARED OAT RUNNER

# ══════════════════════════════════════════════════════════════════════

  

def _run_1d_param_group(

    parameter_values: Dict[str, Tuple[float, ...]],

    title: str,

    scenario: str = "uniform",

    seed: int | None = None,

    seeds: Tuple[int, ...] = DEFAULT_SEEDS,

    recalibrate: bool = False,

) -> Dict:

    """Run independent OFAT sweeps with the same seeds at every value."""

    run_seeds = _resolve_seeds(seed=seed, seeds=seeds)

    frozen_thresholds = None

    if not recalibrate:

        frozen_thresholds = _reference_thresholds(

            scenario=scenario,

        )

  

    print(f"\n{'─'*60}\n{title}  |  scenario={scenario}\n{'─'*60}")

    print(f"  seeds = {run_seeds}")

    if frozen_thresholds is not None:

        print(f"  frozen thresholds = theta_pos={frozen_thresholds[0]:.4g}, "

              f"theta_neg={frozen_thresholds[1]:.4g}")

  

    rows_by_param: Dict[str, List[Dict]] = {}

    raw_rows_by_param: Dict[str, Dict[float, List[Dict]]] = {}

  

    metrics = (

        "pos_pct", "neg_pct", "neglect_pct", "neg_pos_ratio",

        "pos_cfrac", "neg_cfrac", "ngl_cfrac",

        "pos_csteps", "neg_csteps", "ngl_csteps",

        "pos_meanS", "neg_meanS", "ngl_meanS",

    )

  

    for name, values in parameter_values.items():

        values = tuple(values)

        if not values:

            raise ValueError(f"No values supplied for parameter '{name}'")

        if not hasattr(Params(), name):

            raise AttributeError(f"Params has no attribute '{name}'")

  

        default = getattr(Params(), name)

        rows: List[Dict] = []

        raw_rows_by_param[name] = {}

        print(f"  Sweep: {name}  (default={default})")

  

        for value in values:

            seed_rows: List[Dict] = []

            for s in run_seeds:

                p = _make_params(

                    scenario=scenario,

                    seed=s,

                    recalibrate=recalibrate,

                    frozen_thresholds=frozen_thresholds,

                    **{name: value},

                )

                t0 = time.time()

                m = _run_silent(p)

                r = _extract(m)

                r[name] = value

                r["seed"] = s

                seed_rows.append(r)

                print(

                    f"    {name}={value:g}  seed={s:5d}  "

                    f"pos={r['pos_pct']:5.1f}%  neg={r['neg_pct']:5.1f}%  "

                    f"({time.time()-t0:.0f}s)"

                )

  

            row = {name: value}

            for metric in metrics:

                vals = np.asarray([r.get(metric, np.nan) for r in seed_rows], dtype=float)

                finite = vals[np.isfinite(vals)]

                if finite.size:

                    row[metric] = float(np.mean(finite))

                    row[f"{metric}_sd"] = (

                        float(np.std(finite, ddof=1)) if finite.size > 1 else 0.0

                    )

                else:

                    row[metric] = np.nan

                    row[f"{metric}_sd"] = np.nan

  

            rows.append(row)

            raw_rows_by_param[name][value] = seed_rows

            print(

                f"      mean ± SD: pos={row['pos_pct']:5.1f} ± {row['pos_pct_sd']:.1f}%  "

                f"neg={row['neg_pct']:5.1f} ± {row['neg_pct_sd']:.1f}%"

            )

  

        rows_by_param[name] = rows

  

    return {

        "rows": rows_by_param,

        "raw_rows": raw_rows_by_param,

        "parameter_values": {k: list(v) for k, v in parameter_values.items()},

        "scenario": scenario,

        "seeds": list(run_seeds),

        "thresholds_frozen": not recalibrate,

        "frozen_thresholds": frozen_thresholds,

    }

  
  

def _plot_1d_param_group(

    res: Dict,

    title: str,

    filename: str,

    save: bool = True,

    show: bool = False,

    log_x_params: Tuple[str, ...] = (),   # param names that warrant a log x-axis

) -> None:

    """

    Generic plot for an OFAT parameter group.

    Changes from previous version:

    - neg/pos ratio secondary axis removed entirely

    - 4-parameter figures use 2×2 grid instead of 3+1

    - explicit per-panel xlabel prevents cross-row label bleed

    - optional log x-scale for parameters spanning > 1 order of magnitude

    """

    names = list(res["rows"].keys())

    n     = len(names)

  

    # ── Layout: 2×2 for n=4, otherwise ≤3 cols ──────────────────────

    if n == 4:

        ncols, nrows = 2, 2

    elif n <= 3:

        ncols, nrows = n, 1

    else:

        ncols = 3

        nrows = int(np.ceil(n / ncols))

  

    fig, axes = plt.subplots(

        nrows, ncols,

        figsize=(4.8 * ncols, 5.0 * nrows),

        squeeze=False,

    )

    axes_flat = axes.ravel()

    fig.suptitle(title, fontweight="bold", y=1.01)

  

    for ax, name in zip(axes_flat, names):

        rows    = res["rows"][name]

        x       = [r[name]        for r in rows]

        pos     = [r["pos_pct"]   for r in rows]

        neg     = [r["neg_pct"]   for r in rows]

        pos_sd  = [r.get("pos_pct_sd", 0.0) for r in rows]

        neg_sd  = [r.get("neg_pct_sd", 0.0) for r in rows]

  

        # ── Main lines with SD bands ──────────────────────────────────

        ax.plot(x, pos, "o-", color=C_POS, lw=1.5, ms=6, label="Positive sel.")

        ax.plot(x, neg, "s-", color=C_NEG, lw=1.5, ms=6, label="Negative sel.")

        ax.fill_between(x,

                        np.array(pos) - np.array(pos_sd),

                        np.array(pos) + np.array(pos_sd),

                        color=C_POS, alpha=0.12)

        ax.fill_between(x,

                        np.array(neg) - np.array(neg_sd),

                        np.array(neg) + np.array(neg_sd),

                        color=C_NEG, alpha=0.12)

  

        # ── Reference markers ─────────────────────────────────────────

        default = getattr(Params(), name)

        ax.axvline(default, ls=":", color="black", alpha=0.5, label=f"Default ({default})")

  

        # ── Axes ──────────────────────────────────────────────────────

        if name in log_x_params:

            ax.set_xscale("log")

        ax.set_xlabel(name)          # explicit per-panel — prevents cross-row bleed

        ax.set_ylabel("Selection rate (%)")

        ax.set_title(name)

        ax.legend(fontsize=7, loc="best")

  

    # ── Hide surplus panels ───────────────────────────────────────────

    for ax in axes_flat[n:]:

        ax.set_visible(False)

  

    plt.tight_layout(h_pad=3.5, w_pad=2.5)   # h_pad prevents row-label overlap

    if save:

        path = os.path.join(FIG_DIR, filename)

        fig.savefig(path)

        print(f"  → saved: {path}")

    if show:

        plt.show()

    plt.close(fig)

  
  

# ══════════════════════════════════════════════════════════════════════

# SA-2: MECHANOBIOLOGICAL AMPLIFICATION

# Standardized as two independent OFAT sweeps.

# Thresholds are frozen at the baseline calibration.

# ══════════════════════════════════════════════════════════════════════

  

def sa2_mechanobiological_amplification(

    gamma_vals:   Tuple[float, ...] = (0.0, 0.25, 0.50, 0.75, 1.00),

    alpha_C_vals: Tuple[float, ...] = (0.0, 0.25, 0.50, 0.75, 1.00),

    alpha_M_vals: Tuple[float, ...] = (0.0, 0.25, 0.50, 0.75, 1.00),

    scenario: str = "uniform",

    seed: int | None = None,

    seeds: Tuple[int, ...] = DEFAULT_SEEDS,

) -> Dict:

    """Sweep gamma, alpha_C_max and alpha_M_max independently using the shared OAT runner."""

    res = _run_1d_param_group(

        {"gamma": tuple(gamma_vals),

         "alpha_C_max": tuple(alpha_C_vals),

         "alpha_M_max": tuple(alpha_M_vals)},

        "SA-2  Mechanobiological amplification",

        scenario=scenario, seed=seed, seeds=seeds, recalibrate=False,

    )

    return dict(

        gamma_vals=list(gamma_vals),

        gamma_rows=res["rows"]["gamma"],

        alpha_C_vals=list(alpha_C_vals),

        alpha_C_rows=res["rows"]["alpha_C_max"],

        alpha_M_vals=list(alpha_M_vals),

        alpha_M_rows=res["rows"]["alpha_M_max"],

        scenario=scenario,

        seeds=res["seeds"],

        thresholds_frozen=True,

        frozen_thresholds=res["frozen_thresholds"],

    )

  
  

def plot_sa2(

    res: Dict,

    save: bool = True,

    show: bool = False,

) -> None:

  

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

  

    fig.suptitle(

        "SA-2 — Mechanobiological Amplification Parameters",

        fontweight="bold"

    )

  

    sweeps = [

        (

            res["gamma_vals"],

            res["gamma_rows"],

            "gamma",

            "γ (scanning capacity gain)",

            "Scanning capacity amplification",

        ),

        (

            res["alpha_C_vals"],

            res["alpha_M_rows"],

            "alpha_C_max",

            "α_C,max (cytoskeletal amplification gain)",

            "C-channel signal amplification",

        ),

        (

            res["alpha_M_vals"],

            res["alpha_M_rows"],

            "alpha_M_max",

            "α_M,max (mechanotransduction amplification gain)",

            "M-channel signal amplification",

        ),

    ]

  

    for ax, (x_vals, rows, x_key, xlabel, title) in zip(axes, sweeps):

        pos    = [r["pos_pct"] for r in rows]

        neg    = [r["neg_pct"] for r in rows]

        pos_sd = [r.get("pos_pct_sd", 0.0) for r in rows]

        neg_sd = [r.get("neg_pct_sd", 0.0) for r in rows]

  

        ax.plot(x_vals, pos, "o-", color=C_POS, lw=1.5, ms=7, label="Positive sel.")

        ax.plot(x_vals, neg, "s-", color=C_NEG, lw=1.5, ms=7, label="Negative sel.")

        ax.fill_between(x_vals,

                        np.array(pos) - np.array(pos_sd),

                        np.array(pos) + np.array(pos_sd),

                        color=C_POS, alpha=0.12)

        ax.fill_between(x_vals,

                        np.array(neg) - np.array(neg_sd),

                        np.array(neg) + np.array(neg_sd),

                        color=C_NEG, alpha=0.12)

  

        ax.axvline(getattr(Params(), x_key), ls=":", color="black",

                   alpha=0.5, label="Default")

  

        ax.set_xlabel(xlabel)

        ax.set_ylabel("Selection rate (%)")

        ax.set_title(title)

        ax.legend(fontsize=7, loc="best")

  

    plt.tight_layout()

    if save:

        path = os.path.join(

            FIG_DIR,"SA2_mechanobiological_amplification.png"

        )

        fig.savefig(path)

        print(f"  → saved: {path}")

    if show:

        plt.show()

  

    plt.close(fig)

  
  

# ══════════════════════════════════════════════════════════════════════

# SA-3: koff DISTRIBUTION SENSITIVITY

# ══════════════════════════════════════════════════════════════════════

  

def sa3_koff_distribution(

    # These ranges are set relative to the updated reference ligand

    # values: koff_delete=0.057, koff_select=0.36, koff_null=5.0

    # The lognormal mu parameter shifts median affinity of the population.

    # koff_mu = 1.0 gives median koff ~ e^1 = 2.72 s⁻¹ (between select and null)

    # Vary around this to test repertoire assumption sensitivity.

    shift_vals: Tuple[float, ...] = (0.1, 0.2, 0.3, 0.5, 1.0),

    sigma_vals: Tuple[float, ...] = (3.0, 6.0, 9.0),

    scenario: str = "uniform",

    seed: int | None = None,

    seeds: Tuple[int,...] = DEFAULT_SEEDS,

) -> Dict:

    """

    Sweep koff_shift_factor and koff_sigma_factor independently.

    shift_factor < 1 moves population toward stronger binders (left).

    sigma_factor controls repertoire breadth.

    Both parameters now live in Params, so the standard OAT runner

    handles seeding, SD computation and frozen thresholds correctly.

    """

    res = _run_1d_param_group(

        {

            "koff_shift_factor": tuple(shift_vals),

            "koff_sigma_factor": tuple(sigma_vals),

        },

        "SA-3  koff distribution",

        scenario=scenario,

        seed=seed,                                      # ← seed used here

        seeds=seeds,                                    # ← seeds used here

        recalibrate=False,

    )

  

    # ── 2. pull the two result tables out of the runner's return ──────

    shift_rows = res["rows"]["koff_shift_factor"]       # ← DEFINES shift_rows

    sigma_rows = res["rows"]["koff_sigma_factor"]       # ← DEFINES sigma_rows

  

    for r in shift_rows:

        r["koff_mu"] = r["koff_shift_factor"]

    for r in sigma_rows:

        r["koff_mu"] = r["koff_sigma_factor"]

  
  

    return dict(

        shift_vals        = list(shift_vals),

        shift_rows        = shift_rows,

        sigma_vals        = list(sigma_vals),

        sigma_rows        = sigma_rows,

        scenario          = scenario,

        seeds             = res["seeds"],

        thresholds_frozen = True,

        frozen_thresholds = res["frozen_thresholds"],

    )

  

def plot_sa3(res: Dict, save: bool = True, show: bool = False) -> None:

    """

    SA-3: two separate figures, one per swept parameter.

    Combining them on a shared x-axis was incorrect because

    shift_factor and sigma_factor are different quantities with

    different units and different interpretations.

    """

    _plot_sa3_shift(res, save=save, show=show)

    _plot_sa3_sigma(res, save=save, show=show)

  

def _plot_sa3_shift(res: Dict, save: bool = True, show: bool = False) -> None:

    """Panel A: shift_factor sweep — population median location."""

    rows  = res["shift_rows"]

    x_raw = res["shift_vals"]

  

    # Convert shift_factor to resulting median koff (s⁻¹) for biological interpretability

    # median koff = exp(mu_log); mu_log changes with shift_factor.

    # We can compute this from first principles or just label shift_factor directly.

    # Using shift_factor on x-axis is safer; caption explains the conversion.

  

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    fig.suptitle(

        "SA-3A — TCR Affinity Distribution: Population Median (shift_factor)",

        fontweight="bold",

    )

  

    for ax, key, ylabel, fate in zip(

        axes,

        ["pos_pct",               "neg_pct"],

        ["Positive selection (%)", "Negative selection (%)"],

        ["pos",                   "neg"],

    ):

  

        vals = [r[key] for r in rows]

        sds  = [r.get(f"{key}_sd", 0.0) for r in rows]

  

        ax.plot(x_raw, vals, "o-", color=C_POS if fate=="pos" else C_NEG,

                lw=1.5, ms=7)

        ax.fill_between(x_raw,

                        np.array(vals) - np.array(sds),

                        np.array(vals) + np.array(sds),

                        color=C_POS if fate=="pos" else C_NEG, alpha=0.12)

  

        ax.axvline(Params().koff_shift_factor, ls=":", color="black",

                   alpha=0.5, label=f"Default ({Params().koff_shift_factor})")

        ax.set_xlabel("koff_shift_factor  (< 1 → stronger binders, > 1 → weaker)")

        ax.set_ylabel(ylabel)

        ax.set_title(ylabel.split(" (")[0])

        ax.legend(fontsize=7)

  

    plt.tight_layout()

    if save:

        path = os.path.join(FIG_DIR, "SA3A_koff_shift_factor.png")

        fig.savefig(path); print(f"  → saved: {path}")

    if show:

        plt.show()

    plt.close(fig)

  

def _plot_sa3_sigma(res: Dict, save: bool = True, show: bool = False) -> None:

    """Panel B: sigma_factor sweep — repertoire breadth."""

    rows  = res["sigma_rows"]

    x_raw = res["sigma_vals"]

  

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    fig.suptitle(

        "SA-3B — TCR Affinity Distribution: Repertoire Breadth (sigma_factor)",

        fontweight="bold",

    )

  

    for ax, key, ylabel, fate in zip(

        axes,

        ["pos_pct",               "neg_pct"],

        ["Positive selection (%)", "Negative selection (%)"],

        ["pos",                   "neg"],

    ):

        vals = [r[key]              for r in rows]

        sds  = [r.get(f"{key}_sd", 0.0) for r in rows]

  

        ax.plot(x_raw, vals, "o-", color=C_POS if fate=="pos" else C_NEG,

                lw=1.5, ms=7)

        ax.fill_between(x_raw,

                        np.array(vals) - np.array(sds),

                        np.array(vals) + np.array(sds),

                        color=C_POS if fate=="pos" else C_NEG, alpha=0.12)

  

        ax.axvline(Params().koff_sigma_factor, ls=":", color="black",

                   alpha=0.5, label=f"Default ({Params().koff_sigma_factor})")

        ax.set_xlabel("koff_sigma_factor  (larger → wider distribution)")

        ax.set_ylabel(ylabel)

        ax.set_title(ylabel.split(" (")[0])

        ax.legend(fontsize=7)

  

    plt.tight_layout()

    if save:

        path = os.path.join(FIG_DIR, "SA3B_koff_sigma_factor.png")

        fig.savefig(path); print(f"  → saved: {path}")

    if show:

        plt.show()

    plt.close(fig)

  

# ══════════════════════════════════════════════════════════════════════

# SA-4: STROMAL DENSITY

# Standardized as an OFAT sweep using the shared runner.

# Thresholds are frozen at the baseline calibration.

# ══════════════════════════════════════════════════════════════════════

  

def sa4_stromal_density(

    density_vals: Tuple[float, ...] = (0.0005, 0.001, 0.002, 0.005, 0.01),

    scenario: str = "uniform",

    seed: int | None = None,

    seeds: Tuple[int, ...] = DEFAULT_SEEDS,

) -> Dict:

    """Sweep stromal density using the shared OFAT framework."""

    res = _run_1d_param_group(

        {"stromal_density": tuple(density_vals)},

        "SA-4  Stromal density",

        scenario=scenario, seed=seed, seeds=seeds, recalibrate=False,

    )

    return dict(

        density_vals=list(density_vals),

        rows=res["rows"]["stromal_density"],

        scenario=scenario,

        seeds=res["seeds"],

        thresholds_frozen=True,

        frozen_thresholds=res["frozen_thresholds"],

    )

  
  

def plot_sa4(res: Dict, save: bool = True, show: bool = False) -> None:

    rows = res["rows"]; ds = res["density_vals"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    fig.suptitle("SA-4 — Stromal Density Sensitivity", fontweight="bold")

  

    ax = axes[0]

    for key, color, marker, label in [

        ("pos_pct",     C_POS, "o", "Positive sel."),

        ("neg_pct",     C_NEG, "s", "Negative sel."),

        ("neglect_pct", C_NGL, "^", "Neglect"),

    ]:

        vals = [r[key]              for r in rows]

        sds  = [r.get(f"{key}_sd", 0.0) for r in rows]

        ax.plot(ds, vals, f"{marker}-", color=color, lw=1.5, ms=7, label=label)

        ax.fill_between(ds,

                        np.array(vals) - np.array(sds),

                        np.array(vals) + np.array(sds),

                        color=color, alpha=0.12)

  

    ax.axvline(Params().stromal_density, ls=":", color="black",

               alpha=0.5, label=f"Default ({Params().stromal_density})")

    ax.set_xlabel("Stromal density (fraction of grid)")

    ax.set_ylabel("Selection rate (%)")

    ax.set_title("Fate rates vs stromal density")

    ax.legend(fontsize=7)

  

    ax = axes[1]

    for key, color, label in [("pos_cfrac",C_POS,"Positive"),

                              ("neg_cfrac",C_NEG,"Negative"),

                              ("ngl_cfrac",C_NGL,"Neglect")]:

        ax.plot(ds, [r.get(key,np.nan) for r in rows], "o-",

                color=color, lw=1.5, ms=7, label=label)

        sd_key = f"{key}_sd"

        if sd_key in rows[0]:

            vals = np.asarray([r.get(key,np.nan) for r in rows], dtype=float)

            sds = np.asarray([r.get(sd_key,np.nan) for r in rows], dtype=float)

            ax.fill_between(ds, vals-sds, vals+sds, color=color, alpha=0.10)

    ax.set_xlabel("Stromal density")

    ax.set_ylabel("Mean contact fraction")

    ax.set_title("Per-fate contact rate")

    ax.legend(fontsize=7)

  

    plt.tight_layout()

    if save:

        path = os.path.join(FIG_DIR, "SA4_stromal_density.png")

        fig.savefig(path); print(f"  → saved: {path}")

    if show:

        plt.show()

    plt.close(fig)

  
  

# ══════════════════════════════════════════════════════════════════════

# SA-5–9: MECHANISTIC PARAMETER SENSITIVITY

# Standardized OFAT sweeps. Thresholds are frozen at the baseline calibration.

# ══════════════════════════════════════════════════════════════════════

  

def sa5_kinetic_proofreading(

    k_p_vals: Tuple[float, ...] = (1.00, 2.67, 4.5, 10, 21),

    N_kp_vals: Tuple[float, ...] = (4.5, 7, 9, 11, 13),

    k_SHP1_vals: Tuple[float, ...] = (0.00, 0.025, 0.050, 0.075, 0.10),

    scenario: str = "uniform",

    seed: int | None = None,

    seeds: Tuple[int, ...] = DEFAULT_SEEDS,

) -> Dict:

    """Sensitivity of the kinetic-proofreading layer."""

    return _run_1d_param_group(

        {"k_p": tuple(k_p_vals), "N_kp": tuple(N_kp_vals), "k_SHP1": tuple(k_SHP1_vals)},

        "SA-5  Kinetic proofreading",

        scenario=scenario, seed=seed, seeds=seeds, recalibrate=False,

    )

  
  

def plot_sa5(res: Dict, save: bool = True, show: bool = False) -> None:

    _plot_1d_param_group(res, "SA-5 — Kinetic proofreading sensitivity",

                          "SA5_kinetic_proofreading.png", save=save, show=show)

  
  

def sa6_catch_bond(

    kc0_vals: Tuple[float, ...] = (0.25, 0.375, 0.50, 0.625, 0.75),

    ks0_vals: Tuple[float, ...] = (0.25, 0.375, 0.50, 0.625, 0.75),

    Fc_vals: Tuple[float, ...] = (6.0, 9.0, 12.0, 15.0, 18.0),

    Fs_vals: Tuple[float, ...] = (10.0, 15.0, 20.0, 25.0, 30.0),

    scenario: str = "uniform",

    seed: int | None = None,

    seeds: Tuple[int, ...] = DEFAULT_SEEDS,

) -> Dict:

    """Sensitivity of the catch-slip bond layer."""

    return _run_1d_param_group(

        {"kc0": tuple(kc0_vals), "ks0": tuple(ks0_vals),

         "Fc": tuple(Fc_vals), "Fs": tuple(Fs_vals)},

        "SA-6  Catch-slip bond parameters",

        scenario=scenario, seed=seed, seeds=seeds, recalibrate=False,

    )

  
  

def plot_sa6(res: Dict, save: bool = True, show: bool = False) -> None:

    _plot_1d_param_group(res, "SA-6 — Catch-slip bond sensitivity",

                          "SA6_catch_bond.png", save=save, show=show)

  
  

def sa7_force_transmission(

    F_T_vals: Tuple[float, ...] =  (10.0, 15.0, 20.0, 25.0, 30.0),

    K_F_vals: Tuple[float, ...] =  (1.25, 1.875, 2.50, 3.125, 3.75),

    beta_vals: Tuple[float, ...] = (0.0, 0.25, 0.50, 0.75, 1.0),

    scenario: str = "uniform",

    seed: int | None = None,

    seeds: Tuple[int, ...] = DEFAULT_SEEDS,

) -> Dict:

    """Sensitivity of force transmission: F_T, K_F and beta."""

    return _run_1d_param_group(

        {"F_T": tuple(F_T_vals), "K_F": tuple(K_F_vals), "beta": tuple(beta_vals)},

        "SA-7  Force transmission",

        scenario=scenario, seed=seed, seeds=seeds, recalibrate=False,

    )

  
  

def plot_sa7(res: Dict, save: bool = True, show: bool = False) -> None:

    _plot_1d_param_group(res, "SA-7 — Force transmission sensitivity",

                          "SA7_force_transmission.png", save=save, show=show)

  
  

def sa8_migration(

    v_max_vals:  Tuple[float, ...] = (18.0, 25.0, 32.0),

    K_high_vals: Tuple[float, ...] = (2.0, 3.0, 4.0, 5.0, 6.0),

    p_bias_vals: Tuple[float, ...] = (0.0, 0.02, 0.04, 0.08, 0.20),

    scenario: str = "uniform",

    seed: int | None = None,

    seeds: Tuple[int, ...] = DEFAULT_SEEDS,

    ) -> Dict:

    """Sensitivity of the migration layer."""

    return _run_1d_param_group(

        {"v_max": tuple(v_max_vals),

        "K_high": tuple(K_high_vals),

        "p_bias": tuple(p_bias_vals)},

        "SA-8  Migration parameters",

        scenario=scenario, seed=seed, seeds=seeds, recalibrate=False,

    )

  
  

def plot_sa8(res: Dict, save: bool = True, show: bool = False) -> None:

    _plot_1d_param_group(res, "SA-8 — Migration sensitivity",

                          "SA8_migration.png", save=save, show=show)

  
  

def sa9_lambda_X(

    lambda_X_vals: Tuple[float, ...] = (0.005, 0.0075, 0.010, 0.015, 0.020),

    scenario: str = "uniform",

    seed: int | None = None,

    seeds: Tuple[int, ...] = DEFAULT_SEEDS,

) -> Dict:

    """Sensitivity of the cumulative signal decay rate lambda_X."""

    return _run_1d_param_group(

        {"lambda_X": tuple(lambda_X_vals)},

        "SA-9  Signal decay rate",

        scenario=scenario, seed=seed, seeds=seeds, recalibrate=False,

    )

  
  

def plot_sa9(res: Dict, save: bool = True, show: bool = False) -> None:

    rows  = res["rows"]["lambda_X"]

    x_raw = [r["lambda_X"] for r in rows]

    x_scaled = [v * 1e3 for v in x_raw]   # convert to ×10⁻³ min⁻¹ → short numbers

  

    vals_pos = [r["pos_pct"]            for r in rows]

    vals_neg = [r["neg_pct"]            for r in rows]

    sds_pos  = [r.get("pos_pct_sd", 0.) for r in rows]

    sds_neg  = [r.get("neg_pct_sd", 0.) for r in rows]

  

    fig, ax = plt.subplots(figsize=(8, 5))

    fig.suptitle("SA-9 — Signal Decay Rate (λ_X) Sensitivity", fontweight="bold")

  

    ax.plot(x_scaled, vals_pos, "o-", color=C_POS, lw=1.5, ms=7, label="Positive sel.")

    ax.plot(x_scaled, vals_neg, "s-", color=C_NEG, lw=1.5, ms=7, label="Negative sel.")

    ax.fill_between(x_scaled,

                    np.array(vals_pos) - np.array(sds_pos),

                    np.array(vals_pos) + np.array(sds_pos),

                    color=C_POS, alpha=0.12)

    ax.fill_between(x_scaled,

                    np.array(vals_neg) - np.array(sds_neg),

                    np.array(vals_neg) + np.array(sds_neg),

                    color=C_NEG, alpha=0.12)

  

    ax.axvline(Params().lambda_X * 1e3, ls=":", color="black",

               alpha=0.5, label=f"Default ({Params().lambda_X * 1e3:.1f})")

  

    ax.set_xlabel("λ_X  (×10⁻³ min⁻¹)")

    ax.set_ylabel("Selection rate (%)")

    ax.set_title("Cumulative TCR signal decay rate")

    ax.legend(fontsize=7)

  

    plt.tight_layout()

    if save:

        path = os.path.join(FIG_DIR, "SA9_lambda_X.png")

        fig.savefig(path); print(f"  → saved: {path}")

    if show:

        plt.show()

    plt.close(fig)

  
  

# ══════════════════════════════════════════════════════════════════════

# SUMMARY TABLE

# ══════════════════════════════════════════════════════════════════════

  

def print_summary(all_res: Dict) -> None:

    print(f"\n{'═'*80}")

    print("  SENSITIVITY ANALYSIS SUMMARY")

    print(f"  n_initial={_SA_DEFAULTS['n_initial']},"

          f" grid={_SA_DEFAULTS['grid_w']}×{_SA_DEFAULTS['grid_h']},"

          f" T_run={_SA_DEFAULTS['T_run_hr']}h")

    print(f"{'═'*80}")

    # ← header now has a Note column wide enough for SNR annotations

    print(f"  {'Analysis':<38}  {'Δpos':>7}  {'Δneg':>7}  {'Δtotal':>7}  {'Note':<28}")

    print(f"  {'─'*38}  {'─'*7}  {'─'*7}  {'─'*7}  {'─'*28}")

  

    entries = []

  

    # SA-1: noise floor — uses SD directly, not _sensitivity_flag

    if "sa1" in all_res:

        r  = all_res["sa1"]

        dp = np.std([x["pos_pct"] for x in r["rows"]], ddof=1) * 2.0

        dn = np.std([x["neg_pct"] for x in r["rows"]], ddof=1) * 2.0

        entries.append(("SA-1  Seed noise (±2 SD)", dp, dn,

                        "verification — noise floor"))

  

    # SA-2: _sensitivity_flag available (multi-seed via _run_1d_param_group)

    if "sa2" in all_res:

        r = all_res["sa2"]

        dp_g, note_g     = _sensitivity_flag(r["gamma_rows"],  "pos_pct")

        dn_g, _          = _sensitivity_flag(r["gamma_rows"],  "neg_pct")

        dp_a_c, note_a_c = _sensitivity_flag(r["alpha_C_rows"],  "pos_pct")

        dn_a_c, _        = _sensitivity_flag(r["alpha_C_rows"],  "neg_pct")

        dp_a_m, note_a_m = _sensitivity_flag(r["alpha_M_rows"],  "pos_pct")

        dn_a_m, _        = _sensitivity_flag(r["alpha_M_rows"],  "neg_pct")

        entries.append(("SA-2  gamma sweep",       dp_g, dn_g, note_g))

        entries.append(("SA-2  alpha_C_max sweep", dp_a_c, dn_a_c, note_a_c))

        entries.append(("SA-2  alpha_M_max sweep", dp_a_m, dn_a_m, note_a_m))

  

    # SA-3: _sensitivity_flag available once ported to _run_1d_param_group

    #        Falls back to [no SD] note if SA-3 is still the old single-seed version

    if "sa3" in all_res:

        r = all_res["sa3"]

        shift_rows = r.get("shift_rows")

        sigma_rows = r.get("sigma_rows")

        if shift_rows is not None:

            # New version — ported to _run_1d_param_group

            dp_sh, note_sh = _sensitivity_flag(shift_rows, "pos_pct")

            dn_sh, _       = _sensitivity_flag(shift_rows, "neg_pct")

            dp_si, note_si = _sensitivity_flag(sigma_rows, "pos_pct")

            dn_si, _       = _sensitivity_flag(sigma_rows, "neg_pct")

            entries.append(("SA-3  koff shift_factor", dp_sh, dn_sh, note_sh))

            entries.append(("SA-3  koff sigma_factor", dp_si, dn_si, note_si))

        else:

            # Old version — single-seed, no SD

            flat = [x for v in r["sigma_results"].values() for x in v]

            dp, _ = _sensitivity_flag(flat, "pos_pct")

            dn, _ = _sensitivity_flag(flat, "neg_pct")

            entries.append(("SA-3  koff distribution", dp, dn,

                            "[no SD — single seed; port pending]"))

  

    # SA-4: _sensitivity_flag available (multi-seed via _run_1d_param_group)

    if "sa4" in all_res:

        r = all_res["sa4"]

        dp, note = _sensitivity_flag(r["rows"], "pos_pct")

        dn, _    = _sensitivity_flag(r["rows"], "neg_pct")

        entries.append(("SA-4  Stromal density", dp, dn, note))

  

    # SA-5 through SA-9: all use _run_1d_param_group, all have SD data

    for sa_key, sa_label in [

        ("sa5",  "SA-5   Kinetic proofreading"),

        ("sa6",  "SA-6   Catch-slip bond"),

        ("sa7",  "SA-7   Force transmission"),

        ("sa8",  "SA-8   Migration"),

        ("sa9",  "SA-9   lambda_X"),

    ]:

        if sa_key in all_res:

            r    = all_res[sa_key]

            flat = [row for rows in r["rows"].values() for row in rows]

            dp, note = _sensitivity_flag(flat, "pos_pct")

            dn, _    = _sensitivity_flag(flat, "neg_pct")

            entries.append((sa_label, dp, dn, note))

  

    # Sort: SA-1 always first as the noise floor reference,

    # then by combined range descending so the most influential parameters

    # appear at the top.

    sa1_entries  = [e for e in entries if e[0].startswith("SA-1")]

    rest_entries = sorted(

        [e for e in entries if not e[0].startswith("SA-1")],

        key=lambda x: -(x[1] + x[2]),

    )

  

    for name, dp, dn, note in sa1_entries + rest_entries:

        print(f"  {name:<38}  {dp:>7.1f}%  {dn:>7.1f}%"

              f"  {dp+dn:>7.1f}%  {note:<28}")

  

    print(f"{'═'*80}")

    print("  Δ = max − min across the swept range (means across seeds).")

    print("  SNR = Δ / mean within-condition SD across seeds.")

    print("  SNR > 2 → effect exceeds noise floor → REAL.")

    print("  SNR ≤ 2 → indistinguishable from noise → NOISE.")

    print(f"{'═'*80}")

  
  

# ══════════════════════════════════════════════════════════════════════

# MASTER RUNNER

# ══════════════════════════════════════════════════════════════════════

  

def run_all(

        which: Tuple[int, ...] = (1, 2, 3, 4, 5, 6, 7, 8, 9),

        show: bool = False) -> Dict:

    all_res: Dict = {}

    t_total = time.time()

  

    run_map = {

        1: (sa1_seed_reproducibility,            plot_sa1,  "sa1"),

        2: (sa2_mechanobiological_amplification, plot_sa2,  "sa2"),

        3: (sa3_koff_distribution,               plot_sa3,  "sa3"),

        4: (sa4_stromal_density,                 plot_sa4,  "sa4"),

        6: (sa5_kinetic_proofreading,            plot_sa5,  "sa5"),

        7: (sa6_catch_bond,                      plot_sa6,  "sa6"),

        8: (sa7_force_transmission,              plot_sa7,  "sa7"),

        9: (sa8_migration,                       plot_sa8,  "sa8"),

        10:(sa9_lambda_X,                        plot_sa9,  "sa9"),

    }

    for n in sorted(which):

        if n not in run_map:

            print(f"  Warning: SA-{n} not defined. Skipping.")

            continue

        fn, plot_fn, key = run_map[n]

        res = fn()

        plot_fn(res, save=True, show=show)

        all_res[key] = res

  

    print(f"\n  Total wall time: {(time.time()-t_total)/60:.1f} min")

    print_summary(all_res)

    return all_res

  
  

# ══════════════════════════════════════════════════════════════════════

# ENTRY POINT

# ══════════════════════════════════════════════════════════════════════

  

if __name__ == "__main__":

    parser = argparse.ArgumentParser(

        description="Thymic selection model — sensitivity analyses")

    parser.add_argument(

        "--sa", nargs="+", type=int, default=None,

        metavar="N",

        help="Analyses to run, eg. --sa 3 or --sa 3 7 9",

    )

    parser.add_argument(

        "--show",

        action="store_true",

        help="display plots interactively after each completed SA"

    )

    args = parser.parse_args()

    which = (

        tuple(args.sa)

        if args.sa

        else (1, 2, 3, 4, 5, 6, 7, 8, 9)

        )

    run_all(which=which, show=args.show)