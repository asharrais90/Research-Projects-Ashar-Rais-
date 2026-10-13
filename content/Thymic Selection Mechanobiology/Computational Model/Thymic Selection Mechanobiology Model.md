"""

Thymic_selection_mechanobiolgy.py

======================

Hybrid Agent-Based + ODE Model of Thymocyte Selection (Currently at proof of concept level)

  

The main purpose of this model is to investigate whether the stiffness of thymic stroma influences

thymic selection in a consistent and predictable manner. The primary hypothesis regarding the possible

link between stiffness and selection is the mechanotransductional aspect of TCR siganlling, and

thymocytes as whole. This signalling aspect is modelled in 4 layers: a thymocyte's migrational capacity

that influences how much it comes in contact with the antigen presenting cells (APCs), its capacity to scan the APCs,

the state of specific pathways that influence the intensity of downstream TCR signalling and TCR dynamics plus kinetics

that regulate how long a TCR stays bound to pMHC leading to signal production. These layers are dependent upon the state

of cytoskeletal networks which are influenced by the state of extracellular stroma. These modelling decisions lead to

emergence of a short-term (comapred to the residence time in the thymus) "mechanical memory" that is path-dependent and

varies from one thymocyte to another, providing a novel framework to study thymocyte selection.

  

────────────────────

Mesa 3.x  (tested on 3.5.1)

  

Install:  pip install mesa>=3.0 numpy scipy matplotlib

────────────────────

  

"""

  

from __future__ import annotations

  

import warnings

from dataclasses import dataclass

from typing import Dict, Optional, Set, Tuple

  

import numpy as np

from mesa import Agent, DataCollector, Model

from mesa.space import SingleGrid

  

from population_distributions import (

    generate_koff_values,

    generate_zap70_values,

)

  

try:

    from scipy.ndimage import gaussian_filter

    _SCIPY = True

except ImportError:

    _SCIPY = False

    warnings.warn("scipy not found — no smoothing applied.")

  

try:

    import matplotlib.pyplot as plt

    _MPL = True

except ImportError:

    _MPL = False

    warnings.warn("matplotlib not found — plotting disabled.")

  
  

# ══════════════════════════════════════════════════════════════════════

# 1.  PARAMETERS

# ══════════════════════════════════════════════════════════════════════

  

@dataclass

class Params:

    # ── Grid ──────────────────────────────────────────────────────────

    grid_w:       int   = 100

    grid_h:       int   = 100

    grid_spacing: float = 10.0     #μm per grid cell

    cortex_frac:  float = 0.70     #cortex is modelled to be larger than medulla based on imaging studies

    # ── Stiffness field ───────────────────────────────────────────────

    scenario:    str   = "uniform"   # "uniform"|"cortex_stiffer"|"medulla_stiffer"

    S_mean:      float = 9.0        # kPa

    S_std:       float = 8.0         # kPa  (spatial noise)

    S_smooth:    float = 1.5         # Gaussian σ (grid cells)

    S_stiffer:   float = 16.5         # kPa (higher mean for non-uniform scenarios)

    S_softer:    float = 13.5         # kPa (lower mean for non-uniform scenarios)

  

    # ── Stromal cells ─────────────────────────────────────────────────

    stromal_density: float = 0.001

  

    # ── Population dynamics ─────────────────────────────────────

    T_run_hr:      float = 288.0  #12 days  

    n_initial:     float = 1000

  

    # ── koff population distribution ─────────────────────────────────

    # Controls the shape of the log-normal TCR affinity distribution.

    # shift_factor < 1 moves distribution toward stronger binders (left).

    # sigma_factor controls repertoire breadth (current value very wide:

    # sigma_log ≈ 7.6, spanning ~12 decades — see population_distributions.py).

    koff_shift_factor: float = 0.3

    koff_sigma_factor: float = 6.0

    koff_lo:    float = 0.001   # s⁻¹ clip

    koff_hi:    float = 50.0   # s⁻¹ clip

  

    # ── Other permanent property distributions ──────────────────────────────

    zap70_mean: float = 1.0

    zap70_cv:   float = 0.3

    zap70_lo:   float = 0.1

    zap70_hi:   float = 4.0

    Tres_shape: float = 5.0      # Gamma distribution for the residence time of thymocytes;

    Tres_scale: float = 1152.0   # mins; mean = 5760 mins = 4 days & mean ~ 3.9 days; Chosen from pulse-chase studies

  

    # ── Biphasic migration ──────────────────────────────────────

    v_max:  float = 25.0    # μm/min peak speed

    K_low:  float = 0.5     # kPa

    n_low:  float = 1.5

    K_high: float = 4.0     # kPa

    n_high: float = 2.0

    p_bias: float = 0.04    # biased step towards the medulla

  

    # ── Cytoskeletal priming ODE ──────────────────────────────────────

    k_C: float = 0.05        # min⁻¹  adaptation rate  (τ ≈ 20 min)

    K_C: float = 2.5         # kPa   half-max priming stiffness

  

    # ── Mechanotransduction ODE ────────────────────────────

    k_M:     float = 0.03    # min⁻¹  (τ ≈ 33 min)

    K_M:     float = 2.5     # kPa   half-max mechanosensitive stiffness

  

    # ── Signal amplification ────────────────────────────────────

    # α(C,M) = 1 + α_C_max·C + α_M_max·M

    # C-channel: actin-dependent ITAM phosphorylation efficiency

    # M-channel: Piezo1/YAP-TAZ downstream amplification

    alpha_C_max: float = 1.0   # dimensionless (C in [0,1])

    alpha_M_max: float = 0.5   # dimensionless (M in [0,1])

  

    # ── Force transmission ────────────────────────────────────────────

    F_T:   float = 20.0     # pN   saturating actomyosin force

    K_F:   float = 2.5      # kPa  half-max stiffness for force

    beta:  float = 0.5      # C amplifies Feff (bond stabilisation)

  

    # ── Catch-slip bond  ───────────────────

    # G(0) = 1 exactly by construction

    kc0: float = 0.5    # s⁻¹  catch zero-force rate

    Fc:  float = 12.0   # pN   catch characteristic force

    ks0: float = 0.5    # s⁻¹  slip zero-force rate

    Fs:  float = 20.0   # pN   slip characteristic force

  

    # ── Kinetic proofreading ──────────────────────────────────────────

    k_p:    float = 1.0    # s⁻¹  

    k_SHP1: float = 0.05    # s⁻¹  Contribution of SHP‐1 to the dephosphorylation rate (Richez & Khadra, 2026)

    N_kp:   float = 11.0    # proofreading steps;

  

    # ── Scanning capacity  ───────────────────────────────────────

    # k_scan(C) = k_scan_base × (1 + γ·C)

    # Represents how thoroughly a thymocyte scans the pMHC repertoire:

    # microvilli formation, penetration of physical barriers, force-

    # dependent contact formation — all cytoskeleton-dependent.

    k_scan_base: float = 1.0   # base capacity (normalised; dimensionless)

    gamma:       float = 0.5   # priming effect on scanning capacity

  

    # ── Signal ODE ───────────────────────────────────────────────────

    lambda_X: float = 0.01   # min⁻¹  signal decay rate

  

    # ── Selection thresholds  ────────────────────────────────────

    # Computed by calibrate_thresholds(); anchored to TCR 42.12 / H-2Kb data.

    # Set at biochemical baseline (C=0, M=0, S=0 → α=1, G=1, no force).

    theta_pos: float = 0.0   # neglect / positive selection boundary

    theta_neg: float = 0.0   # positive / negative selection boundary

  

    # Reference ligands (Grakoui et al. 1999 | 2B4 TCR binding to MCC-derived peptides presented by I-Eᵏ)

    koff_null:    float = 5.00    # s⁻¹  T102G       (ANERADLIAYLKQAGK) — neglect

    koff_select:  float = 0.36    # s⁻¹  T102S       (ANERADLIAYLKQASK) — positive selection

    koff_delete:  float = 0.057   # s⁻¹  MCC88-103   (ANERADLIAYLKQATK) — negative selection

  

    # ── Numerics ─────────────────────────────────────────────────────

    dt:   float = 1.0    # min per timestep (Euler)

    seed: int   = 42

  
  

# ══════════════════════════════════════════════════════════════════════

# 2.  THRESHOLD CALIBRATION

# ══════════════════════════════════════════════════════════════════════

  

def calibrate_thresholds(p: Params) -> Tuple[float, float]:

    """

    Anchor selection thresholds to 2B4 TCR binding to MCC-derived peptides presented by I-Eᵏ.

    Reference ligands (Grakoui et al. 1999):

       T102G        ANERADLIAYLKQAGK |  koff = 5.00 s⁻¹   → death by neglect

       T102S        ANERADLIAYLKQASK |  koff = 0.36 s⁻¹   → positive selection

       MCC88-103    ANERADLIAYLKQATK |  koff = 0.057 s⁻¹  → negative selection

  

    Reference state: C=0, M=0, S=0 → α=1, G(F=0)=1, k_scan=k_scan_base

    This is the biochemical baseline: no cytoskeletal amplification,

    no force-dependent bond stabilisation, base scanning capacity.

  

    Thresholds are means of adjacent X_ss pairs. After this call, θ_pos and θ_neg

    are FROZEN for all subsequent predictions.  The mechanical model then measures

    whether stiffness helps or hinders cells in crossing these fixed requirements.

  

    """

    kp_Z = p.k_p * p.zap70_mean

  

    def X_ss_baseline(koff_base: float) -> float:

        """X_ss at biochemical baseline — no mechanical context."""

        p_step = kp_Z / (kp_Z + koff_base + p.k_SHP1)   # G=1, koff_eff=koff_base

        KP     = p_step ** p.N_kp

        return p.k_scan_base * 1.0 * KP / p.lambda_X     # α=1, k_scan=base

  

    X_null   = X_ss_baseline(p.koff_null)

    X_select = X_ss_baseline(p.koff_select)

    X_delete = X_ss_baseline(p.koff_delete)

  

    theta_pos = float(np.sqrt(X_null   * X_select))

    theta_neg = float(np.sqrt(X_select * X_delete))

  

    return theta_pos, theta_neg

  
  

# ══════════════════════════════════════════════════════════════════════

# 3.  THYMOCYTE AGENT

# ══════════════════════════════════════════════════════════════════════

  

class Thymocyte(Agent):

    """

    Single thymocyte navigating the thymic stiffness grid.

  

    ODEs evolved per timestep (Euler, Δt = 1 min):

        dC/dt = k_C · [Ceq(S) − C]

        dM/dt = k_M · [Meq(S) − M]                    

        dX/dt = E · k_scan(C) · α(C,M) · KP(Feff(C)) − λX · X

  

    Signal production architecture:

        k_scan(C)   = k_scan_base·(1+γ·C)                   APC surface scanning capacity of a thymocyte

        α(C,M)      = 1 + α_C,max·C + α_M,max·M             Efficiency gain for physical force transmission

        Feff        = F(S)·(1 + β·C)                        Effective force transmitted to intracellular networks

        KP(Feff)    = [kp_Z/(kp_Z+koff·G(Feff)+k_SHP1)]^N   Kinetic proofreading

    """

  

    def __init__(self, model: ThymicEnvironment, koff_base: float, zap70: float) -> None:

        super().__init__(model)      # Mesa 3.x: auto-registers with model.agents

        p   = model.p

        rng = model.rng

  

        # ── Permanent properties (fixed at birth) ─────────────────────

        self.koff_base: float = float(np.clip(

            koff_base, p.koff_lo, p.koff_hi))

        self.ZAP70i: float = float(np.clip(

            zap70, p.zap70_lo, p.zap70_hi))

        self.T_res: float = float(rng.gamma(p.Tres_shape, p.Tres_scale))

  

        # ── Dynamic ODE state ─────────────────────────────────────────

        self.age: float = 0.0

        self.C:   float = 0.0    # cytoskeletal priming      [0, 1]

        self.M:   float = 0.0    # mechanosensitive signal   [0, 1]

        self.X:   float = 0.0    # cumulative TCR signal     [0, ∞)

  

        # ── Fate ──────────────────────────────────────────────────────

        self.fate:  Optional[str] = None

        self.alive: bool          = True

  

        # ── Movement history ──────────────────────────────────────────

        self.history: list = []

        self.contact_steps: int = 0

        self.contact_S_sum: int = 0  #sum of stiffness during contact steps only

  

    # ── Environment queries ───────────────────────────────────────────

  

    def _S(self) -> float:

        x, y = self.pos

        return float(self.model.stiffness[x, y])

  

    def _E(self) -> int:

        """E(t)=1 if any Moore-neighbour is a stromal cell."""

        for nb in self.model.grid.iter_neighborhood(

                self.pos, moore=True, include_center=True):

            if nb in self.model.stromal_pos:

                return 1

        return 0

  

    # ── Mechanical calculations ───────────────────────────────────────

  

    def _F_of_S(self, S: float) -> float:

        p = self.model.p

        return p.F_T * S / (p.K_F + S)

  

    def _G(self, Feff: float) -> float:

        """Catch-slip modulation; G(0)=1 exactly."""

        p = self.model.p

        return (p.kc0 * np.exp(-Feff / p.Fc) +

                p.ks0 * np.exp( Feff / p.Fs)) / (p.kc0 + p.ks0)

  

    def _KP(self, S: float) -> float:

        """

        Kinetic proofreading probability.

        Path of S influencing TCR signalling: C amplifies Feff → reduces koff.

        """

        p    = self.model.p

        Feff = self._F_of_S(S) * (1.0 + p.beta * self.C)

        koff = self.koff_base * self._G(Feff)

        kp_Z = p.k_p * self.ZAP70i

        return (kp_Z / (kp_Z + koff + p.k_SHP1)) ** p.N_kp

  

    # ── ODE right-hand sides ─────────────────────────────────────────

  

    def _dC_dt(self, S: float) -> float:

        p = self.model.p

        return p.k_C * (S / (p.K_C + S) - self.C)

  

    def _dM_dt(self, S: float) -> float:

        """

        Mechanosensitive competence ODE.

        dM/dt = k_M · (M_ss(S) - M),   M_ss(S) = S / (K_M + S)

        Represents integrated Piezo1 / YAP-TAZ / MRTF-SRF activity.

        k_M is a phenonomenological term that encodes how quickly

        downstream mechanosensitive pathways respond to changes in cytoskeletal state.

        """

        p = self.model.p

        return p.k_M * (S / (p.K_M + S) - self.M)

  

    def _dX_dt(self, E: int, S: float) -> float:

        """

        Cumulative TCR signal ODE. Production term incorporates all three mechanical channels:

          k_scan(C)  

          α(C,M)    

          KP(Feff(C))

  

        At C=0, M=0: k_scan=k_scan_base, α=1 → matches calibration baseline.

        Any priming (C>0, M>0) amplifies signal above that baseline.

        This is the mechanobiological effect the model is designed to capture.

        """

        p     = self.model.p

        decay = p.lambda_X * self.X

        if E == 0:

            return -decay

  

        k_scan_val = p.k_scan_base * (1.0 + p.gamma * self.C)           # [C1]

        alpha      = 1.0 + p.alpha_C_max * self.C + p.alpha_M_max * self.M  # [C3]

        KP         = self._KP(S)

        return k_scan_val * alpha * KP - decay

  

    # ── Migration  [F6] ──────────────────────────────────────────────

  

    def _migrate(self, S: float) -> None:

        """

        Biphasic v(S) = vmax·[rising Hill]·[falling Hill].

        Optional bias toward medulla (chemokine gradient proxy).

        """

        p   = self.model.p

        rng = self.model.rng

        S   = max(S, 1e-9)

        rising  = S**p.n_low  / (p.K_low**p.n_low  + S**p.n_low)

        falling = p.K_high**p.n_high  / (p.K_high**p.n_high + S**p.n_high)

  

        v = p.v_max * rising * falling

        n_hops = int((v * p.dt)/p.grid_spacing)            #how many grid cells can be moved across in a timestep

        frac = ((v * p.dt)/p.grid_spacing) - n_hops        #leftover "hopping" capacity

        if rng.random() < frac:

            n_hops += 1

  

        for _ in range(n_hops):

            candidates = list(self.model.grid.iter_neighborhood(

                self.pos, moore=True, include_center=False))

            if not candidates:

                break

            if rng.random() < p.p_bias:

                biased = [c for c in candidates if c[1] > self.pos[1]]

                if biased:

                    candidates = biased

            chosen = candidates[rng.integers(len(candidates))]

            if self.model.grid.is_cell_empty(chosen):

                self.model.grid.move_agent(self, chosen)

            else:

                break

  

    # ── Fate evaluation ───────────────────────────────────────────────

  

    def _check_fate(self) -> None:

        """

        Fate logic against frozen thresholds (θ_pos, θ_neg).

  

        Negative selection: tested every step — immediate deletion if

        X ≥ θ_neg (strong agonist signal → clonal deletion).

  

        End-of-residence:

          X ≥ θ_pos  → positive selection

          else        → death by neglect

        """

        p = self.model.p

        if self.X >= p.theta_neg:

            self.fate  = "negative_selection"

            self.alive = False

            return

        if self.age < self.T_res:

            return

        if self.X >= p.theta_pos:

            self.fate = "positive_selection"

        else:

            self.fate = "death_by_neglect"

        self.alive = False

  

    # ── Main step ─────────────────────────────────────────────────────

  

    def step(self) -> None:

        if not self.alive:

            return

        if not self.history:

            self.history.append((self.age, self.pos[0], self.pos[1]))

        dt = self.model.p.dt

        S  = self._S()

        E  = self._E()

        if E == 1:

            self.contact_steps += 1

            self.contact_S_sum += S

        self._migrate(S)

        # Euler integration

        self.C = float(np.clip(self.C + self._dC_dt(S)    * dt, 0.0, 1.0))

        self.M = float(np.clip(self.M + self._dM_dt(S)    * dt, 0.0, 1.0))

        self.X = max(0.0,      self.X + self._dX_dt(E, S) * dt)

        self.history.append((self.age, self.pos[0], self.pos[1]))

        self._check_fate()

        self.age += dt

  
  

# ══════════════════════════════════════════════════════════════════════

# 4.  THYMIC ENVIRONMENT (MODEL)

# ══════════════════════════════════════════════════════════════════════

  

class ThymicEnvironment(Model):

    """

    200 × 200 spatial grid with stiffness field and stromal cells.

    Mesa 3.x: activation via self.agents.shuffle_do("step").

    """

  

    def __init__(self, params: Optional[Params] = None) -> None:

        super().__init__()

        self.p   = params or Params()

        self.rng = np.random.default_rng(self.p.seed)

        self.t_min: float = 0.0

        self.koff_population = generate_koff_values(self.p.n_initial,self.rng,

                                                    shift_factor=self.p.koff_shift_factor,

                                                    sigma_factor=self.p.koff_sigma_factor,

                                                    )

        self.zap70_population = generate_zap70_values(self.p.n_initial, self.rng, mean=self.p.zap70_mean, cv=self.p.zap70_cv)

  

        self.fate_counts: Dict[str, int] = {

            "positive_selection":   0,

            "negative_selection":   0,

            "death_by_neglect":     0,

        }

        self.grid       = SingleGrid(self.p.grid_w, self.p.grid_h, torus=False)

        self.stiffness  = self._build_stiffness()

        self.stromal_pos: Set[Tuple[int, int]] = self._place_stromal()

        self._introduce_initial()

        self.koff_by_fate: Dict[str, list] = {

            "positive_selection": [],

            "negative_selection": [],

            "death_by_neglect"  : [],

        }

        self.contact_stats: Dict[str, list] = {

            "positive_selection":   [],

            "negative_selection":   [],

            "death_by_neglect":     [],

        }

        self.trajectories: Dict[str, list] = {

            "positive_selection":  [],

            "negative_selection": [],

            "death_by_neglect":   [],

        }

        self.datacollector = DataCollector(model_reporters={

            "t_hr":      lambda m: m.t_min / 60.0,

            "n_alive":   lambda m: sum(1 for a in m.agents

                                       if isinstance(a, Thymocyte) and a.alive),

            "mean_C":    lambda m: _mean(m, "C"),

            "mean_M":    lambda m: _mean(m, "M"),

            "mean_X":    lambda m: _mean(m, "X"),

            "mean_kscan":lambda m: _mean_kscan(m),

            "mean_alpha":lambda m: _mean_alpha(m),

            "pos_sel":   lambda m: m.fate_counts["positive_selection"],

            "neg_sel":   lambda m: m.fate_counts["negative_selection"],

            "neglect":   lambda m: m.fate_counts["death_by_neglect"],

        })

  

    # ── Environment construction ──────────────────────────────────────

  

    def _build_stiffness(self) -> np.ndarray:

        p  = self.p

        W, H = p.grid_w, p.grid_h

        ch   = int(H * p.cortex_frac)

  

        if p.scenario == "uniform":

            base = self.rng.normal(p.S_mean, p.S_std, (W, H))

        elif p.scenario == "cortex_stiffer":

            base = np.empty((W, H))

            base[:, :ch] = self.rng.normal(p.S_stiffer,  p.S_std, (W, ch))

            base[:, ch:] = self.rng.normal(p.S_softer, p.S_std, (W, H - ch))

        elif p.scenario == "medulla_stiffer":

            base = np.empty((W, H))

            base[:, :ch] = self.rng.normal(p.S_softer, p.S_std, (W, ch))

            base[:, ch:] = self.rng.normal(p.S_stiffer,  p.S_std, (W, H - ch))

        else:

            raise ValueError(f"Unknown scenario: {p.scenario!r}")

  

        smoothed = gaussian_filter(base, sigma=p.S_smooth) if _SCIPY else base

        return np.clip(smoothed, 0.001, 50.0).astype(np.float32)

  

    def stiffness_report(model: "ThymicEnvironment",

                     expected_range:  tuple[float, float] | None = None,

                     quartile_range:  tuple[float, float] | None = None,

                     verbose: bool = True) -> dict:

        """

        Summarise the stiffness field of a single run.

  

        Parameters

        ----------

        model          : ThymicEnvironment

        expected_range : (lo, hi) kPa — check whether the *mean* lies in this window

        quartile_range : (lo, hi) kPa — split into 4 equal sub-bins and report

                        the fraction of grid points in each (Q1..Q4)

        verbose        : print a formatted report

  

        Returns

        -------

        dict with point stats, coarse histogram (0–5 / 5–10 / 10–15 / 15+),

        and — when quartile_range is given — per-quartile counts and fractions.

        """

        S = np.asarray(model.stiffness)

        N = S.size

  

        stats = {

            "scenario": model.p.scenario,

            "n_points": int(N),

            "mean":     float(S.mean()),

            "std":      float(S.std()),

            "min":      float(S.min()),

            "max":      float(S.max()),

            "median":   float(np.median(S)),

        }

  

        # ── coarse histogram (fixed buckets) ──────────────────────────────

        counts, _ = np.histogram(S, bins=[0, 5, 10, 15, np.inf])

        for label, f in zip(("frac_0_5", "frac_5_10",

                            "frac_10_15", "frac_15p"), counts / N):

            stats[label] = float(f)

  

        # ── quartile breakdown of a user-specified range ──────────────────

        if quartile_range is not None:

            lo, hi = quartile_range

            if not (hi > lo):

                raise ValueError("quartile_range must satisfy hi > lo")

  

            edges = np.linspace(lo, hi, 5)                # 4 equal-width bins

            sub_counts, _ = np.histogram(S, bins=edges)

            in_range = int(sub_counts.sum())

  

            stats["quartile_range"]  = (float(lo), float(hi))

            stats["quartile_edges"]  = [float(e) for e in edges]

            stats["n_in_range"]      = in_range

            stats["frac_in_range"]   = in_range / N

  

            stats["quartile_counts"]     = [int(c) for c in sub_counts]

            stats["quartile_frac_all"]   = [float(c / N)          for c in sub_counts]

            stats["quartile_frac_range"] = [

                float(c / in_range) if in_range else 0.0 for c in sub_counts

            ]

            stats["n_clipped"] = int((S <= 0.001).sum())

            stats["frac_clipped"] = stats["n_clipped"] / N

  

        # ── printing ──────────────────────────────────────────────────────

        if verbose:

            print(f"[{model.p.scenario:>16}] "

                f"mean={stats['mean']:6.2f}  std={stats['std']:5.2f}  "

                f"min={stats['min']:5.2f}  max={stats['max']:6.2f}  |  "

                f"0–5:{stats['frac_0_5']:5.1%}  "

                f"5–10:{stats['frac_5_10']:5.1%}  "

                f"10–15:{stats['frac_10_15']:5.1%}  "

                f"15+:{stats['frac_15p']:5.1%}  "

                f"n_clipped: {stats['n_clipped']}  "

                f"frac_clipped: {stats['frac_clipped']}")

  

            if expected_range is not None:

                lo, hi = expected_range

                ok = "✓" if lo <= stats["mean"] <= hi else "✗"

                print(f"    {ok} expected mean in [{lo}, {hi}] kPa")

  

            if quartile_range is not None:

                lo, hi = quartile_range

                print(f"    Quartile split of [{lo}, {hi}] kPa  "

                    f"(n_in_range = {stats['n_in_range']} "

                    f"= {stats['frac_in_range']:.1%} of grid):")

                for i, (c, fa, fr) in enumerate(zip(

                        stats["quartile_counts"],

                        stats["quartile_frac_all"],

                        stats["quartile_frac_range"]), start=1):

                    a, b = edges[i-1], edges[i]

                    bar = "█" * int(30 * fa)

                    print(f"      Q{i} [{a:6.3f}, {b:6.3f})  "

                        f"n={c:7d}  "

                        f"{fa:6.2%} of all  |  "

                        f"{fr:6.2%} within  {bar}")

  

        return stats

  

    def _place_stromal(self) -> Set[Tuple[int, int]]:

        p   = self.p

        all_pos = [(x, y) for x in range(p.grid_w) for y in range(p.grid_h)]

        n   = int(len(all_pos) * p.stromal_density)

        idx = self.rng.choice(len(all_pos), size=n, replace=False)

        return {tuple(all_pos[i]) for i in idx}

  

    def plot_stiffness_map(model: ThymicEnvironment, save: str = "") -> None:

        if not _MPL:

            return

        fig, ax = plt.subplots(figsize=(6, 6))

        im = ax.imshow(model.stiffness.T, origin="lower", cmap="viridis",

                    extent=[0, model.p.grid_w, 0, model.p.grid_h])

        if model.stromal_pos:

            xs, ys = zip(*model.stromal_pos)

            ax.scatter(xs, ys, s=4, c="white", edgecolors="black",

                   linewidths=0.3, label="stromal cells")

            ax.legend(fontsize=8, loc="upper right")

        ax.set_xlabel("x (grid units)"); ax.set_ylabel("y (grid units)")

        ax.set_title(f"Stiffness field — {model.p.scenario}")

        fig.colorbar(im, ax=ax, label="Stiffness (kPa)")

        plt.tight_layout()

        if save:

            plt.savefig(save, dpi=150, bbox_inches="tight")

        plt.show()

  

    # ── Agent lifecycle  ─────────────────────────────────────────

  

    def _introduce_initial(self) -> None:

        """

        In this version all the thymocytes are introduced in the cortex at once.

        Although this does not replicate physiological continous migration,

        It works in this version primarily because processes such as crowding effects

        are not being modeled currently and will be introduced in a future version.

        """

        p     = self.p

        ch    = int(p.grid_h * p.cortex_frac)

        empty = [(x, y) for x in range(p.grid_w) for y in range(ch)

                  if self.grid.is_cell_empty((x, y))]

        if not empty:

            return

        n_new = min(p.n_initial , len(empty))

        for j, i in enumerate(

            self.rng.choice(len(empty), size=n_new, replace = False)

        ):

            a = Thymocyte(self,

                          koff_base=self.koff_population[j],

                          zap70=self.zap70_population[j])

            self.grid.place_agent(a, empty[i])

  

    def _reap(self) -> None:

        """Remove dead agents. Mesa 3.x requires both calls."""

        dead = [a for a in self.agents

                if isinstance(a, Thymocyte) and not a.alive]

        for a in dead:

            if a.fate in self.fate_counts:

                self.fate_counts[a.fate] += 1

            if a.fate in self.koff_by_fate:

                self.koff_by_fate[a.fate].append(a.koff_base)

            if a.fate in self.trajectories:

                self.trajectories[a.fate].append(a.history)

            if a.fate in self.contact_stats:                                  

                mean_S_during_contact = (a.contact_S_sum / a.contact_steps  

                                      if a.contact_steps > 0 else 0.0)    

                self.contact_stats[a.fate].append({                        

                    "contact_steps": a.contact_steps,                        

                    "contact_min":   a.contact_steps * self.p.dt,            

                    "contact_frac":  a.contact_steps * self.p.dt / max(a.age, 1e-9),

                    "mean_S":        mean_S_during_contact,                  

                })                                                            

            self.grid.remove_agent(a)   # remove from spatial grid

            a.remove()                   # deregister from model.agents

  

    # ── Model step ───────────────────────────────────────────────────

  

    def step(self) -> None:

        self.t_min += self.p.dt

        self.agents.shuffle_do("step")

        self._reap()

        self.datacollector.collect(self)

  
  

# ── DataCollector helpers ─────────────────────────────────────────────

  

def _mean(m: ThymicEnvironment, attr: str) -> float:

    vals = [getattr(a, attr) for a in m.agents

            if isinstance(a, Thymocyte) and a.alive]

    return float(np.mean(vals)) if vals else 0.0

  

def _mean_kscan(m: ThymicEnvironment) -> float:

    """Population mean scanning capacity k_scan(C)."""

    p    = m.p

    vals = [p.k_scan_base * (1.0 + p.gamma * a.C)

            for a in m.agents if isinstance(a, Thymocyte) and a.alive]

    return float(np.mean(vals)) if vals else p.k_scan_base

  

def _mean_alpha(m: ThymicEnvironment) -> float:

    """Population mean signal amplification α(C,M)."""

    p    = m.p

    vals = [1.0 + p.alpha_C_max * a.C + p.alpha_M_max * a.M

            for a in m.agents if isinstance(a, Thymocyte) and a.alive]

    return float(np.mean(vals)) if vals else 1.0

  
  
  

# ══════════════════════════════════════════════════════════════════════

# 5.  ANALYSIS

# ══════════════════════════════════════════════════════════════════════

  

def print_summary(model: ThymicEnvironment) -> None:

    fc    = model.fate_counts

    total = max(sum(fc.values()), 1)

    p     = model.p

    print(f"\n{'─'*60}")

    print(f"  Scenario : {p.scenario}")

    print(f"  Time     : {model.t_min/60:.1f} h")

    print(f"  θ_pos    = {p.theta_pos:.4e}")

    print(f"  θ_neg    = {p.theta_neg:.4e}")

    print(f"  (calibrated at C=0, M=0, S=0 — biochemical baseline)")

    print(f"{'─'*60}")

    for fate, n in fc.items():

        bar = "█" * int(25 * n / total)

        print(f"  {fate:<28} {n:5d}  {100*n/total:5.1f}%  {bar}")

    print(f"{'─'*60}")

  

def koff_summary(model: ThymicEnvironment) -> None:

    for fate, vals in model.koff_by_fate.items():

        if not vals:

            print(f"{fate:22s}: no resolved cells yet")

            continue

        print(f"{fate:22s} n={len(vals):4d} | "

              f"⟨koff_base⟩={np.mean(vals):.3f} s⁻¹ | "

              f"median={np.median(vals):.3f} | "

              f"std={np.std(vals):.3f}")

  

def plot_stiffness_map(model: ThymicEnvironment, save: str = "") -> None:

    if not _MPL:

        return

  

    p = model.p

    fig, ax = plt.subplots(figsize=(6, 6))

    im = ax.imshow(model.stiffness.T, origin="lower", cmap="RdYlGn_r",

                    extent=[0, model.p.grid_w, 0, model.p.grid_h])

  

    if model.stromal_pos:

        xs, ys = zip(*model.stromal_pos)

        ax.scatter(xs, ys, s=4, c="white", edgecolors="black",

                   linewidths=0.3, label="stromal cells")

        ax.legend(fontsize=8, loc="upper right")

  

    ax.set_xlabel("x (grid units)")

    ax.set_ylabel(f"y (grid units) | cortex limit = {p.grid_h * p.cortex_frac}")

    ax.set_title(f"Stiffness field — {model.p.scenario}")

    fig.colorbar(im, ax=ax, label="Stiffness (kPa)")

  

    plt.tight_layout()

    if save:

        plt.savefig(save, dpi=150, bbox_inches="tight")

    plt.show()

    plt.close(fig)

  

def plot_dynamics(model: ThymicEnvironment, save: str = "") -> None:

    if not _MPL:

        return

    df = model.datacollector.get_model_vars_dataframe()

    t  = df["t_hr"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 7))

    fig.suptitle(f"Thymic dynamics — {model.p.scenario}", fontweight="bold")

  

    axes[0,0].plot(t, df["n_alive"], "#378ADD", lw=1.5)

    axes[0,0].set_xlabel("Time (h)"); axes[0,0].set_ylabel("Alive thymocytes")

    axes[0,0].set_title("Population")

  

    axes[0,1].plot(t, df["mean_X"], "#EF9F27", lw=1.5)

    axes[0,1].set_xlabel("Time (h)"); axes[0,1].set_ylabel("⟨X⟩")

    axes[0,1].set_title("Average cumulative TCR signal of alive thymocytes")

  

    axes[1,0].plot(t, df["pos_sel"], "#378ADD", lw=1.5, label="positive")

    axes[1,0].plot(t, df["neg_sel"], "#E24B4A", lw=1.5, label="negative")

    axes[1,0].plot(t, df["neglect"], "#1D9E75", lw=1.5, label="neglect")

    axes[1,0].set_xlabel("Time (h)"); axes[1,0].set_ylabel("Cumulative count")

    axes[1,0].set_title("Fate outcomes over time")

    axes[1,0].legend(fontsize=8)

  

    fc = model.fate_counts; total = max(sum(fc.values()), 1)

    colors = ["#378ADD", "#E24B4A", "#888780", "#1D9E75"]

    axes[1,1].bar(range(len(fc)), [100*v/total for v in fc.values()],

                  color=colors, width=0.5)

    axes[1,1].set_xticks(range(len(fc)))

    axes[1,1].set_xticklabels([k.replace("_"," ") for k in fc],

                                rotation=18, ha="right", fontsize=7)

    axes[1,1].set_ylabel("%"); axes[1,1].set_title("Final selection outcomes")

  

    plt.tight_layout()

    if save:

        plt.savefig(save, dpi=150, bbox_inches="tight")

    plt.show()

    plt.close(fig)

  

def contact_summary(model: ThymicEnvironment) -> None:

    for fate, records in model.contact_stats.items():

        if not records:

            print(f"{fate:22s}: no resolved cells yet")

            continue

        steps = [r["contact_steps"] for r in records]

        frac  = [r["contact_frac"]  for r in records]

        meanS = [r["mean_S"]        for r in records]

        print(f"{fate:22s} n={len(records):4d} | "

              f"⟨contact steps⟩={np.mean(steps):7.1f} | "

              f"⟨contact frac⟩={np.mean(frac):.3f} | "

              f"⟨S during contact⟩={np.mean(meanS):.2f}")

  

def plot_trajectories(model: ThymicEnvironment, save: str = "",

                       max_per_fate: int = 25) -> None:

    if not _MPL:

        return

  

    fig, axes = plt.subplots(1, 3, figsize=(10, 10), sharex=True, sharey=True)

    fates  = ["positive_selection", "negative_selection",

              "death_by_neglect"]

    colors = ["#378ADD", "#E24B4A", "#1D9E75"]

  

    stromal_xs, stromal_ys = zip(*model.stromal_pos) if model.stromal_pos else ([], [])

  

    for ax, fate, color in zip(axes.flat, fates, colors):

        ax.scatter(stromal_xs, stromal_ys, s=3, c="black",

                   alpha=0.6, zorder=0, label="stromal cells")

  

        trajs = model.trajectories[fate][:max_per_fate]

        for hist in trajs:

            xs = [h[1] for h in hist]

            ys = [h[2] for h in hist]

            ax.plot(xs, ys, color=color, alpha=0.2, lw=0.8, zorder=1)

            ax.scatter(xs[0], ys[0], color=color, s=10, marker="o", zorder=2)

            ax.scatter(xs[-1], ys[-1], color=color, s=10, marker="x", zorder=2)

  

        ax.set_xlim(0, model.p.grid_w); ax.set_ylim(0, model.p.grid_h)

        ax.set_title(f"{fate.replace('_',' ')}  (n={len(model.trajectories[fate])})")

        ax.set_xlabel("Width"); ax.set_ylabel("Height")

        if fate == fates[0]:

            ax.legend(fontsize=7, loc="upper right")

  

    fig.suptitle(f"Cell trajectories by fate — {model.p.scenario}", fontweight="bold")

    plt.tight_layout()

    if save:

        plt.savefig(save, dpi=150, bbox_inches="tight")

    plt.show()

    plt.close(fig)

  
  

# ══════════════════════════════════════════════════════════════════════

# 6.  RUNNER

# ══════════════════════════════════════════════════════════════════════

  

def run(

    scenario: str           = "uniform" or "cortex_stiffer" or "medulla_stiffer" ,

    n_steps:  Optional[int] = None,

    demo:     bool          = True,

    verbose:  bool          = True,

    seed:     int           = 42,

) -> ThymicEnvironment:

    """

    Run the thymic selection simulation.

  

    Parameters

    ----------

    scenario : "uniform" | "cortex_stiffer" | "medulla_stiffer"

    n_steps  : override total steps (default: T_run_hr × 60)

    demo     : if True, use 50×50 grid, 500 cells and lower run time/steps for speed

    verbose  : print progress regularly to see how mechanobiological mechanisms are operating

    """

    p = Params(scenario=scenario, seed=seed)

    if demo:

        p.grid_w = p.grid_h = 50

        p.T_run_hr = 20.0

        p.n_initial = 500

        n_steps = n_steps or int(p.T_run_hr * 60 / p.dt)  

  

    # Calibrate thresholds from OT-I data at biochemical baseline

    p.theta_pos, p.theta_neg = calibrate_thresholds(p)

  

    if verbose:

        print(f"\n{'═'*60}")

        print(f"  Thymic Selection Mechanobiology Model  |   Mesa {__import__('mesa').__version__}")

        print(f"{'═'*60}")

        print(f"  Scenario     : {p.scenario}")

        print(f"  Grid         : {p.grid_w}×{p.grid_h} ({'demo' if demo else 'full'})")

        print(f"  S_mean       = {p.S_mean}")

        print(f"  S_std        = {p.S_std}")

        print(f"  S_smooth     = {p.S_smooth}")

        print(f"  S_stiffer    = {p.S_stiffer}")

        print(f"  S_softer     = {p.S_softer}")

        print(f"  k_scan(C)    : {p.k_scan_base}·(1 + {p.gamma}·C)             ")

        print(f"  α(C,M)       : 1 + {p.alpha_C_max}·C + {p.alpha_M_max}·M           ")

        print(f"  θ_pos        = {p.theta_pos:.4e}                         ")

        print(f"  θ_neg        = {p.theta_neg:.4e}                         ")

        print(f"{'─'*60}")

  

    model   = ThymicEnvironment(p)

    n_steps = n_steps or int(p.T_run_hr * 60 / p.dt)

    report  = max(1, n_steps // 100)  

  

    for step in range(n_steps):

        model.step()

        if verbose and step % report == 0:

            alive = sum(1 for a in model.agents

                        if isinstance(a, Thymocyte) and a.alive)

            fc = model.fate_counts

            kscan = _mean_kscan(model)

            alpha = _mean_alpha(model)

            print(f"  t={model.t_min/60:5.1f}h | alive={alive:4d} | "

                  f"⟨k_scan⟩={kscan:.3f} | ⟨α⟩={alpha:.3f} | "

                  f"pos={fc['positive_selection']:3d} | "

                  f"neg={fc['negative_selection']:3d} | "

                  f"ngl={fc['death_by_neglect']:3d}")

  

    return model

  
  

# ══════════════════════════════════════════════════════════════════════

# 7.  ENTRY POINT

# ══════════════════════════════════════════════════════════════════════

  

if __name__ == "__main__":

    # ── Print calibrated thresholds ───────────────────────────────────

    p_cal = Params()

    tp, tn = calibrate_thresholds(p_cal)

    print(f"\nCalibrated thresholds  (biochemical baseline):")

    print(f"  θ_pos = {tp:.4e}   (neglect / positive selection)")

    print(f"  θ_neg = {tn:.4e}   (positive / negative selection)")

    print(f"  Window width: {np.log10(tn/tp):.2f} log₁₀ decades")

  

    # ── Run ──────────────────────────────────────────────────────

    model = run(scenario="uniform",demo=False,verbose=True)

    print_summary(model)

    koff_summary(model)

    contact_summary(model)

    plot_stiffness_map(model)

    plot_trajectories(model)

    plot_dynamics(model)

    model.stiffness_report(expected_range=(12.0, 18.0), quartile_range=(12.0, 18.0))