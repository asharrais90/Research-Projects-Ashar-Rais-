
### Some relevant data from the literature: 
- Mechanosensing in the immune response (A. Upadhyaya, 2017):
	- The response of T cells to substrate stiffness [19] appears to obey predictions from a simple model derived from active matter theory [92]. The cell exerts an active stress (a) resulting from actin polymerization and myosin contractility, which consumes energy by ATP hydrolysis. Force balance leads to a simple expression for the force, Fss = Fsat ksubs ksubs+kcell , where Fss is the steady state force and Fsat is the saturating force. Fsat = a A ,where A is the cell area, ksubs is the effective spring constant for the substrate (0.1–10 nN/um) and kcell is the stiffness of the cell (∼1 nN/um). Thus, on compliant substrates, the steady state force exerted by the cell is predicted to linearly increase with stiffness and saturate when ksubs  kcell. Assuming a linear force-velocity relationship, the loading rate can be estimated as dF/dt ∼(A x sigma_sat / n) x ksubs, where n = tau (k_subs + k_cell) is the viscous dissipation in the actin gel and A is the area over which the stress is exerted. Using the estimated values for these parameters (tau ∼ 10s, ksubs ∼ kcell), a um2 patch of TCRs would experience loading rates of ∼2-3 pN/sec, which is similar to the force probe experiments [74].
- Cytoskeletal forces during signaling activation in {Jurkat} {T}-cells (Hui et.al., 2014):
	- The response of Jurkat cells to substrate stiffness can be explained using a model derived from active matter theory (Marcq et al., 2011). In a simple one-dimensional version of this model, the extracellular environment can be represented by a linear spring of length lsubs(t) at time t, rest length l0,subs(t), and spring constant ksubs. The cell exerts an active stress (σa) resulting from actin polymerization and myosin contractility, which consumes energy by ATP hydrolysis. Force balance leads to a simple expression for the force, F F k k k eq sat subs subs cell = + , where Feq is the steady-state force; F A sat a = σ, where A is the cell area; ksubs is the effective spring constant for the substrate (0.1–10 nN/μm for the gel stiffnesses considered here); and kcell is the stiffness of the cell (∼1 nN/μm). Thus, on soft substrates, the steady-state force exerted by the cell linearly increases with stiffness and saturates when k k ext c >> ell. Given our measured saturating force of 5 nN and the typical cell area A ≈ 100 μm2, we can estimate the peak active stress, σ ≈ active 50 Pa, which is in agreement with the measured Young’s moduli of Jurkat cells (Rosenbluth et al., 2006). The same model can also be used to estimate the maximal strain rate that can be exerted by these cells, which will pull on TCR–ligand bonds. Assuming a linear force–velocity relationship (appropriate for overdamped systems), the loading rate can be estimated as dF dt A k sat ≈ subs σ η , where η = τ + ( ) k k subs cell is the viscous dissipation in the actin gel and A is the area over which the stress is exerted. Again, using our estimated values for these parameters ( 1 τ ≈ 0 s, ) k k subs ≈ cell , a typical bead (area, ∼1 μm2) would experience a loading rate of ∼2–3 pN/s, which is well in agreement with observed rates experienced by anti-CD3–coated beads in a force sensor apparatus (Husson et al., 2011).





# 1. ABM parameters

|Parameter|Suggested value|Units|Source/justification|
|---|---|---|---|
|Grid size|200 × 200|cells|computational choice|
|Cell diameter|8–10|μm|thymocyte diameter|
|Time step ((\Delta t))|30–60|s|typical migration simulations|
|Maximum residence time|48–96|h|thymocyte positive selection window (mouse estimates)|
|Initial thymocyte number|500–5000|agents|depends on computational cost|

---

# 2. Biphasic migration

## Equation

[  
v(S)=  
v_{opt}  
\frac{(S/K_{low})^{n_{low}}}  
{1+(S/K_{low})^{n_{low}}}  
\frac{1}{1+(S/K_{high})^{n_{high}}}  
]

---

|Parameter|Suggested value|Units|Justification|
|---|---|---|---|
|(v_{opt})|10–15|μm/min|intrathymic migration imaging|
|(K_{low})|0.5–1|kPa|onset of traction increase|
|(K_{high})|10–20|kPa|saturation/stiffening regime|
|(n_{low})|2|—|typical Hill slope|
|(n_{high})|2–3|—|phenomenological|

These values are broadly consistent with mechanobiology migration studies and T-cell traction experiments. ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7041673/?utm_source=chatgpt.com "Multiple actin networks coordinate mechanotransduction at the immunological synapse - PMC"))

---

# 3. Cytoskeletal priming ODE

## Equation

# [  
\frac{dC}{dt}

k_C(C_{eq}(S)-C)  
]

---

## Equilibrium

# [  
C_{eq}(S)

\frac{S^{n_C}}  
{K_C^{n_C}+S^{n_C}}  
]

---

|Parameter|Suggested value|
|---|---|
|(K_C)|1–3 kPa|
|(n_C)|2|
|(k_C)|0.02–0.10 min(^{-1})|

---

Interpretation:

[  
\tau_C=\frac1{k_C}  
]

gives adaptation times

```
10–50 minutes
```

which are very reasonable for actin remodeling and focal adhesion maturation.

Comparable adaptation times have been reported for mechanotransduction and YAP/TAZ nuclear responses. ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7041673/?utm_source=chatgpt.com "Multiple actin networks coordinate mechanotransduction at the immunological synapse - PMC"))

---

# 4. Amplification

## Equation

[  
\alpha(C)=1+\alpha_{max}C  
]

---

|Parameter|Suggested value|
|---|---|
|(\alpha_{max})|0.5–2|

Meaning:

```
no priming → amplification = 1

full priming → amplification = 1.5–3
```

This is phenomenological and should be explored in sensitivity analysis.

---

# 5. Force generation

## Equation

# [  
F(S)

F_{max}  
\frac{S^{n_F}}  
{K_F^{n_F}+S^{n_F}}  
]

---

|Parameter|Suggested value|Units|
|---|---|---|
|(F_{max})|10–20|pN|
|(K_F)|2–5|kPa|
|(n_F)|2||

---

Why?

Single TCRs typically experience

```
5–20 pN
```

during activation.

Catch bonds are observed around

```
~10 pN
```

of applied force. ([Rockefeller University Press](https://rupress.org/jcb/article/225/7/e202601058/282591/T-cell-mechanobiology-How-molecular-forces-shape?utm_source=chatgpt.com "T-cell mechanobiology: How molecular forces shape immune function | Journal of Cell Biology | Rockefeller University Press"))

---

# 6. Cytoskeletal enhancement

## Equation

# [  
F_{eff}

F(S)(1+\beta C)  
]

---

|Parameter|Suggested value|
|---|---|
|(\beta)|0.2–1|

Meaning:

fully primed cells exert

```
20–100% more effective force
```

than unprimed cells.

This is currently a calibration parameter.

---

# 7. Catch-slip bond

## Equation

# [  
k_{off}

k_c^0e^{-F/F_c}  
+  
k_s^0e^{F/F_s}  
]

---

|Parameter|Suggested value|
|---|---|
|(k_c^0)|0.5–2 s(^{-1})|
|(k_s^0)|0.01–0.1 s(^{-1})|
|(F_c)|5–10 pN|
|(F_s)|10–20 pN|

These ranges are consistent with published catch-slip descriptions for TCR–pMHC interactions. ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7041673/?utm_source=chatgpt.com "Multiple actin networks coordinate mechanotransduction at the immunological synapse - PMC"))

---

# 8. Kinetic proofreading

## Equation

# [  
KP

\left(  
\frac{k_p}  
{k_p+k_{off}+k_{shp1}}  
\right)^N  
]

---

|Parameter|Suggested value|
|---|---|
|(N)|3–8|
|(k_p)|1–10 s(^{-1})|
|(k_{shp1})|0.1–1 s(^{-1})|

---

## Why?

Recent analyses suggest that the effective proofreading delay is only a few seconds and that only a modest number of sequential modifications are required. ([eLife](https://elifesciences.org/articles/67092?utm_source=chatgpt.com "The discriminatory power of the T cell receptor | eLife"))

---

# 9. TCR signalling ODE

# [  
\frac{dX}{dt}

## E(t)\alpha(C)KP(F_{eff})

\lambda_X X  
]

---

## Signal decay

|Parameter|Suggested value|
|---|---|
|(\lambda_X)|0.005–0.05 min(^{-1})|

---

Signal half-life

# [  
t_{1/2}

\frac{\ln2}{\lambda_X}  
]

becomes

```
15–140 minutes
```

which is compatible with transient phosphorylation/NFAT activation dynamics.

---

# 10. Mechanosurvival variable (recommended separate from C)

# [  
\frac{dM}{dt}

k_M(M_{eq}(S)-M)  
]

---

|Parameter|Suggested value|
|---|---|
|(k_M)|0.01–0.05 min(^{-1})|
|(K_M)|2–5 kPa|
|(n_M)|2|

These are reasonable estimates for slower transcriptional responses (YAP/TAZ, survival programs).

---

# 11. Selection thresholds

These should **not** be taken from the literature as absolute numbers because (X) is a dimensionless model variable.

Instead normalize

[  
0\le X\le1.  
]

Then use

|Fate|Suggested threshold|
|---|---|
|Death by neglect|(X<0.25)|
|CD8|(0.25\le X<0.45)|
|CD4 (ThPOK-high)|(0.45\le X<0.70)|
|Treg (optional high-intermediate window)|(0.70\le X<0.85)|
|Negative selection|(X\ge0.85)|

These are **calibration thresholds** rather than experimentally measured constants.

---

# My strongest recommendation

I would **not present these as fixed values** in a thesis or paper. Instead, classify them as:

|Category|Parameters|
|---|---|
|**Measured from literature**|(v_{opt}, F_{max}, F_c, F_s, N)|
|**Estimated from experimental ranges**|(K_C, K_F, K_M, k_C, k_M, \lambda_X)|
|**Calibrated by sensitivity analysis**|(\alpha_{max}, \beta), lineage thresholds|

This is exactly how most hybrid ABM–ODE immunology models are parameterized. In fact, I would go one step further and perform a **global sensitivity analysis (e.g., Latin Hypercube Sampling + Sobol indices or Morris screening)** on the calibration parameters. That approach is generally viewed more favorably by reviewers than claiming precise values for quantities that have never been directly measured in thymocytes.