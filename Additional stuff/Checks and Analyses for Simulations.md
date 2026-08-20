
### Claudes take on this:
These analyses fall into several distinct families, each addressing a different type of question about a model's behaviour. The family that stability analysis belongs to — questions about the qualitative structure of a dynamical system — is only one of six or seven independent families that a computational biologist applies routinely. Here is a systematic account of all of them.

---

## Before writing any equations: order-of-magnitude and dimensional analysis

The cheapest analysis in computational biology costs nothing computationally and happens before any code is written. For every quantity in your model, you should be able to estimate its order of magnitude from first principles and check that it is consistent with the biology.

The specific technique is Fermi estimation applied to biological parameters. For the thymocyte model: a thymocyte has diameter ~8 μm and migrates at ~10–25 μm/min, so it crosses its own body length in under a minute. Contact duration with a stromal cell should therefore be on the order of minutes to tens of minutes, not seconds or hours. If your model produced mean contact durations of 0.1 seconds or 3 hours, that would signal a parameter error before you ran a single simulation.

Non-dimensionalisation is the formal version of this. You rescale every variable in your equations by a characteristic scale, producing dimensionless groups that reveal what actually controls the behaviour. For the X accumulation ODE:

```
dX/dt = k_scan · α · KP - λ_X · X
```

Rescaling by `X̃ = X · λ_X / (k_scan · α_ref · KP_ref)` and `t̃ = λ_X · t` gives:

```
dX̃/dt̃ = KP/KP_ref - X̃
```

The entire dynamics now depends only on the ratio `KP/KP_ref` — how far the current productive probability is from the reference. This tells you immediately that what matters is not `k_scan`, `α`, or `KP` individually, but their product relative to `λ_X`. Three parameters collapse into one dimensionless group. The Buckingham Pi theorem formalises this: a system with n variables and k independent physical dimensions has `n - k` independent dimensionless groups that fully characterise the behaviour. Finding these groups before parameterising is standard practice because it tells you exactly how many independent things you need to measure.

---

## Fixed-point and qualitative structure analysis

This is the family the stability analysis belongs to. Several techniques within this family go beyond the linearisation check.

**Nullcline analysis** applies to two-variable systems and is the standard first tool for understanding phase plane geometry. The nullcline of variable X is the curve in (X, Y) space where `dX/dt = 0`. The intersection of the X-nullcline and Y-nullcline gives the fixed point. The shape of the nullclines determines the qualitative structure — whether the fixed point is stable, unstable, a saddle, or whether the system produces limit cycle oscillations. For the two-variable subsystem `(C, X)` in the thymocyte model, you can plot the C-nullcline `C = Ceq(S)` (a horizontal line) and the X-nullcline `X = k_scan · α(C) · KP / λ_X` (a curve in C-X space) and immediately see whether they intersect once or multiple times. Multiple intersections would signal bistability — a regime where the same parameters allow two stable steady states, which is biologically relevant as a switch.

**Phase plane trajectories** generalise this. For a two-variable ODE system, you can draw the vector field `(dX/dt, dC/dt)` as arrows at every point in the phase plane, revealing the direction of flow everywhere. Limit cycles appear as closed orbits in the phase plane. Separatrices — invariant manifolds that divide the phase plane into regions with qualitatively different trajectories — identify which initial conditions lead to which fate. In the thymocyte model, understanding the phase plane geometry of (C, X) under alternating contact (E=1) and non-contact (E=0) conditions would tell you whether transient contacts can push a cell across a selection threshold from which it cannot return.

**Bifurcation analysis** asks how the number and stability of fixed points changes as a parameter varies continuously. A saddle-node bifurcation occurs when two fixed points (one stable, one unstable) collide and annihilate as a parameter crosses a critical value. A Hopf bifurcation is where a stable fixed point becomes unstable and a limit cycle is born. These bifurcation points are where qualitative behaviour changes — selection criteria that are sharp on one side of a parameter value become permissive on the other. Standard tools are AUTO, XPPAUT, or MATCONT for tracing bifurcation diagrams, though the analysis can be done manually for simple systems by asking where the Jacobian eigenvalues change sign.

For the thymocyte model, the biologically interesting bifurcation question is: is there a critical stiffness `S*` above which the model switches from graded selection (fraction of cells selected varies smoothly with stiffness) to threshold selection (most or all cells above a stiffness threshold are selected)? This would appear as a saddle-node bifurcation in the parameter S.

---

## Timescale analysis and model reduction

This family of analyses identifies which parts of a model can be simplified without loss of qualitative behaviour.

**Quasi-steady-state approximation (QSSA)** is used when one process is much faster than everything else. If species Y equilibrates rapidly relative to species X, then Y is approximately always at its steady state conditional on the current X, and you can eliminate the ODE for Y by setting `dY/dt = 0` and solving for `Y_ss(X)`. The Michaelis-Menten enzyme kinetics approximation is the canonical example — the enzyme-substrate complex is assumed to be at quasi-steady state relative to the substrate. In the thymocyte model, the kinetic proofreading cascade was treated this way: because `k_SHP1 = 1 s⁻¹` and the ODE timestep is 60 s, the proofreading cascade equilibrates within the timestep and KP is computed as a static probability rather than integrated dynamically. The validity criterion is `k_SHP1 >> lambda_X / 60` in consistent units — checking this explicitly is the standard sanity check for any QSSA application.

**Timescale separation** is the generalisation. You order all rate constants by magnitude and group them into fast, intermediate, and slow. The fast subsystem can be collapsed to its steady state (QSSA). The slow subsystem evolves on a timescale where the fast variables are always equilibrated. For a system with three timescales (seconds, minutes, hours), this produces a hierarchy of approximations: a fast-collapse approximation, a slow-manifold approximation, and the full system. The Fenichel slow manifold theorem formalises when this reduction is valid. In practice, checking timescale separation means listing all `1/rate_constant` values and verifying they span at least a factor of 10 between adjacent groups.

**Lumping** is a different reduction technique that combines multiple species into aggregated variables when they are kinetically similar. If three intermediate proofreading states have similar equilibration rates, they can be lumped into a single effective state with an effective rate constant derived from the original ones. Exact lumping preserves the dynamics precisely; approximate lumping introduces controlled error. The standard test is to compare the full and lumped model over the parameter range of interest and verify that the key output statistics agree within acceptable tolerance.

---

## Sensitivity analysis

Sensitivity analysis is the family of methods that quantifies which parameters most strongly affect model outputs. It is one of the most important and most commonly underdone analyses in published computational biology models.

**Local sensitivity analysis** computes the partial derivative of an output with respect to each parameter at a specific parameter point. The normalised sensitivity coefficient `S_ij = (∂Y_i/∂p_j) × (p_j/Y_i)` gives the fractional change in output i per fractional change in parameter j, which is interpretable regardless of parameter units. Computing these analytically is ideal; finite difference approximations — `(Y(p_j + δ) - Y(p_j - δ)) / (2δ)` — work when analytic derivatives are unavailable. For the thymocyte model, the local sensitivity of the positive selection fraction to `k_SHP1` tells you how much the selection outcome changes if SHP-1 activity is 10% higher or lower.

The limitation is that local sensitivity is only valid near the point where it is computed. If the model is nonlinear (which this one is, through the N-th power in KP), sensitivities can change dramatically across the parameter space.

**Global sensitivity analysis** samples the parameter space broadly and decomposes output variance across parameters. The two standard approaches are Morris screening and Sobol indices.

Morris screening evaluates the model at a sparse random sample of points and computes an elementary effect `EE_i = (Y(x + Δe_i) - Y(x)) / Δ` for each parameter i, where `e_i` is the unit vector and Δ is a step size. Averaging over many such evaluations gives a mean absolute elementary effect (μ*) and standard deviation (σ). High μ* means the parameter has large average influence. High σ means the parameter's influence is highly nonlinear or interaction-dependent. This is a screening method — cheap, O(k(r+1)) model evaluations for k parameters and r trajectories — that identifies which parameters deserve full variance decomposition.

Sobol variance decomposition computes first-order indices `S_i = Var_i(E_{~i}(Y | x_i)) / Var(Y)` and total-order indices `S_Ti = 1 - Var_{~i}(E_i(Y | x_{~i})) / Var(Y)`. The first-order index measures the fraction of output variance attributable to parameter i alone. The total-order index includes all higher-order interactions involving parameter i. Parameters where `S_Ti >> S_i` have strong interaction effects — they matter primarily through their co-variation with other parameters rather than individually. For the thymocyte model, the difference between F_half_koff and F_half_shp1 sensitivity indices would reveal whether the catch-bond or SHP-1 exclusion channel is the dominant mechanical mechanism.

**Partial rank correlation coefficient (PRCC)** is a sampling-based alternative to Sobol indices that is computationally cheaper but assumes monotonic input-output relationships. You run Latin hypercube sampling across the parameter space, rank-transform all inputs and outputs, and compute the partial correlation between each input and the output after removing the linear effects of all other inputs. PRCC values above 0.5 in absolute value are conventionally considered significant. It is standard practice in epidemiological models because it is fast and interpretable.

---

## Identifiability analysis

Identifiability analysis asks a question that is logically prior to parameter estimation: can the parameters in your model in principle be inferred from your observable outputs, even with infinite data?

**Structural identifiability** is a mathematical property of the model equations and the observation mapping, independent of any specific dataset. A model is structurally identifiable if there is a unique parameter vector consistent with the observable output trajectory. If two different parameter vectors produce identical outputs for all time and all inputs, the model is structurally unidentifiable — you cannot distinguish those parameter values from observations, no matter how much data you collect.

The differential algebra approach (implemented in tools like DAISY, SIAN, or StructuralIdentifiability.jl) computes the input-output relation — the algebraic relationship between observable outputs and model inputs, with parameters appearing as coefficients. If every coefficient is uniquely determined by the input-output relation, the model is globally identifiable. If coefficients appear only as products or sums with other coefficients, those combinations are identifiable but the individual parameters are not.

For the thymocyte model, a specific structural identifiability concern is the pair `(k_scan, KP_ref)` in the threshold calibration. The production term `k_scan × alpha × KP` appears always as a product. If your observable is only the selection fraction at the end of the simulation, you cannot separately identify `k_scan` and the overall scale of KP — only their product matters. This would need to be either fixed by measuring one of them independently, or acknowledged as a non-identifiable combination whose product is the relevant biological quantity.

**Practical identifiability** goes beyond structural identifiability to ask whether parameters can be estimated from a specific finite dataset with specific noise levels. Even a structurally identifiable model may be practically unidentifiable if the data are too sparse or noisy to constrain parameter values. The standard diagnostic is **profile likelihood**: you fix one parameter at a range of values, optimise all other parameters to fit the data at each fixed value, and plot the resulting likelihood profile. A parameter is practically identifiable if the profile likelihood has a clear minimum and confidence interval. A flat profile indicates practical non-identifiability — the data have no information about that parameter. For the thymocyte model, the profile likelihood of `alpha_max` would tell you whether any observable selection outcome statistic could distinguish `alpha_max = 0.5` from `alpha_max = 1.5`.

**Sloppiness analysis** is the eigenvalue perspective on identifiability. The Fisher information matrix (FIM) for a model at a specific parameter point has eigenvalues that span many orders of magnitude in typical biological models — often 10 to 12 orders of magnitude. The eigenvectors corresponding to large eigenvalues are stiff parameter combinations that are well-constrained by data. The eigenvectors corresponding to small eigenvalues are sloppy combinations that can vary widely without affecting the output. Computing the eigenvalue spectrum of the FIM before fitting reveals how many parameter combinations are actually constrained, regardless of how many individual parameters exist. For the thymocyte model, this analysis would tell you that even though you have eight mechanical parameters, perhaps only three or four linear combinations of them are constrained by the observable selection fractions.

---

## Conservation law and constraint checking

Any physically or biologically derived constraint that your model must satisfy is a free consistency check at zero computational cost.

**Conservation laws** in biochemical models include mass conservation (total protein concentration is conserved if there is no synthesis or degradation), moiety conservation (the sum of free ATP and all ATP-containing complexes is constant), and charge conservation in electrochemical models. For every ODE system describing biochemical reactions, you can derive the conservation laws algebraically by finding the left null space of the stoichiometric matrix. If the model violates a conservation law at any timepoint, there is a bug.

For the thymocyte model, the ODE states C and M are bounded between 0 and 1 by their equilibrium expressions — Ceq and Meq are both in [0,1], and the first-order relaxation dynamics ensure C and M stay in [0,1] if they start there. This can be verified analytically: `dC/dt = k_C(Ceq - C)`. At C = 1: `dC/dt = k_C(Ceq - 1) ≤ 0` since Ceq ≤ 1. At C = 0: `dC/dt = k_C × Ceq ≥ 0`. So C cannot leave [0,1]. The `np.clip` in the code is therefore defensive rather than necessary for the nominal parameter set — and knowing that analytically means you can remove it for a production version and trust that boundary violations are signals of numerical errors.

**Probability conservation** is the equivalent for stochastic models. If your model represents a probability distribution over states, the probabilities must sum to 1 at every timepoint. For the kinetic proofreading model viewed as a Markov chain, the probability of being in each state of the cascade must sum to 1. If instead of computing KP as a closed-form expression you were integrating the master equation for the proofreading cascade directly, you would check that `sum(p_i) = 1` at every step as a debugging tool.

**Steady-state flux balance** applies to metabolic models but has analogues in signalling models. At steady state, every internal node of a network must have zero net flux — production equals consumption. For the X variable at steady state with E=1: `k_scan × alpha × KP = lambda_X × X_ss`, which gives the closed-form expression used in threshold calibration. Checking this algebraic relationship against the numerical steady state reached by simulation verifies that the integration is working correctly.

---

## Numerical method validation

These are checks on the computational machinery rather than the biological model.

**Timestep convergence testing** is the most fundamental numerical check. You run the same simulation with several values of Δt (e.g., 1 min, 0.5 min, 0.25 min, 0.1 min) and examine how the output statistics change. If the selection fractions are converging as Δt decreases — the differences between consecutive runs shrink systematically — you have evidence that the Δt you chose is within the convergent regime. If the outputs are still changing substantially at your chosen Δt, you need a smaller step or a higher-order method. This should be reported alongside any quantitative results from the model.

**Grid resolution convergence** is the spatial analogue for ABMs with continuous spatial environments. For the thymocyte model's discrete 200×200 grid, the relevant check is whether results change substantially if you use a 100×100 or 400×400 grid. Because cells are point particles on a lattice here rather than objects with volume, the grid resolution determines the spatial precision of migration paths and contact detection. A result that changes significantly between 100×100 and 200×200 grids but not between 200×200 and 400×400 indicates that 200×200 is sufficient.

**Mass balance checking** at runtime means periodically verifying that the number of agents tracked by `model.agents` equals the number of alive agents on the grid. Any discrepancy indicates a bookkeeping error — an agent that was removed from the grid but not deregistered, or vice versa. This is particularly important in Mesa 3.x because the two-step removal (grid removal and agent deregistration) creates a window where the count can be temporarily inconsistent.

**Identical seed reproducibility** is a basic sanity check: run the simulation twice with identical parameters and seed, and verify that every output is bit-for-bit identical. This confirms that there are no hidden sources of randomness (unset seeds, time-dependent operations, hash randomisation in Python dicts). The check is trivial to implement and immediately catches a class of reproducibility bugs that are otherwise extremely hard to diagnose.

---

## Stochastic analysis

When a model is stochastic — which ABMs inherently are — the statistics of the stochastic output require specific analytical tools beyond what applies to deterministic models.

**Coefficient of variation (CV)** is `σ/μ` — the standard deviation divided by the mean, which normalises variability for comparison across different scales. For the selection fraction across multiple simulation seeds, CV tells you how much of the variability in observed outcomes is due to stochastic fluctuations versus genuine parametric differences. A CV above ~0.1 suggests that your results are substantially seed-dependent and you need more runs for reliable statistics.

**Fano factor** is `σ²/μ` — the variance divided by the mean. For a Poisson process, the Fano factor equals 1. Values above 1 indicate super-Poissonian noise (more variability than a Poisson process), often reflecting positive feedback or burst production. Values below 1 indicate sub-Poissonian noise, often reflecting negative feedback or precise regulatory mechanisms. For the X accumulation process, computing the Fano factor of the signal distribution across cells at a given timepoint tells you whether the signalling variability is consistent with independent TCR encounters (Poisson, Fano ≈ 1) or shows clustering effects.

**Linear noise approximation (LNA)** is a systematic method for propagating stochastic fluctuations through a nonlinear ODE system. You write `y = y_det + ε·η` where `y_det` is the deterministic trajectory, ε is a small noise parameter, and η is a fluctuation. Expanding the master equation to second order in ε produces a linear ODE for the variance of η whose coefficients come from the Jacobian of the deterministic system. This gives you an analytical approximation for the variance of any model variable as a function of time, without running many stochastic simulations. For the thymocyte model, LNA applied to the X variable would predict how the variance in accumulated signal across cells grows over time as a function of `k_scan` and `lambda_X`.

**Ensemble statistics and their convergence** require careful handling. For an ABM with N_run independent runs, the standard error on the mean selection fraction scales as `σ/√N_run`. To achieve a standard error below 1% of the mean, you need approximately `(σ/0.01μ)² = (CV/0.01)²` runs. If CV ≈ 0.05, this requires 25 runs — manageable. If CV ≈ 0.5, this requires 2500 runs — potentially prohibitive. Computing the CV from a pilot run of 10 seeds before committing to a full parameter sweep is therefore standard practice.

---

## Network topology analysis

For models defined on biological networks — signalling pathways, gene regulatory networks, metabolic networks — there is a family of graph-theoretic analyses that apply before any dynamics are simulated.

**Feedback loop enumeration** identifies all positive and negative feedback loops in a network. A positive feedback loop (even number of negative edges) is a candidate mechanism for bistability and switch-like behaviour. A negative feedback loop (odd number of negative edges) is a candidate for adaptation, homeostasis, or oscillation. In the thymocyte model, the question of whether X accumulation has any feedback onto the proofreading machinery — whether productive signalling changes k_p, k_SHP1, or koff — is a network topology question. The current model has no such feedback (by design), and identifying that absence is itself a structural claim that should be stated explicitly.

**Network motif analysis** asks which small subgraph patterns appear more frequently in a network than in a random network of the same degree distribution. The coherent feedforward loop (activation of Z by both X directly and X via Y) produces adaptation — a transient response that returns to baseline even with persistent input. The incoherent feedforward loop produces pulse-like dynamics. In thymic signalling, the relationship between TCR signal, ZAP70 activation, and ITAM phosphorylation is a feedforward structure whose coherence or incoherence determines whether signalling adapts to persistent antigen exposure — a biologically important question for selection threshold stability.

**Structural controllability and observability** are linear algebra properties of the network that determine which states can be driven by inputs (controllability) and which states produce observable outputs (observability). For a model of thymic selection, observability analysis asks: given that you can only measure fate outcomes (positive/negative/neglect/rescued), which combinations of internal state variables (C, M, X, koff_base) are in principle recoverable from those measurements? A state that is unobservable cannot be inferred from the data no matter how many cells you measure.

---

## Model comparison and hypothesis discrimination

These analyses address the question of whether your model is the right one, not just whether it is self-consistent.

**Akaike Information Criterion (AIC)** and **Bayesian Information Criterion (BIC)** compare models with different numbers of parameters fitted to the same data. Both penalise parameter count: `AIC = 2k - 2ln(L)` and `BIC = k·ln(n) - 2ln(L)` where k is the number of parameters, n is the number of data points, and L is the maximum likelihood. The model with lower AIC or BIC is preferred. BIC penalises parameters more heavily than AIC and tends to select simpler models. For the thymocyte model, you would compare the full catch-slip bond model against a simpler version where koff is constant (no force dependence) — if the AIC difference is less than 2, the simpler model is statistically comparable.

**Likelihood ratio tests** provide a formal significance test for nested models. If model B is model A with some parameters fixed to specific values, twice the difference in log-likelihood between them is chi-squared distributed with degrees of freedom equal to the number of fixed parameters. This tells you whether adding the catch-slip bond formulation provides a statistically significant improvement in fit over a constant koff assumption.

**Leave-one-out cross-validation** evaluates predictive performance by training the model on all but one data point, predicting the withheld point, and averaging the prediction error across all withheld points. This directly estimates out-of-sample prediction accuracy rather than in-sample fit. For the thymocyte model, this would involve fitting parameters to selection fractions from one set of experimental conditions and predicting fractions under held-out conditions — which is exactly the kind of test that distinguishes a Ptolemaic from a Keplerian model.

**Bayesian model comparison** via the Bayes factor `B = P(data | M1) / P(data | M2)` compares the total probability of the data under each model, integrating over all parameter values weighted by their prior probability. This naturally penalises complexity without requiring explicit parameter counting, because models with more parameters have their likelihood spread over a larger parameter space and the integral is correspondingly smaller. A Bayes factor above 10 is conventionally considered strong evidence for the favoured model. Computing Bayes factors exactly requires integration over the full parameter space, which is generally intractable and typically approximated by thermodynamic integration or nested sampling.

---

## Robustness analysis

Robustness analysis asks whether a model's qualitative behaviour is maintained under perturbations, which is different from sensitivity analysis (which asks how much output values change) and different from identifiability (which asks whether parameters can be inferred).

**Parameter robustness** asks whether the qualitative outcome — not the exact quantitative value, but the qualitative feature such as "positive selection fraction is higher in cortex_stiffer than uniform scenario" — is preserved across the entire biologically plausible parameter range. A claim is robust if it holds across this entire range, not just at the point estimates. This is directly related to the earlier discussion of Ptolemaic models: a claim that depends on being at exactly the right parameter values is an epicycle; a claim that holds across an order of magnitude in each relevant parameter is a genuine prediction.

**Structural robustness** asks whether the behaviour is preserved under small changes to the model structure — adding or removing a term, changing a Hill function to a Michaelis-Menten function, adding a delay. A model whose predictions are sensitive to these structural assumptions is making claims that are not justified by the available data about model structure. This is particularly relevant for the thymocyte model's dual role of C, which was a structural assumption. If removing one of the two C channels substantially changes the selection fractions, that is evidence that the model's claims about the mechanism are sensitive to a structural assumption rather than robust to it.

**External perturbation analysis** asks what the model predicts when you simulate a specific experimental perturbation — drug treatment, genetic knockout, environmental change. This is the primary way a model generates testable predictions. For the thymocyte model, simulating blebbistatin treatment (which inhibits myosin II and would set beta = 0, eliminating Path A) and Y-27632 treatment (which inhibits ROCK and would reduce F_max) generates specific quantitative predictions about how each drug would shift the selection fractions. These predictions are directly testable in thymic organoid experiments.

---

## Putting it together as a workflow

The standard workflow in computational biology applies these analyses in a specific order that is not arbitrary. You start with dimensional analysis and order-of-magnitude estimation before writing any equations, because these determine whether your proposed model equations are plausible at all. You then analyse the equation structure — fixed points, nullclines, timescale separation — before parameterising, because this tells you whether your model is capable of producing the qualitative behaviour you observe before you spend any effort on parameter values. You then perform structural identifiability analysis before fitting parameters, because there is no point fitting parameters that cannot in principle be determined from your data. You then fit parameters and perform local sensitivity analysis to identify the most important ones for validation experiments. You then perform global sensitivity analysis and robustness analysis to determine which claims are worth making. Finally, you compare models using information criteria or Bayes factors to determine which formulation is best supported by the available evidence.

Each of these analyses is cheap relative to running the full simulation and tells you something about a different aspect of the model. Together they constitute the difference between a model that is a sophisticated interpolation scheme and a model that is a scientific instrument.