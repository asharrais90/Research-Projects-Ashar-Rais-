---
share: true
publish: true
---

## Obtaining k<sub>off</sub> values for the distribution:
Experimental TCR k<sub>off</sub> measurements (given in the table below) were used to fit a log-normal probability distribution. The distribution was parameterized by a location parameter $\mu$<sub>log</sub> and dispersion parameter $\sigma$<sub>log</sub>. A shift-factor was used to translate the distribution along the k<sub>off</sub> axis, while sigma-factor controlled repertoire heterogeneity. Individual thymocytes were then assigned stochastic k<sub>off</sub> values by sampling from this adjusted distribution using a reproducible random seed.

| peptide                    | koff value (s<sup>-1</sup>)          | TCR                  | pMHC                                  | species         | Assay                                                                 | Reference                                                                                                                                               |
| -------------------------- | ------------------------------------ | -------------------- | ------------------------------------- | --------------- | --------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| MCC                        | 0.057                                | 2B4 T-cell hybridoma | class II MHC molecule I-E<sup>k</sup> | mouse           | SPR (with amine-coupled TCRs)                                         | Kinetics of T-cell receptor binding to peptide/I-Ek complexes: Correlation of the dissociation rate with T-cell responsiveness (K. Matsui et.al., 1994) |
| PCC                        | 0.090                                | {same as above}      | {same as above}                       | {same as above} | {same as above}                                                       | {same as above}                                                                                                                                         |
| MCC (102S)                 | 0.1-0.3                              | {same as above}      | {same as above}                       | {same as above} | {same as above}                                                       | {same as above}                                                                                                                                         |
| MCC (agonist)              | 0.063 ± 0.009                        | 2B4 T-cell hybridoma | class II MHC molecule I-E<sup>k</sup> | mouse           | SPR (with cysteine-coupled TCRs)                                      | A TCR Binds to Antagonist Ligands with Lower Affinities and Faster Dissociation Rates Than to Agonists (D. lyons et.al., 1996)                          |
| MCC (T102S) (weak agonist) | 0.36 ± 0.018                         | {same as above}      | {same as above}                       | {same as above} | {same as above}                                                       | {same as above}                                                                                                                                         |
| MCC (T102N) (weak agonist) | 0.44 ± 0.020                         | {same as above}      | {same as above}                       | {same as above} | {same as above}                                                       | {same as above}                                                                                                                                         |
| OVA                        | 0.020 0.001                          | OT-I system          | H-2K<sup>b</sup> (MHC class I)        | mouse           | SPR                                                                   | T-cell receptor affinity and thymocyte positive selection (S.M. Alam et.al., 1996)                                                                      |
| E1                         | 0.068 ± 0.004                        | {same as above}      | {same as above}                       | {same as above} | {same as above}                                                       | {same as above}                                                                                                                                         |
| V-OVA                      | 0.039 ± 0005                         | {same as above}      | {same as above}                       | {same as above} | {same as above}                                                       | {same as above}                                                                                                                                         |
| R4                         | 0.146 ± 0.012                        | {same as above}      | {same as above}                       | {same as above} | {same as above}                                                       | {same as above}                                                                                                                                         |
| K4                         | >0.2                                 | {same as above}      | {same as above}                       | {same as above} | {same as above}                                                       | {same as above}                                                                                                                                         |
| VSV                        | >0.6                                 | {same as above}      | {same as above}                       | {same as above} | {same as above}                                                       | {same as above}                                                                                                                                         |


## Values of Different Parameters Used in the Model
### Values Taken from Literature
- koff<sub>null</sub> = 5.00 s<sup>-1</sup>
- koff<sub>select</sub> = 0.36 s<sup>-1</sup>
- koff<sub>delete</sub>= 0.057 s<sup>-1</sup>
- k<sub>p</sub> = 1.0 s<sup>-1</sup> 
- k<sub>SHP1</sub> = 0.05 s<sup>-1</sup>
- N<sub>kp </sub>= 11

### Estimated from Experimental Ranges
- Stromal density: 0.001 
- zap70<sub>mean</sub>:  1.0
- zap70<sub>cv</sub> = 0.3
- zap70<sub>lo</sub> = 0.1
- zap70<sub>hi </sub>= 4.0
- Residence time distribution: Gamma distribution
- T<sub>res-shape</sub> = 5.0
- T<sub>res-scale</sub> = 1152.0
- v<sub>max</sub> = 25.0
- K<sub>low</sub> = 0.5
- nl<sub>ow</sub> = 1.5
- K<sub>high</sub> = 4.0
- n<sub>high</sub> = 2.0
- p<sub>bias</sub> = 0.04
- k<sub>C</sub> = 0.05
- K<sub>C </sub>= 2.5
- k<sub>M</sub> = 0.03
- K<sub>M </sub>= 2.5
- alpha<sub>C,max</sub> = 1.0
- alpha<sub>M,max</sub> = 0.5
- F<sub>T</sub> = 20.0
- K<sub>F</sub> = 2.5
- beta = 0.5
- k<sub>c0 </sub>= 0.5
- F<sub>c</sub> = 12.0
- k<sub>s0</sub> = 0.5
- F<sub>s</sub> = 25.0
- k<sub>scan,base</sub> = 1.0
- gamma = 0.5
- lambda<sub>X</sub> = 0.01
### Modeler's choice 
- grid width (grid_w): 100
- grid height (grid_h): 100
- grid spacing (grid_spacing): 10
- cortex fraction (cortex-frac): 70%
- Total run time (T<sub>run, hr</sub>): 288 hr = 12 days
- No. of thymocytes introduced (n<sub>initial</sub>) = 1000
- Lower limit of koff: 0.001 s<sup>-1</sup>
- Upper limit of koff: 50 s<sup>-1</sup>
- time step duration (dt): 1 min
- seed: 42


## Sensitivity Analyses
For this model, the method of performing the sensitivity analysis (SA) can be described as a stochastic, replicated one-at-a-time (OAT) analysis, with an additional signal-to-noise screening metric. To elaborate a bit more, this type of test varies one parameter or parameter group while holding the others fixed (OAT part), then repeats each setting over predetermined random seeds which are five seeds in this case (stochastic replication part). The signal-to-noise ratio (SNR) is then a custom descriptive statistic: in the current implementation it is the range of the output across the parameter sweep divided by the mean within-condition SD across seeds. The ABM literature specifically recommends repeated runs with different random seeds because stochasticity can substantially change individual simulation outcomes. 
The current version of the model has roughly 23 tunable parameters across SA-2 through SA-9. A full factorial or even a modest-resolution global design across 23 dimensions is computationally prohibitive at the current per-run cost  OAT is the standard pragmatic choice when per-run cost is high and the number of parameters is large, because its cost scales linearly with the number of parameters rather than exponentially or combinatorialy. For this first generation POC model specifically, the purpose of SAs revolves more around the demonstration of understanding which parameters are the main drivers of the model's behavior and distinguishing signal from noise rather than the production of a definitive, publication-grade characterization of the full parameter space. 
However the limitations of the current approach should also be kept in mind. OAT cannot detect parameter interactions which are cases where the effect of parameter A depends on the value of parameter B. *This is not a hypothetical concern for this model, as a concrete case of exactly this problem can be read about in section of the [Overview, rules and results for the hybrid model on the influence of thymic stromal stiffness upon thymic selection outcomes.](./Overview,%20rules%20and%20results%20for%20the%20hybrid%20model%20on%20the%20influence%20of%20thymic%20stromal%20stiffness%20upon%20thymic%20selection%20outcomes..md). The koff distribution (SA-3) and stromal density (SA-4) interact strongly: the stromal density effect is as large as it is because the koff distribution places most cells below the deletion threshold, making fate contact-limited. If you fixed the koff distribution first (shift-factor > 2.5) and then re-ran SA-4, the density effect would almost certainly shrink substantially, because fate would become affinity-limited rather than contact-limited*. OAT, by design, never tests this joint condition and only shows the stromal density effect at the current, uncorrected koff baseline. Any OAT sweep is implicitly conditional on every other parameter sitting at its current default, and the model's two most dominant default values (koff distribution, stromal density) are both acknowledged to have issues with calibration due to certain reasons (lack of proper data in current literature and mismatch of simplistic representation with real life cell shape and contact properties). This means the entire SA runs from SA-2 to SA-9 should be understood as characterizing sensitivity in the current, known-to-be-suboptimal regime which maybe a useful thing to characterize, but is not the same as characterizing sensitivity in the biologically intended regime. Another important point to keep in mind is that SNR > 2 used in the SA interpretation is not, by itself, a standard statistical significance test, but is better described as an empirical noise-floor or signal-to-noise screening criterion. In light of these limitations, more advanced forms of sensitivity analyses such as Morris screening, Latin Hypercube Sampling (LHS) with Partial Rank Correlation Coefficients (PRCC) and Variance-based global sensitivity analysis (Sobol indices, eFAST) can be performed for future versions of the model whose mechanisms have been refined further.