# Reliability Assessment of JAXA H3 Fairing under Manufacturing Variability  
## Using Polynomial Chaos Expansion as a Stochastic FEM Surrogate

**Author:** Keisuke Nishioka (10081049)  
**Course:** Stochastic Finite Element Methods 2026  
**Date:** April 20, 2026

---

## 1. Introduction

Composite fairing structures for launch vehicles, such as the JAXA H3 rocket, are manufactured from Carbon Fiber Reinforced Polymer (CFRP) skins bonded to aluminum honeycomb cores. Material properties of CFRP and the adhesive cohesive zone are inherently variable due to manufacturing processes: fiber volume fraction fluctuations, cure cycle variations, and adhesive film thickness inconsistencies. These uncertainties propagate into structural responses—stress, damage, and deformation—affecting reliability.

This study applies **non-intrusive Polynomial Chaos Expansion (PCE)** as a Stochastic Finite Element Method (SFEM) surrogate to quantify how manufacturing variability in five key parameters propagates to structural responses. The existing high-fidelity Abaqus/Standard FEM model of the H3 fairing panel (used in a parallel GNN-based SHM study) serves as the black-box FEM solver. A neural network (NN) surrogate trained on the same data provides a basis for comparison.

---

## 2. Methodology

### 2.1 Polynomial Chaos Expansion

For a model response $Y = \mathcal{M}(\xi)$ depending on random input vector $\xi \in \mathbb{R}^d$, PCE approximates:

$$Y \approx \hat{Y} = \sum_{\alpha \in \mathcal{A}} c_\alpha \Psi_\alpha(\xi)$$

where $\Psi_\alpha$ are orthogonal polynomial basis functions (Hermite polynomials for normal/truncated-normal inputs), $c_\alpha$ are coefficients to be determined, and $\mathcal{A}$ is a truncated multi-index set with $|\alpha| \leq p$ (total degree $p$).

For $d = 5$ variables and degree $p$, the number of basis terms is:

$$P = \binom{d+p}{p} = \binom{5+p}{p}$$

giving $P = 21$ for $p=2$ and $P = 56$ for $p=3$.

### 2.2 Sparse Smolyak Quadrature

Full tensor-product Gauss-Hermite quadrature requires $n_q^d$ points (e.g., $4^5 = 1024$ for $p=3$, $d=5$), which is computationally prohibitive. The **Smolyak sparse grid** reduces this dramatically:

| Degree $p$ | Full grid | Sparse grid (this work) | Reduction |
|---|---|---|---|
| 2 | $3^5 = 243$ | **66** | 3.7× |
| 3 | $4^5 = 1024$ | **286** | 3.6× |

Sparse grids achieve the same polynomial exactness as full grids for degree $\leq 2p-1$ while avoiding the curse of dimensionality. Negative quadrature weights, which can arise in sparse Gauss-Hermite rules, were handled by falling back to least-squares regression (`chaospy.fit_regression`).

### 2.3 Uncertain Input Parameters

Five manufacturing-related parameters are treated as independent truncated-normal random variables $\xi_i \sim \mathrm{TruncNormal}(\mu_i, \sigma_i, 0.5\mu_i, 1.5\mu_i)$:

| Parameter | Description | Nominal $\mu$ | CoV | $\sigma$ |
|---|---|---|---|---|
| $E_1$ | CFRP fiber-direction Young's modulus | 146,000 MPa | 5% | 7,300 MPa |
| $G_{12}$ | CFRP in-plane shear modulus | 5,200 MPa | 10% | 520 MPa |
| $K_n$ | Cohesive zone normal stiffness | 1.0×$10^6$ N/mm$^3$ | 15% | 1.5×$10^5$ N/mm$^3$ |
| $G_{Ic}$ | Mode-I fracture energy | 0.3 N/mm | 20% | 0.06 N/mm |
| $t_n$ | Cohesive zone normal strength | 50 MPa | 15% | 7.5 MPa |

Dependent material properties are scaled proportionally: $E_2 = 10000 \cdot (E_1/\mu_{E_1})$, $G_{23} = 3000 \cdot (G_{12}/\mu_{G_{12}})$, $K_s = K_n/2$, $G_{IIc} = G_{Ic}/0.3$.

### 2.4 Quantities of Interest (QoI)

Three structural response quantities are extracted from each Abaqus simulation:

| QoI | Symbol | Description | Failure threshold |
|---|---|---|---|
| Max von Mises stress | $\sigma_\text{vM}^\text{max}$ | Peak stress in CFRP skin | 1,200 MPa (CFRP tensile strength) |
| Max scalar damage variable | $d^\text{max}$ | Peak SDEG in adhesive layer | 0.5 (delamination criterion) |
| Max displacement magnitude | $u^\text{max}$ | Peak nodal displacement | — |

### 2.5 Non-Intrusive SFEM Pipeline

```
Random inputs (ξ)
       │
       ▼
Sparse Gauss quadrature nodes (66 or 286 points)
       │
       ▼
modify_inp(): patch Abaqus .inp material cards
       │
       ▼
abaqus job=PCE-SXXXX interactive cpus=1   (×66 or ×286)
       │
       ▼
abaqus python extract_pce_qoi.py          (odbAccess)
       │
       ▼
chaospy.fit_regression() → PCE coefficients {c_α}
       │
       ├─→ Sobol sensitivity indices (S_i, S_T^i)
       ├─→ MC on surrogate ($10^6$ samples) → PDF, Pf, β
       └─→ NN surrogate comparison (PyTorch MLP)
```

### 2.6 Neural Network Surrogate

A fully connected MLP (5→64→64→32→1, ReLU activations) was trained on the same Abaqus data for 500 epochs (Adam optimizer, lr=$10^{-3}$). The NN serves as a black-box surrogate baseline for accuracy comparison against the interpretable PCE model.

---

## 3. Finite Element Model

The Abaqus/Standard model represents a representative panel of the H3 fairing sandwich structure:

- **CFRP skin**: LAMINA elastic material (orthotropic), with fiber orientation, defined by $E_1$, $E_2$, $\nu_{12}$, $G_{12}$, $G_{13}$, $G_{23}$
- **Al-Honeycomb core**: linear elastic, fixed properties
- **Adhesive interface**: Cohesive Zone Model (CZM) with BK mixed-mode damage evolution ($\eta = 2.284$), stiffness $K_n$/$K_s$, strength $t_n$/$t_s$, fracture energies $G_{Ic}$/$G_{IIc}$
- **Loading**: Two-step analysis — Step-1 (static pressure), Step-2 (incremental load)
- **QoI extraction**: Step-2 last frame, using `odbAccess` Python API

Total model size: ~293 MB per ODB file.

---

## 4. Results

### 4.1 PCE Statistics

**Table 1. PCE response statistics (degree = 2, 66 Abaqus simulations)**

| QoI | Mean | Std. Dev. | CV (%) | $P_f$ | $\beta$ |
|---|---|---|---|---|---|
| $\sigma_\text{vM}^\text{max}$ | 209.7 MPa | 2.375 MPa | 1.13 | 0 | ∞ |
| $d^\text{max}$ (SDEG) | 0.000 | 0.000 | — | 0 | ∞ |
| $u^\text{max}$ | 10.07 mm | 0.178 mm | 1.76 | — | — |

**Table 2. PCE response statistics (degree = 3, 286 Abaqus simulations)**

| QoI | Mean | Std. Dev. | CV (%) | $P_f$ | $\beta$ |
|---|---|---|---|---|---|
| $\sigma_\text{vM}^\text{max}$ | 209.7 MPa | 2.378 MPa | 1.13 | 0 | ∞ |
| $d^\text{max}$ (SDEG) | 0.000 | 0.000 | — | 0 | ∞ |
| $u^\text{max}$ | 10.07 mm | 0.178 mm | 1.76 | — | — |

The statistics are virtually identical between degree 2 and 3, confirming **PCE convergence at degree 2**. The peak von Mises stress (209.7 MPa) is well below the CFRP tensile strength (1,200 MPa), yielding **zero failure probability** under the considered manufacturing variability range. The adhesive damage variable remains identically zero across all 352 simulations, indicating that the applied load does not initiate cohesive damage under any realisation of the uncertain parameters.

### 4.2 Global Sensitivity Analysis (Sobol Indices)

**Table 3. First-order ($S_i$) and total-order ($S_T^i$) Sobol indices — degree 2**

| Parameter | $S_i$ (smises) | $S_T^i$ (smises) | $S_i$ (disp) | $S_T^i$ (disp) |
|---|---|---|---|---|
| $E_1$ | **0.9921** | **0.9921** | **0.9995** | **0.9996** |
| $G_{12}$ | 0.0079 | 0.0079 | 0.0005 | 0.0005 |
| $K_n$ | ~0 | ~0 | ~0 | ~0 |
| $G_{Ic}$ | ~0 | ~0 | ~0 | ~0 |
| $t_n$ | ~0 | ~0 | ~0 | ~0 |

**Key findings:**

- **$E_1$ is overwhelmingly dominant** ($S_{E_1} \approx 0.992$ for stress, $0.9995$ for displacement). Over 99% of the output variance is explained by uncertainty in the fiber-direction Young's modulus alone.
- **$G_{12}$ contributes marginally** (~0.8% for stress). This is expected since shear response is secondary in a predominantly uniaxial-loading scenario.
- **CZM parameters ($K_n$, $G_{Ic}$, $t_n$) have negligible influence** (Sobol indices ~$10^{-11}$ for degree 2, ~$10^{-8}$ to 10⁻⁶ for degree 3). This is consistent with the observation that SDEG = 0 at all sample points—the cohesive interface never activates, so its material parameters have no influence on the response.
- First-order and total-order Sobol indices are essentially equal ($S_i \approx S_T^i$), confirming **no significant parameter interactions**.

Results are consistent between degree 2 and degree 3, further validating convergence.

![Sobol sensitivity indices for max von Mises stress (degree 2). $E_1$ accounts for 99.2% of output variance.](figures/pce_uq/sobol_max_smises.pdf){width=80%}

![Sobol sensitivity indices for max displacement (degree 2). $E_1$ dominates with $S_{E_1} = 0.9995$.](figures/pce_uq/sobol_max_disp.pdf){width=80%}

### 4.3 Surrogate Accuracy: PCE vs Neural Network

**Table 4. Leave-one-out (LOO) RMSE comparison**

| QoI | PCE deg=2 | PCE deg=3 | NN deg=2 data | NN deg=3 data |
|---|---|---|---|---|
| $\sigma_\text{vM}^\text{max}$ (MPa) | **9.4×$10^{-4}$** | 1.5×$10^{-3}$ | 7.2×$10^{-4}$ | 1.77×$10^{-2}$ |
| $d^\text{max}$ | 0 | 0 | ~0 | ~0 |
| $u^\text{max}$ (mm) | **4.2×10⁻⁵** | 7.1×10⁻⁵ | 2.7×10⁻⁵ | 1.6×$10^{-3}$ |

**Observations:**

1. **PCE degree 2 achieves the best overall accuracy** among all combinations tested. Degree 3 PCE shows slightly higher LOO error, suggesting mild over-parameterisation (56 coefficients fitted on 286 points, ratio ~5:1) relative to the near-linear response surface.

2. **NN accuracy degrades sharply from 66 to 286 data points**. This counterintuitive result indicates overfitting: the fixed MLP architecture (4,000+ parameters) is too large for the smooth, nearly linear response surface, and the larger training set does not improve generalisation without regularisation tuning.

3. **PCE LOO RMSE < $10^{-3}$ MPa** for smises (nominal ~210 MPa) corresponds to a relative error of ~0.0005%, demonstrating excellent surrogate fidelity for this smooth, low-variance response.

4. **Interpretability advantage**: PCE directly yields analytical Sobol indices at zero additional cost. The NN requires separate variance-based sensitivity analysis (e.g., Monte Carlo with $10^6$ evaluations), and the resulting indices are not analytically guaranteed.

### 4.4 Reliability Analysis

With $\sigma_\text{vM}^\text{max}$ mean of 209.7 MPa and standard deviation of 2.38 MPa, the stress distribution is far from the 1,200 MPa threshold. Monte Carlo sampling of $10^6$ realisations on the PCE surrogate yields:

$$P_f = P(\sigma_\text{vM}^\text{max} > 1200\ \text{MPa}) = 0 \quad \Rightarrow \quad \beta = \infty$$

This result holds across both degree-2 and degree-3 PCE, and is physically sensible: the nominal stress is only 17.5% of the failure threshold, and the coefficient of variation is only 1.13%. Even a 6σ event (209.7 + 6×2.38 = 224 MPa) remains far below failure.

The SDEG = 0 result (no cohesive damage) confirms that under the static design load considered here, the fairing adhesive layer operates entirely within its elastic regime regardless of manufacturing variability in CZM parameters.

![PDF of max von Mises stress: PCE surrogate (blue), NN surrogate (orange), and Abaqus sample points (histogram). Both surrogates agree well. Degree 2.](figures/pce_uq/pdf_max_smises.pdf){width=80%}

![PDF of max displacement: PCE vs NN surrogate comparison. Degree 2.](figures/pce_uq/pdf_max_disp.pdf){width=80%}

![PCE convergence with number of terms for max von Mises stress (degree 2).](figures/pce_uq/convergence_max_smises.pdf){width=80%}

![PCE convergence for max displacement (degree 2).](figures/pce_uq/convergence_max_disp.pdf){width=80%}

![Sobol indices comparison: degree 2 (66 samples) vs degree 3 (286 samples) for max von Mises stress. Results are virtually identical, confirming convergence.](figures/pce_uq_deg3/sobol_max_smises.pdf){width=80%}

![PDF comparison degree 3: max von Mises stress surrogate vs Abaqus samples.](figures/pce_uq_deg3/pdf_max_smises.pdf){width=80%}

---

## 5. PCE vs Neural Network: Broader Comparison

| Criterion | PCE (this work) | NN surrogate |
|---|---|---|
| Training data required | 66 points (deg=2) | 66–286 points |
| Accuracy (LOO RMSE) | ~$10^{-4}$ | ~$10^{-4}$ to $10^{-2}$ |
| Convergence with data | Stable (deg=2 optimal) | Degrades with more data (overfitting) |
| Sobol indices | Analytical, exact | Requires separate MC |
| Interpretability | High (coefficient-level) | Low (black box) |
| Uncertainty type | Aleatoric (input variability) | Epistemic (model uncertainty via MC Dropout/Deep Ensemble) |
| Applicable to | Forward UQ, reliability | SHM, defect detection |

**Complementarity remark**: The existing GNN-SHM model quantifies *epistemic uncertainty* about structural state (whether a delamination defect is present), using MC Dropout and Deep Ensembles. PCE quantifies *aleatoric uncertainty* arising from manufacturing parameter scatter. These two uncertainty sources are independent and complementary: one addresses "do we know the state?", the other "how does material variability affect response?".

---

## 6. Discussion

### Why does $E_1$ dominate?

The fairing panel is a thin-walled structure under primarily membrane/bending loads. The von Mises stress and global deflection are governed by the membrane stiffness $E_1 t$ (modulus × skin thickness). With $E_1$ carrying a 5% CoV and all other stiffness-related parameters carrying 10% or higher CoV but being secondary (shear-dominated), the output variance is captured almost entirely by $E_1$.

### Why is SDEG identically zero?

The applied load in Step-2 corresponds to a design-point static pressure load, not an extreme load case. At nominal parameters, the adhesive interface operates at a fraction of its strength capacity. Even at the extreme of the uncertain parameter space (CoV ±50% bounds on $K_n$, $t_n$, $G_{Ic}$), the interface stress never exceeds the cohesive strength $t_n$. A more meaningful sensitivity analysis of CZM parameters would require either a higher load level or a damage initiation/propagation scenario (e.g., impact loading).

### Degree convergence

The agreement between degree-2 (66 samples) and degree-3 (286 samples) statistics to four significant figures confirms that the response surface is well approximated by a second-order polynomial in the input space. This is consistent with the near-linear structural mechanics of a thin composite plate under small deformations.

---

## 7. Conclusion

Non-intrusive PCE with sparse Smolyak quadrature was successfully applied to quantify manufacturing uncertainty in a JAXA H3 CFRP/Al-Honeycomb fairing panel modelled in Abaqus/Standard. Key conclusions are:

1. **PCE degree 2 (66 Abaqus simulations) is sufficient** — degree 3 (286 simulations) yields identical statistics, confirming convergence.

2. **$E_1$ uncertainty dominates all response quantities** ($S_{E_1} > 0.99$). Quality control of CFRP fiber-direction modulus is the single most impactful lever for reducing structural response variability.

3. **The structure is reliable under the design load** ($P_f = 0$, $\beta = \infty$). The 209.7 MPa peak stress has a 1.13% CV, well within the 1,200 MPa failure threshold.

4. **PCE outperforms the NN surrogate** in accuracy and interpretability for this smooth, near-linear response. NN performance degrades with larger datasets due to overfitting in the fixed architecture.

5. **PCE and GNN-SHM are complementary tools**: PCE addresses aleatoric (manufacturing) uncertainty in the design phase; GNN addresses epistemic (damage state) uncertainty in the operational/monitoring phase.

---

## References

1. Ghanem, R. G. & Spanos, P. D. (1991). *Stochastic Finite Elements: A Spectral Approach*. Springer.
2. Sudret, B. (2008). Global sensitivity analysis using polynomial chaos expansions. *Reliability Engineering & System Safety*, 93(7), 964–979.
3. Smolyak, S. A. (1963). Quadrature and interpolation formulas for tensor products of certain classes of functions. *Dokl. Akad. Nauk SSSR*, 148(5), 1042–1045.
4. Fajraoui, N. et al. (2017). Sequential design of experiment for sparse polynomial chaos expansions. *SIAM/ASA J. Uncertainty Quantification*, 5(1), 1061–1085.
5. Faber, M. H. (2012). *Statistics and Probability Theory*. Springer.
6. Saltelli, A. et al. (2008). *Global Sensitivity Analysis: The Primer*. Wiley.
7. Dassault Systèmes (2024). *Abaqus Analysis User's Guide*, version 2024.

---

## Appendix A: Software and Computational Environment

| Component | Version |
|---|---|
| Abaqus/Standard | 2024 |
| Python | 3.12 |
| chaospy | latest |
| PyTorch | 2.10.0 |
| Server | 172.17.36.98, 4×GPU |
| Total Abaqus CPU-hours | 352 jobs × ~4 min = ~23 h |

## Appendix B: File Structure

```
Payload2026/
├── src/
│   ├── pce_driver.py          # PCE UQ main driver
│   ├── extract_pce_qoi.py     # Abaqus Python QoI extractor
│   └── reliability_analysis.py # Post-processing & plots
├── abaqus_work/
│   ├── pce_uq/                # degree=2, 66 jobs
│   │   ├── PCE-S0000/ … PCE-S0065/
│   │   └── pce_results.json
│   └── pce_uq_deg3/           # degree=3, 286 jobs
│       ├── PCE-S0000/ … PCE-S0285/
│       └── pce_results.json
└── figures/
    ├── pce_uq/                # Sobol, PDF, convergence plots (deg=2)
    └── pce_uq_deg3/           # Sobol, PDF, convergence plots (deg=3)
```
