# SFEM Final Project — JAXA H3 Fairing UQ

**Course:** Stochastic Finite Element Methods (SoSe 2026), LUH  
**Supervisors:** Prof. Dr. Udo Nackenhorst, Dr. Zhibao Zheng

## Overview

Uncertainty Quantification of the JAXA H3 rocket fairing (CFRP/Al-honeycomb sandwich panel) using Non-Intrusive Polynomial Chaos Expansion (PCE), with three extensions:

| Topic | Method | Key result |
|---|---|---|
| **Forward UQ** | PCE deg. 2/3, Smolyak sparse quadrature | E₁ explains 99.2% of stress variance; Pf ≈ 0 |
| **Topic 6** — Bayesian inverse | Conjugate Gaussian update via PCE likelihood | 5 measurements → 74% posterior std reduction |
| **Topic 9** — KL random field | Anisotropic exponential kernel, KL expansion | Scalar model conservative (σ_eff < σ_scalar) |
| **Topic 13** — Surrogate comparison | PCE vs GP vs NN, N=10–286 | PCE deg.2: RMSE=3×10⁻¹³ MPa (exact recovery) |

## Files

```
report_paper.tex / .pdf          # Main report (33 pages: ~25 main + appendices C/D/E)
slides_sfem.tex                  # Presentation slides (Beamer)
make_report_figs.py              # Generates PCE UQ figures
generate_figures.py              # Additional visualization
make_visuals.py                  # Hero / pipeline figures
topic6_bayesian.py               # Topic 6: Bayesian identification
topic9_kl_expansion.py           # Topic 9: KL expansion of E1(x,y)
topic13_surrogate_comparison.py  # Topic 13: PCE vs GP vs NN convergence
thesis_style.py                  # Matplotlib style (LaTeX fonts, LUH width)
figures/                         # All generated figures
```

## Uncertain Parameters

| Parameter | Mean | CoV | Description |
|---|---|---|---|
| E₁ | 146,000 MPa | 5% | Fiber-direction modulus |
| G₁₂ | 5,200 MPa | 10% | In-plane shear modulus |
| Kₙ | 1,000 N/mm³ | 50% | Cohesive normal stiffness |
| GIc | 0.5 N/mm | 50% | Mode-I fracture energy |
| tₙ | 50 MPa | 50% | Cohesive tensile strength |

## Requirements

```
pip install chaospy numpy scipy matplotlib scikit-learn
```

TeX Live 2025 required for figure rendering (`thesis_style.py`).  
Abaqus ODB files (~293 MB each) are stored on the LUH compute server.
