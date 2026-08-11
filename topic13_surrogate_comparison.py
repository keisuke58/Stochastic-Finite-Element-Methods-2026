#!/usr/bin/env python3
"""
Topic 13: Stochastic Surrogate vs Data-Driven Approach.

Systematically compares PCE, Gaussian Process, and Neural Network surrogates
for the H3 fairing problem as a function of training data size N.
Addresses: "Describing the lack of knowledge due to insufficient
experimental data: stochastic approach vs. data-driven approach."
"""
import numpy as np
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import chaospy as cp
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from thesis_style import use, clean_ax

BLUE, RED, GRAY, GREEN, ORANGE = '#2166ac', '#d6604d', '#888888', '#4dac26', '#f4a582'
OUTDIR = 'figures/topic13_comparison'
os.makedirs(OUTDIR, exist_ok=True)

# ── Ground-truth PCE model (from existing calibration) ─────────────────────
# sigma_vM(xi) ≈ 209.7 + 2.365*xi_E1 + 0.211*xi_G12
# (xi_E1, xi_G12, xi_Kn, xi_GIc, xi_tn) ~ N(0,I_5)
C0   = 209.7
C1   = np.array([2.365, 0.211, 0.0, 0.0, 0.0])
# Small 2nd-order terms inferred from degree convergence (deg2≈deg3)
C11  = 0.008   # xi_E1^2 coefficient (Hermite: (xi^2-1)/sqrt2)


def true_model(xi):
    """True sigma_vM [MPa]. xi: (N,5) standard-normal."""
    return (C0
            + C1[0] * xi[:, 0]
            + C1[1] * xi[:, 1]
            + C11 * (xi[:, 0]**2 - 1) / np.sqrt(2))


# ── Test set ───────────────────────────────────────────────────────────────
rng = np.random.default_rng(42)
N_TEST = 20_000
xi_test = rng.standard_normal((N_TEST, 5))
y_test  = true_model(xi_test)


# ── Helper: fit PCE via chaospy ────────────────────────────────────────────
def fit_pce(xi_train, y_train, degree):
    dist = cp.J(*[cp.Normal(0, 1) for _ in range(5)])
    expansion = cp.generate_expansion(degree, dist)
    approx = cp.fit_regression(expansion, xi_train.T, y_train)
    return approx, dist


def rmse_pce(approx, dist, xi_test, y_test):
    y_pred = approx(*xi_test.T)
    return np.sqrt(np.mean((y_pred - y_test)**2))


# ── Helper: fit GP ────────────────────────────────────────────────────────
def fit_gp(xi_train, y_train):
    d = xi_train.shape[1]
    kernel = ConstantKernel(1.0) * RBF(length_scale=np.ones(d)) + WhiteKernel(1e-3)
    gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=1,
                                   normalize_y=True, random_state=0)
    gpr.fit(xi_train, y_train)
    return gpr


def rmse_gp(gpr, xi_test, y_test):
    y_pred = gpr.predict(xi_test)
    return np.sqrt(np.mean((y_pred - y_test)**2))


# ── Helper: fit NN ─────────────────────────────────────────────────────────
def fit_nn(xi_train, y_train):
    # Normalise y so network learns residuals ≈ 0 rather than ~209.7
    y_mean = float(y_train.mean())
    y_std  = float(y_train.std()) + 1e-8
    y_sc   = (y_train - y_mean) / y_std
    pipe = Pipeline([
        ('scaler', StandardScaler()),
        ('nn', MLPRegressor(hidden_layer_sizes=(64, 64, 32),
                             activation='relu',
                             max_iter=3000,
                             random_state=0,
                             early_stopping=True,
                             validation_fraction=0.15,
                             n_iter_no_change=40)),
    ])
    pipe.fit(xi_train, y_sc)
    pipe._y_mean = y_mean
    pipe._y_std  = y_std
    return pipe


def rmse_nn(pipe, xi_test, y_test):
    y_pred = pipe.predict(xi_test) * pipe._y_std + pipe._y_mean
    return np.sqrt(np.mean((y_pred - y_test)**2))


# ── Convergence study ──────────────────────────────────────────────────────
N_LIST    = [10, 15, 20, 30, 40, 66, 100, 150, 286]
N_REPEAT  = 5   # repeated random seeds

results = {
    'pce1': {n: [] for n in N_LIST},
    'pce2': {n: [] for n in N_LIST},
    'gp':   {n: [] for n in N_LIST},
    'nn':   {n: [] for n in N_LIST},
}

print("Running convergence study...")
for seed in range(N_REPEAT):
    rng_s = np.random.default_rng(seed * 100 + 1)
    xi_pool = rng_s.standard_normal((max(N_LIST) + 50, 5))
    y_pool  = true_model(xi_pool)

    for n in N_LIST:
        xi_tr = xi_pool[:n]
        y_tr  = y_pool[:n]

        # PCE degree 1
        try:
            approx1, dist1 = fit_pce(xi_tr, y_tr, degree=1)
            results['pce1'][n].append(rmse_pce(approx1, dist1, xi_test, y_test))
        except Exception:
            results['pce1'][n].append(np.nan)

        # PCE degree 2
        try:
            approx2, dist2 = fit_pce(xi_tr, y_tr, degree=2)
            results['pce2'][n].append(rmse_pce(approx2, dist2, xi_test, y_test))
        except Exception:
            results['pce2'][n].append(np.nan)

        # GP (skip for large N to save time)
        if n <= 150:
            try:
                gpr = fit_gp(xi_tr, y_tr)
                results['gp'][n].append(rmse_gp(gpr, xi_test, y_test))
            except Exception:
                results['gp'][n].append(np.nan)

        # NN
        if n >= 20:
            try:
                nn = fit_nn(xi_tr, y_tr)
                results['nn'][n].append(rmse_nn(nn, xi_test, y_test))
            except Exception:
                results['nn'][n].append(np.nan)

    print(f"  seed {seed+1}/{N_REPEAT} done")

# Compute median RMSE (robust to occasional outlier fits)
def median_rmse(res_dict):
    meds, lo, hi = {}, {}, {}
    for n, vals in res_dict.items():
        v = np.array([x for x in vals if not np.isnan(x)])
        if len(v) == 0:
            meds[n] = np.nan; lo[n] = np.nan; hi[n] = np.nan
        else:
            meds[n] = np.median(v)
            lo[n]   = np.percentile(v, 25)
            hi[n]   = np.percentile(v, 75)
    return meds, lo, hi

med_pce1, lo_pce1, hi_pce1 = median_rmse(results['pce1'])
med_pce2, lo_pce2, hi_pce2 = median_rmse(results['pce2'])
med_gp,   lo_gp,   hi_gp   = median_rmse(results['gp'])
med_nn,   lo_nn,   hi_nn   = median_rmse(results['nn'])

# ── Fig 1: RMSE convergence curves ────────────────────────────────────────
fig, ax = plt.subplots(figsize=use(1.0, 0.60))
for meds, los, his, col, lbl, ls in [
    (med_pce1, lo_pce1, hi_pce1, GREEN,  'PCE degree 1', '--'),
    (med_pce2, lo_pce2, hi_pce2, BLUE,   'PCE degree 2', '-'),
    (med_gp,   lo_gp,   hi_gp,   RED,    'Gaussian Process (RBF)', '-'),
    (med_nn,   lo_nn,   hi_nn,   ORANGE, 'Neural Network (64-64-32)', '-'),
]:
    N_vals  = [n for n in N_LIST if not np.isnan(meds.get(n, np.nan))]
    med_vals = [meds[n] for n in N_vals]
    lo_vals  = [los[n]  for n in N_vals]
    hi_vals  = [his[n]  for n in N_vals]
    ax.semilogy(N_vals, med_vals, color=col, lw=2.2, ls=ls, marker='o', ms=5, label=lbl)
    ax.fill_between(N_vals, lo_vals, hi_vals, color=col, alpha=0.15)

ax.set_xlabel('Training set size $N$')
ax.set_ylabel('Test RMSE [MPa]')
ax.set_title('Surrogate Accuracy vs Training Data Size\n'
             f'(median $\\pm$ IQR over {N_REPEAT} random seeds)')
ax.legend(fontsize=10)
ax.set_xlim(8, 300)
fig.tight_layout()
fig.savefig(f'{OUTDIR}/t13_rmse_convergence.pdf')
plt.close(fig)
print(f"Saved: {OUTDIR}/t13_rmse_convergence.pdf")

# ── Fig 2: Uncertainty decomposition (epistemic vs aleatory) ───────────────
# GP naturally provides epistemic uncertainty.
# At N=66 (our Abaqus budget), compare:
#   PCE: aleatoric (response std from Sobol) with no epistemic
#   GP:  aleatoric (posterior mean variability) + epistemic (posterior std)
# Use GP posterior std as proxy for epistemic component.

N_CALIB = 66
xi_calib = rng.standard_normal((N_CALIB, 5))
y_calib  = true_model(xi_calib)
gpr66    = fit_gp(xi_calib, y_calib)

N_MC = 50_000
xi_mc = rng.standard_normal((N_MC, 5))
y_pce_mc   = true_model(xi_mc)                 # PCE "truth"
y_gp_mean, y_gp_std = gpr66.predict(xi_mc, return_std=True)

E1_axis = np.linspace(-3.5, 3.5, 300)

fig, axes = plt.subplots(1, 2, figsize=use(1.0, 0.44))

# Left: 1D slice along xi_E1 at xi_-E1 = 0, epistemic band from the 5-D GP.
# Fitting a 1-D GP on xi_E1 alone would fold the G12 aleatory spread (c2 =
# 0.211 MPa) into the WhiteKernel, so its band would be ~30x the true
# epistemic uncertainty and would not be an epistemic band at all.
xi_slice = np.zeros((len(E1_axis), 5))
xi_slice[:, 0] = E1_axis
y1d_mean, y1d_std = gpr66.predict(xi_slice, return_std=True)
stress_E1 = true_model(xi_slice)
axes[0].plot(E1_axis, stress_E1, 'k-', lw=2.5, label='True response at $\\xi_{-E_1}=0$')
axes[0].plot(E1_axis, y1d_mean, color=BLUE, lw=2.0, ls='--', label='GP mean')
axes[0].fill_between(E1_axis, y1d_mean - 1.96*y1d_std, y1d_mean + 1.96*y1d_std,
                      color=BLUE, alpha=0.25, label='GP 95% CI (epistemic)')
axes[0].scatter(xi_calib[:, 0], y_calib, s=15, c=GRAY, alpha=0.6, zorder=5,
                label='Training points (projected)')
axes[0].set_xlabel(r'$\xi_{E_1}$ (standard normal)')
axes[0].set_ylabel(r'$\sigma_\mathrm{vM}^\mathrm{max}$ [MPa]')
axes[0].set_title('Aleatory vs Epistemic Uncertainty\n(1D slice, $N=66$)')
axes[0].legend(fontsize=9)

# Right: full response PDF with epistemic band (marginalise over xi_E1)
x_range = np.linspace(200, 222, 300)
axes[1].fill_between(x_range,
                     norm.pdf(x_range, y_pce_mc.mean(), y_pce_mc.std()),
                     color=GRAY, alpha=0.3, label='Aleatory distribution (PCE)')
axes[1].plot(x_range, norm.pdf(x_range, y_pce_mc.mean(), y_pce_mc.std()),
             color='k', lw=2.0, ls='--')
# epistemic range: perturb mean by ±2*median(epistemic std)
eps_med = np.median(y_gp_std)
for sign, col, lbl in [(+1, RED, 'With epistemic (+)'), (-1, BLUE, 'With epistemic (-)')]:
    mu_shifted = y_pce_mc.mean() + sign * eps_med
    axes[1].plot(x_range, norm.pdf(x_range, mu_shifted, y_pce_mc.std()),
                 color=col, lw=1.8, ls=':', label=lbl)
axes[1].set_xlabel(r'$\sigma_\mathrm{vM}^\mathrm{max}$ [MPa]')
axes[1].set_ylabel('Probability Density')
axes[1].set_title('Response PDF: Aleatory Only vs Aleatory + Epistemic')
axes[1].legend(fontsize=9)

fig.tight_layout()
fig.savefig(f'{OUTDIR}/t13_epistemic_aleatory.pdf')
plt.close(fig)
print(f"Saved: {OUTDIR}/t13_epistemic_aleatory.pdf")

# ── Fig 3: Interpretability comparison — PCE Sobol vs GP SHAP-style ────────
# PCE: analytical Sobol indices
# GP / NN: permutation importance (brute-force sensitivity)
N_PERM = 5_000
xi_perm = rng.standard_normal((N_PERM, 5))
gpr66_5d = fit_gp(xi_calib, y_calib)
nn66     = fit_nn(xi_calib, y_calib)
y_base_gp  = gpr66_5d.predict(xi_perm)
y_base_nn  = nn66.predict(xi_perm)

def perm_importance(y_base, model, xi_perm, rng):
    d = xi_perm.shape[1]
    imp = np.zeros(d)
    for j in range(d):
        xi_perm_j = xi_perm.copy()
        xi_perm_j[:, j] = rng.standard_normal(len(xi_perm))
        y_perm = model.predict(xi_perm_j)
        imp[j] = np.mean((y_base - y_perm)**2)
    return imp / imp.sum()

imp_gp = perm_importance(y_base_gp, gpr66_5d, xi_perm, rng)
imp_nn = perm_importance(y_base_nn, nn66,      xi_perm, rng)
sobol  = np.array([0.9921, 0.0079, 0.0, 0.0, 0.0])

params = [r'$E_1$', r'$G_{12}$', r'$K_n$', r'$G_{Ic}$', r'$t_n$']
x = np.arange(5)
w = 0.25

fig, ax = plt.subplots(figsize=use(1.0, 0.50))
ax.bar(x - w, sobol, w, label='PCE Sobol $S_i$ (analytical)', color=BLUE, zorder=3)
ax.bar(x,     imp_gp, w, label='GP permutation importance',    color=RED,  alpha=0.85, zorder=3)
ax.bar(x + w, imp_nn, w, label='NN permutation importance',    color=GREEN,alpha=0.85, zorder=3)
ax.set_xticks(x)
ax.set_xticklabels(params, fontsize=12)
ax.set_ylabel('Relative Importance')
ax.set_title('Sensitivity Analysis: PCE (Analytical) vs GP / NN (Permutation)')
ax.legend(fontsize=10)
fig.tight_layout()
fig.savefig(f'{OUTDIR}/t13_sensitivity_comparison.pdf')
plt.close(fig)
print(f"Saved: {OUTDIR}/t13_sensitivity_comparison.pdf")

# ── Summary table ──────────────────────────────────────────────────────────
print("\n--- RMSE at N=66 (median over seeds) ---")
for key, meds, lbl in [
    ('pce1', med_pce1, 'PCE deg 1'),
    ('pce2', med_pce2, 'PCE deg 2'),
    ('gp',   med_gp,   'GP       '),
    ('nn',   med_nn,   'NN       '),
]:
    val = meds.get(66, np.nan)
    print(f"  {lbl}: {val:.4e} MPa")

print("\nTopic 13 complete.")
