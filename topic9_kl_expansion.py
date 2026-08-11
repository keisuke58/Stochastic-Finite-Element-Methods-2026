#!/usr/bin/env python3
"""
Topic 9: Spatially varying random field representation of E1 via KL expansion.

Represents the fiber-direction modulus E1(x,y) of the H3 fairing panel as a
Gaussian random field with an anisotropic exponential kernel, discretised via
the Karhunen-Loeve expansion, and quantifies how spatial averaging attenuates
the variability of the effective panel modulus relative to the scalar
(perfectly correlated) model used in the forward study.

Scope: this is an input-side study. The KL realisations are not propagated
through the FE model, so the attenuation bounds the input variability seen by
average-governed responses (e.g. global deflection) but says nothing directly
about the peak von Mises stress, which is a local extreme.

The field is Gaussian while the forward study uses a normal truncated to
[0.5mu, 1.5mu]. For E1 (CoV = 5%) those bounds sit at mu +/- 10 sigma, so the
two laws are numerically indistinguishable here.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from thesis_style import use, clean_ax

BLUE, RED, GRAY, GREEN = '#2166ac', '#d6604d', '#888888', '#4dac26'
OUTDIR = 'figures/topic9_kl'
os.makedirs(OUTDIR, exist_ok=True)

# ── H3 fairing panel geometry ──────────────────────────────────────────────
LX, LY   = 600., 400.   # panel dimensions [mm]
MU_E1    = 146_000.      # mean E1 [MPa]
SIG_E1   = 7_300.        # std  E1 [MPa]   (CoV = 5%)

# ── KL discretisation grid ─────────────────────────────────────────────────
NX, NY   = 20, 15
x_pts    = np.linspace(0, LX, NX)
y_pts    = np.linspace(0, LY, NY)
XX, YY   = np.meshgrid(x_pts, y_pts)
coords   = np.column_stack([XX.ravel(), YY.ravel()])   # (NX*NY, 2)
N_TOTAL  = NX * NY


def cov_matrix(coords, l_cx, l_cy, sigma2=1.0):
    """Anisotropic exponential covariance: C(r) = sigma^2 exp(-|dx|/lx -|dy|/ly)."""
    dx = np.abs(coords[:, None, 0] - coords[None, :, 0])   # (N,N)
    dy = np.abs(coords[:, None, 1] - coords[None, :, 1])
    return sigma2 * np.exp(-dx / l_cx - dy / l_cy)


def kl_decompose(C, n_modes=None):
    """
    KL eigendecomposition: C = Φ Λ Φᵀ, sorted descending.
    Returns eigenvalues λ and eigenvectors φ (columns).
    """
    lam, phi = eigh(C)
    # sort descending
    idx = np.argsort(lam)[::-1]
    lam = lam[idx]
    phi = phi[:, idx]
    # discard negative numerical noise
    lam = np.maximum(lam, 0.0)
    if n_modes is not None:
        lam = lam[:n_modes]
        phi = phi[:, :n_modes]
    return lam, phi


def kl_sample(lam, phi, xi, mu=0.0, sigma=1.0):
    """
    Reconstruct random field realisation given KL coefficients xi ~ N(0,1).
    Returns field values at all grid points.
    """
    return mu + sigma * (phi @ (np.sqrt(lam) * xi))


def explained_variance(lam, n_modes):
    """
    Cumulative fraction of total variance explained by the first n_modes.

    `lam` must be the FULL eigenvalue spectrum: the denominator is the total
    variance. Passing an already-truncated spectrum makes the ratio collapse
    to 1.0 by construction.
    """
    return np.cumsum(lam[:n_modes]) / lam.sum()


# ── Study 1: effect of correlation length ──────────────────────────────────
L_CASES = [(200, 150), (400, 300), (1200, 800)]   # (l_cx, l_cy) in mm
LABELS  = [r'$l_{cx}=200,\,l_{cy}=150$ mm (short)',
           r'$l_{cx}=400,\,l_{cy}=300$ mm (medium)',
           r'$l_{cx}=1200,\,l_{cy}=800$ mm (long)']

# ── Fig 1: Eigenvalue decay for the three correlation lengths ──────────────
fig, axes = plt.subplots(1, 2, figsize=use(1.0, 0.38))
colors_lc = [RED, BLUE, GREEN]

for (lcx, lcy), lbl, col in zip(L_CASES, LABELS, colors_lc):
    C    = cov_matrix(coords, lcx, lcy, sigma2=1.0)
    lam, phi = kl_decompose(C)
    n_show = 60
    axes[0].semilogy(np.arange(1, n_show+1), lam[:n_show] / lam.sum(),
                     color=col, lw=2.0, label=lbl)
    n_modes = np.arange(1, N_TOTAL+1)
    ev = explained_variance(lam, N_TOTAL)
    axes[1].plot(n_modes[:80], ev[:80] * 100, color=col, lw=2.0, label=lbl)

axes[0].set_xlabel('KL mode index $i$')
axes[0].set_ylabel(r'Relative eigenvalue $\lambda_i / \sum\lambda_j$')
axes[0].set_title('Eigenvalue Decay vs Correlation Length')
axes[0].legend(fontsize=9.5)

axes[1].axhline(99, color='k', ls='--', lw=1.2, alpha=0.6)
axes[1].set_xlabel('Number of retained KL modes $M$')
axes[1].set_ylabel('Explained variance [%]')
axes[1].set_title('Explained Variance vs $M$')
axes[1].legend(fontsize=9.5)
axes[1].set_xlim(1, 80)

fig.tight_layout()
fig.savefig(f'{OUTDIR}/kl_eigenvalue_decay.pdf')
plt.close(fig)
print(f"Saved: {OUTDIR}/kl_eigenvalue_decay.pdf")

# ── Fig 2: First 4 KL mode shapes (medium correlation length) ─────────────
lcx, lcy = 400, 300
C = cov_matrix(coords, lcx, lcy, sigma2=SIG_E1**2)
# keep the full spectrum: explained_variance() needs the total variance below
lam_all, phi_all = kl_decompose(C)
lam, phi = lam_all[:4], phi_all[:, :4]

fig, axes = plt.subplots(2, 2, figsize=use(1.0, 0.8))
for i, ax in enumerate(axes.flat):
    mode = phi[:, i].reshape(NY, NX)
    # normalise for display
    vmax = np.abs(mode).max()
    im = ax.contourf(XX, YY, mode, levels=20, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    ax.set_title(f'Mode {i+1}:  $\\lambda_{i+1}$ = {lam[i]/SIG_E1**2:.4f}')
    ax.set_xlabel('$x$ [mm]')
    ax.set_ylabel('$y$ [mm]')
    ax.set_aspect('equal')
    fig.colorbar(im, ax=ax, fraction=0.04)

ev4 = explained_variance(lam_all, 4)[-1] * 100
fig.suptitle(
    rf'First 4 KL Mode Shapes of $E_1(x,y)$  '
    rf'($l_{{cx}}=400,\,l_{{cy}}=300$ mm; {ev4:.1f}\% variance)',
    fontsize=13)
fig.tight_layout()
fig.savefig(f'{OUTDIR}/kl_mode_shapes.pdf')
plt.close(fig)
print(f"Saved: {OUTDIR}/kl_mode_shapes.pdf")

# ── Fig 3: Sample realisations of E1(x,y) field ───────────────────────────
rng = np.random.default_rng(7)
N_MODES_USE = 10   # truncate at 10 modes for panel (medium lc)
lcx, lcy = 400, 300
C = cov_matrix(coords, lcx, lcy, sigma2=SIG_E1**2)
lam_full, phi_full = kl_decompose(C)
lam10, phi10 = lam_full[:N_MODES_USE], phi_full[:, :N_MODES_USE]
ev10 = explained_variance(lam_full, N_MODES_USE)[-1] * 100

fig, axes = plt.subplots(2, 3, figsize=use(1.0, 0.72))
for ax in axes.flat:
    xi = rng.standard_normal(N_MODES_USE)
    field = kl_sample(lam10, phi10, xi, mu=MU_E1, sigma=1.0)
    field = field.reshape(NY, NX)
    vmin, vmax = MU_E1 - 3*SIG_E1, MU_E1 + 3*SIG_E1
    im = ax.contourf(XX, YY, field / 1e3, levels=20, cmap='coolwarm',
                     vmin=vmin/1e3, vmax=vmax/1e3)
    ax.set_xlabel('$x$ [mm]')
    ax.set_ylabel('$y$ [mm]')
    ax.set_aspect('equal')
    fig.colorbar(im, ax=ax, label='$E_1$ [GPa]')

fig.suptitle(
    rf'Sample Realisations of $E_1(x,y)$ via KL Expansion '
    rf'($M={N_MODES_USE}$ modes, {ev10:.1f}\% variance; $l_{{cx}}=400$ mm)',
    fontsize=12)
fig.tight_layout()
fig.savefig(f'{OUTDIR}/kl_sample_fields.pdf')
plt.close(fig)
print(f"Saved: {OUTDIR}/kl_sample_fields.pdf")

# ── Fig 4: Impact of M on effective E1 variance (homogenisation argument) ─
# Effective E1 for a panel = spatial average ∫ E1(x,y) dA / A
# Var[E1_eff] = (1/A^2) ∫∫ C(x1,x2) dx1 dx2
# Numerically:
def variance_effective(lam, phi, NX, NY):
    """Var of spatial mean from KL truncation."""
    # spatial average weight: w_i = (1/N) * sum_k phi_{ki}
    w = phi.mean(axis=0)   # (M,)
    return np.sum(lam * w**2)

M_list  = np.arange(1, 81)
var_eff_short  = []
var_eff_medium = []
var_eff_long   = []

for lc_pair, lst in zip([(200, 150), (400, 300), (1200, 800)],
                         [var_eff_short, var_eff_medium, var_eff_long]):
    lcx_, lcy_ = lc_pair
    C_ = cov_matrix(coords, lcx_, lcy_, sigma2=SIG_E1**2)
    lam_, phi_ = kl_decompose(C_)
    for M in M_list:
        lst.append(np.sqrt(variance_effective(lam_[:M], phi_[:, :M], NX, NY)))

scalar_std = SIG_E1   # std of scalar (homogeneous) E1

fig, ax = plt.subplots(figsize=use(1.0, 0.55))
ax.axhline(scalar_std / 1e3, color='k', ls='--', lw=1.5, alpha=0.8,
           label='Scalar E1 std (current model, homogeneous)')
for lst, col, lbl in zip([var_eff_short, var_eff_medium, var_eff_long],
                          [RED, BLUE, GREEN], LABELS):
    ax.plot(M_list, np.array(lst) / 1e3, color=col, lw=2.0, label=lbl)
ax.set_xlabel('Truncation $M$ (number of KL modes)')
ax.set_ylabel(r'Std dev of $\overline{E_1}$ [GPa]')
ax.set_title(r'Spatial Mean $\overline{E_1}$ Variability vs KL Truncation')
ax.legend(fontsize=9.5)
fig.tight_layout()
fig.savefig(f'{OUTDIR}/kl_effective_variance.pdf')
plt.close(fig)
print(f"Saved: {OUTDIR}/kl_effective_variance.pdf")

# ── Summary stats ──────────────────────────────────────────────────────────
print("\n--- Modes needed for 99% variance ---")
for lc_pair, lbl in zip([(200, 150), (400, 300), (1200, 800)], LABELS):
    lcx_, lcy_ = lc_pair
    C_ = cov_matrix(coords, lcx_, lcy_)
    lam_, _ = kl_decompose(C_)
    ev = explained_variance(lam_, N_TOTAL)
    M99 = int(np.searchsorted(ev, 0.99)) + 1
    print(f"  {lbl[:30]:30s}: M={M99} modes")

print("\nTopic 9 complete.")
