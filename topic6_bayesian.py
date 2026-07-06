#!/usr/bin/env python3
"""
Topic 6: Bayesian Parameter Identification via PCE surrogate.

Uses the analytical PCE reconstruction of the H3 fairing forward model
(no new Abaqus runs required) to infer the posterior distribution of E1
from synthetic displacement/stress measurements.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm, truncnorm
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from thesis_style import use, clean_ax

# ── Style ───────────────────────────────────────────────────────────────────
BLUE, RED, GRAY, GREEN, ORANGE = '#2166ac', '#d6604d', '#888888', '#4dac26', '#f4a582'
OUTDIR = 'figures/topic6_bayesian'
os.makedirs(OUTDIR, exist_ok=True)

# ── Analytical PCE reconstruction ──────────────────────────────────────────
# From existing results:
#   mean_stress = 209.7 MPa, std_stress = 2.375 MPa
#   S_E1 = 0.9921, S_G12 = 0.0079  (first-order Sobol indices)
#   E1: mu=146000 MPa, sigma=7300 MPa  (CoV=5%)
#   G12: mu=5200 MPa, sigma=520 MPa   (CoV=10%)
#
# In standard-normal space (xi_E1, xi_G12 ~ N(0,1)):
#   sigma_vM(xi) ≈ c0 + c1*xi_E1 + c2*xi_G12
# where:
#   Var = c1^2 + c2^2 = 2.375^2 = 5.6406
#   S_E1 = c1^2 / (c1^2+c2^2) = 0.9921  → c1 = sqrt(0.9921*5.6406) = 2.365
#   S_G12 = c2^2 / (c1^2+c2^2) = 0.0079 → c2 = sqrt(0.0079*5.6406) = 0.2111

MU_STRESS  = 209.7   # MPa
STD_STRESS = 2.375   # MPa
C0         = MU_STRESS
C1_E1      = np.sqrt(0.9921 * STD_STRESS**2)   # = 2.365
C2_G12     = np.sqrt(0.0079 * STD_STRESS**2)   # = 0.211

MU_E1, SIG_E1     = 146_000., 7_300.   # MPa
MU_G12, SIG_G12   = 5_200.,   520.    # MPa


def pce_forward(xi_E1, xi_G12=0.0):
    """PCE surrogate: sigma_vM [MPa] given standard-normal inputs."""
    return C0 + C1_E1 * xi_E1 + C2_G12 * xi_G12


def xi_from_E1(E1):
    return (E1 - MU_E1) / SIG_E1


def E1_from_xi(xi):
    return MU_E1 + SIG_E1 * xi


# ── Bayesian setup ─────────────────────────────────────────────────────────
# Prior: xi_E1 ~ TruncNormal(0, 1, lo, hi)  (physical bounds ±50% of mu)
XI_LO = (0.5 * MU_E1 - MU_E1) / SIG_E1   # = -0.5*mu/sig = -10.0
XI_HI = (1.5 * MU_E1 - MU_E1) / SIG_E1   # = +10.0
prior_a = (XI_LO - 0) / 1.
prior_b = (XI_HI - 0) / 1.
prior_dist = truncnorm(a=prior_a, b=prior_b, loc=0, scale=1)

# Synthetic experiment: "true" E1 is 2% above nominal
E1_TRUE  = 153_000.   # MPa  (representative of a good but above-spec batch)
XI_TRUE  = xi_from_E1(E1_TRUE)   # = 0.9589
Y_TRUE   = pce_forward(XI_TRUE)  # ≈ 212.0 MPa

# Measurement noise standard deviation (digital strain gauge accuracy)
SIGMA_MEAS = 2.0   # MPa

rng = np.random.default_rng(0)

# ── Analytical Bayesian update (conjugate Gaussian) ────────────────────────
# Observed stress (k repeated measurements from the same panel batch)
def bayesian_update_analytical(y_obs_mean, k=1, sigma_meas=SIGMA_MEAS):
    """
    Conjugate Gaussian-Gaussian posterior for xi_E1.
    Prior: xi_E1 ~ N(0,1)
    Likelihood: y | xi_E1 ~ N(C0 + C1*xi_E1, sigma_meas^2 / k)
    """
    sigma_eff2 = sigma_meas**2 / k  # effective noise with k measurements
    h = C1_E1

    # Posterior parameters in xi-space
    sig_post2  = 1.0 / (1.0 / 1.0 + h**2 / sigma_eff2)
    mu_post    = sig_post2 * (h * (y_obs_mean - C0) / sigma_eff2)
    return mu_post, np.sqrt(sig_post2)


# Generate single noisy observation
y_obs = Y_TRUE + rng.normal(0, SIGMA_MEAS)
print(f"True E1     : {E1_TRUE:,.0f} MPa  (xi = {XI_TRUE:.4f})")
print(f"True stress : {Y_TRUE:.2f} MPa")
print(f"Observed    : {y_obs:.2f} MPa  (noise = {y_obs-Y_TRUE:.2f} MPa)")

mu_post_xi, sig_post_xi = bayesian_update_analytical(y_obs, k=1)
mu_post_E1  = E1_from_xi(mu_post_xi)
sig_post_E1 = SIG_E1 * sig_post_xi
print(f"\n--- Posterior (k=1) ---")
print(f"xi_E1 posterior: N({mu_post_xi:.4f}, {sig_post_xi:.4f}²)")
print(f"E1 posterior   : N({mu_post_E1:,.0f}, {sig_post_E1:,.0f}²) MPa")

# ── Fig 1: Prior vs Posterior of E1 for k = 1, 5, 20 ─────────────────────
E1_plot = np.linspace(MU_E1 - 4*SIG_E1, MU_E1 + 4*SIG_E1, 800)
xi_plot = xi_from_E1(E1_plot)
prior_pdf = prior_dist.pdf(xi_plot) / SIG_E1   # Jacobian for E1 axis

fig, ax = plt.subplots(figsize=use(1.0, 0.55))
ax.plot(E1_plot / 1e3, prior_pdf * 1e3, color=GRAY, lw=2.5,
        ls='--', label='Prior $p(E_1)$')
ax.axvline(E1_TRUE / 1e3, color='k', lw=1.5, ls=':', alpha=0.6, label=f'True $E_1$ = {E1_TRUE/1e3:.0f} GPa')

colors_k = [BLUE, RED, GREEN]
for k, col in zip([1, 5, 20], colors_k):
    mu_xi, sig_xi = bayesian_update_analytical(y_obs, k=k)
    mu_E1  = E1_from_xi(mu_xi)
    sig_E1 = SIG_E1 * sig_xi
    pdf = norm.pdf(E1_plot, mu_E1, sig_E1)
    ax.plot(E1_plot / 1e3, pdf * 1e3, color=col, lw=2.2,
            label=f'Posterior $k={k}$ (mean={mu_E1/1e3:.1f} GPa)')

ax.set_xlabel(r'$E_1$ [GPa]')
ax.set_ylabel('Probability Density [GPa$^{-1}$]')
ax.set_title('Bayesian Identification of $E_1$: Prior vs Posterior')
ax.legend(fontsize=10)
fig.tight_layout()
fig.savefig(f'{OUTDIR}/bayes_prior_posterior_E1.pdf')
plt.close(fig)
print(f"Saved: {OUTDIR}/bayes_prior_posterior_E1.pdf")

# ── Fig 2: Likelihood surface and posterior in (E1, G12) space ────────────
# Full 2D: (xi_E1, xi_G12) ~ N(0,I), observe y_obs
# Likelihood: p(y_obs | xi_E1, xi_G12) = N(pce(xi_E1,xi_G12), sigma_meas^2)
xi_E1_grid  = np.linspace(-3.5, 3.5, 300)
xi_G12_grid = np.linspace(-3.5, 3.5, 300)
XE, XG = np.meshgrid(xi_E1_grid, xi_G12_grid)
Y_pred = pce_forward(XE, XG)
log_lik  = -0.5 * ((y_obs - Y_pred) / SIGMA_MEAS)**2
log_prior = -0.5 * (XE**2 + XG**2)
log_post = log_lik + log_prior
post = np.exp(log_post - log_post.max())

fig, axes = plt.subplots(1, 2, figsize=use(1.0, 0.42))
# — Likelihood —
cm1 = axes[0].contourf(xi_E1_grid, xi_G12_grid, np.exp(log_lik), 20, cmap='YlOrRd')
axes[0].axvline(XI_TRUE, color='k', ls=':', lw=1.5, label=f'True $\\xi_{{E_1}}$={XI_TRUE:.2f}')
fig.colorbar(cm1, ax=axes[0])
axes[0].set_xlabel(r'$\xi_{E_1}$')
axes[0].set_ylabel(r'$\xi_{G_{12}}$')
axes[0].set_title('Likelihood $p(y_\\mathrm{obs}|\\xi)$')
axes[0].legend(fontsize=10)
# — Posterior —
cm2 = axes[1].contourf(xi_E1_grid, xi_G12_grid, post, 20, cmap='Blues')
axes[1].axvline(XI_TRUE, color='k', ls=':', lw=1.5, label=f'True $\\xi_{{E_1}}$={XI_TRUE:.2f}')
axes[1].axhline(0.0,     color=RED, ls=':', lw=1.5, label=r'True $\xi_{G_{12}}$=0')
fig.colorbar(cm2, ax=axes[1])
axes[1].set_xlabel(r'$\xi_{E_1}$')
axes[1].set_ylabel(r'$\xi_{G_{12}}$')
axes[1].set_title('Unnormalised Posterior $p(\\xi|y_\\mathrm{obs})$')
axes[1].legend(fontsize=10)
fig.suptitle(r'Bayesian Identification: Likelihood and Posterior in $(\xi_{E_1},\,\xi_{G_{12}})$ Space',
             fontsize=13)
fig.tight_layout()
fig.savefig(f'{OUTDIR}/bayes_likelihood_posterior_2d.pdf')
plt.close(fig)
print(f"Saved: {OUTDIR}/bayes_likelihood_posterior_2d.pdf")

# ── Fig 3: Posterior stress distribution (before vs after identification) ──
N_MC = 500_000
# Prior samples
xi_E1_prior  = rng.standard_normal(N_MC)
xi_G12_prior = rng.standard_normal(N_MC)
stress_prior = pce_forward(xi_E1_prior, xi_G12_prior)

# Posterior samples (k=5 measurements, analytical)
mu_xi5, sig_xi5 = bayesian_update_analytical(y_obs, k=5)
xi_E1_post5 = rng.normal(mu_xi5, sig_xi5, N_MC)
stress_post5 = pce_forward(xi_E1_post5, xi_G12_prior)

mu_xi20, sig_xi20 = bayesian_update_analytical(y_obs, k=20)
xi_E1_post20 = rng.normal(mu_xi20, sig_xi20, N_MC)
stress_post20 = pce_forward(xi_E1_post20, xi_G12_prior)

x_range = np.linspace(200, 222, 400)

fig, ax = plt.subplots(figsize=use(1.0, 0.55))
for samples, col, lbl in [
    (stress_prior,  GRAY,   'Prior propagation'),
    (stress_post5,  BLUE,   'Posterior propagation ($k=5$)'),
    (stress_post20, RED,    'Posterior propagation ($k=20$)'),
]:
    mu, sig = samples.mean(), samples.std()
    ax.plot(x_range, norm.pdf(x_range, mu, sig), color=col, lw=2.2, label=f'{lbl}')
ax.axvline(Y_TRUE, color='k', ls=':', lw=1.5, label=f'True stress = {Y_TRUE:.1f} MPa')
ax.set_xlabel(r'Max von Mises Stress $\sigma_\mathrm{vM}^\mathrm{max}$ [MPa]')
ax.set_ylabel('Probability Density')
ax.set_title('Propagated Response Uncertainty: Before vs After Bayesian Identification')
ax.legend(fontsize=10)
fig.tight_layout()
fig.savefig(f'{OUTDIR}/bayes_propagated_stress.pdf')
plt.close(fig)
print(f"Saved: {OUTDIR}/bayes_propagated_stress.pdf")

# ── Fig 4: Posterior mean and 95% CI vs number of measurements ─────────────
K_list  = np.array([1, 2, 3, 5, 8, 10, 15, 20, 30, 50])
mu_E1s  = []
lo95_E1s = []
hi95_E1s = []
for k in K_list:
    mu_xi_, sig_xi_ = bayesian_update_analytical(y_obs, k=k)
    mu_ = E1_from_xi(mu_xi_)
    sig_ = SIG_E1 * sig_xi_
    mu_E1s.append(mu_)
    lo95_E1s.append(mu_ - 1.96 * sig_)
    hi95_E1s.append(mu_ + 1.96 * sig_)

mu_E1s   = np.array(mu_E1s)   / 1e3
lo95_E1s = np.array(lo95_E1s) / 1e3
hi95_E1s = np.array(hi95_E1s) / 1e3

fig, ax = plt.subplots(figsize=use(1.0, 0.55))
ax.fill_between(K_list, lo95_E1s, hi95_E1s, color=BLUE, alpha=0.25, label='95% credible interval')
ax.plot(K_list, mu_E1s,                   color=BLUE,  lw=2.2, marker='o', ms=5, label='Posterior mean')
ax.axhline(E1_TRUE / 1e3, color='k', ls=':', lw=1.5, label=f'True $E_1$ = {E1_TRUE/1e3:.1f} GPa')
ax.axhline((MU_E1 - 1.96*SIG_E1) / 1e3,  color=GRAY, ls='--', lw=1.2, alpha=0.6)
ax.axhline((MU_E1 + 1.96*SIG_E1) / 1e3,  color=GRAY, ls='--', lw=1.2, alpha=0.6,
           label='Prior 95% interval')
ax.set_xlabel('Number of Stress Measurements $k$')
ax.set_ylabel(r'$E_1$ [GPa]')
ax.set_title('Posterior Convergence: $E_1$ Estimate vs Number of Measurements')
ax.legend(fontsize=10)
fig.tight_layout()
fig.savefig(f'{OUTDIR}/bayes_convergence_k.pdf')
plt.close(fig)
print(f"Saved: {OUTDIR}/bayes_convergence_k.pdf")

print("\nTopic 6 complete.")
