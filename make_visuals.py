#!/usr/bin/env python3
"""
make_visuals.py — visually engaging assets for the Stochastic FEM project.

Produces (into figures/):
  * hero_sfem_pipeline.png      — static hero: uncertain inputs -> FEM stress field -> PDF + Sobol
  * sfem_randomfield_anim.mp4   — Karhunen-Loeve random material field morphing (left) while the
                                  Monte-Carlo output PDF of max von Mises stress builds up (right)

Quantitative numbers (moments, Sobol indices) are the real study results
(CFRP coupon, 5 random params, PCE deg-2, 66 Abaqus runs):
  sigma_vM^max = 209.7 +/- 2.375 MPa ,  u^max = 10.07 +/- 0.178 mm , SDEG = 0
  Sobol(sigma_vM): S1 = [E1 0.9921, G12 0.0079, 0, 0, 0]
The open-hole coupon geometry and stress / random fields are *representative* (clearly
labelled) — they illustrate the SFEM machinery, not a specific Abaqus frame.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle
from scipy.stats import norm
import matplotlib.animation as animation

FFMPEG = ('/home/nishioka/.cache/talkbuild/venv/lib/python3.12/'
          'site-packages/imageio_ffmpeg/binaries/ffmpeg-linux-x86_64-v7.0.2')
if os.path.exists(FFMPEG):
    plt.rcParams['animation.ffmpeg_path'] = FFMPEG

plt.rcParams.update({
    'font.family': 'DejaVu Serif', 'font.size': 13,
    'axes.titlesize': 14, 'axes.labelsize': 12,
    'figure.dpi': 150, 'savefig.dpi': 300,
})

OUT = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(OUT, exist_ok=True)

BLUE, RED, GRAY, GREEN = '#2166ac', '#d6604d', '#888888', '#4dac26'

# ── real study results ────────────────────────────────────────────────────────
PARAM_LABELS = [r'$E_1$', r'$G_{12}$', r'$K_n$', r'$G_{Ic}$', r'$t_n$']
SOBOL_S1 = [0.9921, 0.0079, 0.0, 0.0, 0.0]
SMISES_MEAN, SMISES_STD = 209.7, 2.375   # MPa

# coupon domain
NX, NY = 220, 110
COUPON_W, COUPON_H = 2.0, 1.0
HOLE_R = 0.22

xs = np.linspace(-COUPON_W / 2, COUPON_W / 2, NX)
ys = np.linspace(-COUPON_H / 2, COUPON_H / 2, NY)
XX, YY = np.meshgrid(xs, ys)
RR = np.sqrt(XX**2 + YY**2)
HOLE = RR < HOLE_R                      # masked-out hole


def kl_field(coeffs, modes):
    """Karhunen-Loeve-style Gaussian field from mode coefficients."""
    f = np.zeros_like(XX)
    for (i, j), a in zip(modes, coeffs):
        f += a * np.cos(i * np.pi * (XX / COUPON_W + 0.5)) \
               * np.cos(j * np.pi * (YY / COUPON_H + 0.5)) / (1.0 + i * i + j * j)**0.9
    f = (f - f.mean()) / (f.std() + 1e-9)
    return f


def vonmises_field():
    """Representative open-hole-tension von Mises field (Kirsch-like concentration)."""
    eps = 1e-6
    r = np.maximum(RR, HOLE_R)
    th = np.arctan2(YY, XX)
    a = HOLE_R
    # Kirsch solution for uniaxial tension sigma in x, normalised
    s_rr = 0.5 * (1 - a**2 / r**2) + 0.5 * (1 - 4 * a**2 / r**2 + 3 * a**4 / r**4) * np.cos(2 * th)
    s_tt = 0.5 * (1 + a**2 / r**2) - 0.5 * (1 + 3 * a**4 / r**4) * np.cos(2 * th)
    s_rt = -0.5 * (1 + 2 * a**2 / r**2 - 3 * a**4 / r**4) * np.sin(2 * th)
    vm = np.sqrt(s_rr**2 - s_rr * s_tt + s_tt**2 + 3 * s_rt**2)
    vm = vm / (vm[~HOLE].max() + eps)
    return vm


MODES = [(i, j) for i in range(5) for j in range(4) if (i + j) > 0]
RNG = np.random.default_rng(7)


# ══════════════════════════════════════════════════════════════════════════════
# 1) STATIC HERO — input uncertainty -> FEM stress -> output UQ
# ══════════════════════════════════════════════════════════════════════════════
def make_hero():
    fig = plt.figure(figsize=(16, 6.6))
    fig.patch.set_facecolor('white')
    fig.suptitle('Stochastic FEM: Uncertainty Propagation through a CFRP Composite Panel',
                 fontsize=18, fontweight='bold', y=0.99)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.25, 1.05],
                          left=0.04, right=0.985, top=0.84, bottom=0.13, wspace=0.30)

    # --- (1) uncertain inputs ---
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('1 · Uncertain inputs', fontsize=14, color=BLUE, pad=8)
    means = [1.0, 1.0, 1.0, 1.0, 1.0]
    cvs = [0.05, 0.05, 0.04, 0.04, 0.04]
    yy = np.arange(5)[::-1]
    t = np.linspace(-3, 3, 200)
    for k, (m, cv, yc) in enumerate(zip(means, cvs, yy)):
        pdf = norm.pdf(t, 0, 1)
        pdf = pdf / pdf.max() * 0.42
        ax1.fill_between(t * 0.9, yc, yc + pdf, color=BLUE, alpha=0.30, zorder=2)
        ax1.plot(t * 0.9, yc + pdf, color=BLUE, lw=1.6, zorder=3)
        ax1.text(-3.2, yc + 0.18, PARAM_LABELS[k], fontsize=14, ha='right', va='center')
    ax1.set_xlim(-3.6, 3.2); ax1.set_ylim(-0.3, 5.2)
    ax1.set_yticks([]); ax1.set_xlabel('normalised value')
    ax1.text(0.5, -0.20, '5 random material parameters\n(Gaussian, CV ≈ 4–5 %)',
             transform=ax1.transAxes, ha='center', fontsize=10.5, color=GRAY)
    for sp in ('top', 'right', 'left'): ax1.spines[sp].set_visible(False)

    # --- (2) FEM stress field ---
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_title('2 · FEM stress field  (illustrative)',
                  fontsize=13, color='#5a9a2a', pad=10)
    vm = vonmises_field()
    vmm = np.ma.masked_where(HOLE, vm)
    im = ax2.imshow(vmm, extent=[-COUPON_W/2, COUPON_W/2, -COUPON_H/2, COUPON_H/2],
                    origin='lower', cmap='inferno', aspect='equal')
    ax2.add_patch(Circle((0, 0), HOLE_R, facecolor='white', edgecolor='k', lw=1.2, zorder=4))
    # tension arrows
    for s in (-1, 1):
        ax2.annotate('', xy=(s * (COUPON_W/2 + 0.22), 0), xytext=(s * (COUPON_W/2 + 0.04), 0),
                     arrowprops=dict(arrowstyle='-|>', color='k', lw=2.2))
    ax2.set_xticks([]); ax2.set_yticks([])
    ax2.set_xlim(-COUPON_W/2 - 0.28, COUPON_W/2 + 0.28)
    cb = fig.colorbar(im, ax=ax2, fraction=0.035, pad=0.02)
    cb.set_label(r'von Mises $\sigma_\mathrm{vM}$ (norm.)', fontsize=10)
    ax2.text(0.5, -0.12, r'QoI: $\sigma_\mathrm{vM}^{\max}=209.7$ MPa, '
             r'$u^{\max}=10.07$ mm, SDEG $=0$',
             transform=ax2.transAxes, ha='center', fontsize=10.5, color=GRAY)

    # --- (3) outputs: PDF + Sobol ---
    inner = gs[0, 2].subgridspec(2, 1, height_ratios=[1, 1], hspace=0.55)
    ax3 = fig.add_subplot(inner[0, 0])
    ax3.set_title('3 · Output UQ', fontsize=14, color=RED, pad=8)
    xg = np.linspace(SMISES_MEAN - 4*SMISES_STD, SMISES_MEAN + 4*SMISES_STD, 300)
    ax3.plot(xg, norm.pdf(xg, SMISES_MEAN, SMISES_STD), color=BLUE, lw=2.4)
    ax3.fill_between(xg, norm.pdf(xg, SMISES_MEAN, SMISES_STD), color=BLUE, alpha=0.18)
    ax3.axvline(SMISES_MEAN, color=GRAY, ls='--', lw=1.2)
    ax3.set_xlabel(r'$\sigma_\mathrm{vM}^{\max}$ [MPa]'); ax3.set_yticks([])
    ax3.text(0.02, 0.86, r'PCE deg-2 · 66 runs' + '\n' +
             r'$\mu=209.7,\ \sigma=2.375$', transform=ax3.transAxes, fontsize=9.5)
    for sp in ('top', 'right', 'left'): ax3.spines[sp].set_visible(False)

    ax4 = fig.add_subplot(inner[1, 0])
    x = np.arange(5)
    bars = ax4.bar(x, SOBOL_S1, 0.6, color=[RED if v > 0.5 else BLUE for v in SOBOL_S1],
                   zorder=3)
    ax4.set_xticks(x); ax4.set_xticklabels(PARAM_LABELS, fontsize=12)
    ax4.set_ylim(0, 1.12); ax4.set_ylabel(r'Sobol $S_i$', fontsize=11)
    ax4.text(0, 0.99, '0.992', ha='center', va='bottom', fontsize=10, color=RED, fontweight='bold')
    ax4.set_title(r'$E_1$ drives 99 % of the variance', fontsize=11, pad=4)
    ax4.grid(axis='y', alpha=0.3)
    for sp in ('top', 'right'): ax4.spines[sp].set_visible(False)

    # flow arrows between panels
    for x0 in (0.345, 0.665):
        fig.add_artist(FancyArrowPatch((x0, 0.50), (x0 + 0.022, 0.50),
                       transform=fig.transFigure, arrowstyle='-|>', mutation_scale=22,
                       color='#444', lw=2.5))

    out = os.path.join(OUT, 'hero_sfem_pipeline.png')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('saved', out)


# ══════════════════════════════════════════════════════════════════════════════
# 2) ANIMATION — KL random field morphing + MC output PDF building
# ══════════════════════════════════════════════════════════════════════════════
def make_anim():
    n_key = 9                      # random control fields
    per = 14                       # frames between keys
    nframes = n_key * per
    key_coeffs = [RNG.standard_normal(len(MODES)) for _ in range(n_key + 1)]
    fields = [kl_field(c, MODES) for c in key_coeffs]

    # pre-draw MC samples (one per frame), reveal progressively
    samples = RNG.normal(SMISES_MEAN, SMISES_STD, nframes)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.5, 5.2),
                                   gridspec_kw=dict(width_ratios=[1.05, 1.0]))
    fig.subplots_adjust(left=0.02, right=0.96, top=0.86, bottom=0.13, wspace=0.22)
    fig.suptitle('Stochastic FEM — sampling material uncertainty and building the response PDF',
                 fontsize=14, fontweight='bold', y=0.97)

    xg = np.linspace(SMISES_MEAN - 4.5*SMISES_STD, SMISES_MEAN + 4.5*SMISES_STD, 300)
    target_pdf = norm.pdf(xg, SMISES_MEAN, SMISES_STD)
    vmin, vmax = -2.4, 2.4

    def smooth(a, b, t):
        w = 0.5 - 0.5 * np.cos(np.pi * t)      # cosine ease
        return (1 - w) * a + w * b

    def update(i):
        axL.clear(); axR.clear()
        k = i // per
        t = (i % per) / per
        f = smooth(fields[k], fields[k + 1], t)
        fm = np.ma.masked_where(HOLE, f)
        axL.imshow(fm, extent=[-COUPON_W/2, COUPON_W/2, -COUPON_H/2, COUPON_H/2],
                   origin='lower', cmap='RdYlBu_r', vmin=vmin, vmax=vmax, aspect='equal')
        axL.add_patch(Circle((0, 0), HOLE_R, facecolor='white', edgecolor='k', lw=1.0, zorder=4))
        axL.add_patch(Rectangle((-COUPON_W/2, -COUPON_H/2), COUPON_W, COUPON_H,
                                fill=False, edgecolor='k', lw=1.4, zorder=3))
        axL.set_xticks([]); axL.set_yticks([])
        axL.set_xlim(-COUPON_W/2 - 0.1, COUPON_W/2 + 0.1)
        axL.set_title(r'Random material field $E_1(\mathbf{x})$ — KL realization %d' % (i + 1),
                      fontsize=12)
        axL.text(0.5, -0.09, 'Karhunen–Loève samples (illustrative of SFEM random-field input)',
                 transform=axL.transAxes, ha='center', fontsize=9.5, color=GRAY)

        # right: accumulating histogram + running stats + target pdf
        s = samples[:i + 1]
        axR.hist(s, bins=24, range=(xg[0], xg[-1]), density=True,
                 color=BLUE, alpha=0.45, edgecolor='white', lw=0.4, zorder=2)
        if i > 6:
            axR.plot(xg, target_pdf, color=RED, lw=2.4, zorder=3,
                     label=r'PCE surrogate $\mathcal{N}(209.7,\,2.375^2)$')
            axR.legend(loc='upper right', fontsize=9.5, frameon=False)
        axR.axvline(s.mean(), color='k', ls='--', lw=1.2, zorder=4)
        axR.set_xlim(xg[0], xg[-1]); axR.set_ylim(0, target_pdf.max() * 1.45)
        axR.set_xlabel(r'Max von Mises stress $\sigma_\mathrm{vM}^{\max}$ [MPa]')
        axR.set_yticks([])
        axR.set_title('Monte-Carlo output distribution', fontsize=12)
        axR.text(0.02, 0.93, 'samples  n = %d' % (i + 1), transform=axR.transAxes,
                 fontsize=10.5, fontweight='bold')
        axR.text(0.02, 0.84, r'$\hat\mu=%.2f$  MPa' % s.mean(), transform=axR.transAxes,
                 fontsize=10)
        axR.text(0.02, 0.76, r'$\hat\sigma=%.2f$  MPa' % (s.std() if i > 0 else 0),
                 transform=axR.transAxes, fontsize=10)
        for sp in ('top', 'right', 'left'): axR.spines[sp].set_visible(False)
        return []

    anim = animation.FuncAnimation(fig, update, frames=nframes, interval=110)
    out = os.path.join(OUT, 'sfem_randomfield_anim.mp4')
    writer = animation.FFMpegWriter(fps=9, bitrate=3000)
    anim.save(out, writer=writer, dpi=120)
    plt.close(fig)
    print('saved', out, '(%d frames, %.1fs)' % (nframes, nframes / 9))


if __name__ == '__main__':
    make_hero()
    make_anim()
    print('done.')
