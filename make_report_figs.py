#!/usr/bin/env python3
"""
make_report_figs.py — extra honest, on-topic figures for the SFEM report.

All figures are faithful to THIS study's model (H3 fairing CFRP/Al-honeycomb sandwich
panel with a cohesive-zone adhesive, 5 uncertain params, PCE surrogate). They are
generated from the documented description/results — no figures are borrowed from the
GNN guided-wave work (different physics/model).

Outputs (figures/):
  * h3_fairing_context.png   — rocket + fairing structural context (schematic)
  * fem_panel_model.png      — sandwich-panel FEM model schematic with the 5 params located
  * reliability_margin.png   — limit-state view: response PDF vs failure limit, margin & beta
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch, FancyBboxPatch, Circle
from matplotlib.lines import Line2D
from scipy.stats import norm

plt.rcParams.update({
    'font.family': 'DejaVu Serif', 'font.size': 12,
    'axes.titlesize': 13, 'axes.labelsize': 12,
    'figure.dpi': 150, 'savefig.dpi': 300,
})
OUT = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(OUT, exist_ok=True)

BLUE, RED, GRAY, GREEN, ORANGE = '#2166ac', '#c62828', '#888888', '#2e7d32', '#e08020'
CFRP = '#2b3a55'
CORE = '#f0c33c'


# ── 1) rocket + fairing structural context ────────────────────────────────────
def make_structure_context():
    fig, (axr, axp) = plt.subplots(1, 2, figsize=(11, 5.4),
                                   gridspec_kw=dict(width_ratios=[1.0, 1.25]))
    fig.suptitle('Application: the H3 launch-vehicle payload fairing',
                 fontsize=15, fontweight='bold', y=0.98)

    # --- left: H3 rocket silhouette ---
    axr.set_xlim(-3.2, 3.2); axr.set_ylim(-0.3, 11.5); axr.set_aspect('equal'); axr.axis('off')
    # core stage
    axr.add_patch(Rectangle((-0.7, 0.2), 1.4, 7.0, facecolor='#dfe6ef', edgecolor='k', lw=1.3))
    # fairing (nose) — highlighted
    fairing = Polygon([(-0.7, 7.2), (0.7, 7.2), (0.55, 9.2), (0.0, 10.6), (-0.55, 9.2)],
                      closed=True, facecolor=BLUE, edgecolor='k', lw=1.5, alpha=0.9, zorder=3)
    axr.add_patch(fairing)
    axr.text(0.0, 11.15, 'payload fairing', ha='center', va='bottom', color=BLUE,
             fontsize=11, fontweight='bold', zorder=4)
    # two solid boosters
    for s in (-1, 1):
        axr.add_patch(Rectangle((s * 1.45 - 0.32, 0.6), 0.64, 5.4,
                                facecolor='#c3ccd8', edgecolor='k', lw=1.1))
        axr.add_patch(Polygon([(s * 1.45 - 0.32, 6.0), (s * 1.45 + 0.32, 6.0), (s * 1.45, 6.9)],
                              closed=True, facecolor='#c3ccd8', edgecolor='k', lw=1.1))
        axr.add_patch(Polygon([(s * 1.45 - 0.32, 0.6), (s * 1.45 + 0.32, 0.6),
                               (s * 1.45 + 0.55, 0.0), (s * 1.45 - 0.55, 0.0)],
                              closed=True, facecolor='#9aa6b5', edgecolor='k', lw=1.0))
    # core nozzle
    axr.add_patch(Polygon([(-0.7, 0.2), (0.7, 0.2), (1.0, -0.3), (-1.0, -0.3)],
                          closed=True, facecolor='#9aa6b5', edgecolor='k', lw=1.0))
    axr.text(0.0, -0.05, 'H3', ha='center', fontsize=9, color=GRAY)
    # callout from fairing to the right panel
    axr.annotate('', xy=(3.0, 8.4), xytext=(0.6, 8.4),
                 arrowprops=dict(arrowstyle='-|>', color=BLUE, lw=2))

    # --- right: a fairing wall patch -> sandwich stack ---
    axp.set_xlim(0, 10); axp.set_ylim(0, 10); axp.axis('off')
    axp.text(5, 9.4, 'A representative wall panel is a bonded sandwich',
             ha='center', fontsize=12)
    # sandwich stack (exploded)
    x0, w = 2.2, 5.6
    layers = [(6.6, 1.0, CFRP, 'CFRP skin  (outer)'),
              (5.0, 1.4, CORE, 'Al-honeycomb core'),
              (3.4, 1.0, CFRP, 'CFRP skin  (inner)')]
    for (y, h, col, lab) in layers:
        axp.add_patch(Rectangle((x0, y), w, h, facecolor=col, edgecolor='k', lw=1.3))
        axp.text(x0 + w + 0.25, y + h / 2, lab, va='center', fontsize=11)
    # honeycomb hint in the core band
    for xx in np.linspace(x0 + 0.4, x0 + w - 0.4, 9):
        axp.add_patch(Polygon([(xx, 5.15), (xx + 0.28, 5.4), (xx + 0.28, 5.9),
                               (xx, 6.15), (xx - 0.28, 5.9), (xx - 0.28, 5.4)],
                              closed=True, fill=False, edgecolor='#9a7a05', lw=0.8))
    # adhesive CZM interfaces
    for y in (6.6, 6.45, 4.85, 4.7):
        pass
    for y in (6.55, 4.78):
        axp.plot([x0, x0 + w], [y, y], color=RED, lw=2.0, ls=(0, (4, 2)), zorder=5)
    axp.text(x0 + w + 0.25, 6.55, 'adhesive (CZM)', va='center', fontsize=10, color=RED)
    axp.text(x0 + w + 0.25, 4.78, 'adhesive (CZM)', va='center', fontsize=10, color=RED)
    axp.text(5, 2.4, 'Manufacturing variability in the skin and the bond\n'
             'propagates to stress, damage and deflection.',
             ha='center', fontsize=10.5, color=GRAY)

    out = os.path.join(OUT, 'h3_fairing_context.png')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig); print('saved', out)


# ── 2) sandwich-panel FEM model schematic ─────────────────────────────────────
def make_fem_schematic():
    fig, ax = plt.subplots(figsize=(11, 6.0))
    fig.suptitle('Finite-element model: bonded sandwich panel with cohesive interfaces',
                 fontsize=15, fontweight='bold', y=0.97)
    ax.set_xlim(0, 14); ax.set_ylim(0, 8.4); ax.axis('off')

    x0, w = 1.3, 8.0
    skin_t, core_t = 0.9, 2.0
    y_inner = 1.7
    y_core = y_inner + skin_t
    y_outer = y_core + core_t
    y_top = y_outer + skin_t

    def mesh(x, y, w, h, nx, ny, fc):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor='none', zorder=1))
        for i in range(nx + 1):
            ax.plot([x + i * w / nx] * 2, [y, y + h], color='k', lw=0.4, alpha=0.45, zorder=2)
        for j in range(ny + 1):
            ax.plot([x, x + w], [y + j * h / ny] * 2, color='k', lw=0.4, alpha=0.45, zorder=2)
        ax.add_patch(Rectangle((x, y), w, h, fill=False, edgecolor='k', lw=1.4, zorder=3))

    mesh(x0, y_inner, w, skin_t, 24, 2, '#9fb0cc')          # inner skin
    mesh(x0, y_core, w, core_t, 24, 3, CORE)                # core
    mesh(x0, y_outer, w, skin_t, 24, 2, '#9fb0cc')          # outer skin
    # CZM interfaces
    for y in (y_core, y_outer):
        ax.plot([x0, x0 + w], [y, y], color=RED, lw=2.6, ls=(0, (5, 2)), zorder=5)
    ax.text(x0 + w + 0.2, y_core, 'cohesive interface', va='center', color=RED, fontsize=10)

    # layer labels (left)
    ax.text(x0 - 0.2, y_inner + skin_t / 2, 'CFRP skin', ha='right', va='center', fontsize=11)
    ax.text(x0 - 0.2, y_core + core_t / 2, 'Al-honeycomb\ncore', ha='right', va='center', fontsize=11)
    ax.text(x0 - 0.2, y_outer + skin_t / 2, 'CFRP skin', ha='right', va='center', fontsize=11)

    # pressure load on top
    for xx in np.linspace(x0 + 0.4, x0 + w - 0.4, 11):
        ax.annotate('', xy=(xx, y_top), xytext=(xx, y_top + 0.9),
                    arrowprops=dict(arrowstyle='-|>', color=BLUE, lw=1.8))
    ax.text(x0 + w / 2, y_top + 1.15, 'pressure load (Step 1) + incremental load (Step 2)',
            ha='center', fontsize=11, color=BLUE)

    # fixed BC (hatched ground) at both bottom edges
    for xx in (x0, x0 + w):
        ax.add_patch(Rectangle((xx - 0.45, y_inner - 0.65), 0.9, 0.5,
                               facecolor='none', edgecolor='k', hatch='////', lw=1.0))
        ax.plot([xx, xx], [y_inner - 0.15, y_inner], color='k', lw=1.2)
    ax.text(x0, y_inner - 0.9, 'clamped', ha='center', fontsize=9.5)
    ax.text(x0 + w, y_inner - 0.9, 'clamped', ha='center', fontsize=9.5)

    # parameter callouts
    cw = 2.7
    cx = x0 + w + 1.5
    def callout(x, y, txt, col):
        ax.add_patch(FancyBboxPatch((x, y), cw, 0.62, boxstyle='round,pad=0.05',
                                    facecolor=col, edgecolor='k', lw=1.0, alpha=0.92, zorder=8))
        ax.text(x + cw / 2, y + 0.31, txt, ha='center', va='center', fontsize=11,
                color='white', fontweight='bold', zorder=9)
    callout(cx, y_outer + 0.1, r'$E_1,\ G_{12}$', BLUE)
    ax.annotate('', xy=(x0 + w - 0.2, y_outer + skin_t / 2), xytext=(cx, y_outer + 0.4),
                arrowprops=dict(arrowstyle='-|>', color=BLUE, lw=1.6))
    callout(cx, y_core + 0.7, r'$K_n,\ G_{Ic},\ t_n$', RED)
    ax.annotate('', xy=(x0 + w, y_core), xytext=(cx, y_core + 1.0),
                arrowprops=dict(arrowstyle='-|>', color=RED, lw=1.6))

    ax.text(x0 + w / 2, 0.35,
            'Skin stiffness ($E_1,G_{12}$) and the adhesive cohesive law ($K_n,G_{Ic},t_n$)'
            ' are the 5 uncertain inputs.', ha='center', fontsize=10.5, color=GRAY)

    out = os.path.join(OUT, 'fem_panel_model.png')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig); print('saved', out)


# ── 3) reliability limit-state figure ─────────────────────────────────────────
def make_reliability():
    mean, std = 209.7, 2.375          # max von Mises (MPa), PCE deg-2
    limit = 1200.0                    # failure threshold (MPa)
    beta = (limit - mean) / std

    fig, ax = plt.subplots(figsize=(11, 4.8))
    fig.suptitle('Reliability: the response sits far below the failure limit',
                 fontsize=15, fontweight='bold', y=0.99)

    xg = np.linspace(mean - 5 * std, mean + 5 * std, 400)
    pdf = norm.pdf(xg, mean, std)
    ax.fill_between(xg, pdf, color=BLUE, alpha=0.25, zorder=2)
    ax.plot(xg, pdf, color=BLUE, lw=2.4, zorder=3,
            label=r'response PDF  $\sigma_\mathrm{vM}^{\max}\sim\mathcal{N}(209.7,\,2.375^2)$')
    ax.axvline(mean, color=BLUE, ls='--', lw=1.2, zorder=3)
    ax.axvline(limit, color=RED, lw=2.6, zorder=4, label=f'failure limit = {limit:.0f} MPa')

    # safety-margin arrow
    ax.annotate('', xy=(limit, pdf.max() * 0.45), xytext=(mean, pdf.max() * 0.45),
                arrowprops=dict(arrowstyle='<|-|>', color=GREEN, lw=2))
    ax.text((mean + limit) / 2, pdf.max() * 0.52,
            f'safety margin  ≈ {limit - mean:.0f} MPa', ha='center', color=GREEN, fontsize=12)

    ax.text(limit, pdf.max() * 0.92,
            r'$\beta=\dfrac{\mu_{\mathrm{limit}}-\mu_R}{\sigma_R}\approx %.0f$' % beta + '\n'
            r'$P_f \approx 0$', ha='right', va='top', color=RED, fontsize=12)
    ax.text(mean, pdf.max() * 1.02, 'mean 209.7 MPa', ha='center', fontsize=10, color=BLUE)

    ax.set_xlim(mean - 6 * std, limit + 90)
    ax.set_ylim(0, pdf.max() * 1.18)
    ax.set_xlabel(r'Max von Mises stress $\sigma_\mathrm{vM}^{\max}$ [MPa]')
    ax.set_yticks([])
    ax.legend(loc='center right', frameon=False, fontsize=10.5)
    for sp in ('top', 'right', 'left'): ax.spines[sp].set_visible(False)
    ax.text(0.5, -0.20, 'A degree-2 PCE surrogate places the entire stress distribution '
            'far below the limit; cohesive damage never initiates (SDEG $=0$).',
            transform=ax.transAxes, ha='center', fontsize=10.5, color=GRAY)

    out = os.path.join(OUT, 'reliability_margin.png')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig); print('saved', out)


if __name__ == '__main__':
    make_structure_context()
    make_fem_schematic()
    make_reliability()
    print('done.')
