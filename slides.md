---
title: "Reliability Assessment of JAXA H3 Fairing under Manufacturing Variability"
subtitle: "PCE as a Stochastic FEM Surrogate"
author: "Keisuke Nishioka (10081049)"
date: "Stochastic Finite Element Methods 2026"
institute: ""
theme: "Madrid"
colortheme: "beaver"
fontsize: 10pt
aspectratio: 169
header-includes:
  - \usepackage{booktabs}
  - \usepackage{graphicx}
  - \setbeamertemplate{navigation symbols}{}
  - \newcommand{\vM}{\sigma_\text{vM}^\text{max}}
---

# Motivation

## Why Uncertainty Quantification for Rocket Fairings?

\begin{columns}
\begin{column}{0.55\textwidth}
\textbf{JAXA H3 Rocket Fairing}
\begin{itemize}
  \item CFRP skin + Al-Honeycomb core sandwich structure
  \item Manufacturing scatter is unavoidable:
  \begin{itemize}
    \item Fiber volume fraction variation
    \item Cure cycle fluctuations
    \item Adhesive film thickness variation
  \end{itemize}
  \item How much does this affect structural reliability?
\end{itemize}

\vspace{0.5em}
\textbf{Approach:} Non-intrusive PCE as SFEM surrogate\\
\small{(existing Abaqus model as black-box FEM solver)}
\end{column}
\begin{column}{0.45\textwidth}
\textbf{Uncertain parameters (5):}
\vspace{0.3em}

\small
\begin{tabular}{lcc}
\toprule
Parameter & CoV \\
\midrule
$E_1$ (CFRP modulus) & 5\% \\
$G_{12}$ (shear mod.) & 10\% \\
$K_n$ (CZM stiffness) & 15\% \\
$G_{Ic}$ (fracture energy) & 20\% \\
$t_n$ (CZM strength) & 15\% \\
\bottomrule
\end{tabular}
\end{column}
\end{columns}

---

# Methodology

## Polynomial Chaos Expansion (PCE)

**PCE approximates the model response as:**
$$Y \approx \hat{Y} = \sum_{\alpha \in \mathcal{A}} c_\alpha \Psi_\alpha(\xi)$$

- $\Psi_\alpha$: Hermite polynomial basis (orthogonal w.r.t. Gaussian measure)
- $c_\alpha$: coefficients fitted by least-squares regression
- Non-intrusive: Abaqus treated as a **black box**

\vspace{0.5em}

**Sparse Smolyak quadrature** avoids the curse of dimensionality:

\begin{center}
\begin{tabular}{cccc}
\toprule
Degree $p$ & Full grid & \textbf{Sparse grid} & Reduction \\
\midrule
2 & $3^5 = 243$ & \textbf{66} & 3.7$\times$ \\
3 & $4^5 = 1024$ & \textbf{286} & 3.6$\times$ \\
\bottomrule
\end{tabular}
\end{center}

\vspace{0.3em}
\small{Negative weights $\Rightarrow$ fallback to \texttt{fit\_regression} (chaospy)}

---

## Non-Intrusive SFEM Pipeline

\begin{center}
\texttt{Random inputs $\xi$}\\
$\downarrow$\\
\texttt{Sparse Gauss nodes (66 or 286 points)}\\
$\downarrow$\\
\texttt{Patch Abaqus .inp material cards}\\
$\downarrow$\\
\texttt{abaqus job=PCE-SXXXX interactive} $\times$ 66 (or 286)\\
$\downarrow$\\
\texttt{abaqus python extract\_pce\_qoi.py (odbAccess)}\\
$\downarrow$\\
\texttt{chaospy.fit\_regression()} $\rightarrow$ PCE coefficients\\
$\downarrow$\\
\texttt{Sobol indices / MC ($10^6$) / NN comparison}
\end{center}

---

# Results

## PCE Statistics: Convergence at Degree 2

\begin{center}
\begin{tabular}{lccccc}
\toprule
QoI & Degree & Mean & Std & CV & $P_f$ \\
\midrule
$\vM$ & 2 (66 pts) & 209.7 MPa & 2.38 MPa & 1.13\% & 0 \\
$\vM$ & 3 (286 pts) & 209.7 MPa & 2.38 MPa & 1.13\% & 0 \\
\midrule
$u^\text{max}$ & 2 (66 pts) & 10.07 mm & 0.178 mm & 1.76\% & -- \\
$u^\text{max}$ & 3 (286 pts) & 10.07 mm & 0.178 mm & 1.76\% & -- \\
\bottomrule
\end{tabular}
\end{center}

\vspace{0.5em}
\begin{block}{Key result}
Statistics identical between degree 2 and 3 $\Rightarrow$ \textbf{converged at degree 2}\\
$P_f = 0$: peak stress 209.7 MPa $\ll$ 1,200 MPa failure threshold
\end{block}

---

## Sobol Sensitivity Indices

\begin{columns}
\begin{column}{0.5\textwidth}
\includegraphics[width=\textwidth]{figures/pce_uq/sobol_max_smises.pdf}
\end{column}
\begin{column}{0.5\textwidth}
\includegraphics[width=\textwidth]{figures/pce_uq/sobol_max_disp.pdf}
\end{column}
\end{columns}

\vspace{0.3em}
\begin{block}{Finding}
$E_1$ dominates \textbf{over 99\%} of output variance for all QoIs.\\
CZM parameters ($K_n$, $G_{Ic}$, $t_n$): negligible ($\approx 0$) — adhesive never damages.\\
No parameter interactions ($S_i \approx S_T^i$).
\end{block}

---

## PDF Comparison: PCE vs NN Surrogate

\begin{columns}
\begin{column}{0.5\textwidth}
\includegraphics[width=\textwidth]{figures/pce_uq/pdf_max_smises.pdf}
\end{column}
\begin{column}{0.5\textwidth}
\includegraphics[width=\textwidth]{figures/pce_uq/pdf_max_disp.pdf}
\end{column}
\end{columns}

\vspace{0.3em}
Both PCE and NN surrogates agree with Abaqus samples.\\
PCE: analytical PDF via moment propagation. NN: MC sampling required.

---

## PCE Convergence

\begin{columns}
\begin{column}{0.5\textwidth}
\includegraphics[width=\textwidth]{figures/pce_uq/convergence_max_smises.pdf}
\end{column}
\begin{column}{0.5\textwidth}
\includegraphics[width=\textwidth]{figures/pce_uq/convergence_max_disp.pdf}
\end{column}
\end{columns}

\vspace{0.3em}
Mean and variance converge rapidly. Degree 2 is sufficient.

---

## Degree 2 vs Degree 3: LOO Accuracy

\begin{center}
\begin{tabular}{lcccc}
\toprule
QoI & PCE deg=2 & PCE deg=3 & NN (66 pts) & NN (286 pts) \\
\midrule
$\vM$ (MPa) & \textbf{9.4e-4} & 1.5e-3 & 7.2e-4 & 1.8e-2 \\
$u^\text{max}$ (mm) & \textbf{4.2e-5} & 7.1e-5 & 2.7e-5 & 1.6e-3 \\
\bottomrule
\end{tabular}
\end{center}

\vspace{0.5em}
\begin{alertblock}{Counterintuitive NN result}
NN accuracy \textbf{degrades} with more data (66 $\to$ 286 points).\\
Fixed MLP architecture overfits the near-linear response surface.
\end{alertblock}

\begin{block}{PCE advantage}
Degree 2 PCE is optimal. Degree 3 shows mild over-parameterisation.
\end{block}

---

# Discussion

## Physical Interpretation

\textbf{Why does $E_1$ dominate ($S_{E_1} > 0.99$)?}

\begin{itemize}
  \item Fairing panel: thin-walled under membrane/bending loads
  \item Response governed by membrane stiffness $E_1 \cdot t$
  \item $E_1$ CoV = 5\% is relatively small, but it is the \textit{only} stiffness-path parameter
  \item $\Rightarrow$ Quality control of fiber-direction modulus is the key manufacturing lever
\end{itemize}

\vspace{0.5em}
\textbf{Why is SDEG identically zero?}

\begin{itemize}
  \item Design-point static load — adhesive operates far below cohesive strength $t_n$
  \item Even at extreme CoV bounds, interface stress never reaches damage threshold
  \item CZM sensitivity requires higher load or impact scenario
\end{itemize}

---

## PCE vs NN: Complementary Tools

\begin{center}
\small
\begin{tabular}{lll}
\toprule
Criterion & \textbf{PCE (this work)} & NN surrogate \\
\midrule
Data required & 66 pts (deg=2) & 66--286 pts \\
LOO accuracy & $\sim 10^{-4}$ & $10^{-4}$ to $10^{-2}$ \\
Sobol indices & Analytical, exact & Requires separate MC \\
Interpretability & High & Low (black box) \\
Uncertainty type & \textbf{Aleatoric} (manufacturing) & \textbf{Epistemic} (model) \\
Use case & Forward UQ, reliability & SHM, defect detection \\
\bottomrule
\end{tabular}
\end{center}

\vspace{0.3em}
\begin{block}{Complementarity with existing GNN-SHM work}
\begin{itemize}
  \item GNN-SHM: ``Is there a defect?'' (epistemic uncertainty via MC Dropout)
  \item PCE-SFEM: ``How does manufacturing scatter affect response?'' (aleatoric)
  \item Independent and complementary — design phase vs. monitoring phase
\end{itemize}
\end{block}

---

# Conclusion

## Summary

\begin{enumerate}
  \item \textbf{PCE degree 2 (66 Abaqus runs) is sufficient} — degree 3 (286 runs) gives identical statistics
  \vspace{0.3em}
  \item \textbf{$E_1$ uncertainty dominates all QoIs} ($S_{E_1} > 0.99$)\\
  $\Rightarrow$ CFRP fiber-direction modulus is the key quality control target
  \vspace{0.3em}
  \item \textbf{Structure is reliable} under design load: $P_f = 0$, $\beta = \infty$\\
  (209.7 MPa peak stress vs. 1,200 MPa threshold, CV = 1.13\%)
  \vspace{0.3em}
  \item \textbf{PCE outperforms NN} for this smooth, near-linear response\\
  — better accuracy, analytical Sobol indices, no overfitting
  \vspace{0.3em}
  \item \textbf{PCE and GNN-SHM are complementary}: aleatoric vs. epistemic uncertainty
\end{enumerate}

\vspace{0.5em}
\begin{block}{Future work}
Higher load cases (impact, acoustic fatigue) to activate CZM damage and reveal CZM parameter sensitivity
\end{block}

---

## References

\small
1. Ghanem \& Spanos (1991). \textit{Stochastic Finite Elements: A Spectral Approach}. Springer.
2. Sudret (2008). Global sensitivity analysis using PCE. \textit{Rel. Eng. \& Sys. Safety}, 93(7).
3. Smolyak (1963). Quadrature and interpolation formulas. \textit{Dokl. Akad. Nauk SSSR}.
4. Saltelli et al. (2008). \textit{Global Sensitivity Analysis: The Primer}. Wiley.
5. Dassault Systèmes (2024). \textit{Abaqus Analysis User's Guide}, v2024.

\vspace{1em}
\begin{center}
\Large Thank you for your attention.\\
\vspace{0.5em}
\normalsize Questions?
\end{center}
