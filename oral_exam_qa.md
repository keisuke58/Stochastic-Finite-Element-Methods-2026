# Stochastic FEM (#08) — Oral Exam Q&A Prep

Examiners: Dr. Z. Zheng / Prof. U. Nackenhorst. Format: report submitted, then oral defense.
Study: non-intrusive PCE surrogate of an H3 fairing CFRP/Al-honeycomb sandwich panel; 5 uncertain
material params; QoIs = max von Mises stress, max SDEG, max displacement.

> **試験概要（和訳）**：審査=Zheng博士／Nackenhorst教授。形式=レポート提出後に口頭ディフェンス。
> 研究内容=H3フェアリングのCFRP/アルミハニカム・サンドイッチパネルに対する**非侵入型PCE**サロゲート。
> 不確かな材料パラメータ5個。評価量(QoI)=最大ミーゼス応力・最大SDEG（損傷変数）・最大変位。

Key numbers to have on the tip of your tongue:
- d = 5 inputs, degree p = 2 → **P = C(d+p, p) = C(7,2) = 21** PCE terms.
- Smolyak sparse grid: **66** Abaqus runs (deg 2), **286** (deg 3).
- σ_vM^max = **209.7 ± 2.375 MPa** (CoV 1.13%); u^max = 10.07 ± 0.178 mm; SDEG ≡ 0.
- Sobol: **S_E1 = 0.992** (E1 explains 99% of variance); G12 ≈ 0.008; rest ≈ 0.
- Failure limit 1200 MPa → margin ≈ 416σ → **P_f ≈ 0, β → ∞**.

> **即答すべき数値（和訳）**：
> - 入力 d=5、次数 p=2 → **項数 P=C(7,2)=21**。
> - Smolyakスパースグリッド：**66回**のAbaqus実行（deg2）、**286回**（deg3）。
> - 最大ミーゼス応力 σ_vM=**209.7±2.375 MPa**（変動係数CoV 1.13%）／最大変位 u=10.07±0.178 mm／SDEG≡0。
> - Sobol指数：**S_E1=0.992**（E1が分散の99%を説明）／G12≈0.008／他≈0。
> - 破壊限界1200 MPa → 余裕≈416σ → **破壊確率 P_f≈0、信頼性指標 β→∞**。

---

## A. PCE fundamentals

**Q1. What is polynomial chaos and why Hermite polynomials here?**
PCE expands the response as `Y ≈ Σ_α c_α Ψ_α(ξ)` in polynomials orthogonal w.r.t. the input
probability measure. My inputs are Gaussian (truncated), so the matching orthogonal basis is
**Hermite** polynomials — this is the Wiener–Askey scheme: Gaussian↔Hermite, uniform↔Legendre, etc.
Matching the basis to the measure gives the fastest (spectral) convergence.
*(What they test: do you know the basis is dictated by the input distribution, not arbitrary.)*

> **【和訳】Q1. 多項式カオス(PCE)とは何か。なぜここでHermite多項式か？**
> PCEは応答を `Y ≈ Σ_α c_α Ψ_α(ξ)` と、**入力の確率測度に関して直交する多項式**で展開する手法。
> 本研究の入力はガウス分布（打ち切り）なので、対応する直交基底は**Hermite多項式**になる。これは
> Wiener–Askeyスキーム（ガウス↔Hermite、一様↔Legendre 等）に従う。基底を測度に合わせると
> **収束が最速（スペクトル収束）**になる。
> *（試験の狙い：基底は入力分布で決まる必然であって任意ではない、と理解しているか。）*

**Q2. Why does orthogonality matter — what do you get for free?**
Because `E[Ψ_α Ψ_β] = 0` for α≠β:
- **Mean** = `c_0` (the zeroth coefficient).
- **Variance** = `Σ_{α≠0} c_α² E[Ψ_α²]`.
- **Sobol indices** follow analytically by grouping the squared coefficients per variable — **no extra
model runs**. That's the headline advantage over Monte Carlo or a black-box NN.

> **【和訳】Q2. 直交性がなぜ重要か。何が「タダで」得られるか？**
> α≠β で `E[Ψ_α Ψ_β]=0` が成り立つため：
> - **平均** = `c_0`（0次の係数そのもの）。
> - **分散** = `Σ_{α≠0} c_α² E[Ψ_α²]`。
> - **Sobol指数**も、係数の二乗を変数ごとにまとめるだけで解析的に出る → **追加のモデル実行ゼロ**。
> これがモンテカルロや黒箱NNに対する最大の利点。

**Q3. How many terms, and how does that scale?**
Total terms `P = C(d+p, p)`. For d=5, p=2 → 21; p=3 → 56. It grows polynomially in p but
combinatorially in d — that's the curse of dimensionality the sparse grid fights.

> **【和訳】Q3. 項数はいくつで、どうスケールするか？**
> 総項数 `P=C(d+p,p)`。d=5, p=2 → 21、p=3 → 56。次数pに対しては多項式的だが、次元dに対しては
> **組合せ的に爆発**する。これが次元の呪いで、スパースグリッドが対抗する相手。

**Q4. Non-intrusive vs intrusive PCE — which did you use and why?**
**Non-intrusive**: I treat Abaqus as a black box, evaluate it at quadrature nodes, and fit/project the
coefficients. Intrusive (stochastic Galerkin) would require rewriting the FE equations — impossible
with a commercial solver. Non-intrusive reuses the existing high-fidelity model unchanged.

> **【和訳】Q4. 非侵入型 vs 侵入型PCE — どちらを使ったか、なぜ？**
> **非侵入型**を使用。Abaqusを黒箱とみなし、求積点で評価して係数を回帰/射影で求める。
> 侵入型（確率Galerkin）はFE方程式そのものの書き換えが必要で、商用ソルバーでは不可能。
> 非侵入型なら既存の高忠実度モデルを**無改変のまま**再利用できる。

---

## B. Smolyak / quadrature

**Q5. Why Smolyak sparse quadrature instead of a full tensor grid?**
A full tensor grid needs `(p+1)^d` points — for 5 dims this explodes. Smolyak combines low-order
1-D rules to keep accuracy for the total-degree polynomial space while using far fewer points (66 vs
hundreds–thousands). It exploits that I only need to integrate total-degree-≤p polynomials, not full
tensor-degree ones.

> **【和訳】Q5. なぜフルテンソルグリッドでなくSmolyakスパース求積か？**
> フルテンソルグリッドは `(p+1)^d` 点を要し、5次元では爆発する。Smolyakは低次の1次元則を組み合わせ、
> **総次数（total-degree）多項式空間**の精度を保ちつつ、はるかに少ない点で済む（66点 対 数百〜数千点）。
> 「総次数≤pの多項式だけ積分できればよく、テンソル次数まで要らない」性質を利用している。

**Q6. 66 points for degree 2 — how do you fit 21 coefficients from 66 runs?**
Over-determined, which is good: I can fit by **regression** (least squares) or by **quadrature
projection** when the Smolyak weights are usable. The extra points give robustness and let me do
leave-one-out cross-validation for the surrogate error.

> **【和訳】Q6. deg2で66点なのに、なぜ21係数を66回の実行から決められるのか？**
> **優決定系（方程式が未知数より多い）**なので、むしろ好都合。**回帰（最小二乗）**で当てるか、Smolyak
> の重みが使えるときは**求積射影**で求める。余分な点は頑健性を与え、サロゲート誤差の
> **Leave-One-Out交差検証（LOO-CV）**も可能にする。

---

## C. Degree / convergence

**Q7. How do you know degree 2 is enough?**
Two checks: (1) **moment convergence** — mean and std are unchanged to four significant figures
between degree 2 (66 runs) and degree 3 (286 runs); (2) **LOO-CV error** is already minimal at
degree 2 and does **not** improve at degree 3 — it slightly worsens, a sign of over-parameterization
(fitting 56 coefficients to 286 points adds variance in redundant high-order terms). The response is
nearly linear in the inputs, so a low order suffices.
*(Trap: don't say "higher degree is always better." Defend the model-order choice.)*

> **【和訳】Q7. なぜ次数2で十分と言えるのか？**
> 2つの確認：(1)**モーメント収束** — 平均と標準偏差がdeg2（66回）とdeg3（286回）で**有効数字4桁まで一致**。
> (2)**LOO-CV誤差**がすでにdeg2で最小で、deg3にしても**改善しない**（むしろ微増）。これは過剰パラメータ化の
> サイン（286点に56係数を当てると冗長な高次項に分散が乗る）。応答が入力にほぼ線形なので低次で足りる。
> *（罠：「高次ほど良い」と言ってはダメ。モデル次数の選択を論理で守ること。）*

---

## D. Sensitivity — the E1 story

**Q8. Why does E1 dominate 99% of the variance? Give a first-order argument.**
A first-order (linear) sensitivity gives variance contribution `≈ (∂Y/∂x_i)² σ_i²`. Two reasons E1 wins:
1. **Largest absolute CoV-weighted effect**: the panel response under bending/membrane load is
governed by the **fiber-direction stiffness E1** of the skins; the QoI is the peak skin stress.
2. The cohesive params (Kn, GIc, tn) have `∂Y/∂x_i ≈ 0` here because **the bond never damages**
(SDEG≡0) — they don't move any output, so they cannot contribute variance regardless of their CoV.
So effectively the 5-D problem collapses to ~1-D in E1.

> **【和訳】Q8. なぜE1が分散の99%を支配するのか。一次の議論で説明せよ。**
> 一次（線形）感度では分散寄与は `≈ (∂Y/∂x_i)² σ_i²`。E1が勝つ理由は2つ：
> 1. **CoV重み付き効果が最大**：曲げ/膜荷重下のパネル応答は皮膚の**繊維方向剛性E1**が支配し、QoIは皮膚の
>    ピーク応力だから。
> 2. 凝集要素パラメータ(Kn, GIc, tn)はここで `∂Y/∂x_i≈0`。理由は**接着が一切損傷しない**（SDEG≡0）から。
>    出力を一切動かさない以上、CoVがいくら大きくても分散に寄与できない。
> よって実質的に5次元問題はE1の**ほぼ1次元**に縮退する。

**Q9. First-order vs total Sobol — were they different?**
No — S_i ≈ S_T^i for all inputs. That means **interaction effects are negligible** and the response is
essentially additive (consistent with the near-linear, degree-2-sufficient finding).

> **【和訳】Q9. 一次Sobol指数と全Sobol指数は違ったか？**
> いいえ。全入力で S_i ≈ S_T^i。これは**交互作用効果が無視できる**＝応答がほぼ**加法的**を意味する
> （ほぼ線形・deg2で十分という結果と整合）。

---

## E. Reliability

**Q10. You report P_f = 0 and β = ∞. Is that physical?**
No — read it as **"below the resolution of the analysis," not a literal guarantee.** Two limits:
1. Monte Carlo on the surrogate (10^6 samples) can only resolve `P_f` down to ~10^-6.
2. More fundamentally, my inputs are **truncated at ±50%**, so the surrogate is only valid in that
range; genuine tail risk lives outside it where the polynomial must not be extrapolated.
The honest statement is: the mean sits ~416σ below the 1200 MPa limit, so failure is far outside the
modeled range.

> **【和訳】Q10. P_f=0、β=∞ と報告したが、それは物理的か？**
> いいえ。**「解析の分解能以下」であって文字通りの保証ではない**と読むべき。限界は2つ：
> 1. サロゲート上のモンテカルロ（10^6 サンプル）では `P_f` は ~10^-6 までしか分解できない。
> 2. より本質的に、入力を**±50%で打ち切っている**ので、サロゲートはその範囲でのみ有効。真の裾(tail)リスクは
>    その外側にあり、多項式を外挿してはならない。
> 正直な言い方：平均は1200 MPa限界の約416σ下にあり、破壊はモデル化範囲のはるか外側、ということ。

**Q11. Definition of β and link to P_f?**
`β = -Φ^{-1}(P_f)` (reliability index = number of standard deviations to the limit state for a Gaussian
margin). Large β ↔ tiny P_f.

> **【和訳】Q11. βの定義とP_fとの関係は？**
> `β = -Φ^{-1}(P_f)`（信頼性指標＝ガウス余裕において限界状態までの標準偏差の個数）。
> βが大きい ↔ P_fが極小。

---

## F. The SDEG = 0 result (a favorite probe)

**Q12. SDEG is identically zero. What does that imply?**
At this static design load the adhesive stays fully elastic — **no cohesive damage initiates** in any of
the 352 runs. Methodologically this means the three cohesive parameters (Kn, GIc, tn) are
**structurally unidentifiable**: since they never change any output, no sampling can assign them a
non-zero Sobol index. The nominally 5-D UQ problem is effectively 1-D (E1) at this load.

> **【和訳】Q12. SDEGが恒等的にゼロ。これは何を意味するか？**
> この静的設計荷重では接着剤は完全に弾性のまま — **352回の全実行で凝集損傷が一切開始しない**。
> 手法的には、3つの凝集パラメータ(Kn, GIc, tn)が**構造的に同定不能**であることを意味する。出力を一切
> 変えない以上、どんなサンプリングでも非ゼロのSobol指数を割り当てられない。名目上5次元のUQ問題は、
> この荷重では実質1次元（E1）になる。

**Q13. So was including them a mistake?**
No — it's a **finding**, not a bug: it tells the designer the bond is not load-limiting here, and it
flags that to *study the adhesive you must drive it into damage* — a higher load or a damage-onset
load case. That's exactly the future work.

> **【和訳】Q13. では、それらを含めたのは間違いだったのか？**
> いいえ。バグではなく**知見**。設計者に「ここでは接着が荷重制限要因ではない」と伝え、さらに
> 「**接着を研究したいなら損傷まで追い込む必要がある**（より高い荷重、または損傷開始荷重ケース）」
> ことを示す。それがまさに今後の課題（future work）。

---

## G. Input modeling

**Q14. Why truncated normal, not plain Gaussian or lognormal?**
Stiffnesses and fracture energies are **strictly positive**; an unbounded Gaussian assigns probability
to non-physical negative draws. I truncate at [0.5μ, 1.5μ] because flight hardware is screened against
incoming-inspection limits, so extreme outliers aren't representative. Lognormal is an alternative for
positivity; truncated normal keeps the Hermite-PCE machinery clean and is justified by the inspection
bounds.

> **【和訳】Q14. なぜ素のガウスや対数正規でなく、打ち切り正規分布か？**
> 剛性や破壊エネルギーは**厳密に正値**。非有界なガウスでは非物理的な負値に確率を割り当ててしまう。
> [0.5μ, 1.5μ] で打ち切るのは、フライト品は受入検査の上下限でスクリーニングされ、極端な外れ値は
> 代表的でないから。正値性なら対数正規も代替案だが、打ち切り正規はHermite-PCEの仕組みを綺麗に保て、
> 検査限界という根拠もある。

**Q15. Are the inputs independent? What if they're correlated?**
I assume independence (standard for a first study, and the Hermite tensor basis assumes it).
Correlation would require either decorrelation (e.g., Nataf/Rosenblatt transform) or a basis adapted to
the joint measure. A known limitation.

> **【和訳】Q15. 入力は独立か？相関があったらどうなる？**
> 独立を仮定（最初の研究としては標準で、Hermiteのテンソル基底も独立を前提とする）。
> 相関がある場合は、**Nataf/Rosenblatt変換**などで無相関化するか、同時測度に適合した基底が必要になる。
> 既知の限界点。

---

## H. CZM / FE model

**Q16. Explain the cohesive zone model and the BK criterion.**
The adhesive interface is modeled with a traction–separation law: linear elastic up to strength `t_n`,
then softening governed by fracture energy `G_Ic`. Mixed-mode damage evolution uses the
**Benzeggagh–Kenane** law with exponent η = 2.284, which interpolates the critical energy between
mode I and mode II. Stiffness `K_n` is the initial interface penalty stiffness.

> **【和訳】Q16. 凝集域モデル(CZM)とBK基準を説明せよ。**
> 接着界面は**牽引力–分離(traction–separation)則**でモデル化：強度 `t_n` までは線形弾性、その後は破壊
> エネルギー `G_Ic` が支配する軟化。混合モードの損傷発展には指数 η=2.284 の**Benzeggagh–Kenane(BK)則**を
> 使い、モードIとモードIIの間で臨界エネルギーを補間する。`K_n` は界面の初期ペナルティ剛性。

**Q17. What are the QoIs and their thresholds?**
Peak von Mises stress in the CFRP skin (limit 1200 MPa), peak SDEG in the adhesive (limit 0.5), and
peak displacement. Extracted from the last frame of the two-step (pressure + incremental) analysis via
`odbAccess`.

> **【和訳】Q17. QoIとその閾値は？**
> CFRP皮膚のピークミーゼス応力（限界1200 MPa）、接着剤のピークSDEG（限界0.5）、ピーク変位。
> 2ステップ解析（圧力＋増分）の最終フレームから `odbAccess` で抽出する。

---

## I. PCE vs NN

**Q18. Your NN did worse. Does PCE beat NN in general?**
No — and I'm careful to say so. PCE wins **here** because it encodes a strong, *correct* prior: the
response is a low-order polynomial. With only 66–286 samples on a smooth, near-linear map, that prior
extracts accurate statistics and gives Sobol indices for free. The NN (thousands of weights) overfits
on so few samples — its LOO error grows an order of magnitude from 66 to 286 points. It's a
**match-between-method-and-problem** statement, not a universal ranking.

> **【和訳】Q18. NNの方が悪かった。PCEは一般にNNに勝るのか？**
> いいえ — そしてそう言い切らないよう注意している。PCEが**ここで**勝つのは、強くて*正しい*事前知識
> 「応答は低次多項式」を埋め込んでいるから。滑らかでほぼ線形な写像に対し66〜286サンプルしかなくても、
> その事前知識で正確な統計量を引き出し、Sobol指数もタダで得る。NN（数千の重み）はこれほど少ない
> サンプルでは過剰適合し、LOO誤差は66点→286点で一桁悪化した。これは**手法と問題の相性**の話であって、
> 普遍的な優劣ではない。

---

## J. Random variables vs random fields (lecture-linked probe)

**Q19. The course covered random *fields* (KL expansion). You used 5 random *variables* — why?**
Spatially-varying material fields would be the next step (KL expansion to reduce a field to a few
modes). I modeled each property as a single random variable, i.e., a perfectly correlated field, which
is appropriate when the dominant manufacturing variability is **batch/coupon-level** (a panel cured
together has one effective E1) rather than within-panel spatial fluctuation. Adding a random field for
E1 is a natural extension and would raise the effective dimension (handled again by sparse PCE on the
KL modes).

> **【和訳】Q19. 講義は確率*場*（KL展開）を扱った。あなたは5つの確率*変数*を使った — なぜ？**
> 空間的に変動する材料場は次のステップ（KL展開で場を少数モードに縮約）。本研究では各物性を単一の
> 確率変数、すなわち**完全相関の場**としてモデル化した。これは支配的な製造ばらつきがパネル内の空間変動
> ではなく**バッチ/クーポン単位**である場合に妥当（同時硬化したパネルは実効E1が1つ）。E1に確率場を
> 加えるのは自然な拡張で、実効次元が上がる（これも**KLモードに対するスパースPCE**で扱える）。

---

## K. Limitations / future work (close strong)

**Q20. Biggest limitations?**
1. Single static load case; SDEG=0 ⇒ cohesive params unidentifiable — need a damage-driving load.
2. Inputs truncated ±50% ⇒ no genuine tail / extreme-value statistics.
3. Independence + single-RV-per-property (no spatial random field yet).
4. The honeycomb core is a homogenized orthotropic continuum, not resolved cells.

> **【和訳】Q20. 最大の限界は？**
> 1. 単一の静的荷重ケース。SDEG=0 ⇒ 凝集パラメータが同定不能 — 損傷を駆動する荷重が必要。
> 2. 入力を±50%で打ち切り ⇒ 真の裾/極値統計が得られない。
> 3. 独立仮定＋物性ごとに単一確率変数（空間的な確率場はまだ無し）。
> 4. ハニカムコアは均質化した直交異方性連続体で、セルを解像していない。

**Q20b. What's the future-work roadmap? (staged)**
The single unifying weakness is that the response is *too benign* — linear, SDEG≡0, P_f trivially 0,
cohesive params unidentifiable. The roadmap below progressively makes failure *reachable*, which is
what gives the analysis teeth:
1. **CZM damage onset (static).** Reuse the existing cohesive model but apply a load that drives
   SDEG > 0 — this finally makes K_n, G_Ic, t_n move the output, so they earn non-zero Sobol indices.
   Cheapest entry point.
2. **Hashin progressive damage → ultimate failure-load distribution.** Add in-ply Hashin damage
   (fiber/matrix × tension/compression) so the QoI becomes the *failure load itself*, now a random
   variable. Turns the trivial P_f = 0 into a *meaningful* reliability problem; highest leverage.
3. **Adaptive sampling for the tail (subset / importance sampling).** Crucial coupling: once (1)/(2)
   introduce a failure bifurcation, the response is **non-smooth**, so the low-order-polynomial PCE
   assumption breaks. Either go to adaptive / multi-element PCE, or fall back to sampling — subset
   simulation / importance sampling — to resolve the small-P_f tail. So sampling returns *as a
   necessary consequence* of adding damage, not as a separate axis.
4. **Fatigue delamination life (VCCT + Paris).** Largest scope: VCCT for energy-release-rate G at the
   crack front, Paris-type da/dN vs ΔG under cyclic load → delamination *life*. A paper on its own;
   do last.
Plus the input-modeling fixes throughout: correlated CFRP constants (Nataf), KL random fields, and
PCE-based tolerance allocation (tighten E1 QC — the highest-leverage action).

> **【和訳】Q20b. 今後のロードマップは？（段階的）**
> 全限界を貫く一つの弱点は「応答が良性すぎる」こと — 線形・SDEG≡0・P_f が自明に0・凝集パラメータ
> 同定不能。以下のロードマップは**破壊を到達可能にしていく**もので、それが解析に意味（牙）を与える：
> 1. **CZM損傷開始（静的）**：既存の凝集モデルのまま、SDEG>0 を起こす荷重を与える。これで初めて
>    K_n, G_Ic, t_n が出力を動かし、非ゼロのSobol指数を得る。最も安い入口。
> 2. **Hashin進行性損傷 → 最終破壊荷重の分布**：層内Hashin損傷（繊維/母材×引張/圧縮）を入れ、QoIを
>    **破壊荷重そのもの**＝確率変数にする。自明な P_f=0 を**意味のある信頼性問題**に変える。最大レバレッジ。
> 3. **裾のための適応サンプリング（subset / 重点サンプリング）**：重要な連結点 — (1)(2)が破壊分岐を
>    入れると応答は**非平滑**になり、低次多項式というPCEの前提が壊れる。適応PCE/multi-element PCEに
>    上げるか、サンプリング（subset simulation・重点サンプリング）に戻して小さい P_f の裾を解く。
>    サンプリングは別軸でなく、**損傷を入れた帰結として必然的に**復活する。
> 4. **疲労による層間剥離寿命（VCCT + Paris）**：最大スコープ。VCCTでき裂先端のエネルギー解放率Gを取り、
>    繰返し荷重下で da/dN vs ΔG（Paris則）→ 剥離**寿命**。それ単体で論文1本。最後にやる。
> 並行して入力モデリングの改善：CFRP定数の相関（Nataf）、KL確率場、PCEベースの公差配分
> （E1の品質管理を締める＝最もレバレッジの高い施策）。

---

## 日本語・要点（直前確認）

- **PCE**：入力がガウス→Hermite基底（Wiener–Askey）。直交性で 平均=c0、分散=Σc_α²E[Ψ²]、**Sobolは係数から解析的（追加計算ゼロ）**。
- **項数** P=C(d+p,p)=C(7,2)=**21**。Smolyakで**66点**（deg2）/286（deg3）。フルテンソル(p+1)^dの回避。
- **非侵入型**：Abaqusは黒箱。回帰orquadratureで係数。
- **deg2十分**：モーメント4桁一致＋LOO最小（deg3は過剰適合で悪化）。応答ほぼ線形。
- **E1支配0.992**：一次感度≈(∂Y/∂x)²σ²。皮膚応力は繊維方向剛性E1支配。CZMは損傷ゼロ→∂Y/∂≈0で寄与不能。S_i≈S_T（交互作用なし＝加法的）。
- **Pf=0/β=∞**：「解析の分解能以下」と言う。MC 10^-6 限界＋±50%打ち切りで外挿不可。平均は限界の~416σ下。
- **SDEG≡0**：接着は弾性のまま→Kn,GIc,tnは**構造的に同定不能**。バグでなく知見（高荷重/損傷開始ケースが必要）。
- **truncated normal**：剛性・破壊エネは正値＋検査上下限。lognormalも可。独立仮定は限界（相関ならNataf）。
- **PCE>NN**は一般論でなく「問題に合致した事前知識」のため。NNは少標本で過剰適合。
- **罠回避**：「高次ほど良い」と言わない／Pf=0を文字通り保証と言わない／PCEがNNに普遍的に勝るとは言わない。
