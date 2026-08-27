# 確率有限要素法 — JAXA H3 ロケットフェアリングの不確実性定量化（UQ）

[English](README.md) | **日本語**

![Course](https://img.shields.io/badge/Course-Stochastic%20FEM%20(SoSe%202026)-blue)
![Institution](https://img.shields.io/badge/LUH-IBNM-004b7c)
![Solver](https://img.shields.io/badge/FEM-Abaqus%2FStandard-orange)
![UQ](https://img.shields.io/badge/UQ-Non--intrusive%20PCE-green)
![Python](https://img.shields.io/badge/Python-chaospy%20%7C%20scikit--learn-3776ab)

> **非介入型多項式カオス展開（Non-intrusive PCE）を実在の航空宇宙構造に適用:**
> JAXA H3 ロケット・ペイロードフェアリングの CFRP/アルミハニカム・サンドイッチパネルを対象に、
> 順方向 UQ・ベイズ逆問題・Karhunen–Loève 確率場・サロゲートモデル比較の4本立てで解析。

![SFEM パイプライン](figures/hero_sfem_pipeline.png)

---

## ハイライト

- **Abaqus シミュレーション 352 回**を非介入的にオーケストレーション（Smolyak 疎格子: 次数2で66点、次数3で286点）。ソルバーの改修は一切不要。
- **支配パラメータはただ一つ:** 繊維方向弾性率 E₁ が応力分散の **99.2 %**、変位分散の **99.95 %** を説明（1次 Sobol 指数）。
- **全352回で破損ゼロ**（SDEG = 0）: P_f ≈ 0、信頼性指標 β → ∞。フェアリングの保守的設計余裕を定量的に裏付け。
- **次数収束を確認:** PCE 次数2と3で平均・分散が有効数字4桁一致。Leave-one-out RMSE ≈ 10⁻⁴ MPa。
- 順方向 UQ に加えて**3つの発展課題**を実施: ベイズ更新（Topic 6）、KL 確率場モデリング（Topic 9）、PCE / GP / NN サロゲート比較（Topic 13）。

## 主要結果一覧

| 解析 | 手法 | 主要結果 |
|---|---|---|
| **順方向 UQ** | PCE 次数 2/3 + Smolyak 疎求積 | 最大 von Mises 応力に対し S₁(E₁) = 0.992、P_f ≈ 0、β → ∞ |
| **Topic 6 — ベイズ逆問題** | PCE サロゲート尤度による共役ガウス更新 | ノイズ付き計測10点で E₁ の事後標準偏差が 7,300 → 1,885 MPa（**74 %** 削減） |
| **Topic 9 — 非ガウス確率場** | KL 展開 + 対数正規変換場 | 空間平均化で実効パネル剛性の標準偏差が 7,300 → **約4,700 MPa** に低下。平均支配の応答についてスカラーモデルは保守側（入力側の結果） |
| **Topic 13 — サロゲート比較** | PCE vs ガウス過程 vs ニューラルネット、N = 10–286 | N=66 で GP 7.3×10⁻⁴ vs NN 3.6×10⁻¹ MPa（約500倍差）。PCE 次数2は構成上厳密 — 参照応答がその模型クラスに属するため |

<p align="center">
  <img src="figures/fairing_odb_mises.png" width="49%" alt="フェアリングパネルの von Mises 応力場">
  <img src="figures/fairing_overview3d.png" width="49%" alt="H3 フェアリングモデルの3D概観">
</p>

## 解析モデル

**構造:** JAXA H3 ペイロードフェアリングの CFRP 表皮／アルミハニカムコア・サンドイッチパネル。代表的な上昇時荷重条件を想定し、**Abaqus/Standard** で解析（表皮–コア界面は凝集力要素でモデル化）。

**不確実パラメータ**（独立な5変数）:

| パラメータ | 平均 | 変動係数 | 説明 |
|---|---|---|---|
| E₁ | 146 000 MPa | 5 % | 繊維方向弾性率（CFRP） |
| G₁₂ | 5 200 MPa | 10 % | 面内せん断弾性率 |
| Kₙ | 1.0 × 10⁶ N/mm³ | 15 % | 凝集力法線剛性 |
| G_Ic | 0.3 N/mm | 20 % | Mode-I 破壊エネルギー |
| tₙ | 50 MPa | 15 % | 凝集力引張強度 |

**注目量（QoI）:** 最大 von Mises 応力、最大変位、最大凝集損傷（SDEG）。

## リポジトリ構成

```
├── report_paper.tex / .pdf          # 本レポート（elsarticle テンプレート、約20ページ）
├── References.bib                   # レポートの参考文献
├── slides_sfem.tex / .pdf           # 発表スライド（Beamer）
├── summary_ja.tex / .pdf            # 日本語要約
│
├── pce_driver.py                    # PCE オーケストレーション: 疎格子設計 → Abaqus ジョブ
├── extract_pce_qoi.py               # Abaqus ODB からの QoI 抽出
├── reliability_analysis.py         # 破損確率・信頼性指標の評価
├── topic6_bayesian.py               # Topic 6: E1 のベイズ同定
├── topic9_kl_expansion.py           # Topic 9: 確率場 E1(x,y) の KL 展開
├── topic13_surrogate_comparison.py  # Topic 13: PCE / GP / NN 収束比較
│
├── make_report_figs.py              # レポート用 UQ 図の生成
├── generate_figures.py              # 追加の可視化
├── make_visuals.py                  # ヒーロー図・パイプライン図
├── thesis_style.py                  # Matplotlib スタイル（LaTeX フォント、LUH 段組幅）
├── figures/                         # 生成した全図（アニメーション含む）
│
├── SFEM_2026#*.pdf                  # 講義資料（Intro 〜 Collocation）
├── SFEM_Exercise_*.pdf              # 演習問題・解答
├── stochastic_fem_cheatsheet.md     # 理論チートシート
└── oral_exam_qa.md                  # 口頭試問 Q&A 対策
```

## ポスト処理の再現方法

Abaqus の ODB ファイル（各約293 MB）は LUH の計算サーバー上にあります。抽出済み QoI テーブル以降の処理はすべてローカルで実行できます:

```bash
pip install chaospy numpy scipy matplotlib scikit-learn

python topic6_bayesian.py              # ベイズ更新の図
python topic9_kl_expansion.py          # KL 固有値分解とサンプル場
python topic13_surrogate_comparison.py # サロゲート収束比較
python make_report_figs.py             # Sobol / PDF / 収束プロット
```

LaTeX フォントによる図のレンダリング（`thesis_style.py`）と `report_paper.tex` のビルドには TeX Live 2025 が必要です。

## 手法の要約

応答 u(ξ) を標準化した確率入力 ξ の多変量直交多項式で展開し、評価点を **Smolyak 疎求積格子**上に配置します。疎な Gauss–Hermite 則は一般に負の重みを生じるため、係数は求積射影ではなく**最小二乗回帰**で当てています（66点に対し未知数21で優決定。Leave-one-out 検証も可能）。FE ソルバーの*入力ファイルと出力データベース*にしか触れないため、完全に**非介入型** — Abaqus をブラックボックスとして利用できます。統計量（平均・分散）、大域感度（Sobol 指数）、出力の確率密度は PCE 係数から*解析的に*得られ、同じサロゲートが Topic 6 のベイズ尤度評価をほぼゼロコストで駆動します。

## コース情報

| | |
|---|---|
| **講義** | Stochastic Finite Element Methods（2026年夏学期） |
| **機関** | ハノーファー大学（Leibniz Universität Hannover）— IBNM |
| **指導教員** | Prof. Dr.-Ing. Udo Nackenhorst, Dr. Zhibao Zheng |
| **評価** | セメスタープロジェクト + レポート + 口頭発表（5 ECTS、成績評価あり） |

## ライセンス / 免責

大学の課題として作成したものです。構造モデルは公開情報に基づく教育目的の簡易表現であり、JAXA の非公開データには**一切**基づいていません。
