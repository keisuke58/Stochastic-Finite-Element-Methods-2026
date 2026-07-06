# #08 Stochastic FEM — 戦略ノート

> 🟢 **ほぼ完了** — レポート・スライド・コード全て済み。口頭発表の日程調整のみ。
> 担当: Dr. Zhibao Zheng (IBNM)

## 成績評価
| 種別 | 形式 | 評価 | ECTS |
|------|------|------|------|
| PL | Semester Project + Report (max 20p) + Oral Presentation | 成績あり（graded） | 5 |

筆記試験なし。**プロジェクトの質が全て。**

## ✅ 現在の状況（2026-05-26時点）

| 成果物 | 状況 |
|--------|------|
| プロジェクトテーマ | ✅ **PCE UQ for JAXA H3 CFRP/Al-Honeycomb Fairing** |
| コード実装 | ✅ `pce_driver.py`, `extract_pce_qoi.py`, `reliability_analysis.py` |
| Abaqusシミュレーション | ✅ 352回（deg=2: 66回, deg=3: 286回） |
| 数値結果 | ✅ S₁(E₁)=0.992 (stress), 0.9995 (disp), Pf=0, β=∞ |
| 図・グラフ | ✅ Sobol / PDF / Convergence（deg=2 & deg=3） |
| **Report** | ✅ **`report.md` + `report.pdf` 完成**（~2500 words） |
| **Slides** | ✅ **`slides.md` + `slides.pdf` 完成** |
| 口頭発表 | ❌ **日程未調整** |

**GitHub**: `keisuke58/Stochastic-Finite-Element-Methods-2026`（パブリック）

## ⚡ 残り作業（〜8/31）

- [ ] Dr. Zheng に口頭発表の日程を相談（メール送付）
- [ ] Report の最終確認・PDF化（Overleafで仕上げ）
- [ ] 発表スライドのリハーサル（15〜20分）
- [ ] Q&A対策（下記参照）

## Q&Aに備える重要ポイント

```
Q: なぜNon-intrusiveなのか？
A: Abaqusのソースコードを変更せずに使えるから。Intrusive PCEはソルバー改修が必要で現実的でない。

Q: E₁だけが支配的なのはなぜ？
A: CFRP は繊維方向の剛性が圧倒的。繊維方向弾性率E₁のばらつきが直接応力分散を支配。

Q: Pf=0, β=∞ の意味は？
A: 352回全シミュレーションでSDEG=0（破損なし）。破損限界から十分遠い。保守設計の裏付け。

Q: deg=2 vs deg=3 の違いは？
A: 統計量（平均・分散）が4桁一致。PCEの次数収束確認済み。deg=2で十分。

Q: モンテカルロと比べた利点は？
A: 66点で同精度（MC比）。計算コストが大幅削減。Smolyak疎格子で次元の呪い緩和。
```

## プロジェクト概要メモ（発表前に再確認）

```
不確実パラメータ（5個）:
  E₁ (繊維方向弾性率), G₁₂ (面内せん断剛性),
  Kn (接触剛性), GIc (Mode-I破壊靱性), tn (引張強度)

手法:
  Non-intrusive PCE + Smolyak sparse quadrature
  deg=2: 66 points, deg=3: 286 points
  Abaqus/Standard + chaospy

主要結果:
  ・E₁がほぼ全分散を説明（第一Sobol指数S≈0.99）
  ・応力/変位ともに正規分布に近い
  ・LOO RMSE ≈ 10⁻⁴ MPa（deg=2でも高精度）
  ・Pf=0, 設計余裕大
```

## 演習について（注意）

**Exercise 1–4はプロジェクト提出物ではない。** 演習は授業の課題（SL相当）であり、Semester Projectは別に必要。
ただし今回は keisuke58 リポジトリのオリジナル研究がそのまま提出物になるため問題なし。
