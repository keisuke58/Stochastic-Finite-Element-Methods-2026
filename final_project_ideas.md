# SFEM ファイナルプロジェクト

作成: 2026-04-13

---

## 決定テーマ

### 製造ばらつきを考慮した CFRP/Al-Honeycomb フェアリングの PCE サロゲートと信頼性評価
#### — Abaqus × chaospy × PyTorch による非侵入型確率的有限要素解析 —

**対象構造:** JAXA H3 ロケット CFRP/Al-Honeycomb フェアリング（既存モデル流用）
**言語:** Python 3 + Abaqus Python (embedded)
**ソルバー:** Abaqus 2024（静解析 Step-2: 熱＋機械荷重 Max-Q 条件）

---

## 既存プロジェクトとの関係

| 既存研究 (Payload2026) | 本プロジェクト (SFEM課題) |
|----------------------|------------------------|
| GNN による欠陥検出 (SHM) | PCE サロゲートによる UQ |
| MC Dropout / Deep Ensemble（GNN の認識論的不確かさ） | PCE + Sobol 指数（入力パラメータ感度） |
| 欠陥位置・サイズのランダム化 | **材料特性・CZM パラメータのランダム化** |
| `uncertainty.py` | `pce_driver.py` + `reliability_analysis.py` |
| `manufacturing_variability.py` | **これを PCE の確率入力として定式化** |

---

## 問題設定

### 物理モデル

Abaqus 静解析（2 ステップ）:
- **Step-1**: 熱荷重のみ（CTE ミスマッチ）
- **Step-2**: 熱＋機械荷重（Max-Q 条件）
  - 差圧: 0.005 MPa
  - 打上加速度: 3g（29,430 mm/s²）

### 不確かな入力変数（5変数）

| 記号 | 物理量 | 公称値 | CoV | 分布 | PCE 基底 |
|------|--------|--------|-----|------|---------|
| E1 | CFRP 繊維方向ヤング率 [MPa] | 160,000 | 5% | 切断正規 | Hermite |
| G12 | CFRP 面内せん断弾性係数 [MPa] | 5,000 | 10% | 切断正規 | Hermite |
| Kn | 接着剤 CZM 法線剛性 [MPa/mm] | 100,000 | 15% | 切断正規 | Hermite |
| GIc | Mode-I 破壊エネルギー [N/mm] | 0.3 | 20% | 切断正規 | Hermite |
| tn | Mode-I 引張強度 [MPa] | 50 | 15% | 切断正規 | Hermite |

**派生変数（名目比を保持）:**
- E2 = 10,000 × (E1/160,000), G13 = G12, G23 = 3,000 × (G12/5,000)
- Ks = Kn/2, ts = 0.8 × tn, GIIc = GIc / 0.3

### 関心量（QoI）

| QoI | 記号 | 意味 | 破壊判定閾値 |
|-----|------|------|------------|
| 最大 von Mises 応力 | σ_max [MPa] | CFRP 面板の応力集中 | 1,200 MPa（引張強度/SF=2） |
| 最大損傷変数 | SDEG_max [-] | 接着層デボンディング度 | 0.5（部分デボンディング） |
| 最大変位 | u_max [mm] | フェアリング変形量 | — |

---

## 手法・パイプライン

```
manufacturing_variability.py の CoV 定義
            ↓
  build_joint_distribution() [chaospy]
  TruncNormal × 5変数
            ↓
  generate_quadrature(degree=3)
  → 35 点（Gauss 求積点）
            ↓
  modify_inp() で .inp マテリアルカード書き換え
            ↓
  abaqus job=PCE-S000N interactive   × 35回
            ↓
  abaqus python extract_pce_qoi.py
  → max_smises, max_sdeg, max_disp
            ↓
  fit_quadrature() → PCE サロゲート
            ↓
  ┌─────────────────┬──────────────────┐
  │ PCE 解析的統計量  │  NN サロゲート比較 │
  │ E[Q], Std[Q]    │  (同 35 点で学習)  │
  │ Sobol 指数       │  RMSE LOO 比較    │
  └─────────────────┴──────────────────┘
            ↓
  MC on surrogate (N=100,000)
  → P(σ_max > 1200 MPa), P(SDEG > 0.5)
  → 信頼性指標 β = -Φ⁻¹(Pf)
```

### 必要な Abaqus 実行回数

| 手法 | Abaqus 実行数 |
|------|-------------|
| Gauss 求積 (degree=3, 5変数) | **35 回** |
| 回帰ベース PCE (Halton) | ~70 回 |
| MC 参照解 | ~500 回（収束確認用） |

---

## ファイル構成

```
Payload2026/src/
├── pce_driver.py          ← NEW: PCE メインドライバー
├── extract_pce_qoi.py     ← NEW: Abaqus Python QoI 抽出
├── reliability_analysis.py ← NEW: 信頼性解析 + 図生成
├── manufacturing_variability.py  ← 既存: CoV 定義（流用）
├── material_properties.py        ← 既存: 公称値（流用）
└── uncertainty.py                ← 既存: GNN UQ（比較対象）
```

### 実行手順（サーバー上）

```bash
cd /home/nishioka/Payload2026

# 1. chaospy インストール確認
pip install chaospy

# 2. ドライラン（Abaqus 未実行、サンプル点のみ確認）
python3 src/pce_driver.py \
    --template abaqus_work/batch_s12_100/Job-S12-D001/Job-S12-D001.inp \
    --workdir  abaqus_work/pce_uq \
    --degree 3 --rule gaussian --dry_run

# 3. 本番実行（35 回 Abaqus）
python3 src/pce_driver.py \
    --template abaqus_work/batch_s12_100/Job-S12-D001/Job-S12-D001.inp \
    --workdir  abaqus_work/pce_uq \
    --degree 3 --rule gaussian --cpus 8

# 4. 信頼性解析・図生成
python3 src/reliability_analysis.py \
    --results abaqus_work/pce_uq/pce_results.json \
    --outdir  figures/pce_uq
```

---

## 期待される成果物

### 図
1. **Sobol 感度指数バーチャート** — どの製造ばらつきが σ_max / SDEG に最も影響するか
2. **PDF 比較** — MC参照 vs PCE代理モデル vs NN代理モデル
3. **収束プロット** — MC サンプル数 vs 平均・標準偏差の収束
4. **信頼性指標** — Pf と β の比較表

### 考察ポイント
1. **支配パラメータ**: E1 vs GIc vs tn — Sobol 指数で定量化
2. **PCE vs NN**: 35 点という少サンプルでどちらが精度高いか（RMSE LOO）
3. **製造ばらつきの影響**: CV=20% の GIc がデボンディング確率に与える影響
4. **既存 GNN（認識論的不確かさ）との対比**: 入力不確かさ vs モデル不確かさ

---

## Action Items
- [x] pce_driver.py 作成
- [x] extract_pce_qoi.py 作成
- [x] reliability_analysis.py 作成
- [ ] サーバーで dry_run 実行・サンプル点確認
- [ ] degree=2 で小規模テスト（15点）→ 結果確認
- [ ] degree=3 本番実行（35点）
- [ ] reliability_analysis.py で図生成
- [ ] レポート執筆
