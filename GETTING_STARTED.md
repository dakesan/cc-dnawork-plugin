# Getting Started with cc-dnawork-plugin

このガイドでは、cc-dnawork-plugin をインストールして、DNA 関連のバイオインフォマティクス分析を開始する方法を説明します。

## Prerequisites（前提条件）

- **Claude Code**: Latest version installed
- **Python**: 3.9 以上（3.12 推奨）
- **uv**: Python パッケージマネージャー

## Installation（インストール）

### Step 1: Claude Code のセットアップ

まだ Claude Code がインストールされていない場合：

**macOS:**
```bash
curl -fsSL https://claude.ai/install.sh | bash
```

**Windows (PowerShell):**
```powershell
irm https://claude.ai/install.ps1 | iex
```

### Step 2: uv のインストール

uv は Python パッケージ管理ツールです。

**macOS / Linux:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows (PowerShell):**
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

**インストール確認:**
```bash
uv --version
```

### Step 3: プラグインのインストール

Claude Code で以下を実行：

```bash
/plugin marketplace add dakesan/cc-dnawork-plugin
```

### Step 4: プラグインの有効化

1. Claude Code を開く
2. `/plugin list` を実行
3. **dnawork-skills** を選択
4. **Install** をクリック

これで完了です！

## Quick Start（クイックスタート）

### 例1: DNA配列解析

```
Analyze this DNA sequence for:
1. CpG islands and regulatory elements
2. ORF (Open Reading Frames) prediction
3. GC content and sequence statistics
4. Search for homologous sequences in UniProt database
5. Create visualizations of the results

Use BioPython, gget, and UniProt database integration.
```

### 例2: RNA-seq データ解析

```
Process this RNA-seq dataset:
1. Load data with Scanpy
2. Quality control and filtering
3. Normalize the counts
4. Identify cell types
5. Differential expression analysis with PyDESeq2
6. Visualize results with Seaborn
7. Map DEGs to KEGG pathways
```

### 例3: バリアント注釈付け

```
Analyze this VCF file:
1. Parse with pysam
2. Annotate variants with Ensembl database
3. Check pathogenicity in ClinVar
4. Find cancer-related mutations in COSMIC
5. Search PubMed for disease associations
6. Generate a clinical report
```

### 例4: 分子ドッキング

```
Perform molecular docking:
1. Query ChEMBL for candidate compounds
2. Prepare PDB structure with BioPython
3. Analyze drug-likeness with RDKit and MedChem
4. Perform docking with DiffDock
5. Visualize binding poses
6. Create a summary table
```

## Available Skills（利用可能なスキル）

### DNA配列解析（6個）
- biopython, pysam, scikit-bio, bioservices, gget, gtars

### 単一細胞・RNA-seq（4個）
- scanpy, cellxgene-census, pydeseq2, arboreto

### ゲノミクスツール（4個）
- etetoolkit, deeptools, geniml, esm

### 化学・分子設計（6個）
- rdkit, datamol, deepchem, diffdock, medchem, molfeat

### ゲノムデータベース（14個）
- Ensembl, NCBI Gene, UniProt, PDB, PubMed, ClinVar, COSMIC, ChEMBL, PubChem, ZINC, DrugBank, KEGG, Reactome, STRING

### 可視化・解析（4個）
- matplotlib, seaborn, plotly, networkx

### その他（26個）
- 科学コミュニケーション、ラボ統合、ドキュメント処理など

完全なリストは [README.md](README.md) を参照してください。

## Common Workflows（一般的なワークフロー）

### ワークフロー1: 遺伝子同定とアノテーション

```
1. BioPython で FASTA ファイルをパース
2. gget で遺伝子情報を取得
3. Ensembl で遺伝子構造を確認
4. UniProt でタンパク質機能を調査
5. STRING で相互作用パートナーを検索
6. KEGG で生化学経路をマッピング
7. Seaborn で可視化
```

### ワークフロー2: 多層オミクス解析

```
1. Scanpy で単一細胞データを読込
2. PyDESeq2 で発現差異を検出
3. Arboreto で遺伝子制御ネットワークを推定
4. HMDB で代謝物をマッピング
5. STRING でタンパク質相互作用を統合
6. NetworkX でネットワークトポロジーを解析
7. Matplotlib で複合可視化を作成
```

### ワークフロー3: 薬物発見パイプライン

```
1. ChEMBL で既知阻害剤を検索
2. RDKit で構造-活性相関を分析
3. Datamol で新規類似体を設計
4. DiffDock でバーチャルドッキング実施
5. DeepChem で ADMET 特性を予測
6. MedChem で薬物適性を評価
7. 結果をレポートに集約
```

### ワークフロー4: 臨床バリアント解釈

```
1. pysam で VCF ファイルをパース
2. Ensembl で機能的影響を予測
3. ClinVar で病原性を確認
4. COSMIC で癌関連変異を検索
5. ClinPGx で薬理ゲノミクスを確認
6. PubMed で関連論文を検索
7. ClinicalTrials.gov でマッチング試験を見つける
8. 臨床レポートを生成
```

## Tips & Best Practices（ヒントとベストプラクティス）

### 1. 大量データの処理

Scanpy、Polars、Dask などの大規模データ処理ツールを使用：

```
Use Scanpy for large single-cell datasets (>100,000 cells)
Use Polars for fast tabular data processing
Use Dask for out-of-core computation
```

### 2. 可視化の最適化

出版品質の図を作成：

```
Use matplotlib for publication-ready static figures
Use Seaborn for statistical visualizations
Use Plotly for interactive visualizations
Use NetworkX + matplotlib for network visualization
```

### 3. API rate limits の処理

```
Implement caching for database queries
Use batch requests when available
Check API documentation for rate limits
```

### 4. ファイル形式の選択

- **FASTA**: DNA/RNA/タンパク質配列
- **GenBank**: アノテーション付き配列
- **SAM/BAM**: アライメント結果
- **VCF**: バリアント情報
- **H5AD**: 単一細胞データ（Scanpy）
- **CSV/TSV**: 表形式データ

## Troubleshooting（トラブルシューティング）

### 問題: スキルが見つからない

```bash
# プラグイン一覧を確認
/plugin list

# プラグインを再インストール
/plugin marketplace remove dakesan/cc-dnawork-plugin
/plugin marketplace add dakesan/cc-dnawork-plugin
```

### 問題: 依存関係が足りない

```bash
# パッケージをインストール
uv pip install package-name

# または特定のバージョン
uv pip install package-name==1.0.0

# インストール確認
uv pip list | grep package-name
```

### 問題: メモリ不足

大規模データセットの処理：

```python
# Scanpy で分割読込
adata = sc.read_h5ad('file.h5ad', first_n_obs=10000)

# Polars で遅延読込
df = pl.scan_csv('large_file.csv').collect()

# Dask で分散処理
import dask.dataframe as dd
df = dd.read_csv('large_file.csv')
```

### 問題: ネットワーク接続

オフラインでも使用可能なスキル：
- BioPython, pysam, scikit-bio (ローカルファイル処理)
- matplotlib, seaborn, plotly (可視化)
- scikit-learn (機械学習)

オンライン必須のスキル：
- すべてのデータベーススキル
- 文献検索スキル

## Getting Help（ヘルプを得る）

1. **個別スキルのドキュメント**: `scientific-skills/[skill-name]/SKILL.md`
2. **K-Dense Scientific Skills**: https://github.com/K-Dense-AI/claude-scientific-skills
3. **Claude Code Documentation**: https://docs.claude.com/
4. **個別パッケージのドキュメント**: 各スキルの references/ フォルダ

## Next Steps（次のステップ）

- [README.md](README.md) で全スキル一覧を確認
- 個別のスキルドキュメントを探索
- サンプルデータで試験的に分析を実行
- 自分の研究データで実際の分析を開始

---

Happy bioinformatics! 🧬
