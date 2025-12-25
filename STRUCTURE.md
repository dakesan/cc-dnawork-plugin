# Project Structure

cc-dnawork-plugin のプロジェクト構造と概要です。

## Directory Layout

```
cc-dnawork-plugin/
├── .claude-plugin/
│   └── marketplace.json          # Claude Code plugin configuration
├── scientific-skills/
│   ├── sequence-io/              # 本番スキル（整理済み）
│   │   ├── SKILL.md
│   │   └── references/
│   ├── inbox/                    # 未整理スキル（68個）
│   │   ├── biopython/
│   │   ├── pysam/
│   │   ├── scanpy/
│   │   └── ...
│   └── archived/                 # アーカイブ（統合元など）
├── README.md                      # プロジェクト概要（英語）
├── GETTING_STARTED.md            # 使用開始ガイド（日本語）
├── STRUCTURE.md                  # このファイル
├── LICENSE                       # MIT ライセンス
└── .gitignore                    # Git除外ファイル
```

## Skill Organization Strategy

スキルは責任範囲を明確にするため、段階的に整理中です。

| ディレクトリ | 状態 | 説明 |
|-------------|------|------|
| `sequence-io/` | 本番 | FASTA/GenBank/FASTQ の読み書き |
| `inbox/` | 未整理 | K-Dense からの移行スキル（68個） |
| `archived/` | アーカイブ | 統合済みの旧スキル |

## Files

### Configuration Files

**`.claude-plugin/marketplace.json`**
- Claude Code plugin configuration
- 68 個のスキルパスを定義
- Plugin メタデータを含む

### Documentation

**`README.md`**
- プロジェクト概要（英語）
- スキル一覧とカテゴリ
- インストール方法
- 使用例とワークフロー

**`GETTING_STARTED.md`**
- クイックスタートガイド（日本語）
- インストール詳細手順
- よくあるワークフロー
- トラブルシューティング

**`STRUCTURE.md`**
- このファイル
- プロジェクト構成の説明

**`LICENSE`**
- MIT ライセンス全文

**`.gitignore`**
- Git から除外するファイル定義

## Production Skills（本番スキル）

### sequence-io

FASTA/GenBank/FASTQ ファイルの読み書きと配列操作に特化。

| 機能 | ツール |
|------|--------|
| ファイル読み書き | Bio.SeqIO |
| 配列操作 | Bio.Seq |
| インデックスアクセス | pysam.FastaFile + faidx |
| 配列統計 | Bio.SeqUtils (GC%, Tm, MW) |

**参照ファイル:**
- `references/biopython_seqio.md` - Bio.Seq, Bio.SeqIO
- `references/faidx.md` - pysam による高速アクセス
- `references/formats.md` - フォーマット仕様
- `references/utilities.md` - 統計計算

---

## Inbox Skills（未整理スキル）

以下は `inbox/` 内のスキル一覧です。順次整理予定。

### 1. DNA Sequence Analysis（DNA配列解析）
- `biopython/` - Comprehensive molecular biology toolkit (→ sequence-io に統合予定)
- `pysam/` - SAM/BAM/VCF file processing
- `scikit-bio/` - Biological sequence operations
- `bioservices/` - Biological web services
- `gget/` - Genome information retrieval
- `gtars/` - Genomic tools and resources

**Use cases:**
- FASTA/GenBank ファイル処理
- 配列比較とアライメント
- ORF 予測
- 配列統計計算

### 2. Single-Cell & RNA-seq（単一細胞・RNA-seq）
- `scanpy/` - Single-cell RNA-seq analysis
- `cellxgene-census/` - Large-scale cell integration
- `pydeseq2/` - Differential expression
- `arboreto/` - Gene regulatory networks

**Use cases:**
- 単一細胞データの品質管理
- 細胞タイプ同定
- 発現差異検出
- 遺伝子制御ネットワーク推定

### 3. Genomics Tools（ゲノミクスツール）
- `etetoolkit/` - Phylogenetic analysis
- `deeptools/` - Genomic signal processing
- `geniml/` - Machine learning for genomics
- `esm/` - Protein language models

**Use cases:**
- 系統樹解析
- ゲノム領域の可視化
- 機械学習による特性予測
- タンパク質機能予測

### 4. Chemistry & Molecular Design（化学・分子設計）
- `rdkit/` - Cheminformatics
- `datamol/` - Molecular data processing
- `deepchem/` - Deep learning for chemistry
- `diffdock/` - Molecular docking
- `medchem/` - Drug-likeness assessment
- `molfeat/` - Molecular feature computation

**Use cases:**
- 分子構造操作
- 薬物類似性予測
- バーチャルスクリーニング
- リード最適化

### 5. Genomic Databases（ゲノムデータベース）
14 個のデータベーススキル：
- `alphafold-database/` - Protein structures
- `ensembl-database/` - Genome annotations
- `gene-database/` - NCBI Gene
- `uniprot-database/` - Protein information
- `pdb-database/` - Protein 3D structures
- `pubmed-database/` - Literature search
- `clinvar-database/` - Variant pathogenicity
- `cosmic-database/` - Cancer mutations
- `chembl-database/` - Bioactive compounds
- `pubchem-database/` - Chemical structures
- `zinc-database/` - Virtual screening library
- `kegg-database/` - Metabolic pathways
- `reactome-database/` - Pathway database
- `string-database/` - Protein interactions

**Use cases:**
- 遺伝子情報取得
- バリアント注釈付け
- パスウェイマッピング
- 文献検索

### 6. Visualization & Analysis（可視化・解析）
- `matplotlib/` - Publication-quality plots
- `seaborn/` - Statistical visualization
- `plotly/` - Interactive visualization
- `networkx/` - Network analysis

**Use cases:**
- グラフ作成
- 統計可視化
- インタラクティブプロット
- ネットワーク可視化

### 7. Scientific Communication（科学コミュニケーション）
10 個のスキル：
- `literature-review/` - Literature synthesis
- `scientific-writing/` - Research writing
- `scientific-visualization/` - Figure creation
- `citation-management/` - Reference management
- `research-lookup/` - Research discovery
- `scientific-brainstorming/` - Idea generation
- `hypothesis-generation/` - Hypothesis development
- `clinical-decision-support/` - Clinical workflows
- `clinical-reports/` - Report generation
- `market-research-reports/` - Market analysis

**Use cases:**
- 論文作成サポート
- 仮説生成
- 文献レビュー
- レポート生成

### 8. Laboratory Integration（ラボ統合）
- `benchling-integration/` - Lab workflow platform
- `dnanexus-integration/` - Cloud genomics
- `latchbio-integration/` - Bioinformatics platform
- `omero-integration/` - Microscopy data
- `opentrons-integration/` - Liquid handling
- `protocolsio-integration/` - Protocol repository
- `labarchive-integration/` - Electronic lab notebook

**Use cases:**
- ラボワークフロー自動化
- クラウドコンピューティング
- データ管理
- 機器制御

### 9. Supporting Skills（補助スキル）
- `biomni/` - Multi-omics integration
- `denario/` - Multi-omics workflow
- `generate-image/` - AI image generation
- `research-grants/` - Grant writing
- `market-research-reports/` - Market research
- `document-skills/` - Document processing (DOCX, PDF, PPTX, XLSX)
- `get-available-resources/` - Resource discovery

## Skill Metadata

各スキルディレクトリの構成：

```
skill-name/
├── SKILL.md              # スキル説明とドキュメント
├── references/           # (optional) 参考資料
│   ├── package_docs.md   # パッケージドキュメント
│   ├── examples.md       # 使用例
│   └── ...
├── scripts/              # (optional) 実行可能なスクリプト
│   └── script.py         # Python スクリプト
└── assets/               # (optional) テンプレートやリソース
    └── template.py       # テンプレートファイル
```

### SKILL.md の構造

```yaml
---
name: skill-name
description: "Simple one-line description"
---

# Skill Title

## Overview
Comprehensive description

## When to Use This Skill
When to trigger this skill

## Core Capabilities
- Capability 1
- Capability 2

## Installation and Setup
Installation instructions

## Key Examples
Practical usage examples
```

## Design Principles

### Single Responsibility

各スキルは単一の責任範囲を持つべき。

**悪い例（旧 biopython）:**
- 配列 I/O + BLAST + 構造解析 + 系統解析 + ... → 3,730 行

**良い例（sequence-io）:**
- FASTA/GenBank/FASTQ の読み書きのみ → 1,267 行

### Progressive Disclosure

1. **Metadata** (~100 words) - 常にコンテキストに
2. **SKILL.md** (<5k words) - トリガー時に読み込み
3. **references/** - 必要に応じて参照

## Integration with K-Dense Scientific Skills

このプロジェクトは K-Dense-AI/claude-scientific-skills をベースに、責任範囲を明確化したバージョンです。

### アプローチの違い

| 項目 | K-Dense | cc-dnawork |
|------|--------|-----------|
| スキル数 | 125+ | 1 本番 + 68 inbox |
| 設計方針 | ツール単位 | 責任単位 |
| SKILL.md サイズ | 大きい | 小さい（<150行） |

## How to Use

### インストール後

1. **スキルの確認**
   ```bash
   /plugin list
   ```

2. **スキルの探索**
   - `scientific-skills/[skill-name]/SKILL.md` を開く
   - `references/` フォルダで詳細を確認

3. **実行**
   ```
   Use [skill-name] to [task description]
   ```

## Version Information

- **Version**: 2.0.0 (restructured)
- **Base**: K-Dense Scientific Skills v2.10.1
- **Python**: 3.9+
- **License**: MIT

## Contributing

改善提案やバグ報告：

1. Issue を作成
2. Pull request を提交

## Related Projects

- [K-Dense Scientific Skills](https://github.com/K-Dense-AI/claude-scientific-skills)
- [Claude Code](https://github.com/anthropics/claude-code)
- [Anthropic SDK](https://github.com/anthropics/anthropic-sdk-python)

---

For detailed information, see [README.md](README.md) and [GETTING_STARTED.md](GETTING_STARTED.md)
