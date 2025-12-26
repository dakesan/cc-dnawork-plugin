# Project Structure

cc-dnawork-plugin のプロジェクト構造と概要です。

## Directory Layout

```
cc-dnawork-plugin/
├── .claude-plugin/
│   └── marketplace.json          # Claude Code plugin configuration
├── scientific-skills/
│   ├── sequence-io/              # 本番スキル
│   ├── blast-search/             # 本番スキル
│   ├── bam-toolkit/              # 本番スキル
│   ├── vcf-toolkit/              # 本番スキル
│   ├── igv-integration/          # 本番スキル
│   └── cosmic-toolkit/           # 本番スキル
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
| `blast-search/` | 本番 | NCBI BLAST 配列類似性検索 |
| `blat-api-searching/` | 本番 | BLAT ゲノムマッピング |
| `bam-toolkit/` | 本番 | BAM/SAM/CRAM アライメントファイル操作 |
| `vcf-toolkit/` | 本番 | VCF/BCF バリアントファイル操作 |
| `igv-integration/` | 本番 | IGV 自動スナップショット生成 |
| `cosmic-toolkit/` | 本番 | COSMIC がん遺伝子データベース照合 |

## Files

### Configuration Files

**`.claude-plugin/marketplace.json`**
- Claude Code plugin configuration
- 64 個のスキルパスを定義
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

### blast-search

NCBI BLAST による配列類似性検索。BioPython 実装。

| 機能 | ツール |
|------|--------|
| Web BLAST | Bio.Blast.NCBIWWW (qblast) |
| 結果パース | Bio.Blast.NCBIXML |
| JSON 出力 | scripts/run_blast_biopython.py |

**スクリプト:**
- `scripts/run_blast_biopython.py` - BioPython qblast + JSON 出力

### blat-api-searching

UCSC BLAT API によるゲノムマッピング。高速配列アライメント。

| 機能 | ツール |
|------|--------|
| ローカル BLAT | pxblat (C 拡張) |
| JSON 出力 | scripts/run_local_blat.py |
| BLAST URL API | scripts/run_blat_url.py |

**スクリプト:**
- `scripts/run_local_blat.py` - pxblat によるローカル実行
- `scripts/run_blat_url.py` - UCSC BLAT URL API

### vcf-toolkit

VCF/BCF バリアントファイルの統計情報計算、フィルタリング、JSON 出力。WGS/WES 解析結果の確認と品質管理。

| 機能 | ツール |
|------|--------|
| 統計情報計算 | scripts/vcf_stats.py |
| VCF フィルタリング | scripts/filter_vcf.py |
| JSON エクスポート | scripts/inspect_vcf.py |
| 基盤ライブラリ | pysam |

**スクリプト:**
- `scripts/vcf_stats.py` - バリアント統計計算（品質、深度、AF）
- `scripts/filter_vcf.py` - VCF フィルタリング（VCF 出力）
- `scripts/inspect_vcf.py` - バリアント抽出（JSON 出力）

### bam-toolkit

BAM/SAM/CRAM アライメントファイルの解析。リード抽出、indel 検出、カバレッジ統計計算。WGS/WES 解析結果の確認と品質管理。

| 機能 | ツール |
|------|--------|
| リード抽出 | scripts/extract_reads.py |
| Indel 検出 | scripts/extract_indels.py |
| カバレッジ統計 | scripts/calculate_coverage.py |
| 基盤ライブラリ | pysam |

**スクリプト:**
- `scripts/extract_reads.py` - 特定領域のリード抽出（BAM/JSON 出力）
- `scripts/extract_indels.py` - 挿入・欠失の抽出と集計
- `scripts/calculate_coverage.py` - カバレッジ統計計算（mean, median）

### igv-integration

IGV (Integrative Genomics Viewer) の自動化。複数 BAM ファイルのバッチ可視化とスナップショット生成。

| 機能 | ツール |
|------|--------|
| スナップショット生成 | scripts/generate_igv_snapshots.py |
| バッチモード実行 | IGV batch script |
| 複数領域処理 | BED ファイル入力 |
| 基盤ツール | IGV, Java |

**スクリプト:**
- `scripts/generate_igv_snapshots.py` - IGV batch script 生成と実行、PNG スナップショット出力

### cosmic-toolkit

COSMIC Cancer Gene Census データベース照合。がん遺伝子の確認と注釈付け。

| 機能 | ツール |
|------|--------|
| がん遺伝子クエリ | scripts/query_cosmic_genes.py |
| Cancer Gene Census 照合 | CSV 動的読み込み |
| JSON 出力 | 全フィールド保持 |
| 基盤ライブラリ | pandas |

**スクリプト:**
- `scripts/query_cosmic_genes.py` - 遺伝子名から Cancer Gene Census 情報を取得（JSON 出力）

**データ:**
- `data/cancer_gene_census.csv` - COSMIC から手動ダウンロード（ユーザー提供）

---

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
| スキル数 | 125+ → 7 | 7 本番（厳選） |
| 設計方針 | ツール単位 | 責任単位 |
| SKILL.md サイズ | 大きい | 小さい（<400行） |

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
