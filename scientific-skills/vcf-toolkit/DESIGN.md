# VCF Toolkit Skill Design

## Purpose

VCF/BCF ファイルの操作に特化したスキル。bcftools を活用したフィルタリング、小規模データの JSON 変換、統計情報の取得を提供。

## Scope（責任範囲）

### 含めるもの ✅

- **VCF/BCF フィルタリング** - bcftools を使った高速フィルタリング
- **JSON 変換** - 小規模 VCF（≤100 エントリ）を JSON に
- **VCF 統計** - バリアント数、品質分布などの基本統計

### 除外するもの ❌

- **BAM/アラインメント操作** - pysam スキルに任せる
- **大規模データの JSON 化** - メモリ効率のため制限
- **高度なアノテーション** - 専用ツール（SnpEff など）に任せる

## Design Principles

1. **bcftools 活用** - 標準ツールを使い、車輪の再発明を避ける
2. **typer で柔軟な引数** - フィルタ条件を CLI で指定可能
3. **JSON 出力** - BLAT/BLAST と同様、他ツールと連携しやすい形式
4. **データ量制限** - 100 エントリ超えたらエラー（JSON 変換時）

### JSON 形式を選択した理由

タブ区切り形式と比較検討した結果、JSON 形式を採用：

**トークン効率**：
- タブ区切り（VCF 形式）の方がトークン数は少ない（約1/3）
- しかし、100エントリ制限があるためトークン爆発は防げる

**構造化データの必要性**：
- INFO フィールドは複数のキー=値ペア（DP=50;AF=0.5;AC=1）
- サンプルデータは階層構造（各サンプル × 複数フィールド）
- ALT は配列（複数の alternate allele）
- タブ区切りでは文字列として保存され、下流でパース必要

**AI による処理**：
- JSON は self-descriptive（列の意味が明確）
- Claude による解析時にエラーが少ない
- プログラムで処理しやすい

**結論**: トークン効率よりも構造化データの利点を優先

## Scripts 仕様

### inspect_vcf.py - VCF フィルタリング & JSON 出力

**目的**: VCF をフィルタリングして全エントリを JSON で出力

**設計原則**:
- VCF の全カラムを漏れなく JSON 化
- デフォルト 100 エントリ制限（`--force` で無視可能）
- pysam の VariantFile で読み込み
- フィルタ条件は Python で適用（小規模データ想定）

**引数**:
```
# 必須
--vcf PATH              入力 VCF ファイル
--chrom TEXT            染色体指定（例: chr1）※--region と排他
--region TEXT           領域指定（例: chr1:10000-20000）※--chrom と排他

# フィルタ条件（全てオプショナル）
--min-qual FLOAT        最小品質スコア（QUAL >= X）
--min-dp INT            最小深度（INFO/DP >= X）
--min-af FLOAT          最小アレル頻度（INFO/AF >= X）
--max-af FLOAT          最大アレル頻度（INFO/AF <= X）
--pass-only             FILTER=PASS のみ（デフォルト: True）
--all-filters           全フィルタを含む（--pass-only を無効化）

# 制限
--max-variants INT      最大バリアント数（デフォルト: 100）
--force                 エントリ数制限を無視（フラグ）

# 出力
--output PATH           JSON 出力パス（オプション、未指定時は標準出力）
```

**注意**:
- `--chrom` または `--region` のどちらか一方は必須
- `--pass-only` はデフォルト有効（全フィルタを含む場合は `--all-filters` を指定）

**使用例**:
```bash
# Test 1: chr1 の PASS フィルタのみ（デフォルト）
uv run python scripts/inspect_vcf.py \
  --vcf test-data/sample.vcf.gz \
  --chrom chr1 \
  --output test1_chr1_pass.json
# → 9 variants

# Test 2: chr1 の全フィルタを含む
uv run python scripts/inspect_vcf.py \
  --vcf test-data/sample.vcf.gz \
  --chrom chr1 \
  --all-filters \
  --output test2_chr1_all.json
# → 10 variants (LowQual も含む)

# Test 3: chr1 で品質スコア >= 90
uv run python scripts/inspect_vcf.py \
  --vcf test-data/sample.vcf.gz \
  --chrom chr1 \
  --min-qual 90 \
  --output test3_chr1_qual90.json
# → 7 variants

# Test 4: 領域指定（chr1:10000-14000）
uv run python scripts/inspect_vcf.py \
  --vcf test-data/sample.vcf.gz \
  --region chr1:10000-14000 \
  --output test4_region.json
# → 6 variants

# 制限を無視（大量出力）
uv run python scripts/inspect_vcf.py \
  --vcf huge.vcf.gz \
  --chrom chr1 \
  --force \
  --output large.json
```

**出力形式（JSON）**:
```json
{
  "num_variants": 9,
  "samples": ["sample1", "sample2"],
  "variants": [
    {
      "chrom": "chr1",
      "pos": 10177,
      "id": "rs367896724",
      "ref": "A",
      "alts": ["AC"],
      "qual": 100.0,
      "filter": ["PASS"],
      "info": {
        "DP": 50,
        "AF": [0.5],
        "AC": [1]
      },
      "samples": {
        "sample1": {"GT": [0, 1], "DP": 25, "GQ": 99},
        "sample2": {"GT": [0, 0], "DP": 25, "GQ": 99}
      }
    }
  ]
}
```

**注意**: pysam によってパースされた形式で出力される：
- `GT` は配列形式：`[0, 1]`（元の VCF: `0/1`）
- `AF` などの INFO フィールドも配列形式：`[0.5]`（元の VCF: `AF=0.5`）
- 全てのフィールドは適切な型に変換される（文字列 → 数値/配列）

**エラーハンドリング**:
```bash
# ケース 1: 制限超過（1,234 エントリ）
$ python scripts/inspect_vcf.py --vcf huge.vcf --output out.json

Error: VCF contains 1,234 variants after filtering (limit: 100).

Suggestions:
  - Apply more restrictive filters: --min-qual, --min-dp, --pass-only
  - Specify a genomic region: --region chr1:1000-2000
  - Override limit with --force (warning: may produce very large JSON)
  - Use bcftools directly for large-scale processing

Current filter conditions:
  (no filters applied)

# ケース 2: --force で制限を無視
$ python scripts/inspect_vcf.py --vcf huge.vcf --force --output out.json
⚠ Warning: Exporting 1,234 variants (limit bypassed with --force)
✓ Successfully exported 1,234 variants to out.json (15.2 MB)

# ケース 3: 成功（45 エントリ）
$ python scripts/inspect_vcf.py --vcf input.vcf --min-qual 30 --output out.json
✓ Successfully exported 45 variants to out.json
```

**実装の流れ**:
```python
1. VCF を開く（pysam.VariantFile）
2. フィルタ条件を適用しながらバリアントを収集
3. エントリ数をカウント
4. 制限チェック（--force がない場合）
5. JSON 形式に変換
6. 出力（ファイルまたは標準出力）
```

## SKILL.md 構成案

```markdown
---
name: vcf-toolkit
description: "VCF/BCF variant file operations using bcftools. Filter variants by quality/depth/frequency, convert small VCFs to JSON, and calculate statistics. Use when working with variant call results from WGS/WES pipelines."
---

# VCF Toolkit

VCF/BCF ファイル操作ツールキット。bcftools を活用した高速フィルタリング、JSON 変換、統計情報取得。

## Quick Start

### Install
```bash
uv pip install pysam typer
```

### Filter variants
```bash
# 高品質バリアントのみ
python scripts/filter_vcf.py \
  --vcf input.vcf \
  --min-qual 30 \
  --min-dp 10 \
  --pass-only \
  --output filtered.vcf
```

### Convert to JSON (small VCF only)
```bash
python scripts/read_vcf.py --vcf filtered.vcf --output result.json
```

## Scripts

### filter_vcf.py - VCF Filtering

bcftools を使った VCF フィルタリング。品質、深度、アレル頻度などで条件指定。

**Options:**
- `--min-qual`: 最小品質スコア
- `--min-dp`: 最小深度
- `--min-af/--max-af`: アレル頻度範囲
- `--pass-only`: FILTER=PASS のみ
- `--region`: 領域指定

[詳細な使用例]

### read_vcf.py - JSON Conversion

小規模 VCF（≤100 エントリ）を JSON に変換。

**制限事項:**
- 100 エントリを超える場合はエラー停止
- 大規模データは filter_vcf.py で絞り込んでから使用

[詳細な使用例]

## Best Practices

1. **大規模 VCF は先にフィルタリング** - read_vcf.py の前に filter_vcf.py を使う
2. **bcftools を直接使う選択肢も** - 複雑な条件は bcftools コマンドが柔軟
3. **インデックスを作成** - `bcftools index input.vcf.gz` で高速化

## bcftools との使い分け

| 用途 | vcf-toolkit | bcftools 直接 |
|------|-------------|---------------|
| 基本フィルタリング | filter_vcf.py | どちらでも可 |
| 複雑な条件 | - | bcftools |
| JSON 出力 | read_vcf.py | - |
| 大規模処理 | - | bcftools |
```

## References 構成案

### Option A: 最小限の references

```
references/
└── bcftools_filters.md   # bcftools フィルタ式のリファレンス（150行以内）
```

内容：
- INFO/FORMAT フィールドの指定方法
- よく使うフィルタパターン
- 複雑な条件の例

### Option B: References なし

- SKILL.md に統合
- bcftools 公式ドキュメントへのリンク

## 既存スキルとの棲み分け

| スキル | 責任範囲 |
|--------|----------|
| **vcf-toolkit**（新規） | VCF/BCF 操作、フィルタリング、JSON 変換 |
| **pysam**（簡素化） | BAM/SAM/CRAM 操作、カバレッジ計算 |
| **sequence-io** | FASTA/FASTQ 操作 |
| **blast-search** | BLAST 検索 |
| **blat-api-searching** | BLAT アラインメント |

## Implementation Plan

1. ✅ 設計書作成（このファイル）
2. ✅ Scripts 実装
   - ✅ inspect_vcf.py（filter_vcf.py と read_vcf.py を統合）
   - ❌ vcf_stats.py（スコープ外、必要に応じて将来実装）
3. ✅ SKILL.md 作成
4. ❌ References 作成（不要、SKILL.md で十分）
5. ✅ テスト（実際の VCF ファイルで）
   - ✅ Test 1: chr1 PASS only → 9 variants
   - ✅ Test 2: chr1 all filters → 10 variants
   - ✅ Test 3: chr1 min-qual 90 → 7 variants
   - ✅ Test 4: region chr1:10000-14000 → 6 variants
6. ⏳ 本番環境へ移行（コミット済み、パッケージング待ち）

## Questions / Decisions Needed

1. ✅ vcf-toolkit という名前で良いか？
   - **決定**: vcf-toolkit で確定
2. ✅ inspect_vcf.py の引数は十分か？
   - **決定**: --chrom or --region 必須、--pass-only デフォルト有効で確定
3. ✅ References は必要か？
   - **決定**: 不要、SKILL.md で十分
4. ✅ vcf_stats.py は実装するか？
   - **決定**: スコープ外、必要に応じて将来実装
5. ✅ 出力形式は JSON か TSV か？
   - **決定**: JSON 形式で確定（階層構造の表現が必須）
   - タブ区切りはトークン効率は良いが、構造化データの利点を優先

---

**Note**: 設計完了、実装・テスト済み、本番環境移行待ち
