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

# 出力
--output PATH           JSON 出力パス（オプション、未指定時は標準出力）

# フィルタ条件
--min-qual FLOAT        最小品質スコア（QUAL >= X）
--min-dp INT            最小深度（INFO/DP >= X）
--min-af FLOAT          最小アレル頻度（INFO/AF >= X）
--max-af FLOAT          最大アレル頻度（INFO/AF <= X）
--pass-only             FILTER=PASS のみ（フラグ）
--region TEXT           領域指定（例: chr1:1000-2000）

# 制限
--max-variants INT      最大バリアント数（デフォルト: 100）
--force                 エントリ数制限を無視（フラグ）
```

**使用例**:
```bash
# 基本：フィルタなしで全エントリを JSON 出力
python scripts/inspect_vcf.py --vcf input.vcf --output all.json

# フィルタリング：高品質バリアントのみ
python scripts/inspect_vcf.py \
  --vcf input.vcf \
  --min-qual 30 \
  --min-dp 10 \
  --pass-only \
  --output filtered.json

# 領域指定
python scripts/inspect_vcf.py \
  --vcf input.vcf \
  --region chr1:1000000-2000000 \
  --output region.json

# 制限を無視（大量出力）
python scripts/inspect_vcf.py \
  --vcf huge.vcf \
  --force \
  --output large.json
```

**出力形式（JSON）**:
```json
{
  "num_variants": 45,
  "samples": ["sample1", "sample2"],
  "variants": [
    {
      "chrom": "chr1",
      "pos": 12345,
      "id": "rs123456",
      "ref": "A",
      "alts": ["G"],
      "qual": 100.0,
      "filter": ["PASS"],
      "info": {
        "DP": 50,
        "AF": [0.5],
        "AC": [25]
      },
      "samples": {
        "sample1": {"GT": "0/1", "DP": 25, "GQ": 99},
        "sample2": {"GT": "0/0", "DP": 25, "GQ": 99}
      }
    }
  ]
}
```

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
2. ⏳ Scripts 実装
   - [ ] filter_vcf.py（優先度: 高）
   - [ ] read_vcf.py（優先度: 高）
   - [ ] vcf_stats.py（優先度: 低、後回し可）
3. ⏳ SKILL.md 作成
4. ⏳ References 作成（必要に応じて）
5. ⏳ テスト（実際の VCF ファイルで）
6. ⏳ 本番環境へ移行

## Questions / Decisions Needed

1. ✅ vcf-toolkit という名前で良いか？（または vcf-tools, vcf-ops など）
2. ✅ filter_vcf.py の引数は十分か？
3. ⏳ References は必要か？
4. ⏳ vcf_stats.py は実装するか（後回しでも可）？

---

**Note**: この設計書はユーザーと協議しながら更新する
