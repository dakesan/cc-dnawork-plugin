---
name: vcf-toolkit
description: "VCF/BCF variant file toolkit. Filter and export variants as JSON from WGS/WES results. Use when you need to inspect variants from a specific chromosome or region, apply quality filters, and get structured JSON output for downstream analysis."
---

# VCF Toolkit

VCF/BCF ファイルから特定の chromosome や region のバリアントを抽出し、JSON 形式で出力するツールキット。WGS/WES 解析結果の確認に最適。

## Quick Start

### Install

```bash
uv pip install pysam typer
```

### Basic Usage

```bash
# chr1 の PASS バリアントを JSON 出力（≤100 エントリ）
python scripts/inspect_vcf.py --vcf variants.vcf --chrom chr1 --output chr1.json

# 特定領域を指定
python scripts/inspect_vcf.py --vcf variants.vcf --region chr1:1000000-2000000 --output region.json

# フィルタ条件を追加
python scripts/inspect_vcf.py \
  --vcf variants.vcf \
  --chrom chr1 \
  --min-qual 30 \
  --min-dp 10 \
  --output filtered.json
```

## Scripts

### inspect_vcf.py - VCF Inspection & JSON Export

VCF ファイルから chromosome または region を指定してバリアントを抽出し、JSON 形式で出力。

#### 必須引数

- `--vcf PATH` - 入力 VCF ファイル
- `--chrom TEXT` または `--region TEXT` - どちらか必須
  - `--chrom`: chromosome 全体（例: `chr1`）
  - `--region`: 特定領域（例: `chr1:1000000-2000000`）

#### オプション引数

**出力：**
- `--output PATH` - JSON 出力パス（未指定時は標準出力）

**フィルタ条件：**
- `--min-qual FLOAT` - 最小品質スコア（QUAL >= X）
- `--min-dp INT` - 最小深度（INFO/DP >= X）
- `--min-af FLOAT` - 最小アレル頻度（INFO/AF >= X）
- `--max-af FLOAT` - 最大アレル頻度（INFO/AF <= X）
- `--pass-only` / `--all-filters` - PASS のみ（デフォルト）/ 全フィルタ含む

**制限：**
- `--max-variants INT` - 最大バリアント数（デフォルト: 100）
- `--force` - エントリ数制限を無視（大量 JSON 出力を許可）

#### 出力形式（JSON）

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

## 使用例

### 例 1: 基本的な使い方

```bash
# chr1 の PASS バリアント（デフォルト ≤100 エントリ）
python scripts/inspect_vcf.py --vcf variants.vcf --chrom chr1 --output chr1.json
```

### 例 2: 高品質バリアントのみ

```bash
# QUAL >= 30 かつ DP >= 10 のバリアント
python scripts/inspect_vcf.py \
  --vcf variants.vcf \
  --chrom chr1 \
  --min-qual 30 \
  --min-dp 10 \
  --output high_quality.json
```

### 例 3: レアバリアント検索

```bash
# アレル頻度 <= 0.01 のレアバリアント
python scripts/inspect_vcf.py \
  --vcf variants.vcf \
  --chrom chr1 \
  --max-af 0.01 \
  --output rare_variants.json
```

### 例 4: 特定領域の詳細確認

```bash
# 興味のある領域（例: 遺伝子座）
python scripts/inspect_vcf.py \
  --vcf variants.vcf \
  --region chr17:41196312-41277500 \
  --output brca1_region.json
```

### 例 5: 全フィルタを含める

```bash
# PASS 以外のバリアントも含める
python scripts/inspect_vcf.py \
  --vcf variants.vcf \
  --chrom chr1 \
  --all-filters \
  --output all_variants.json
```

### 例 6: 制限を無視（大量出力）

```bash
# 100 エントリ超えても出力（注意：ファイルが巨大になる可能性）
python scripts/inspect_vcf.py \
  --vcf variants.vcf \
  --chrom chr1 \
  --force \
  --output chr1_all.json
```

## エラー処理

### エラー 1: エントリ数超過

```bash
$ python scripts/inspect_vcf.py --vcf huge.vcf --chrom chr1 --output out.json

Error: VCF contains 1,234+ variants after filtering (limit: 100).

Suggestions:
  - Apply more restrictive filters: --min-qual, --min-dp, --pass-only
  - Specify a genomic region: --region chr1:1000-2000
  - Override limit with --force (warning: may produce very large JSON)
  - Use bcftools directly for large-scale processing

Current filter conditions:
  --chrom chr1 --pass-only
```

**解決策：**
- より厳しいフィルタを適用（`--min-qual 30`, `--min-dp 10` など）
- 領域を狭める（`--region chr1:1000000-1100000`）
- `--force` で制限を無視

### エラー 2: chromosome/region 未指定

```bash
$ python scripts/inspect_vcf.py --vcf variants.vcf --output out.json

Error: Either --chrom or --region must be specified.
```

**解決策：**
- `--chrom chr1` または `--region chr1:1000-2000` を追加

## Best Practices

### 1. chromosome 指定は必須

全 VCF をそのまま JSON 化するのは非効率。必ず chromosome または region を指定する。

```bash
# ❌ Bad: 全 VCF は扱えない
python scripts/inspect_vcf.py --vcf variants.vcf

# ✅ Good: chromosome を指定
python scripts/inspect_vcf.py --vcf variants.vcf --chrom chr1
```

### 2. フィルタを活用

デフォルトで PASS のみだが、さらに品質・深度でフィルタすると効率的。

```bash
# ✅ Good: 高品質バリアントに絞る
python scripts/inspect_vcf.py \
  --vcf variants.vcf \
  --chrom chr1 \
  --min-qual 30 \
  --min-dp 10
```

### 3. 100 エントリ制限を意識

JSON 出力は小規模データ向け。大量のバリアントは bcftools で処理。

```bash
# 大規模データは bcftools で前処理
bcftools view -i 'QUAL>=30 && DP>=10' -r chr1:1000000-2000000 variants.vcf > filtered.vcf

# その後 JSON 化
python scripts/inspect_vcf.py --vcf filtered.vcf --chrom chr1 --output filtered.json
```

### 4. --force は慎重に

`--force` で制限を無視できるが、数千エントリの JSON は数 MB〜数十 MB になる可能性がある。

## bcftools との使い分け

| 用途 | vcf-toolkit | bcftools |
|------|-------------|----------|
| 小規模データの JSON 化 | ✅ inspect_vcf.py | - |
| 大規模フィルタリング | - | ✅ bcftools view |
| 複雑な条件 | - | ✅ bcftools |
| VCF → VCF 変換 | - | ✅ bcftools |
| 統計情報 | - | ✅ bcftools stats |

**推奨ワークフロー：**
1. bcftools で大規模フィルタリング
2. inspect_vcf.py で JSON 化して詳細確認
3. 下流解析（Python, R など）で JSON を活用

## 関連スキル

- **pysam** - BAM/CRAM アラインメントファイル操作
- **sequence-io** - FASTA/配列ファイル操作
- **blast-search** - BLAST 相同性検索
- **blat-api-searching** - BLAT ゲノムマッピング

## Troubleshooting

### Q: VCF が大きすぎてエラーが出る

A: より狭い領域を指定するか、bcftools で前処理してください。

```bash
# 領域を狭める
python scripts/inspect_vcf.py --vcf variants.vcf --region chr1:1000000-1100000

# または bcftools で前処理
bcftools view -i 'QUAL>=50' variants.vcf | python scripts/inspect_vcf.py --vcf - --chrom chr1
```

### Q: インデックスエラーが出る

A: VCF ファイルにインデックスが必要です。

```bash
# bgzip で圧縮
bgzip variants.vcf

# tabix でインデックス作成
tabix -p vcf variants.vcf.gz

# インデックス付き VCF を使用
python scripts/inspect_vcf.py --vcf variants.vcf.gz --chrom chr1
```

### Q: PASS 以外のバリアントも見たい

A: `--all-filters` フラグを使用してください。

```bash
python scripts/inspect_vcf.py --vcf variants.vcf --chrom chr1 --all-filters
```
