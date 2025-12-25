# Local BLAT setup

ローカルで BLAT を実行する前提のメモです。

## Reference directory

- 2bit 参照は `~/.local/share/blat` に置く。
- `hg38` は UCSC の 2bit をダウンロードする。
- `CHM13` は UCSC の 2bit をダウンロードする。

```
https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit
https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.2bit
```

## Required tools

- `blat`
- `curl`

## Utility script

`scripts/run_blat_local.py` は以下を行う:

- `hg38`: UCSC から `hg38.2bit` をダウンロード。
- `CHM13`: UCSC から `hs1.2bit` をダウンロードし、`CHM13.2bit` として保存。

### Examples

```bash
python scripts/run_blat_local.py run --reference hg38 --fasta path/to/query.fasta
python scripts/run_blat_local.py run --reference CHM13 --sequence ACTG...
```

## Output

- 出力は PSL 形式。
- `--output` を指定しない場合は stdout に PSL を出力。
