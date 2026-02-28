# 星カタログ生成ガイド

最終更新: 2026-03-01

`zstarview` の星データを Hipparcos（必須）+ Tycho-2（任意）から再生成する手順です。

## 1. 生成スクリプト

- スクリプト:
  - `src/zstarview/data/stars/generate_star_catalog.py`
- 入力:
  - Hipparcos: `hip_main.dat` または `hip_main.dat.zip`
  - Tycho-2: 正規化済み CSV（任意）または I/259 ディレクトリ（`tyc2.dat.*.gz`）
  - IAU 名称 CSV（任意。HIP への名前付け）

## 2. Hipparcos のみで生成

```bash
uv run -p .venv/bin/python src/zstarview/data/stars/generate_star_catalog.py
```

既定では以下を出力します。

- `src/zstarview/data/stars.csv`（統合版, `vmag <= 9`）
- `src/zstarview/data/stars/stars_base.csv`（`<= 6`）
- `src/zstarview/data/stars/stars_extra7.csv`（`6 < ... <= 7`）
- `src/zstarview/data/stars/stars_extra8.csv`（`7 < ... <= 8`）
- `src/zstarview/data/stars/stars_extra9.csv`（`8 < ... <= 9`）

## 3. Tycho-2 を追加して生成

### 3.1 正規化CSVを使う場合

```bash
uv run -p .venv/bin/python src/zstarview/data/stars/generate_star_catalog.py \
  --tycho-csv /path/to/tycho2_normalized.csv
```

### 3.2 I/259生データを使う場合

`src/zstarview/data/I-259/` が存在すれば、既定で自動読込されます。
明示する場合は以下です。

```bash
uv run -p .venv/bin/python src/zstarview/data/stars/generate_star_catalog.py \
  --tycho-i259-dir src/zstarview/data/I-259
```

## 4. 主なパラメータ

- `--max-vmag`: 生成上限等級（既定: `9.0`）
- `--tycho-csv`: Tycho-2 正規化CSV入力
- `--tycho-i259-dir`: Tycho-2 I/259 ディレクトリ入力（`tyc2.dat.*.gz`）
- `--hip-priority-vmag`: この等級以下は Hipparcos を優先（既定: `6.0`）
- `--dup-sep-arcsec`: 重複判定の角距離しきい値（既定: `5.0`）
- `--dup-mag-diff`: 重複判定の等級差しきい値（既定: `0.75`）
- `--output`: 統合出力先 CSV
- `--output-dir`: 分割出力先ディレクトリ

## 5. 検証観点

- 生成ログで件数・欠損率を確認する
- `vmag <= 6` の件数/見え方が大きく崩れていないことを確認する
- `vmag <= 9` で天の川の連続性が改善していることを確認する
