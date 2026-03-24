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

補足:

- 生成結果は分割CSVのみです（`stars_base.csv`, `stars_extra7/8/9/10.csv`, `stars_extra_faint.csv`）。
- 旧形式の統合 `src/zstarview/data/stars.csv` は現行ランタイムでは使用しません。
- リポジトリへ同梱する現行データは `--max-vmag 10.5` を前提とし、`stars_extra_faint.csv` は `10 < vmag <= 10.5` の範囲を保持します。

既定では以下を出力します。

- `src/zstarview/data/stars/stars_base.csv`（`<= 6`）
- `src/zstarview/data/stars/stars_extra7.csv`（`6 < ... <= 7`）
- `src/zstarview/data/stars/stars_extra8.csv`（`7 < ... <= 8`）
- `src/zstarview/data/stars/stars_extra9.csv`（`8 < ... <= 9`）
- `src/zstarview/data/stars/stars_extra10.csv`（`9 < ... <= 10`）
- `src/zstarview/data/stars/stars_extra_faint.csv`（`10 < ... <= 10.5` for the committed dataset; generator output remains `> 10` when `--max-vmag` is larger）

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

- `--max-vmag`: 生成上限等級（既定: `9.0`。同梱データ更新時の推奨値: `10.5`）
- `--tycho-csv`: Tycho-2 正規化CSV入力
- `--tycho-i259-dir`: Tycho-2 I/259 ディレクトリ入力（`tyc2.dat.*.gz`）
- `--hip-priority-vmag`: この等級以下は Hipparcos を優先（既定: `6.0`）
- `--dup-sep-arcsec`: 重複判定の角距離しきい値（既定: `5.0`）
- `--dup-mag-diff`: 重複判定の等級差しきい値（既定: `0.75`）
- `--output-dir`: 分割出力先ディレクトリ

## 5. 検証観点

- 生成ログで件数・欠損率を確認する
- `vmag <= 6` の件数/見え方が大きく崩れていないことを確認する
- `vmag <= 9` および `<= 10` で目的に応じた密度になっていることを確認する
- 同梱データ更新時は `10 < vmag <= 10.5` の星が `stars_extra_faint.csv` に出力されることを確認する

## 6. 互換運用メモ

- 古い開発環境で `src/zstarview/data/stars.csv` が残っていても、現行コードは `src/zstarview/data/stars/stars_base.csv` を起点に分割読み込みします。
- データ更新時は統合CSVを作らず、分割CSV一式を更新してください。
