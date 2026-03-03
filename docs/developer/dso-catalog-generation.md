# DSO カタログ生成ガイド

最終更新: 2026-03-04

`zstarview` の銀河・星団・星雲データ（DSO）を OpenNGC (`pyongc`) から生成する手順です。

## 1. 前提

- 開発依存をインストール済みであること:

```bash
uv pip install -e ".[dev]"
```

- `pyongc` は開発時のみ使用し、ランタイム依存には含めません。

## 2. 生成スクリプト

- スクリプト:
  - `src/zstarview/data/dso/generate_dso_catalog.py`
- 既定出力:
  - `src/zstarview/data/dso.csv`

## 3. 実行例

### 3.1 既定（Messier/NGC/IC）

```bash
uv run -p .venv/bin/python src/zstarview/data/dso/generate_dso_catalog.py
```

既定では `Type` が次の3種類のみ出力されます。

- `galaxy`
- `open_cluster`
- `globular_cluster`

また、既定で `Vmag <= 12.0` のみを出力します。

### 3.2 対象カタログと等級上限を指定

```bash
uv run -p .venv/bin/python src/zstarview/data/dso/generate_dso_catalog.py \
  --include messier,ngc,ic \
  --types all \
  --max-vmag 14 \
  --output src/zstarview/data/dso.csv
```

## 4. 出力カラム

- `Id`
- `Name`
- `Type`
- `RAh`
- `Dec`
- `Vmag`
- `MajorArcmin`
- `MinorArcmin`
- `PAdeg`
- `SourceCatalog`

## 5. 検証観点

- 生成件数がゼロでないこと
- `Id` が重複していないこと
- `RAh`, `Dec` が空でない行が大半であること
- `--include` / `--max-vmag` の変更で件数が意図通りに変化すること
- `--types all` を付けると、非銀河/非星団の種別も含めて出力されること
