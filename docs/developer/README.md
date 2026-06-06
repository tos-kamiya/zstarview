# Developer Docs

このディレクトリは、開発者向けのセットアップ手順・運用メモを置く場所です。

## Contents

- `docs/developer/uv-workflow.md`
  - `uv` 環境でのセットアップ、実行、テスト、型チェック、ビルド手順
- `docs/developer/python-version-test-matrix.md`
  - CPython 3.10 / 3.11 / 3.12 / 3.13 / 3.14 で `pytest` / `mypy` を確認する手順
  - 自動実行スクリプト `scripts/run_python_matrix.sh` の使い方
- `docs/developer/cloud-snapshot-script.md`
  - 都市指定でリアルタイム雲画像（PNG）を保存する開発用スクリプトの使い方
- `docs/developer/star-catalog-generation.md`
  - Hipparcos + Tycho-2 を使った星カタログ再生成手順
- `docs/developer/dso-catalog-generation.md`
  - OpenNGC (`pyongc`) を使った DSO カタログ生成手順
- `docs/developer/viewpoint-dataset-generation.md`
  - tower / mountain viewpoint dataset の生成フローと運用方針
- `docs/developer/urban-outline-shape-audit-script.md`
  - derived building tiles に含まれる小さな footprint の分布を調べる監査スクリプト
  - 1m ビンの CSV ヒストグラム出力を含む
- `dev-samples/water_tile_overview.py`
  - sea mask tile の解像度別グリッド、`0`/`1`/`.tif` の分布、近傍プローブ、バンド別の実際の open tile 数を表示する確認用スクリプト
  - `--observer-ground-m 8651` のように高地の地表高さを渡して、観測地点まわりでどの tile がどう見えているかをざっくり確認するために使う
- `dev-samples/render_water_footprints_svg.py`
  - water overlay の cache JSON や footprint JSON を読み込み、各フットプリントを lon/lat 平面上の SVG として可視化するスクリプト
  - 簡約化前後の cache を見比べて、リング形状や削減結果を確認する用途に使う
  - `--show-labels` や `--ele-only` で表示を絞り込める
  - 実行例: `uv run -p .venv/bin/python dev-samples/render_water_footprints_svg.py --input ~/.cache/zstarview/water_overlay/<file>.json`

## Notes

- 新しい設計判断や運用フローを追加したら、ここに短いノートを増やしてください。
- 仕様/設計の本体は `docs/specification.md`, `docs/design.md`, `docs/design/*.md` を更新してください。
- 開発環境セットアップ時は `docs/developer/uv-workflow.md` の `uv pip install -e ".[dev]"` を使用してください。
