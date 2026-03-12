# Developer Docs

このディレクトリは、開発者向けのセットアップ手順・運用メモを置く場所です。

## Contents

- `docs/developer/uv-workflow.md`
  - `uv` 環境でのセットアップ、実行、テスト、型チェック、ビルド手順
- `docs/developer/python-version-test-matrix.md`
  - CPython 3.10 / 3.11 / 3.12 / 3.13 で `pytest` / `mypy` を確認する手順
  - 自動実行スクリプト `scripts/run_python_matrix.sh` の使い方
- `docs/developer/cloud-snapshot-script.md`
  - 都市指定でリアルタイム雲画像（PNG）を保存する開発用スクリプトの使い方
- `docs/developer/star-catalog-generation.md`
  - Hipparcos + Tycho-2 を使った星カタログ再生成手順
- `docs/developer/dso-catalog-generation.md`
  - OpenNGC (`pyongc`) を使った DSO カタログ生成手順
- `docs/developer/viewpoint-dataset-generation.md`
  - tower / mountain viewpoint dataset の生成フローと運用方針
- `docs/developer/urban-skyline-generation-ja_JP.md`
  - PLATEAU CityGML から都市スカイラインを生成する手順
- `docs/developer/urban-debug-layer-derived-format-ja_JP.md`
  - PLATEAU 由来の都市デバッグレイヤ用中間タイルデータ形式と同梱方針

## Notes

- 新しい設計判断や運用フローを追加したら、ここに短いノートを増やしてください。
- 仕様/設計の本体は `docs/specification.md`, `docs/design.md` を更新してください。
- 開発環境セットアップ時は `docs/developer/uv-workflow.md` の `uv pip install -e ".[dev]"` を使用してください。
