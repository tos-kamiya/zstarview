# Developer Docs

このディレクトリは、開発者向けのセットアップ手順・運用メモを置く場所です。

## Contents

- `docs/developer/uv-workflow.md`
  - `uv` 環境でのセットアップ、実行、テスト、型チェック、ビルド手順
- `docs/developer/cloud-snapshot-script.md`
  - 都市指定でリアルタイム雲画像（PNG）を保存する開発用スクリプトの使い方
- `docs/developer/star-catalog-generation.md`
  - Hipparcos + Tycho-2 を使った星カタログ再生成手順
- `docs/developer/dso-catalog-generation.md`
  - OpenNGC (`pyongc`) を使った DSO カタログ生成手順

## Notes

- 新しい設計判断や運用フローを追加したら、ここに短いノートを増やしてください。
- 仕様/設計の本体は `docs/specification.md`, `docs/design.md` を更新してください。
- 開発環境セットアップ時は `docs/developer/uv-workflow.md` の `uv pip install -e ".[dev]"` を使用してください。
