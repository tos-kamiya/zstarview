# Cloud Snapshot Script

`scripts/save_cloud_snapshot.py` は、都市名または座標を指定して
「現在時刻の雲ディスク画像」を PNG として保存するための開発用スクリプトです。

## 用途

- GUI を起動せずに、雲画像取得・保存だけを単発で実行したいとき
- 地点を変えながらクラウドデータ取得結果を比較したいとき
- 失敗時の衛星選択や時刻情報をメタデータとして記録したいとき

## 実行例

```bash
uv run -p .venv/bin/python scripts/save_cloud_snapshot.py Tokyo
uv run -p .venv/bin/python scripts/save_cloud_snapshot.py JP/Tokyo -o tokyo_cloud.png
uv run -p .venv/bin/python scripts/save_cloud_snapshot.py "35.68;139.76" --meta tokyo_cloud.json
```

## 主な引数

- `location`:
  - 都市名（例: `Tokyo`）
  - `CC/都市名`（例: `JP/Tokyo`）
  - geonameid（数値）
  - 座標 `lat;lon`（例: `35.68;139.76`）
- `-o, --output`: 出力PNGパス（未指定時は自動命名）
- `--meta`: メタデータJSONの出力先
- `--alt`, `--az`: 視線方向（既定: `90`, `180`）
- `--radius-px`: 画像半径ピクセル（既定: `256`）
- `--sat-priority`: 優先衛星（例: `AUTO`, `HIMAWARI,G18,G16`）
- `--search-back-minutes`: 過去データ探索幅（分）

## 出力内容

- PNG画像（LAベースの雲ディスク）
- 標準出力に実行結果JSON
- `--meta` 指定時は同内容をJSONファイルにも保存

JSONには以下が含まれます。

- 解決した地点情報（lat/lon, geonameid）
- 使用衛星・プロダクト
- データ時刻（UTC）
- 元データファイルパス
- 出力PNGパス

## 補足

- 都市名解決は `cities1000.txt` を使い、複数候補時は人口上位を採用します。
- ネットワーク/S3取得に失敗した場合は、エラー終了（終了コード1または2）します。
