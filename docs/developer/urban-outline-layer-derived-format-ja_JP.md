# 都市アウトラインレイヤ中間データ形式

最終更新: 2026-03-12

## 1. この文書の位置づけ

この文書は、任意地点向けの都市アウトラインレイヤを高速に生成するための、中間データ形式を定義する。  
ここでいう中間データは、PLATEAU CityGML の `udx/bldg/*.gml` をそのまま実行時に読む代わりに、必要最小限の建物情報だけを保持した軽量タイルデータである。

この文書で扱うのは次である。

- なぜ中間データが必要か
- タイル単位で保持する理由
- 推奨 JSON フォーマット
- 保持すべき項目と削除すべき項目
- アプリ同梱対象都市

## 2. 背景

現在の CityGML 直読方式は、実行時に次を毎回行う。

- 多数の `*.gml` を開く
- XML をパースする
- `measuredHeight` / `storeysAboveGround` を読む
- `lod0RoofEdge` を取り出す

このコストが大きく、任意地点でのオンデマンド生成には不向きである。  
そのため、実行時に不要な情報を事前に落とした軽量中間データを作る。

## 3. 基本方針

- 元の `udx/bldg/*.gml` ごとに 1 ファイルを作る
- タイルの境界構造は PLATEAU と同じまま維持する
- 建物は `40m` 以上だけ残す
- 高さは `measuredHeight` を優先し、欠損時のみ `storeysAboveGround * 3.5m` を使う
- 実行時に必要な情報だけ持つ

## 4. 形式

### 4.1 推奨形式

- UTF-8 JSON
- ファイル拡張子は `*.json` または `*.json.gz`
- まずは可読性重視で JSON を採用する

`json.gz` にするかどうかは、同梱サイズとロード時間のバランスで後で決めてよい。  
まずは非圧縮 JSON で仕様を固める。

### 4.2 1 タイル 1 ファイル

PLATEAU の `bldg/*.gml` 1 枚に対して、中間 JSON も 1 枚作る。

例:

```text
raw-data/13100_tokyo23-ku_2020_citygml_4_2_op/udx/bldg/53393599_bldg_6697_2_op.gml
```

に対して

```text
src/zstarview/data/plateau_derived/13100_tokyo23/bldg/53393599_bldg_6697_2_op.json
```

を作る。

## 5. JSON スキーマ案

```json
{
  "schema_version": 1,
  "source": {
    "format": "PLATEAU-CityGML",
    "city_code": "13100",
    "path": "13100_tokyo23-ku_2020_citygml_4_2_op/udx/bldg/53393599_bldg_6697_2_op.gml"
  },
  "tile": {
    "id": "53393599_bldg_6697_2_op",
    "bbox": {
      "min_lat": 35.6500,
      "min_lon": 139.7400,
      "max_lat": 35.6600,
      "max_lon": 139.7500
    }
  },
  "filters": {
    "min_height_m": 40.0,
    "storey_height_m": 3.5
  },
  "buildings": [
    {
      "id": "bldg_xxx",
      "height_m": 82.5,
      "height_source": "measuredHeight",
      "bbox": {
        "min_lat": 35.6550,
        "min_lon": 139.7440,
        "max_lat": 35.6558,
        "max_lon": 139.7449
      },
      "rings": [
        [
          [139.7440, 35.6550],
          [139.7449, 35.6550],
          [139.7449, 35.6558],
          [139.7440, 35.6558],
          [139.7440, 35.6550]
        ]
      ]
    }
  ]
}
```

## 6. 各項目の意味

### 6.1 ルート

- `schema_version`
  - 中間データ仕様のバージョン
- `source`
  - 元の PLATEAU CityGML を特定する情報
- `tile`
  - タイル自体の識別子と bbox
- `filters`
  - 生成時に適用した高さ閾値など
- `buildings`
  - このタイル内で残した建物一覧

### 6.2 建物

- `id`
  - 元 GML 由来の建物 id
- `height_m`
  - skyline / urban outline layer 用に確定した建物高さ
- `height_source`
  - `measuredHeight`
  - `storeysAboveGround*3.5`
- `bbox`
  - 建物 footprint の bbox
- `rings`
  - 外周 ring 群
  - 各座標は `[lon, lat]`
  - 閉路にして持つ

## 7. 保持する項目

最低限必要なのは次だけである。

- タイル bbox
- 建物 id
- 建物 bbox
- `height_m`
- `height_source`
- 外周 `rings`

## 8. 削除する項目

実行時に不要なので、中間データには入れない。

- CityGML XML 全体
- 材質、用途、属性群
- 建物名
- ground elevation
- 3D 面情報
- 内周 ring
- `40m` 未満の建物

## 9. なぜタイル単位か

タイル単位にしておくと、実行時は次の流れにできる。

1. 観測地点から近いタイル bbox だけ読む
2. そのタイル内の建物 bbox でさらに粗く絞る
3. 必要な建物だけ投影する

これにより、全都市全建物を毎回ロードしなくて済む。

## 10. 同梱対象都市

アプリ同梱の derived dataset は、当面次の 3 都市に限定する。

- 東京23区
- 京都市
- 大阪市

理由は次の通り。

- viewpoint からの見え方に一定の効果がある
- 実際に使う可能性が高い
- 同梱サイズを抑えたい

現時点では、横浜市や川崎市の derived dataset は同梱対象にしない。  
必要なら後から追加する。

## 11. 想定ディレクトリ

同梱データの配置先案は次である。

```text
src/zstarview/data/plateau_derived/
  13100_tokyo23/
    bldg/*.json
  26100_kyoto/
    bldg/*.json
  27100_osaka/
    bldg/*.json
```

`README` やコード上では、これを「bundled derived PLATEAU tiles」として扱う。

## 12. 実行時キャッシュとの関係

この中間データは「実行時キャッシュ」ではなく「事前生成済み入力」である。  
任意地点向け urban outline layer の本番設計では、さらにその上で

- 観測地点
- 観測者高さ
- 半径
- 最低建物高さ
- 生成ロジックバージョン

をキーにした描画用キャッシュを持ってよい。

## 13. tile_index.json

各 `bldg/` ディレクトリには、個々のタイル JSON を読む前の索引として
`tile_index.json` を置く。

主な役割は次である。

- 都市全体 bbox で、関係ない都市ディレクトリを先に落とす
- tile bbox で、必要なタイル JSON だけを選ぶ
- index と実ファイルの不整合を検出しやすくする

### 13.1 フォーマット

```json
{
  "schema_version": 1,
  "city_code": "13100",
  "city_name": "tokyo23",
  "generated_at": "2026-03-12T12:34:56Z",
  "bbox": {
    "min_lat": 35.52,
    "min_lon": 139.56,
    "max_lat": 35.84,
    "max_lon": 139.93
  },
  "tile_count": 672,
  "min_height_m": 40.0,
  "tiles": [
    {
      "id": "53393599_bldg_6697_2_op",
      "path": "53393599_bldg_6697_2_op.json",
      "bbox": {
        "min_lat": 35.6500,
        "min_lon": 139.7400,
        "max_lat": 35.6600,
        "max_lon": 139.7500
      },
      "building_count": 42
    }
  ]
}
```

ルート `bbox` は都市データ全体の bbox、`tiles[].bbox` は各タイルの bbox である。

### 13.2 生成スクリプト

現在の index 生成スクリプトは次である。

```text
src/zstarview/data/build_plateau_tile_index.py
```

例:

```bash
.venv/bin/python src/zstarview/data/build_plateau_tile_index.py \
  --derived-dir src/zstarview/data/plateau_derived/13100_tokyo23/bldg
```

これにより、`src/zstarview/data/plateau_derived/13100_tokyo23/bldg/tile_index.json`
が生成される。

## 14. 今後の作業

次段階で必要なのは次である。

- bundled 3 都市分の derived dataset と tile index の作成
- ユーザー追加の derived dataset を同じ規約で検出する実装

## 15. 変換スクリプト

現在の変換スクリプトは次である。

```text
src/zstarview/data/build_plateau_derived_tiles.py
```

このスクリプトは、`udx/bldg/*.gml` を 1 タイルずつ軽量 JSON に変換する。  
高さは `measuredHeight` を優先し、欠損時のみ `storeysAboveGround * 3.5m` を使い、`40m` 未満は落とす。

例:

```bash
.venv/bin/python src/zstarview/data/build_plateau_derived_tiles.py \
  --citygml-dir raw-data/13100_tokyo23-ku_2020_citygml_4_2_op/udx/bldg \
  --workers 8 \
  --output-dir src/zstarview/data/plateau_derived/13100_tokyo23/bldg
```

`--workers` を `2` 以上にすると、タイル変換はマルチプロセスで実行し、JSON の書き込みだけを親プロセスで行う。

同梱対象 3 都市の生成コマンドは次である。

```bash
.venv/bin/python src/zstarview/data/build_plateau_derived_tiles.py \
  --citygml-dir raw-data/13100_tokyo23-ku_2020_citygml_4_2_op/udx/bldg \
  --workers 8 \
  --output-dir src/zstarview/data/plateau_derived/13100_tokyo23/bldg
```

```bash
.venv/bin/python src/zstarview/data/build_plateau_derived_tiles.py \
  --citygml-dir raw-data/26100_kyoto-shi_city_2023_citygml_1_op/udx/bldg \
  --workers 8 \
  --output-dir src/zstarview/data/plateau_derived/26100_kyoto/bldg
```

```bash
.venv/bin/python src/zstarview/data/build_plateau_derived_tiles.py \
  --citygml-dir raw-data/27100_osaka-shi_city_2024_citygml_1_op/udx/bldg \
  --workers 8 \
  --output-dir src/zstarview/data/plateau_derived/27100_osaka/bldg
```

各都市について、derived tile 生成後に tile index も作る。

```bash
.venv/bin/python src/zstarview/data/build_plateau_tile_index.py \
  --derived-dir src/zstarview/data/plateau_derived/13100_tokyo23/bldg
```

```bash
.venv/bin/python src/zstarview/data/build_plateau_tile_index.py \
  --derived-dir src/zstarview/data/plateau_derived/26100_kyoto/bldg
```

```bash
.venv/bin/python src/zstarview/data/build_plateau_tile_index.py \
  --derived-dir src/zstarview/data/plateau_derived/27100_osaka/bldg
```

## 16. urban outline layer の one-off 出力

この節のコマンドは検証用・共有用の one-off 出力である。  
必要なら、derived tile から viewpoint 単位の urban outline を JSON として書き出せる。  
ただしアプリ本体はこの集約 JSON を使わず、derived tile を直接読む。

```bash
.venv/bin/python src/zstarview/data/urban_debug_layer_from_citygml.py \
  --derived-dir src/zstarview/data/plateau_derived/13100_tokyo23/bldg \
  --all-covered-towers \
  --output-json src/zstarview/data/viewpoints/urban_debug_layer/tokyo23_urban_debug_layer.json
```

```bash
.venv/bin/python src/zstarview/data/urban_debug_layer_from_citygml.py \
  --derived-dir src/zstarview/data/plateau_derived/26100_kyoto/bldg \
  --all-covered-towers \
  --output-json src/zstarview/data/viewpoints/urban_debug_layer/kyoto_urban_outline_layer.json
```

```bash
.venv/bin/python src/zstarview/data/urban_debug_layer_from_citygml.py \
  --derived-dir src/zstarview/data/plateau_derived/27100_osaka/bldg \
  --all-covered-towers \
  --output-json src/zstarview/data/viewpoints/urban_debug_layer/osaka_urban_outline_layer.json
```

スクリプト名は過去互換のためそのままだが、入力は `--citygml-dir` ではなく `--derived-dir` である。

## 17. アプリ本体での利用

アプリ本体は urban outline 集約 JSON を読まない。  
`src/zstarview/data/plateau_derived/*/bldg` を走査し、現在地点の近傍で該当する derived tile が見つかる場合は、その場で derived tile を読んで urban outline layer を生成する。

タイル選択は `tile_index.json` がある場合はそれを優先し、まず都市全体 bbox で都市を落とし、次に `tiles[].bbox` で必要タイルだけを選ぶ。  
`tile_index.json` が無い古い derived dataset に対してだけ、`bldg/*.json` を直接走査する fallback を使う。

既定の同梱データは東京23区・京都市・大阪市を想定している。  
同じ形式の derived dataset と `tile_index.json` が追加されれば、同じ仕組みで扱える。
