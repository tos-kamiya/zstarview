# 都市デバッグレイヤ中間データ形式

最終更新: 2026-03-12

## 1. この文書の位置づけ

この文書は、任意地点向けの都市デバッグレイヤを高速に生成するための、中間データ形式を定義する。  
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
  - skyline / debug layer 用に確定した建物高さ
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
任意地点向け debug layer の本番設計では、さらにその上で

- 観測地点
- 観測者高さ
- 半径
- 最低建物高さ
- 生成ロジックバージョン

をキーにした描画用キャッシュを持ってよい。

## 13. 今後の作業

次段階で必要なのは次である。

- CityGML -> derived JSON 変換スクリプトの追加
- タイル bbox index の生成
- bundled 3 都市分の derived dataset の作成
- 任意地点 debug layer 生成時に CityGML 直読ではなく derived dataset を使う実装
