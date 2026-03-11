# 都市スカイライン生成手順メモ

最終更新: 2026-03-11

## 1. この文書の位置づけ

この文書は、`zstarview` 用の都市スカイライン JSON を作成するための補助手順を記録する appendix である。  
利用者向けの機能仕様は `docs/specification.md`、内部設計は `docs/design.md`、時系列の変更履歴は `docs/implementation-history.md` を参照する。

ここでは主に次を扱う。

- PLATEAU CityGML の取得
- `raw-data/` への配置
- CityGML から建物 GeoJSON を作る手順
- `src/zstarview/data/urban_skyline_demo.py` の実行方法
- `src/zstarview/data/urban_skyline_from_citygml.py` の実行方法
- Tokyo Skytree で実際に通した最小例と東京23区全体実行例

## 2. 前提

現時点の都市スカイライン生成は、次の前提で動く。

- 観測地点は `src/zstarview/data/viewpoints/tower_viewpoints.json` の bundled tower viewpoint を使う
- 観測者高さは tower viewpoint の `height_m` をそのまま使う
- 建物入力は GeoJSON FeatureCollection とする
- 各 Feature は Polygon または MultiPolygon とし、建物高さを property に持つ
- 出力は `id -> { name, profile }` 形式の JSON と確認用 PNG とする

## 3. ダウンロード

### 3.1 入手元

PLATEAU の対象都市データは公式ポータルから取得する。

- https://www.mlit.go.jp/plateau/open-data/
- https://front.geospatial.jp/plateau_portal_site/

### 3.2 取得対象

都市スカイライン生成で必要なのは、各都市の CityGML package に含まれる建物レイヤである。
具体的には展開後の次のディレクトリを使う。

```text
.../udx/bldg/*.gml
```

### 3.3 このリポジトリで確認済みの ZIP

`raw-data/` に次の ZIP を配置済みで確認した。

- `raw-data/13100_tokyo23-ku_2020_citygml_4_2_op.zip`
- `raw-data/14100_yokohama-shi_city_2024_citygml_2_op.zip`
- `raw-data/26100_kyoto-shi_city_2023_citygml_1_op.zip`
- `raw-data/27100_osaka-shi_city_2024_citygml_1_op.zip`

## 4. 展開

ZIP は `raw-data/` 直下または任意の作業ディレクトリへ展開してよい。  
今回の確認では `raw-data/` 直下へ展開した。

例:

```bash
unzip raw-data/13100_tokyo23-ku_2020_citygml_4_2_op.zip -d raw-data
```

展開後の建物 GML は次のように配置される。

```text
raw-data/13100_tokyo23-ku_2020_citygml_4_2_op/udx/bldg/*.gml
```

## 5. 都市スカイライン生成スクリプト

生成スクリプト本体は次に置いてある。

```text
src/zstarview/data/urban_skyline_demo.py
```

このスクリプトは次を行う。

- bundled tower viewpoint を読み込む
- GeoJSON の建物 footprint と高さを読む
- 方位ごとに最大仰角を計算する
- 確認用 PNG を出す
- 集約 JSON を `id -> { name, profile }` 形式で出す

また、PLATEAU CityGML の `udx/bldg` ディレクトリを直接読んでタイルごとに skyline を作り、最後に `max` 合成する CLI も追加した。

```text
src/zstarview/data/urban_skyline_from_citygml.py
```

このスクリプトは次を行う。

- `udx/bldg/*.gml` を走査する
- 各 GML の `gml:Envelope` を読んで観測地点周辺タイルを選ぶ
- 各タイルから `lod0RoofEdge` と `measuredHeight` を読む
- タイルごとに複数半径の cumulative skyline を計算する
- 方位ごとに `max(tile_altitude)` を取って全体 skyline を合成する
- 確認用 PNG と集約 JSON を出す

既定の出力先は次である。

```text
src/zstarview/data/viewpoints/urban_skyline/
```

## 6. Tokyo Skytree で実際に通した最小手順

### 6.1 Skytree を含む建物タイルの特定

Tokyo Skytree の bundled viewpoint 座標は次である。

- latitude: `35.710055555`
- longitude: `139.810722222`

東京23区の CityGML から、この座標を含む `bldg` タイルを調べたところ、次の GML が該当した。

```text
raw-data/13100_tokyo23-ku_2020_citygml_4_2_op/udx/bldg/53394654_bldg_6697_2_op.gml
```

### 6.2 CityGML から最小 GeoJSON を作る

現時点ではリポジトリ内に汎用の CityGML -> GeoJSON 変換スクリプトはまだ入れていない。  
そのため、Tokyo Skytree の検証では `lod0RoofEdge` と `measuredHeight` だけを抜く最小スクリプトを 1 回実行した。

```bash
.venv/bin/python - <<'PY'
from pathlib import Path
import json
import xml.etree.ElementTree as ET

src = Path('raw-data/13100_tokyo23-ku_2020_citygml_4_2_op/udx/bldg/53394654_bldg_6697_2_op.gml')
dst = Path('raw-data/derived/tokyo_skytree_53394654_buildings.geojson')
dst.parent.mkdir(parents=True, exist_ok=True)

ns = {
    'bldg': 'http://www.opengis.net/citygml/building/2.0',
    'gml': 'http://www.opengis.net/gml',
}

root = ET.parse(src).getroot()
features = []
for index, building in enumerate(root.findall('.//bldg:Building', ns)):
    height_text = building.findtext('bldg:measuredHeight', default='', namespaces=ns)
    try:
        height = float(height_text)
    except ValueError:
        continue

    pos = building.find('.//bldg:lod0RoofEdge//gml:posList', ns)
    if pos is None or not pos.text:
        continue

    values = [float(v) for v in pos.text.split()]
    if len(values) < 12:
        continue

    coords = []
    for i in range(0, len(values), 3):
        lat, lon = values[i], values[i + 1]
        point = [lon, lat]
        if not coords or coords[-1] != point:
            coords.append(point)

    if len(coords) < 4:
        continue
    if coords[0] != coords[-1]:
        coords.append(coords[0])

    features.append({
        'type': 'Feature',
        'properties': {
            'id': building.attrib.get('{http://www.opengis.net/gml}id', f'bldg-{index}'),
            'measuredHeight': height,
        },
        'geometry': {
            'type': 'Polygon',
            'coordinates': [coords],
        },
    })

dst.write_text(
    json.dumps({'type': 'FeatureCollection', 'features': features}, ensure_ascii=False),
    encoding='utf-8',
)
print(dst)
print('features', len(features))
PY
```

この実行では次の GeoJSON が生成された。

```text
raw-data/derived/tokyo_skytree_53394654_buildings.geojson
```

生成件数は `2864` building features だった。

### 6.3 都市スカイライン生成

生成した GeoJSON を `urban_skyline_demo.py` に渡して実行した。

```bash
.venv/bin/python src/zstarview/data/urban_skyline_demo.py \
  --buildings raw-data/derived/tokyo_skytree_53394654_buildings.geojson \
  --tower "Tokyo Skytree" \
  --output-dir src/zstarview/data/viewpoints/urban_skyline \
  --write-json
```

この実行で次が生成された。

- `src/zstarview/data/viewpoints/urban_skyline/tokyo_skytree_urban.png`
- `src/zstarview/data/viewpoints/urban_skyline/urban_skyline_profiles.json`

コンソール出力は次だった。

```text
[ok] Tokyo Skytree: src/zstarview/data/viewpoints/urban_skyline/tokyo_skytree_urban.png  buildings=2864/659  peak=-27.19deg@306.5
[ok] skyline-json: src/zstarview/data/viewpoints/urban_skyline/urban_skyline_profiles.json
```

## 7. CityGML ディレクトリ直読での実行

### 7.1 方針

GeoJSON へ一旦まとめて変換する方法でも動くが、東京23区全体のように対象範囲が広い場合は中間 GeoJSON が大きくなる。  
そのため、現在は `udx/bldg` ディレクトリを直接受ける CLI を追加し、各タイルの skyline を後で `max` 合成する方式も使えるようにした。

この方式では、観測地点の緯度経度から各タイルの `gml:Envelope` までの距離を見て候補タイルを選び、各タイル内の建物だけで skyline を作る。  
最終結果は各方位について次で合成する。

```python
global_alt = max(tile_alt_1, tile_alt_2, ...)
```

### 7.2 CLI

東京23区の `udx/bldg` ディレクトリを直接指定して Tokyo Tower を計算するコマンドは次である。

```bash
.venv/bin/python src/zstarview/data/urban_skyline_from_citygml.py \
  --citygml-dir raw-data/13100_tokyo23-ku_2020_citygml_4_2_op/udx/bldg \
  --tower "Tokyo Tower" \
  --radius-km 50 \
  --workers 8 \
  --output-dir src/zstarview/data/viewpoints/urban_skyline \
  --write-json
```

`--radius-km 50` は、Tokyo Tower 観測地点から見て東京23区の建物タイルを実質的に全て含める意図で使っている。  
このスクリプトは「全タイルを無条件で読む」のではなく、「観測地点から半径以内にかかるタイルを選ぶ」方式である。
`--workers` を 2 以上にすると、選ばれた GML タイルをプロセス並列で処理する。

現在の既定 skyline 半径は次である。

- `2 km`
- `3 km`
- `4.5 km`
- `6.75 km`
- `10.125 km`
- `15.1875 km`
- `22.78125 km`

また、既定の方位刻みは `0.1°` である。  
初期の `0.5°` より細かくしたのは、遠距離の細い構造物が bin に乗らず見えづらくなるリスクを下げるためである。

各半径は「その距離までの累積最大値」ではなく、その半径から外側へ `90 m` 幅の帯でスキャンした最大遮蔽仰角として扱う。  
つまり `2 km` レイヤは `2000 m - 2090 m` の帯、`3 km` レイヤは `3000 m - 3090 m` の帯を見ている。

選ばれたタイル名を出したい場合は次を使う。

```bash
.venv/bin/python src/zstarview/data/urban_skyline_from_citygml.py \
  --citygml-dir raw-data/13100_tokyo23-ku_2020_citygml_4_2_op/udx/bldg \
  --tower "Tokyo Tower" \
  --radius-km 50 \
  --workers 8 \
  --output-dir src/zstarview/data/viewpoints/urban_skyline \
  --write-json \
  --print-selected-tiles
```

### 7.3 東京23区全体での実行結果

上記コマンドを実行した結果は次だった。

```text
[tile] 53395731_bldg_6697_2_op.gml  buildings=88/8  peak=-3.33deg@42.0
[ok] Tokyo Skytree: src/zstarview/data/viewpoints/urban_skyline/tokyo_skytree_urban.png  tiles=672  buildings=1768062/48675  peak=-1.40deg@235.5
[ok] skyline-json: src/zstarview/data/viewpoints/urban_skyline/urban_skyline_profiles.json
```

重要な値は次である。

- 対象タイル数: `672`
- 建物数: `1,768,062`
- skyline に寄与した建物数: `48,675`
- 最大遮蔽仰角: `-1.40°`
- 最大遮蔽方位: `235.5°`

最大遮蔽仰角が負なのは、Tokyo Skytree の bundled viewpoint では観測者高さを `634 m` としているためである。  
東京タワーのような遠方高層建物は見えるが、Skytree 観測点より十分低いため、仰角としては負になる。

## 8. 出力 JSON 形式

現在の集約 JSON は `tower id` をキーにし、値に `name` と複数半径ぶんの `profiles` を持つ。

```json
{
  "wikidata:Q57965": {
    "name": "Tokyo Skytree",
    "profiles": [
      {
        "radius_km": 0.1,
        "profile": [
          {"az": 0.0, "alt": -12.3},
          {"az": 0.1, "alt": -12.6}
        ]
      },
      {
        "radius_km": 20.0,
        "profile": [
          {"az": 0.0, "alt": -1.9},
          {"az": 0.1, "alt": -2.1}
        ]
      }
    ]
  }
}
```

各 `profile` は、その `radius_km` から外側 `90 m` 幅の帯で取った方位ごとの最大遮蔽仰角である。  
アプリ側ではこれらを近距離から遠距離まで複数本の白線として重ね描きする。

現在は近距離帯を `2 km` 以上に絞っている。  
高所観測地点ではそれより近い帯が画面下端付近に張り付きやすく、スカイライン確認用の情報量よりノイズが勝ちやすいためである。

## 9. 注意点

- 現在の最小変換スクリプトは `lod0RoofEdge` と `measuredHeight` のみを使う簡易版である
- 建物 ground elevation はまだ使っていない
- そのため現時点の結果は「都市スカイライン生成の試作確認」用であり、最終版ではない
- Tokyo Skytree のように観測者高度が高い地点では、周辺建物は大半が負の仰角になる
- 全国タワーに対しては、都市ごとに対応する CityGML または GeoJSON を用意して個別に回す必要がある
- `urban_skyline_from_citygml.py` は半径内のタイルを選ぶ方式であり、都道府県や市区町村全域を必ずしも無条件には読まない

## 10. 今後の整理候補

- CityGML から建物 GeoJSON を作る正式な前処理スクリプトを `src/zstarview/data/` に追加する
- 複数都市をまとめて処理する driver script を追加する
- `urban_skyline_profiles.json` をアプリ本体が直接読み込めるようにする
- terrain horizon と urban skyline の合成ロジックを本体側へ追加する
