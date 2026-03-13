# Viewpoint Dataset Generation

この文書は、`zstarview` の bundled viewpoint データセット生成フローを開発者向けにまとめたものである。  
対象は主に次の 2 系統である。

- `src/zstarview/data/viewpoints/tower_viewpoints.json`
- `src/zstarview/data/viewpoints/mountain_viewpoints.json`

都市アウトラインレイヤ用の PLATEAU 派生タイルについては、別文書 [urban-outline-layer-derived-format-ja_JP.md](/home/toshihiro/playground/zstarview/docs/developer/urban-outline-layer-derived-format-ja_JP.md) を参照する。

## 1. 位置づけ

アプリ本体は bundled JSON を読むだけで、viewpoint データを実行時に外部取得しない。  
そのため viewpoint 更新は「オフライン生成して JSON を同梱し直す」運用である。

関連する主な実行時モジュール:

- `src/zstarview/tower_viewpoints.py`
- `src/zstarview/mountain_viewpoints.py`
- `src/zstarview/startup.py`

関連する生成物:

- `src/zstarview/data/viewpoints/tower_viewpoints.json`
- `src/zstarview/data/viewpoints/mountain_viewpoints.json`

作業用の退避名として、マージ前のオリジナルを
`src/zstarview/data/viewpoints/tower_viewpoints.original.json`
に残す運用を取ることがある。

## 2. Tower Viewpoints

### 2.1 現状のデータソース

`tower_viewpoints.json` は Wikidata ベースの生成物である。  
現行 JSON ヘッダにも次が残っている。

- `source: wikidata`
- `source_query_result: dev-samples/wikidata_tower_query_raw_result_2026-03-08.json`

### 2.2 現状の選定ルール

現行 dataset は歴史的には「展望塔 dataset」として始まったが、現在は同じ `tower_viewpoints.json` に「高所観測地点」も混在させる方針である。  
主な正規化ルールは次の通り。

- `deduplicate_by: qid`
- `min_height_m: 100.0`
- `require_observation_tower: true`
- `observation_tower_qid: Q1440300`
- `height_rule: max`
- `coord_rule: first`

このため、旧クエリだけでは次のような展望フロア付き高層ビルは落ちやすかった。

- あべのハルカス
- 六本木ヒルズ
- 超高層複合ビルの展望台

現在はユーザー向けの dataset を分けず、同じ `tower_viewpoints.json` に追加する。  
その代わり schema 側で `viewpoint_height_m` を持てるようにし、建物全高 `height_m` と観測基準点の高さを分離できるようにする。

tower dataset の高さ関連キーは次の意味で使い分ける。

- `height_m`
  - 構造物そのものの高さ
- `viewpoint_height_m`
  - 地表からの観測基準点の高さ
- CLI `--observer-height-m`
  - 観測基準点から観測者の目線までの高さ

この bundled dataset は viewpoint 候補の一覧であって、各 viewpoint に対応する都市スカイライン前計算データが常に存在するとは限らない。  
例えば `Kobe Port Tower` は bundled viewpoint としては採用しているが、2026-03-12 時点では神戸市の PLATEAU Open Data が無いため、PLATEAU 由来 skyline JSON は未生成である。

### 2.3 生成スクリプト

テストから見る限り、旧 tower 側の生成スクリプトは `dev-samples/build_wikidata_tower_viewpoints.py` を前提としていた。  
現在は追加候補を bundled dataset に混ぜるための補助スクリプトとして `src/zstarview/data/build_wikidata_viewpoints.py` も使う。

関連テスト:

- `tests/test_wikidata_tower_viewpoints.py`

現在のリポジトリでは、この生成スクリプト自体は `dev-samples/` 配下の運用前提で、正式な `docs/developer/` 文書は未整備だった。

### 2.4 実際に使った WDQS クエリ

tower 側は、`https://query.wikidata.org/` で広めに候補を取り、後段のローカル正規化で重複除去とノイズ除外を行う方針だった。  
調査メモは次に残っている。

- `dev-samples/wikidata_tower_query_investigation_2026-03-08.md`
- `dev-samples/wikidata_tower_query_raw_result_2026-03-08.json`

実際に raw result を作るのに使ったクエリはこれである。

```sparql
SELECT ?item ?itemLabel ?class ?classLabel ?coord ?height WHERE {
  VALUES ?class {
    wd:Q1440300   # observation tower
    wd:Q11166728  # television tower
    wd:Q1068623   # transmitter mast
    wd:Q1798641   # communication tower
  }

  ?item wdt:P31 ?instance .
  ?instance wdt:P279* ?class .

  ?item wdt:P625 ?coord .
  ?item wdt:P2048 ?height .

  FILTER(?height >= 100)

  SERVICE wikibase:label { bd:serviceParam wikibase:language "ja,en". }
}
ORDER BY DESC(?height)
```

このクエリは 1 item あたり複数行を返しうる。  
理由は、1 item が複数 class に当たりうることと、height 値が複数ある場合があるためである。  
そのため、query 側で無理に 1 行化せず、`build_wikidata_tower_viewpoints.py` で次を行う設計になっていた。

- `qid` 単位で畳み込み
- `height_m` は最大値採用
- `classes` / `class_qids` は集合化
- `Q1440300` を持たないものは既定では除外
- fire tower / lookout tower 系は名前で除外

### 2.5 高所観測地点候補を足すためのクエリ

展望塔だけでは不足するため、追加候補の発見には別クエリも使う。  
実運用では 1 本の重いクエリではなく、用途別に分けた軽めのクエリを使う。

#### query-1: 展望塔 / 展望台候補

```sparql
SELECT ?item ?itemLabel ?coord ?height WHERE {
  VALUES ?class {
    wd:Q1440300   # observation tower
    wd:Q177305    # observation deck
  }
  ?item wdt:P31/wdt:P279* ?class .
  ?item wdt:P17 wd:Q17 .
  ?item wdt:P625 ?coord .
  OPTIONAL { ?item wdt:P2048 ?height . }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "ja,en". }
}
ORDER BY DESC(?height)
```

#### query-2: 展望用途付き高層建築候補

```sparql
SELECT ?item ?itemLabel ?coord ?height WHERE {
  VALUES ?class {
    wd:Q11303     # skyscraper
    wd:Q41176     # building
  }
  ?item wdt:P31/wdt:P279* ?class .
  ?item wdt:P17 wd:Q17 .
  ?item wdt:P625 ?coord .
  ?item wdt:P366 wd:Q177305 .   # has use = observation deck
  OPTIONAL { ?item wdt:P2048 ?height . }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "ja,en". }
}
ORDER BY DESC(?height)
```

#### query-3: ラベルなしの軽量版

```sparql
SELECT ?item ?coord ?height WHERE {
  VALUES ?class {
    wd:Q1440300
    wd:Q177305
  }
  ?item wdt:P31/wdt:P279* ?class .
  ?item wdt:P17 wd:Q17 .
  ?item wdt:P625 ?coord .
  OPTIONAL { ?item wdt:P2048 ?height . }
}
ORDER BY DESC(?height)
```

`query-2.json` は、あべのハルカスや横浜ランドマークタワーのような「タワーではない高所観測地点」を拾う用途に向いている。  
そのままでは重複や false positive があるため、`src/zstarview/data/build_wikidata_viewpoints.py` で bundled dataset にマージする前に整理する。

## 3. Mountain Viewpoints

### 3.1 現状のデータソース

`mountain_viewpoints.json` は「Wikipedia で候補選定し、Wikidata で正規化した dataset」という扱いである。

README / design / implementation history の記述を総合すると、流れは次の通り。

1. Wikipedia 起点で候補を集める
2. Wikidata metadata を引く
3. review 用 JSON を作る
4. curated seed を作る
5. 最終 `mountain_viewpoints.json` を生成する

### 3.2 関連スクリプト

テストから見えている mountain 側の生成スクリプトは次である。

- `dev-samples/build_wikidata_mountain_candidates.py`
- `dev-samples/build_curated_mountain_seed.py`
- `dev-samples/build_curated_mountain_viewpoints.py`

関連テスト:

- `tests/test_wikidata_mountain_candidates.py`
- `tests/test_curated_mountain_seed.py`
- `tests/test_wikidata_mountain_viewpoints.py`

### 3.3 現状のメタデータ

`mountain_viewpoints.json` には次のような由来情報が含まれる。

- `source_seed`
- `coordinate_rule`
- `elevation_rule`

実運用では、最終 JSON だけでなく途中 seed も保存しておく前提で考えるのがよい。

### 3.4 実際に使った WDQS クエリ

mountain 側は、最終 dataset を WDQS だけで自動生成するのではなく、候補発見とメタデータ補強に WDQS を使う方針だった。  
クエリのテンプレートと investigation 結果は次に残っている。

- `dev-samples/wikidata_mountain_query_investigation_TEMPLATE.md`
- `dev-samples/wikidata_mountain_query_investigation_result_2026-03-09.json`

テンプレート上、実際に使う前提だったクエリは少なくとも次の 3 系統である。

#### Broad mountain discovery

```sparql
SELECT ?item ?itemLabel ?coord ?elevation ?country ?countryLabel WHERE {
  ?item wdt:P31 ?instance .
  ?instance wdt:P279* wd:Q8502 .   # mountain

  ?item wdt:P625 ?coord .
  OPTIONAL { ?item wdt:P2044 ?elevation . }
  OPTIONAL { ?item wdt:P17 ?country . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "ja,en". }
}
ORDER BY DESC(?elevation) ?itemLabel
```

#### Highest point of each country

これは実際に `build_wikidata_mountain_candidates.py` へ流し込む raw result を作る起点として使われたものに近い。

```sparql
SELECT ?country ?countryLabel ?item ?itemLabel ?coord ?elevation WHERE {
  ?country wdt:P31 wd:Q3624078 .   # sovereign state
  ?country wdt:P610 ?item .        # highest point

  ?item wdt:P31 ?instance .
  ?instance wdt:P279* wd:Q8502 .   # mountain

  ?item wdt:P625 ?coord .
  OPTIONAL { ?item wdt:P2044 ?elevation . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "ja,en". }
}
ORDER BY ?countryLabel
```

#### Representative peaks by mountain range

```sparql
SELECT ?range ?rangeLabel ?item ?itemLabel ?coord ?elevation WHERE {
  VALUES ?range {
    wd:Q513   # Himalayas
    wd:Q5463  # Andes
    wd:Q5462  # Alps
  }

  ?item wdt:P4552 ?range .         # mountain range
  ?item wdt:P31 ?instance .
  ?instance wdt:P279* wd:Q8502 .   # mountain
  ?item wdt:P625 ?coord .
  OPTIONAL { ?item wdt:P2044 ?elevation . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "ja,en". }
}
ORDER BY ?rangeLabel DESC(?elevation)
```

mountain 側の実務上のポイントは、WDQS を discovery に留めて、最終採用は curated seed で人手選定することである。

## 4. 共通の運用方針

### 4.1 実行時と生成時を分ける

アプリ実行時は生成処理を行わない。  
生成は開発時のオフライン作業として扱う。

### 4.2 JSON をソース・オブ・トゥルースにする

実行時は bundled JSON のみを参照する。  
CLI の一覧取得や `--show-viewpoint-json` も、ローカル JSON を読むだけで外部アクセスしない。

### 4.3 ASCII フォールバック

`tower_viewpoints.py` / `mountain_viewpoints.py` は、表示名とは別に `ascii_name` を扱う。  
一覧表示や名前解決ではこれが重要なので、生成時の別名クリーニング規則を安定させる必要がある。

### 4.4 IDs を安定キーとして使う

viewpoint の安定キーは表示名ではなく `id` / `qid` で扱う。  
特に都市アウトラインレイヤのような付随データは `name` キーではなく `id` キーにぶら下げる。

## 5. Urban Skyline との関係

都市アウトラインレイヤ用の派生タイル生成は viewpoint dataset の後段処理である。  
流れとしては次の順が前提になる。

1. bundled viewpoint dataset を確定する
2. 対象 viewpoint を選ぶ
3. PLATEAU などの都市建物データから skyline を前計算する
4. `id -> { name, profiles }` 形式で保存する

この手順は [urban-outline-layer-derived-format-ja_JP.md](/home/toshihiro/playground/zstarview/docs/developer/urban-outline-layer-derived-format-ja_JP.md) に分離してある。

## 6. 今後の整理候補

- tower / mountain / 高層ビル展望台を含む `viewpoint` 生成スクリプト群を `src/zstarview/data/` へ整理する
- `dev-samples/` 依存の生成フローを `docs/developer/` の正式文書に寄せる
- `viewpoint_height_m` を実測値や展望フロア高で補正する curated seed を追加する
- `tower_viewpoints.json` というファイル名が実態に合わなくなったら、将来 `viewpoints_high.json` などへの改名を検討する
