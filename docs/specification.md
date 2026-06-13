# zstarview 仕様書

最終更新: 2026-06-09

この文書は、`zstarview` の機能仕様を利用者視点でまとめた正本である。
README より詳細に、何ができるか、どう振る舞うか、どのような制約があるかを記述する。
内部モジュール構成や実装手順は `docs/design.md` と `docs/design/*.md` に分離し、過去の実装履歴はアーカイブとして `docs/implementation-archive.md` に残す。

## 1. この文書の位置づけ

本書は、`zstarview` が利用者に対して何を提供し、どのように振る舞うかをまとめた機能仕様である。
内部モジュール構成やスレッドモデルは `docs/design.md` と `docs/design/*.md`、過去の実装履歴は `docs/implementation-archive.md`、日々の判断記録は `dev-notes/session-YYYY-MM-DD.md` に分離する。

## 2. 製品概要

`zstarview` は、指定した地点と時刻における空の見え方をデスクトップ上で表示するアプリケーションである。

主な目的は次の通り。

- 地点ごとの星空を視覚的に把握できること
- 昼間、薄明、曇天でも空の状態を把握しやすいこと
- 現在時刻だけでなく任意時刻の空を確認できること
- 位置指定、表示方向変更、重ね合わせ表示を軽快に行えること

## 3. 想定利用者と主な用途

### 3.1 想定利用者

- 任意の都市や座標から見た星空を確認したいユーザー
- タワーや展望地点からの見え方を確認したいユーザー
- 観望前に雲や地形も含めて空の状況を把握したいユーザー
- 特定の星、惑星、アステリズム、DSO の位置関係を見たいユーザー

### 3.2 代表的な用途

1. 都市名を指定して現在の星空を表示する。
2. 緯度経度を直接指定して、その地点の空を表示する。
3. タワー名または山名を指定して、そのビューポイントからの見え方を表示する。
4. 視線方向を変えて、天頂から地平線付近まで任意方向を観察する。
5. 指定日時の空を表示して、過去または未来の見え方を確認する。
6. 雲表示や地形地平線を補助情報として重ねて観察する。
7. 星名検索や代表天体ジャンプで目的の天体へ素早く移動する。
8. 準惑星や小天体を検索して、その方向へ視界中心を合わせる。

雲表示は観望補助用の近似表示であり、利用可能な衛星時刻、観測地点、部分欠損に応じて見え方が変わってよい。
雲の状態表示は、衛星ソースの取得中は `Downloading`、取得済みソースからの画像生成中は `Projecting` としてよい。
視線変更中は、旧視線条件で生成された雲画像を暫定表示せず、一時的に非表示としてよい。

## 4. 仕様範囲

### 4.1 対象機能

- 地点解決
  - 都市名
  - 国コード付き都市名
  - GeoNames ID
  - タワー名
  - Wikidata ID
  - 緯度経度
- 時刻指定
  - 現在時刻
  - 相対時間シフト
  - 絶対日時指定
- 星空表示
  - 恒星
  - 太陽、月、惑星
  - 地平線、天の赤道、黄道
  - 天の極マーカー
  - 地平線下の地面ティント
  - 地平線下の地球裏面ガイド
  - 赤緯に基づく never-rises 領域
- 補助表示
  - DSO
  - アステリズム
  - 水面レイヤー
  - 航空機オーバーレイ
  - 雲オーバーレイ
  - 台風・サイクロン補助レイヤー
  - 地形地平線
  - 夜間光オーバーレイ
- GUI 操作
  - キーボード操作
  - ハンバーガーメニュー
  - 通常ウィンドウのメニューバー
  - ホバー情報
  - テーマごとの文字色、文字アウトライン、背景表現
- CLI 起動
- GNOME 向け `.desktop` ファイル生成

### 4.2 対象外

- 学術用途の厳密保証
- 小惑星や人工衛星などの網羅表示
- 複数ユーザー間での設定共有
- オフラインでの雲データ提供保証
- 衛星雲データの完全性、時刻厳密一致、全球一様な解像度保証
- 夜間光データの同梱保証

## 5. 起動と入力

### 5.1 起動形式

`zstarview [options] [location]` 形式で起動する。
通常の GUI 起動では `location` を前回値から引き継いでよい。
初回起動時の既定地点は `Tokyo` である。

### 5.2 地点入力

`location` には次の形式を許容する。

- 都市名
  - 例: `Tokyo`
- 国コード付き都市名
  - 例: `JP/Tokyo`
- GeoNames ID
  - 数値のみ
- タワー名
  - 例: `Tokyo Skytree`
- 明示タワー名
  - 例: `t/Tokyo Skytree`
- 山名
  - 例: `Mount Fuji`
- 明示山名
  - 例: `m/Mount Fuji`
- Wikidata ID
  - 例: `wikidata:Q57965`
- 緯度経度
  - 例: `35.68;139.76`
  - 例: `N35.68;E139.76`
  - 例: `@35.68,139.76`
  - 例: `https://www.google.com/maps/@35.68,139.76,17z`
  - 例: `maps.google.com/maps/@35.68,139.76`
  - 例: `https://www.google.com/maps/place/...!3d35.68!4d139.76...`

Google Maps URL は、`!3dLAT!4dLON` があればその座標を優先し、なければ `@LAT,LON` を用いてよい。
`lat;lon`、`@lat,lon`、Google Maps URL は、緯度 `[-90, 90]`、経度 `[-180, 180]` の範囲チェックを行う。

### 5.3 地点解決の優先順位

- `--place` が指定されている場合は、Nominatim を使った明示的なオンライン検索を優先する。
- `t/NAME` は tower viewpoint としてのみ解決する。
- `m/NAME` は mountain viewpoint としてのみ解決する。
- そのほかの名前入力は、GeoNames、内蔵 viewpoint dataset、緯度経度解析の順序で解決してよい。
- `auto` 指定時は IP アドレスから現在地を推定する。

### 5.4 `--place` 検索

- `--place QUERY`
  - OpenStreetMap Nominatim を使って自由入力文字列を検索する。
  - 位置引数 `location` とは相互排他とする。
  - 対話選択は行わず、正規化後の先頭候補を採用して通常起動する。
- `--place-countrycode CODE`
  - `countrycodes` パラメータとして検索対象国を ISO 3166-1 alpha-2 で制限する。
- `--place-lang LANG`
  - Nominatim への `Accept-Language` として送る。
  - 既定値は `en` とする。

`--place` は 1 回だけ検索を行い、候補一覧を反復照会する経路は持たない。
候補が複数ある場合は、起動時の対話選択を行わず、先頭候補を採用する。
HTTP エラー、レート制限、通信失敗、JSON 解析失敗、0 件結果は起動中断としてよい。

### 5.5 `auto` による現在地取得

`location=auto` の場合は `ip-api.com/json` を使って現在地を推定する。

- 解決結果として緯度、経度、タイムゾーン、都市名、国名を取得してよい。
- 失敗時は起動を中断し、原因を利用者へ表示する。
- この経路は `ip-api.com` の利用条件に従う。

### 5.6 追加の地点補正

- `--use-building-top`
  - 都市名、`--place`、緯度経度、Google Maps URL で解決した地点について、近傍の建物頂部を観測基準にしてよい。
  - tower / mountain viewpoint には適用しない。
  - 候補が複数ある場合は、その建物グループの最大 `height_m` を使ってよい。
  - 候補建物が見つからない場合は地表基準を使う。
- `--height-add-m`
  - 観測基準の上に足す追加高さを指定する。
  - 既定値は `1.7m` とする。
- `--observer-height-m`
  - `--height-add-m` の互換オプションとして扱う。

都市名または緯度経度入力時は、観測基準を地表とし、実効観測点は「地表 + 追加高さ」とする。
山名入力時は山頂ビューポイント、タワー名入力時はタワービューポイントを基準にする。
Google Maps URL に高度らしき数値が含まれていても、追加高さとしては使わない。

### 5.7 時刻仕様

- 現在時刻
  - 実時間表示を使う。
- 相対時間シフト
  - `--hours` または `--days` によるオフセットを適用する。
- 絶対日時表示
  - `--datetime` により日時を明示指定する。

タイムゾーンは、地点から解決した値を使い、`--timezone` が指定された場合は最優先で上書きする。
解決できない場合のみ UTC を使う。

### 5.8 重要な描画入力オプション

以下のオプションは、表示スケール、見える範囲、外観、補助レイヤーの有無を制御する。

- `--edge-fov-deg`
  - ウィンドウ端に対応する角距離を指定する。
  - 投影スケールを変える。
  - 既定値は `95`。
- `--content-fov-deg`
  - 描画対象として保持する角距離上限を指定する。
  - 投影スケールは変えない。
  - 既定値は `110`。
- `-Z`, `--view-center-az`
  - 視線中心の方位を指定する。
- `-A`, `--view-center-alt`
  - 視線中心の仰角を指定する。
- `--sky-opacity`
  - sky disc の色強度を直接制御する。
- `--cloud-opacity`
  - 雲レイヤーの最終合成強度を制御する。
  - 既定値は `0.07` としてよい。
- `--terrain-horizon-opacity`
  - 地形地平線の表示強度を制御する。
  - 既定値は `0.003` としてよい。
- `--earth-guide-opacity`
  - 地平線下の地球裏面ガイドの表示強度を制御する。
  - 既定値は `0.028` としてよい。
- `--night-light-opacity`
  - 夜間光オーバーレイの表示強度を制御する。
  - 既定値は `0.022` としてよい。
- `--water-surface-opacity`
  - 水面レイヤーの表示強度を制御する。
  - 既定値は `0.4` としてよい。
- `--urban-outline-opacity`
  - 都市アウトラインの表示強度を制御する。
- `--overlay-font-size`
  - キャンバス上のラベルや HUD のベース文字サイズを制御する。
- `--visibility-boost`
  - 補助レイヤーや小さい図形の見え方を持ち上げる倍率として扱う。
- `--bright-bodies {outline,fill}`
  - 明るい天体の描画モードを指定する。
  - `outline` では明るい星、太陽、月、惑星の輪郭を強調してよい。
- `--theme`
  - `night`、`day`、`white`、`black`、`transparent`、`transparent-10` から `transparent-90` までの 10 刻みを受け付ける。
  - `transparent` は `transparent-40` の別名として扱ってよい。
- `--window-frame`
  - `frameless` または `window` を受け付ける。

### 5.9 補助レイヤーの有効無効

- `--show-guidelines-initial`
  - 地平線、天の赤道、黄道、方位ラベル、天頂マーカー、never-rises 領域の初期表示を指定する。
- `--cloud-opacity 0`
  - 雲表示を無効化する。
- `--terrain-horizon-opacity 0`
  - 地形地平線を無効化する。
- `--earth-guide-opacity 0`
  - Earth guide を無効化する。
- `--night-light-opacity 0`
  - 夜間光を無効化する。
- `--water-surface-opacity 0`
  - 水面表示を無効化する。
- `--urban-outline-opacity 0`
  - 都市アウトラインを無効化する。

0 を指定して無効化したレイヤーは、そのセッションでは GUI から再有効化できないものがある。

## 6. 表示仕様

### 6.1 表示の基本

画面上には次を表示する。

- 恒星
- 太陽
- 月
- 惑星
- 地平線
- 天の赤道
- 黄道
- 方位ラベル
- 天頂マーカー
- DSO
- アステリズム
- 航空機オーバーレイ
- 雲オーバーレイ
- 地形地平線
- 都市アウトライン
- 地点名、地点要約、時刻、ステータス情報

表示は円形の sky disc を基本とする。
視線中心は方位と高度で指定し、高度は天頂が `90`、地平線が `0` とする。
`--edge-fov-deg` は投影スケール、`--content-fov-deg` は保持する描画対象の角距離上限を決める。
`--content-fov-deg` の外側には、打ち切りを急に見せないための短いフェード帯を置いてよい。

### 6.2 テーマと配色

テーマは、少なくとも次を変えてよい。

- 通常テキスト
- ステータス行テキスト
- ウィンドウ背景
- ウィンドウ枠
- sky disc の不透明度
- スプラッシュスクリーンの配色

`day` と `white` は明るい配色、`night` と `black` は暗い配色、`transparent` 系は半透明の背景を持つテーマとして扱ってよい。
`transparent` 系では、黒寄りの半透明背景と低い sky disc opacity を使ってよい。

### 6.3 Sky Guides

Sky Guides とは、幾何学的地平線、天の赤道、黄道、方位ラベル、天頂マーカー、天の極マーカーをまとめた表示群である。

- 起動時の既定値は表示ありとしてよい。
- GUI 実行中はメニューから切り替えられてよい。
- CLI では `--show-guidelines-initial true|false` で初期状態を指定してよい。
- never-rises 領域と天の極マーカーはこのトグルに従ってよい。

### 6.4 星と天体

- 恒星は等級に応じて明るさとサイズを変えて表示する。
- 色のある恒星は色指数を反映した見た目を持ってよい。
- `vmag < 2.0` の恒星は、通常モードではダイヤ形の強調を重ねてよい。
- 太陽、月、惑星は恒星とは別スタイルで識別しやすく描画する。
- 月は位相を反映して表示する。
- `Enlarge Moon` が有効なときは、月を通常時の見た目半径に対して拡大表示してよい。
- 旧来の「暗い星だけを特別扱いで背景境界外へ出さない」動作は行わない。

### 6.5 地平線と地面

- 地平線上と地平線下が視覚的に区別されること。
- 地平線下の領域には地面色のティントを適用する。
- 地形地平線が有効な場合は地形地平線を境界として使う。
- 地形地平線が無効な場合は幾何学的地平線を使う。
- 地面ティント、Earth guide、never-rises 領域は、用途に応じて別の補助色として扱ってよい。

### 6.6 DSO とアステリズム

- DSO は名称付きの主要天体に限定してよい。
- アステリズムは IAU 星座境界ではなく一般的な線パターンとして表示する。
- どちらも通常は淡い補助表示とし、ホバー時は強調してよい。
- 表示有無は GUI から切り替えられる。

### 6.7 押下中の簡易表示

- マウスを背景上で押下し、まだドラッグやサイズ変更に移行していない間は、表示を一時的に簡略化してよい。
- この押下中の簡易表示は、クリックやメニュー操作ではなく、押下継続中だけ有効な内部状態として扱ってよい。
- 押下中の簡易表示を正式な表示モードとして扱ってよい。
- 押下中は次を非表示にしてよい。
  - 地点名、地点要約、時刻、ステータス情報
  - DSO
  - アステリズム
  - 雲オーバーレイ
  - 夜間光オーバーレイ
  - 航空機オーバーレイ
  - 人工衛星
  - 台風・サイクロン
  - 都市アウトライン
  - 水面レイヤー
  - sky disc
  - Earth guide
  - 副稜線
- 押下中でも、主稜線は fast-mode 相当の経路で、細い線として残してよい。
- 押下中でも、星、太陽、月、惑星、地平線、天の赤道、黄道、方位ラベル、天頂マーカー、天の極マーカー、および hover によるラベル表示は残してよい。
- 押下が終わったら、表示はすぐ通常状態へ戻してよい。

### 6.8 雲オーバーレイ

雲オーバーレイは、利用可能な衛星データに基づく近似表示として扱う。

- `--cloud-opacity 0` はそのセッション中の雲表示を無効化する。
- 取得済みソースがあれば、視点変更では再取得せず投影だけをやり直してよい。
- 取得中は `Downloading`、投影中は `Projecting` と表示してよい。
- 雲データが壊れている、欠損している、または再取得中である場合は、既存画像を維持してよい。

#### 6.8.1 実験的 Geo-satellite 経路

- `Geo-satellite` は、既存の GOES/Himawari 雲オーバーレイとは別の実験的な Europe 補助経路である。
- 対象範囲は、おおむね緯度 `32N` 〜 `73N`、経度 `15W` 〜 `35E` としてよい。
- 観測者がこの責務領域内にいる場合は、雲表示のデータ源として Geo-satellite のみを使ってよい。
- 責務領域の外では、Geo-satellite を無理に使わず、既存の雲表示を使ってよい。
- データ源は MET Norway Geo-Satellite API を使ってよい。
- 画像種別は既定で `infrared` としてよい。

#### 6.8.2 取得失敗時

- 取得失敗時は、既存の雲画像を維持してよい。
- 壊れた画像や未定義の代替画像を表示してはならない。

### 6.9 地形地平線

- DEM から推定した地形地平線を補助線として表示できる。
- 距離帯ごとの副稜線を表示してよい。
- 地形地平線は地点依存であり、時刻依存ではない。
- `opacity == 0` の場合、取得、計算、描画を行わない。
- キャッシュ済み DEM が古い場合でも、再取得失敗時は既存キャッシュを警告付きで使ってよい。

### 6.10 都市アウトライン

- 観測地点周辺の建物輪郭を白線で表示する。
- 通常近距離レイヤーと遠距離スカイスクレーパー補助レイヤーを分けて扱ってよい。
- 取得や再取得はバックグラウンドで行ってよい。
- `opacity == 0` のセッションでは取得しなくてよい。
- 高密度な輪郭は視認性を優先して間引いてよい。

### 6.11 水面レイヤー

- 海と inland water を別ソースとして扱ってよい。
- 水面レイヤーは `Water Surface` と表示してよい。
- inland water は、`ele` / `water_level` がなければ DEM から推定した高さを使ってよい。
- 外部 API へ送る水面取得リクエストは、`zstarview/1.30.6 (+water-overlay)` を付けてよい。
- サービスごとに識別が必要な場合は、`zstarview/1.30.6 (+service)` 形式の接尾辞を足してよい。
- `Terrain Horizon` が無効な間は水面レイヤーも表示しなくてよい。
- `W` で on/off してよい。

### 6.12 航空機オーバーレイ

- 観測地点周辺の民間航空機を補助表示として重畳してよい。
- 視覚的な位置関係の把握を目的とし、管制用途の厳密表示ではない。
- `0.0` で無効化しているセッションでは取得も cache 読込も省略してよい。
- 取得失敗時でも、古いスナップショットが残っていれば継続利用してよい。

### 6.13 人工衛星レイヤー

- 観測地点から見える人工衛星を補助表示として重畳してよい。
- `--search` による検索対象とは別に、GUI では継続表示の対象として扱ってよい。
- 地平線下や詳細条件による可視/不可視の判定を持ってよい。

### 6.14 夜間光オーバーレイ

- 夜間光 GeoTIFF から方位ごとの glow band を生成してよい。
- 副稜線レイヤーを使う場合、配列の `0` 番は最初の副稜線として扱ってよい。
- 表示は地形地平線や副稜線の補助として扱ってよい。
- sky glow は、空ディスクの地平線付近の色と night lights の基準色を混ぜた補助色で描画してよい。夜間などで空ディスクがほぼ黒になっても、sky glow が極端に暗くなりすぎないようにしてよい。
- sky glow の不透明度は street-light 側と独立に制御してよい。

### 6.15 台風・サイクロン補助レイヤー

- 公開 ArcGIS FeatureServer から台風・サイクロンの補助レイヤーを取得してよい。
- 観測補助として扱い、厳密な気象利用を保証しない。

### 6.16 外部 API と `User-Agent`

`zstarview` は、外部 HTTP API を利用するとき、識別可能な `User-Agent` を送ってよい。  
現行のアプリケーション版は `1.30.6` であり、この節の `zstarview/1.30.6` はその現行版を表す。  
将来バージョンを更新する場合は、`zstarview/<current-version>` の基底部分だけを差し替えればよい。  
以下は、現行実装で使っている主要な外部 API と `User-Agent` の対応である。

| 外部 API | 用途 | 現行 `User-Agent` |
| --- | --- | --- |
| Overpass API (`https://overpass-api.de/api/interpreter`) | 水面レイヤーの取得 | `zstarview/1.30.6 (+water-overlay)` |
| OpenStreetMap Nominatim (`https://nominatim.openstreetmap.org/search`) | `--place` 検索 | `zstarview/1.30.6 (+nominatim)` |
| NASA Earth at Night page (`https://science.nasa.gov/earth/earth-observatory/earth-at-night/maps/`) | 夜間光データの取得元ページ参照 | `zstarview/1.30.6 (+night-lights)` |
| Overture release catalog (`https://stac.overturemaps.org/catalog.json`) | Overture 更新確認 | `zstarview/1.30.6 (+overture-release)` |
| MET Norway Geo-Satellite API (`https://api.met.no/weatherapi/geosatellite/1.4/`) | 雲オーバーレイの元画像取得 | `zstarview/1.30.6 (+geosatellite)` |
| ArcGIS FeatureServer (`https://services9.arcgis.com/RHVPKKiFTONKtxq3/arcgis/rest/services/Active_Hurricanes_v1/FeatureServer`) | 台風・サイクロン補助レイヤー | `zstarview/1.30.6 (+tropical-cyclone)` |
| ip-api.com (`http://ip-api.com/json`) | `auto` の現在地推定 | `zstarview/1.30.6 (+ip-api)` |
| OpenSky API (`https://opensky-network.org/api/states/all`) | 航空機オーバーレイ | `zstarview/1.30.6 (+opensky)` |
| CelesTrak (`https://celestrak.org/NORAD/elements/gp.php`) | 衛星 OMM 取得 | `zstarview/1.30.6 (+satellites-celestrak)` |
| JPL Horizons (`https://ssd.jpl.nasa.gov/api/horizons_lookup.api`, `https://ssd.jpl.nasa.gov/api/horizons.api`) | 衛星・小天体の位置取得 | `zstarview/1.30.6 (+satellites-horizons)` |
| WhereTheISS.at (`https://api.wheretheiss.at/v1`) | ISS TLE 取得 | `zstarview/1.30.6 (+satellites-wheretheiss)` |
| Copernicus DEM (`https://copernicus-dem-90m.s3.eu-central-1.amazonaws.com/`) | DEM タイル取得 | `zstarview/1.30.6 (+copernicus-dem)` |
| AWS S3 (`https://*.s3.amazonaws.com/`) | 衛星用の S3 バケット一覧・オブジェクト取得 | `zstarview/1.30.6 (+s3)` |
| Skyfield loader (`https://naif.jpl.nasa.gov/pub/`, `https://ssd.jpl.nasa.gov/`) | ephemeris などの Skyfield 取得 | `zstarview/1.30.6 (+skyfield-loader)` |

## 7. GUI 操作

### 7.1 基本操作

- キーボード操作で視線移動、天体ジャンプ、表示切り替えを行ってよい。
- フレームレスウィンドウでは独自外枠、右上メニュー、右下のサイズ変更グリップを表示してよい。
- 通常ウィンドウでは OS 標準のタイトルバーと枠を使ってよい。

### 7.2 メニュー

- `Sky Guides`
- `Earth Guide`
- `Terrain Horizon`
- `Water Surface`
- `Night Lights`
- `Asterisms`
- `DSO`
- `Aircraft`
- `Clouds`
- `Urban Outline`

などの切り替えを GUI メニューから行ってよい。

### 7.3 起動前設定ダイアログ

`zstarview-gui` は、起動前設定ダイアログを先に開いてから GUI 本体を起動する専用エントリポイントとして扱ってよい。

- タブは `Location`、`View`、`Time`、`Stars`、`Overlays`、`General` の順でよい。
- `Overlays` は `Sky`、`Clouds`、`Tropical Cyclone`、`Aircraft and Satellites`、`Ground and Guides`、`Urban Outline` に分けてよい。
- `City` 欄には複数行入力欄を用いてよい。
- `Auto Search` ボタンは現在地自動取得の結果を反映してよい。
- `Search ...` ボタンは専用の place search dialog を開いてよい。
- `View` タブには視線中心と FOV 系の値を置いてよい。
- `Time` タブでは `Current time`、`Relative time`、`Absolute time` を切り替えてよい。
- `Reset` は前回起動値を既定値へ戻すために使ってよい。
- `Cancel` は GUI 起動を中止してよい。

### 7.4 GUI 検索

GUI では `--search QUERY` に相当する検索を 1 件だけ実行してよい。
検索結果が 1 件の場合は、その対象を中心にして画像を生成してよい。
0 件または複数件の場合は、`--list` がない限り候補一覧を出して終了してよい。

### 7.5 ホバーとハイライト

- ホバー時に天体名やラベルを表示してよい。
- `jump highlight` は短時間の強調表示として扱ってよい。
- label は候補を集約してから、HUD 直前にまとめて描画してよい。

### 7.6 ビューポート操作中の簡易描画

視線変更またはウィンドウリサイズ直後は、短時間の簡易描画モードに入ってよい。
簡易描画では、明るい星のみを即時再投影して表示し、重い補助レイヤーは後段更新でもよい。

## 8. CLI 派生コマンド

### 8.1 `zstarview-export-image`

`zstarview-export-image` は、1 枚の画像を生成して終了する headless CLI である。

- GUI ウィンドウは常駐しない。
- 地点、時刻、視線、テーマ、各レイヤー opacity、`--enlarge-moon`、`--bright-bodies` などを受け付けてよい。
- `--output` は必須としてよい。
- `--image-size` で出力サイズを指定してよい。
- `--layer-timeout-seconds` で各レイヤーの待ち時間上限を指定してよい。
- `--allow-partial-data` により部分データ保存を許可してよい。
- `--include-direction-grid` で静的な方向グリッドを重ねてよい。
- `--sixel` で端末へ SIXEL 出力してよい。

既定では、必要なレイヤーが揃うまで保存してはならない。
`--output -` は PNG bytes を stdout へ流す用途とし、`--sixel` と併用してはならない。
`--sixel` 指定時は、事前に `img2sixel` の存在確認と SIXEL 対応端末の確認を行ってよい。

### 8.2 dataset 参照専用 CLI

次の即時終了オプションは、GUI を起動せずローカルデータだけを参照してよい。

- `--list-viewpoints KIND`
- `--list-viewpoint-names KIND`
- `--show-viewpoint-json NAME`

これらは tower / mountain の内蔵 dataset のみを対象とし、GeoNames 読込や設定保存を行わない。

## 9. 設定保持

- `zstarview` と `zstarview-export-image` は `~/.config/zstarview/config.json` を使い、前回地点と window geometry を引き継いでよい。
- `zstarview-gui` は `~/.config/zstarview/gui-launch-profile.json` を別に使ってよい。
- `zstarview` と `zstarview-export-image` は、`zstarview-gui` 専用の保存ファイルを読んで既定値を上書きしてはならない。
- 前回地点は文字列形式だけでなく、Nominatim 由来の構造化オブジェクトを保存してよい。

## 10. エラーと状態表示

### 10.1 起動失敗

- 地点解決失敗
- `--place` 検索失敗
- 入力不正
- 初期化失敗
- 画像保存失敗

これらは非 0 終了としてよい。
原因は terminal や splash に表示してよい。

### 10.2 継続可能な失敗

- 雲取得失敗時でも星空表示は継続する。
- 地形地平線取得失敗時でも本体表示は継続する。
- 航空機データ取得失敗時でも本体表示は継続する。
- 一時的なネットワーク不調で全体が異常終了しないこと。

### 10.3 ステータス表示

- 雲、地形地平線、航空機、都市アウトラインの取得状態や失敗状態を表示してよい。
- 都市アウトラインは、必要に応じて base レイヤーと skyscraper レイヤーの件数を `base+skyscraper` で表示してよい。
- ステータス表示は画面操作を妨げない補助情報として提供する。

## 11. 非機能仕様

### 11.1 応答性

- UI は操作に対して即応すること。
- 重い処理は利用者操作を長時間ブロックしないこと。

### 11.2 耐障害性

- 外部データ取得失敗が発生しても、本体表示機能はできるだけ維持すること。

### 11.3 保守性

- 利用者向け仕様、内部設計、実装履歴を別文書として維持すること。

### 11.4 出力の文字種

- terminal、console、log、CLI help、exception message、subprocess の stdout/stderr に出る文字列は、ASCII-only を原則とする。
- UI 専用文字列は、画面表示のために必要な場合だけ非 ASCII を許容してよい。
- 外部データ由来の非 ASCII を出力に含める必要がある場合は、要約または escape して出してよい。

## 12. 既知の制約

- 雲表示は外部の衛星データ供給状況に依存する。
- 航空機表示は OpenSky Network の提供状況、認証条件、利用枠に依存する。
- 航空機表示は観測値と短時間予想に基づく近似表示であり、実際の飛行位置と一致しない場合がある。
- オフライン環境では雲表示は利用できない。
- 高密度な表示設定では計算量と描画負荷が増える。
- 表示精度は利用ライブラリ、元データ品質、地点解決結果に依存する。

## 13. 補足

- 詳細な CLI オプション一覧は `README.md` を正本とし、本書では動作仕様を優先している。
- 開発用の投影検証アセットは `docs/design.md` 側で扱う。
