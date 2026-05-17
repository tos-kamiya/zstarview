# 観測地点と時刻

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `-p`, `--place QUERY` | OpenStreetMap Nominatim で地名・駅名・施設名を検索し、最上位候補を観測地点として使います。位置引数 `location` とは併用できません。 | |
| `--place-countrycode CODE` | `--place` の検索対象国を ISO 3166-1 alpha-2 形式の国コード（例: `jp`）で制限します。 | |
| `--place-lang LANG` | `--place` の検索結果に対して Nominatim へ送る `Accept-Language` です。 | `en` |
| `--timezone TZ` | 解決された観測地点のタイムゾーンを上書きして、`--datetime` の既定タイムゾーンと画面表示に使います。`JST`、`Asia/Tokyo`、`UTC+9` などを受け付けます。 | |
| `-H`, `--hours HOURS` | 現在時刻に加算する時間数を指定します。※1 | `0` |
| `-D`, `--days DAYS` | 現在時刻に加算する日数を指定します。※1 | `0` |
| `--datetime "YYYY-MM-DD HH[:MM[:SS]] [TZ]"` | 絶対的な日時を指定します。時刻は `HH`、`HH:MM`、`HH:MM:SS` のいずれでも指定でき、タイムゾーン省略時は UTC 扱いです。※1 | |
| `-Z`, `--view-center-az VIEW_CENTER_AZ` | 表示中心の方位角を指定します（度数または方位記号）。 | `180` |
| `-A`, `--view-center-alt VIEW_CENTER_ALT` | 表示中心の高度角を指定します（90=天頂、0=地平線）。 | `90` |
| `--edge-fov-deg DEGREES` | ウィンドウ端に対応する投影スケールを指定します。`96` なら、ウィンドウ端は視線中心から `96°` に対応します。 | `96` |
| `--content-fov-deg DEGREES` | 全レイヤー共通の overscan 視野角を指定します。ウィンドウ端は引き続き視線中心から `90°` の位置に対応し、`90` を超える値では空・雲・背景などがウィンドウ外へはみ出して描画され、四隅の空白を減らせます。許容範囲は `90`〜`127` です。 | `115` |
| `--height-add-m METERS` | 観測基準の上に追加する高さをメートルで指定します。既定の追加高さ `1.7` を置き換えます。タワーや山のビューポイント自体の高さ・標高とは別に扱われます。 | `1.7` |
| `--use-building-top` | 実験中。都市名、`--place`、直接座標、対応 Google Maps URL で解決した地点について、解決地点の約 5m 以内に建物が見つかった場合、その建物の最も高い頂部を観測基準として使います。タワー/山ビューポイントには適用しません。 | off |

注意: `--observer-height-m` は `--height-add-m` の互換オプションとして残ります。
`--use-building-top` を有効にすると、観測基準は近傍の建物頂部に切り替わり、その上に `--height-add-m` が加算されます。

#### `--place` の挙動

`--place` は、通常の offline-first な `location` 引数とは別の、明示的な online resolver 経路です。

- アプリは Nominatim へ 1 回だけ検索リクエストを送り、候補を importance 順で扱います。
- 起動時に対話選択は行わず、最上位候補を自動で使います。
- 候補が複数ある場合も最上位候補を使い、候補一覧はターミナルへ出力します。
- GUI の地点表示には、Nominatim が返した長い地名をそのまま表示します。
- 採用した結果は設定へ保存され、次回起動時は Nominatim へ再問い合わせせず再利用します。
- `--use-building-top` は `--place` と併用できます。このモードは実験中で、ウィンドウ表示前に近傍建物データを解決するため起動が遅くなることがあります。

#### タワー名入力

Wikidata 由来の同梱タワー/展望地点データから起動することもできます。

* 例:
  * `Tokyo Skytree`
  * `t/Tokyo Skytree`（タワー明示指定）
  * `Tsutenkaku`（`Tsūtenkaku` の ASCII 代替表記）
  * `Tokyo Tower`
  * `wikidata:Q57965`
* タワー名を使った場合、そのタワーの構造物高またはビューポイント高が基準観測点として使われます。
* `--height-add-m` は、その基準観測点に加える高さだけを指定します（既定値 `1.7m`）。タワー自体の高さは置き換えません。
* ダイアクリティカルマーク付き名称についても ASCII の代替表記で解決できます。
* 画面上の地点情報では、タワービューポイントに対して `Height: ground ..., building ..., add ...` を表示することがあります。

例:

```bash
zstarview "Tokyo Skytree"
zstarview "Tokyo Tower" --height-add-m 150
```

#### 山名入力

同梱の山ビューポイントデータから起動することもできます。

* 例:
  * `Mount Fuji`
  * `m/Mount Fuji`（山明示指定）
  * `Aconcagua`
  * `Snezka`（`Sněžka` の ASCII 代替表記）
  * `wikidata:Q39231`
* 山名を使った場合、山頂ビューポイントが基準観測点として使われます。
* `--height-add-m` は、その山頂ビューポイントに加える高さだけを指定します（既定値 `1.7m`）。
* 山名についても ASCII の代替表記で解決できます。
* 画面上の地点情報では、山ビューポイントに対して `Height: ground ..., add ...` を表示することがあります。

例:

```bash
zstarview "Mount Fuji"
zstarview "Snezka"
```

#### 日時指定オプションについて

`--datetime "YYYY-MM-DD HH[:MM[:SS]] [TZ]"` で絶対的な日時を指定できます。
時刻部分は「時」だけ、「時:分」、「時:分:秒」のいずれも使用可能です。
タイムゾーン（TZ）を省略した場合は、解決された観測地点のタイムゾーンを使います。`--timezone TZ` を指定した場合はそちらを優先します。

タイムゾーンは以下のいずれかの形式で指定できます:

* よく使われる略称（JST, UTC, GMT, KST, HKT, AWST, ACST, AEST, NZST, NZDT, MSK, EAT）
* IANA タイムゾーン名（例: `Asia/Tokyo`, `Europe/Moscow`）
* UTC オフセット（例: `UTC+9`, `UTC-07:30`）

例:

```bash
zstarview --datetime "2025-08-17 21:00:00 JST" Tokyo
zstarview --datetime "2025-09-12 9" Tokyo         # 9時ちょうど
zstarview --datetime "2025-09-12 09:00" Tokyo     # 9:00
zstarview --datetime "2025-09-12 9:0:0 JST" Tokyo # 9:00:00 JST
```

#### 直接座標入力

都市名の代わりに、直接座標を指定できます。

* 形式:
  * `緯度;経度`（セミコロン区切り）
  * `@緯度,経度`
  * `maps.google.com/` または `www.google.com/maps/` で始まる対応 Google Maps 共有 URL
* 例:
  * `35.68;139.76`
  * `N35.68;E139.76`
  * `-35.68;139.76`
  * `S35.68;W139.76`
  * `@35.68,139.76`
  * `https://www.google.com/maps/@35.68,139.76,17z`
  * `maps.google.com/maps/@35.68,139.76`
  * `https://www.google.com/maps/place/...!3d35.68!4d139.76...`
* 緯度は -90〜90、経度は -180〜180 の範囲でなければなりません。
* `N/S/E/W` の方向記号を使えます（数値の負号がある場合はそちらを優先）。
* 対応する Google Maps URL は、現在広く観測される共有リンク形式として `maps.google.com/` または `www.google.com/maps/` で始まるものを受け付けます。`https://` は省略できます。
* Google Maps URL では、`!3dLAT!4dLON` があればその座標を優先し、なければ `@LAT,LON` を使います。
* zoom、擬似的な高度、heading、pitch などの後続 URL 情報は無視します。
* 直接座標で起動した場合、タイムゾーンは `--place` と同様に解決された地点座標から補完します。`--timezone TZ` を指定した場合はそちらを優先します。
* 追加高さは `--height-add-m` で指定します。Google Maps URL 内の高度らしき値は使いません。
* `--use-building-top` は直接座標入力とも併用できます。この実験中モードでは、解決地点の近傍建物を探し、見つかった場合はその頂部を観測基準として使います。

例:

```bash
zstarview "35.68;139.76"
zstarview "@35.68,139.76"
zstarview "www.google.com/maps/@35.68,139.76,17z"
zstarview "N35.68;E139.76" --datetime "2025-09-12 21 JST"
zstarview "35.68;139.76" --height-add-m 120
```
