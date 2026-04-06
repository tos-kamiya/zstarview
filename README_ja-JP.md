# zstarview 🌌

雲があっても、太陽が出ていても、満天の星空を。

**Zenith Star View** は、選んだ場所の空を表示するデスクトップ向けのスカイビューアです。

指定した場所と時刻の天球に、恒星、太陽、月、惑星、DSO、アステリズムを表示します。
必要に応じて、地形地平線、都市アウトライン、近傍の航空機、人工衛星も重ねて表示できます。
観測地点は都市名やビューポイント名、緯度経度、オンライン地名検索、Google Maps の URL などで指定できます。

**特徴:**

- **太陽系天体**: 太陽・月・主要惑星に対応しています。小惑星（アステロイド）は未対応です。
- **DSO 表示**: 名前付きの DSO（銀河/散開星団/球状星団）を薄い青系の領域として表します。
- **アステリズム表示**: （IAU の正式な星座境界ではなく通称のパターンとしての）星座（アステリズム）を暗い線で常時表示します。アステリズムに含まれる恒星にマウスホバーすると、そのアステリズムを明るく強調してラベルを表示します。同一の恒星にが複数のアステリズムに含まれる場合は 3 秒ごとに切り替えます。
- **衛星雲画像**: リアルタイムに Himawari/GOES 衛星のデータをダウンロードし、縞模様（ハッチ）の重ね描きとして表示します。衛星データが部分的な場合は欠損領域を薄い黄色で示します。[部分カバー時の黄色い欠損表示の例](docs/images/screenshot5.png) も参照してください。
- **航空機オーバーレイ**: OpenSky の近傍航空機を、予想移動方向付きの紫系ポリラインとして表示できます。
- **人工衛星オーバーレイ**: ISS を、惑星レイヤーと航空機レイヤーの間に小さな紫色のクロスマーカーとして表示できます。
- **都市アウトライン表示**: 現在の観測地点に対して、主要な建物屋根線を白い都市アウトラインとして表示します。高層建築が多い一部の都市では、半径 60km 以内の遠距離スカイスクレーパーも追加で表示されます。
- **地形地平線と地面塗り**: Copernicus DEM データをダウンロードして、地形地平線オーバーレイを表示します。観測者の地点に沿った、薄い黄土色がかった地形線を表示します。地形地平線（地形地平線を表示しない場合には水平線）より下は、向きの把握を助けるため地面色で塗り分けます。
- **ガイド表示**: 昇らない領域の赤い表示、地平線まわりの方位ラベル、天頂マーカーなどのガイドを重ねて表示します。
- **柔軟な場所指定**: CLIの引数により観測者の地点を、都市名、タワー名、山名、緯度経度、対応する Google Maps 座標 URL、または Nominatim を使った地名・駅名検索により指定できます。
- **表示中心の調整**: CLIのオプション `-A`（高度）/`-Z`（方位）、あるいは矢印キーで表示中心を調整でき、視線の変更中やウィンドウのサイズ変更中は一時的に簡易描画モードへ切り替えて応答性を保ちます。
- **端末向け画像出力**: `zstarview-export-image` により、ヘッドレスで空を描画してファイルへ保存したり、sixel 対応端末へ直接表示したりできます。
- **Python 対応**: CPython 3.10, 3.11, 3.12, 3.13 で継続的にテストしています。

## スクリーンショット

1枚目の画像は、アステリズム表示と地形地平線の例を示しています。
2枚目の画像は、航空機オーバレイと昇らない領域を示しています。
3枚目の画像は、`-V10.5 -s4.5` でより高密度に星を描画した例です。
4枚目の画像は、`zstarview-export-image` による sixel 端末出力の例です。

  <p align="center">
    <img src="docs/images/screenshot1.png" alt="アステリズム表示と地形地平線の例を示すスクリーンショット" width="49%" />
    <img src="docs/images/screenshot4.png" alt="航空機オーバレイと昇らない領域を示すスクリーンショット" width="49%" />
  </p>

  <p align="center">
    <img src="docs/images/screenshot3.png" alt="-V10.5 -s4.5 で高密度な星空を描画したスクリーンショット" width="49%" />
    <img src="docs/images/screenshot6.png" alt="zstarview-export-image による sixel 端末出力のスクリーンショット" width="49%" />
  </p>

注意: 等級上限を大きくすると描画時間も増えます。[等級上限オプションについて](#about-magnitude-limit) も参照してください。

都市アウトライン表示の都市別スクリーンショット例:

<table>
  <tr>
    <td align="center" width="25%"><img src="docs/images/screenshot-tokyotower.png" alt="東京タワー付近（東京）" width="100%" /></td>
    <td align="center" width="25%"><img src="docs/images/screenshot-downtowndubai.png" alt="ダウンタウン・ドバイ" width="100%" /></td>
    <td align="center" width="25%"><img src="docs/images/screenshot-marinabay.png" alt="マリーナベイ（シンガポール）" width="100%" /></td>
    <td align="center" width="25%"><img src="docs/images/screenshot-sydney.png" alt="サーキュラー・キー（シドニー）" width="100%" /></td>
  </tr>
  <tr>
    <td align="center"><sub>東京タワー付近（東京）</sub></td>
    <td align="center"><sub>ダウンタウン・ドバイ</sub></td>
    <td align="center"><sub>マリーナベイ（シンガポール）</sub></td>
    <td align="center"><sub>サーキュラー・キー（シドニー）</sub></td>
  </tr>
</table>

## インストール方法（推奨：`pipx`）

**前提条件:**

都市アウトライン表示を使うには、`overturemaps` CLI も別途インストールしてください。
インストール方法は <https://pypi.org/project/overturemaps/> を参照してください。

次で確認できます。

```bash
overturemaps --help
```

Zstarview 自体は、[`pipx`](https://pypa.github.io/pipx/) を使ってインストールする想定です。

> 注記: 主に検証している環境は Linux x86_64 です。cloud-disc 系は
> `pyresample` 依存を削除したため、以前の Windows Arm64
> インストール阻害要因は解消されています。

```bash
pipx install zstarview
```

アップグレード:

```bash
pipx upgrade zstarview
```

> 注記: X11 ライブラリやネットワークが細い場合の回避策などは、下のトラブルシューティングを参照してください。

## 使い方

```bash
zstarview [options] [location]
```

> 注記（Ubuntu/Wayland, GNOME）: ターミナル起動時にタスクバーのアイコンが表示されない場合は、後述の [ツール](#ツール) 内の [`.desktop` ランチャーの生成（GNOME専用）](#desktop-ランチャーの生成gnome専用) を実行してください。

よく使う起動例:

```bash
zstarview Tokyo
zstarview "Tokyo Skytree"
zstarview "35.68;139.76" --datetime "2025-09-12 21 JST"
zstarview --place "Matsue Station" --place-countrycode jp
```

CLI では、場所・時刻・描画設定を細かく指定できます。

<details>
  <summary>CLI リファレンス</summary>

### CLI

#### 引数

| 引数 | 説明 | デフォルト |
| :--- | :--- | :--- |
| `location` | 表示する都市名、タワー名、山名、明示指定の `t/NAME`・`m/NAME`、または `"<lat>;<lon>"`、`"@<lat>,<lon>"`、対応する Google Maps URL などの直接座標形式を指定できます。例: `Tokyo`, `Tokyo Skytree`, `t/Tokyo Skytree`, `Mount Fuji`, `m/Mount Fuji`, `35.68;139.76`, `N35.68;E139.76`, `@35.68,139.76`, `www.google.com/maps/@35.68,139.76,17z`, `www.google.com/maps/place/...!3d35.68!4d139.76...`。省略時は前回起動時の location を使い、初回は `Tokyo` になります。 | 前回の location（初回は `Tokyo`） |

#### 観測地点と時刻

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
| `--content-fov-deg DEGREES` | 全レイヤー共通の overscan 視野角を指定します。ウィンドウ端は引き続き視線中心から `90°` の位置に対応し、`90` を超える値では空・雲・背景などがウィンドウ外へはみ出して描画され、四隅の空白を減らせます。許容範囲は `90`〜`127` です。 | `100` |
| `--observer-height-m METERS` | 観測地点の基準面から見た観測者の目線高さをメートルで指定します。既定の目線高さ `1.7` を置き換えます。タワーや山のビューポイント自体の高さ・標高とは別に扱われます。 | `1.7` |
| `--use-building-top` | 実験中。都市名、`--place`、直接座標、対応 Google Maps URL で解決した地点について、解決地点の約 5m 以内に建物が見つかった場合、その建物の最も高い頂部を観測基準として使います。タワー/山ビューポイントには適用しません。 | off |

#### 星空と天体

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `--sky-opacity SKY_OPACITY` | 空の色ディスクの不透明度を指定します（0.0〜1.0）。0.0 で描画を無効化します。 | `0.15` |
| `-m`, `--enlarge-moon` | 月を 5 倍に拡大して表示します。 | |
| `-s`, `--star-base-radius STAR_BASE_RADIUS` | 2 等星の基本サイズを指定します。 | `4.0` |
| `-w`, `--expected-render-width EXPECTED_RENDER_WIDTH` | 恒星をフル解像度で描画する想定ウィンドウ幅を指定します。天球幅がこの値を超える場合、恒星レイヤーは平方根スケーリングで描画します。 | `600` |
| `-V`, `--vmag-limit V_MAG_LIMIT` | 表示する恒星の等級上限を指定します。 | `6.0` |
| `--vmag-brightness-multiplier MULTIPLIER` | 等級 1 段階あたりの光量変化倍率（`1.58`〜`2.512`、デフォルト `2.5`。Pogson の定義は `2.512`）を指定します。※3 | `2.5` |
| `-i`, `--sky-update-interval SECONDS` | 星空を更新する時間間隔（秒）を指定します。 | `60` |
| `--show-dso-initial true\|false` | 起動時に DSO を表示するかを指定します。 | 自動（カタログがあれば表示） |
| `--show-asterisms-initial true\|false` | 起動時にアステリウムを表示するかを指定します。 | `show` |
| `--show-observation-info-initial true\|false` | 起動時に観測情報ブロックを表示するかを指定します。 | `show` |

#### オーバーレイ

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `-c`, `--cloud-opacity CLOUD_OPACITY` | 雲の不透明度を指定します（0.0〜1.0）。0.0 で描画を無効化します。※2 | `0.07` |
| `--cloud-stripe MODE[,COUNT[,WIDTH]]` | 雲ストライプの方式を指定します。`width` は雲量に応じて見かけの線幅を変え、`alpha` は線幅を固定したまま線の alpha を変えます。`width` は `width,50,0.85`、`alpha` は `alpha,50,0.2` に展開されます。count または width を `0` にすると雲描画を無効化します。 | `width,50,0.85` |
| `--cloud-missing-tint-opacity OPACITY` | 雲欠損領域を示す黄色の濃さを指定します（0.0〜1.0）。 | `0.176` |
| `-a`, `--aircraft-opacity OPACITY` | 航空機オーバーレイの不透明度を指定します（0.0〜1.0）。0.0 で、その起動中の航空機問い合わせと描画を無効化します。 | `0.5` |
| `--satellite-opacity OPACITY` | 人工衛星オーバーレイの不透明度を指定します（0.0〜1.0）。0.0 で、その起動中の軌道要素取得と描画を無効化します。 | `0.5` |
| `--show-guidelines-initial true\|false` | 起動時にガイドライン表示を有効にするかを指定します。対象は幾何学的地平線、天の赤道、黄道、方位ラベル、天頂マーカーです。 | `show` |
| `--terrain-horizon-opacity OPACITY` | 地形地平線ポリラインの不透明度を指定します（0.0〜1.0）。0.0 で DEM ダウンロード・地形地平線計算・描画を無効化します。※4 | `0.05` |
| `--ground-tint-opacity OPACITY` | 幾何学的地平線または地形地平線より下の地面色塗りの強さを指定します（0.0〜1.0）。 | `0.1` |
| `--urban-outline-opacity OPACITY` | 都市アウトライン重ね表示の不透明度を指定します（0.0〜1.0）。0.0 でその起動中は表示を無効化します。 | `0.2` |
| `-r`, `--urban-outline-radius-km RADIUS_KM` | 観測地点からこの半径内の建物を都市アウトラインとして取得・描画します。この値はキャッシュキーにも含まれます。 | `2.5` |
| `--urban-outline-skyscraper-radius-km RADIUS_KM` | 遠距離スカイスクレーパー補助レイヤーの外側半径です。`0` を指定するとその起動中は skyscraper tile 探索を無効化します。それ以外の値は `--urban-outline-radius-km` 以上でなければなりません。 | `60.0` |
| `-b`, `--urban-outline-min-building-height-m METERS` | この高さ未満の建物を都市アウトライン取得時に除外します。この値はキャッシュキーにも含まれます。 | `0.0` |
| `--urban-outline-feature-type {both,building}` | 都市アウトライン用の Overture キャッシュモードを指定します。`both` は `building` と `building_part` を組み合わせ、part がある場合はそちらを優先します。 | `both` |

#### 一般

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `-h`, `--help` | ヘルプメッセージを表示して終了します。 | |
| `--window-geometry restore\|X,Y,W,H` | 初期ウィンドウ位置と大きさを指定します。`restore` で前回終了時の位置/サイズを復元し、`X,Y,W,H` で整数値を直接指定できます。Wayland ではウィンドウ位置の復元は利用できません（サイズ復元は有効です）。 | |
| `--window-frame {frameless,window}` | ウィンドウ装飾モードを選びます。`frameless` は従来の枠なし表示、`window` は OS 標準のタイトルバーと枠を使います。 | `frameless` |
| `-t`, `--theme {night,day,white,black}` | 背景と星の見え方のテーマを指定します。 | `night` |
| `--clear-long-lived-cache` | トラブルシュート用オプションです。起動前に長寿命の DEM / 都市アウトラインキャッシュを削除します。3 日以内に再度使うと起動を拒否し、再実行可能日時を表示します。 | |

※1 `--hours`、`--days`、`--datetime` でリアルタイムではない星空を表示した場合、雲、航空機、人工衛星は描画されません。

※2 雲の描画は気象衛星（**Himawari** / **NOAA GOES**）の赤外線データを公開 S3 バケットから取得して行います。ネットワーク関連の注意や回避策は「トラブルシューティング」を参照してください。

※3 最も明るい等級差の倍率は、古典的な Pogson 値 \(100^{1/5}\approx2.512\) を超えられません。

※4 地形地平線表示は初回利用時に Copernicus DEM タイルをダウンロードし、以後はローカルキャッシュを再利用します。有効時はディスク内の地面/空の塗り分け境界にも地形プロファイルを使います。

※5 `--place` は公開の OpenStreetMap Nominatim 検索サービスを使います。User-Agent と Accept-Language を付けて 1 回だけ検索リクエストを送ります。高頻度利用や自動化で使う場合は、Nominatim の利用ポリシーを確認してください。

#### `--place` の挙動

`--place` は、通常の offline-first な `location` 引数とは別の、明示的な online resolver 経路です。

- アプリは Nominatim へ 1 回だけ検索リクエストを送り、候補を importance 順で扱います。
- 起動時に対話選択は行わず、最上位候補を自動で使います。
- 候補が複数ある場合も最上位候補を使い、候補一覧はターミナルへ出力します。
- GUI の地点表示には、Nominatim が返した長い地名をそのまま表示します。
- 採用した結果は設定へ保存され、次回起動時は Nominatim へ再問い合わせせず再利用します。
- `--use-building-top` は `--place` と併用できます。このモードは実験中で、ウィンドウ表示前に近傍建物データを解決するため起動が遅くなることがあります。

#### 起動時のオーバーレイ表示設定

起動時の初期表示をメニュー操作なしで切り替えるには、次を使います。

```bash
# DSO は非表示、アステリウムは表示で起動
zstarview --show-dso-initial false --show-asterisms-initial true Tokyo
```

#### 表示中心オプションについて

`-Z`（方位角）と `-A`（高度角）のオプションで、画面の表示中心を指定できます。

デフォルトでは `-Z 180`（南向き）、`-A 90`（天頂）です。画面下が南、画面左が東で、天頂を見上げたような円形の表示になります。

例えば、`-Z 90`（東向き）および `-A 25`（高度 25°、地平線から 25° 上）を指定すると、東の空を見上げる視野が生成されます。  

方位角は度数または方位記号（大小区別なし）で指定できます。
例: `-Z E`, `-Z ne`, `-Z SSW`（202.5°）。
（対応表: 0=N, 90=E, 180=S, 270=W。N, NNE, NE, ENE, E, ESE, SE, SSE, S, SSW, SW, WSW, W, WNW, NW, NNW を受け付けます）

<a id="about-magnitude-limit"></a>

#### 等級上限オプションについて

`-V 等級` で、指定した等級までの明るさの星を描画します。
デフォルトは `-V 6.0` です。現在の同梱星表では最大 `-V 10.5` までサポートしており、その場合の描画候補となる恒星は約 536,000 個です。
この値を大きくすると描画時間も増えます。

#### テーマプリセットについて

`--theme` を使うと、背景とコントラストの見え方を切り替えられます。

* `night`: 標準の暗色テーマ
* `black`: より黒く不透明な背景
* `day`: 明るい昼空寄りの背景表現
* `white`: 最も明るい淡色テーマ

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

#### 直接座標指定

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
* 観測者の目線高さは `--observer-height-m` だけで指定します。Google Maps URL 内の高度らしき値は使いません。
* `--use-building-top` は直接座標入力とも併用できます。この実験中モードでは、解決地点の近傍建物を探し、見つかった場合はその頂部を観測基準として使います。

例:

```bash
zstarview "35.68;139.76"
zstarview "@35.68,139.76"
zstarview "www.google.com/maps/@35.68,139.76,17z"
zstarview "N35.68;E139.76" --datetime "2025-09-12 21 JST"
zstarview "35.68;139.76" --observer-height-m 120
```

#### タワー名入力

Wikidata 由来の同梱タワー/展望地点データから起動することもできます。

* 例:
  * `Tokyo Skytree`
  * `t/Tokyo Skytree`（タワー明示指定）
  * `Tsutenkaku`（`Tsūtenkaku` の ASCII 代替表記）
  * `Tokyo Tower`
  * `wikidata:Q57965`
* タワー名を使った場合、そのタワーの構造物高またはビューポイント高が基準観測点として使われます。
* `--observer-height-m` は、その基準観測点からの観測者の目線高さだけを置き換えます（既定値 `1.7m`）。タワー自体の高さは置き換えません。
* ダイアクリティカルマーク付き名称についても ASCII の代替表記で解決できます。
* 画面上の地点情報では、タワービューポイントに対して `Tower height ... m` を表示することがあります。`--observer-height-m` を明示指定した場合は、別行で `Observer height ... m` を表示することがあります。

例:

```bash
zstarview "Tokyo Skytree"
zstarview "Tokyo Tower" --observer-height-m 150
```

山名の同梱ビューポイントデータから起動することもできます。

* 例:
  * `Mount Fuji`
  * `m/Mount Fuji`（山明示指定）
  * `Aconcagua`
  * `Snezka`（`Sněžka` の ASCII 代替表記）
  * `wikidata:Q39231`
* 山名を使った場合、山頂ビューポイントが基準観測点として使われます。
* `--observer-height-m` は、その山頂ビューポイントからの観測者の目線高さだけを置き換えます（既定値 `1.7m`）。
* 山名についても ASCII の代替表記で解決できます。
* 画面上の地点情報では、山ビューポイントに対して `Elevation ... m` を表示することがあります。`--observer-height-m` を明示指定した場合は、別行で `Observer height ... m` を表示することがあります。

例:

```bash
zstarview "Mount Fuji"
zstarview "Snezka"
```

`--datetime` のタイムゾーン指定例:

- IANA ゾーン名: `--datetime "2025-09-12 21 Asia/Tokyo"`
- UTC オフセット: `--datetime "2025-09-12 21 UTC+9"`

### 対応しているアステリウム

ここでの「アステリウム」は、見かけ上の星の並びを結んだ通称パターンです。  
**IAU（国際天文学連合）が定義する正式な 88 星座の境界とは別概念**です。

- 冬: `Winter Triangle`（冬の大三角）, `Orion's Belt`（オリオンの三ツ星）, `Winter Hexagon`（冬のダイヤモンド）, `Southern Cross`（南十字）, `Southern Pointers`, `Diamond Cross`, `False Cross`
- 春: `Big Dipper`（北斗七星）, `Little Dipper`, `Spring Triangle`, `Arc to Arcturus`, `Leo Sickle`, `Southern Triangle`
- 夏: `Summer Triangle`（夏の大三角）, `Northern Cross`（北十字）, `Teapot`, `Keystone`
- 秋: `Great Square of Pegasus`（ペガススの四辺形）, `Circlet of Pisces`, `Water Jar of Aquarius`, `Cassiopeia W`, `House of Cepheus`, `Job's Coffin`

</details>

<details>
  <summary>ツール</summary>

## ツール

### `zstarview-export-image`

GUI を起動せずに画像ファイルを 1 枚生成して終了する場合は、`zstarview-export-image` を使います。

```bash
zstarview-export-image Matsue -o matsue.png
```

`zstarview-export-image` は、通常 GUI 左上に出る場所・時刻・視線方向・vmag limit の要約を、描画後に `stderr` へ出力します。`--sixel` の場合は端末画像を出す直前に出力します。

トラブルシュートや手動キャッシュ削除のために、描画せず cache root だけを表示することもできます。

```bash
zstarview-export-image --print-cache-dir
```

`--clear-long-lived-cache` のクールダウンを回避したい場合は、まず `zstarview-export-image --print-cache-dir` で cache root を確認し、その配下の次のサブディレクトリを起動前に手動削除してください。

- `copernicus-dem`
- `overture_buildings`
- `overture_skyscrapers`

### ビューポイントデータ参照ツール

GUI を起動せずに、同梱タワー/展望地点データと山頂ビューポイントデータを参照できます。

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `-h`, `--help` | ヘルプメッセージを表示して終了します。 | |
| `--list-viewpoints {t,m}` | 同梱タワー (`t`) または山 (`m`) の主表示名を出力して終了します。各行は `t/NAME` または `m/NAME` 形式で、利用可能な場合は ASCII 代替名を優先します。 | |
| `--list-viewpoint-names {t,m}` | 同梱タワー (`t`) または山 (`m`) の名前を、多言語名と ASCII 代替名込みで一覧出力して終了します。各行は `t/NAME` または `m/NAME` 形式です。 | |
| `--show-viewpoint-json NAME` | 指定名で同梱ビューポイントを解決し、利用可能な場合は `ascii_name` を含む JSON メタデータを出力して終了します。`t/` または `m/` を付けると対象 kind を明示できます。 | |

```bash
zstarview --list-viewpoints t
zstarview --list-viewpoint-names t
zstarview --show-viewpoint-json "t/Tokyo Skytree"
zstarview --list-viewpoints m
zstarview --show-viewpoint-json "m/Mount Fuji"
```

これらのオプションは相互排他で、`location` 引数や時刻・描画オプションとは併用できません。
`--list-viewpoints` では、利用可能な場合は ASCII 代替名を優先表示します。
`--list-viewpoint-names` では、元の綴りと ASCII 代替綴りの両方を含みます。
prefix なしの `--show-viewpoint-json` で tower と mountain の両方に完全一致した場合は、`t/...` / `m/...` 候補を列挙して曖昧一致エラーにします。

### `.desktop` ランチャーの生成（GNOME専用）

GNOME 系デスクトップ環境（Ubuntu Dock や DockToPanel を含む）では、
タスクバーに正しいアイコンを表示するために `.desktop` ファイルが必要です。

本アプリにはこれを生成する補助コマンドが付属しています。

```bash
# カレントディレクトリに zstarview.desktop を作成
zstarview-make-desktop-file

# ~/.local/share/applications にインストール
zstarview-make-desktop-file --write
```

* `--write` を付けない場合は、カレントディレクトリに `zstarview.desktop` が作成されます。
* `--write` を付けると `~/.local/share/applications` に書き込み、デスクトップデータベースに登録します。

> **注:** このランチャー機能は GNOME 系環境専用です。  
> 他のデスクトップ環境では不要、または正しく動作しない場合があります。

</details>

GUI では、キーボード操作とメニュー操作で視点移動、検索、各種オーバーレイ切り替えを行えます。

<details>
  <summary>GUI 操作</summary>

### GUI

#### キー操作

* **← / →**: 視線の方位を ±5° 回転
* **↑ / ↓**: 視線の高度を ±5° 変更（0°..90° にクランプ）
  方向キー入力が続く間と最後の入力から約 0.7 秒の間は、ビューポート操作用の簡易描画モードになります。この間は `Vmag <= 4.0` の恒星、天の赤道、黄道、地平線、地形地平線、方位ラベル、天頂マーカーのみを表示し、惑星、全星等の星空、空ディスク、雲、DSO、アステリウム、都市アウトラインは一時的に非表示になります。
* **M**: 月の 5 倍表示をトグル
* **D**: DSO 重ね表示の表示/非表示を切り替え
* **A**: アステリウム重ね表示の表示/非表示を切り替え
* **G**: ガイドライン表示の表示/非表示を切り替え
* **S**: 空ディスク表示をグラデーションとフラットディスクで切り替え
* **C**: 雲の重ね表示の表示/非表示を切り替え
* **P**: 航空機オーバーレイの表示/非表示を切り替え
* **I**: 人工衛星オーバーレイの表示/非表示を切り替え
* **T**: 地形地平線の重ね表示の表示/非表示を切り替え
* **U**: 都市アウトラインの重ね表示の表示/非表示を切り替え
* **Ctrl+J**: Jump to Named Star を開く
* **Ctrl+F**: Search Stars and Asterisms を開く
* **F11**: フルスクリーン表示の切り替え
* **ESC**: フルスクリーンから復帰
* **Q**: 終了

#### メニュー操作

ハンバーガーメニュー（`☰`）から次を利用できます。

* **Jump to Named Star...**: 代表的な固有名星（`Vmag <= 2.0`）を北天 / 赤道付近 / 南天で選んで、視点中心をその星へ移動します。
* **Search Stars and Asterisms...**: 固有名付き恒星、対応アステリウム、ISS を横断検索し、選択した対象へ移動します。
* **Search Places...**: OpenStreetMap Nominatim を使う別ダイアログを開き、地名・駅名・施設名の候補から選んだ地表地点の方向へ視点中心を移動します。
* **Enlarge Moon**: 月の 5 倍表示を切り替えます。
* **DSO**: DSO の重ね表示の表示/非表示を切り替えます。
* **Asterisms**: アステリウムの重ね表示の表示/非表示を切り替えます（有効時は暗い線を常時表示し、構成星にホバーすると該当アステリウムを明るく強調してラベルを表示します）。
* **Guidelines**: 幾何学的地平線、天の赤道、黄道、方位ラベル、天頂マーカーの表示/非表示を切り替えます。
* **Sky Color Disc**: 空ディスク表示を、空色グラデーション表示とフラットな暗色ディスク表示で切り替えます。
* **Clouds**: リアルタイム雲の重ね表示の表示/非表示を切り替えます。
* **Aircraft**: OpenSky ベースの航空機オーバーレイの表示/非表示を切り替えます。CLI で `-a 0` / `--aircraft-opacity 0` を指定して起動した場合、その起動中はメニューから再有効化できません。
* **Satellites**: ISS の人工衛星オーバーレイの表示/非表示を切り替えます。CLI で `--satellite-opacity 0` を指定して起動した場合、その起動中はメニューから再有効化できません。
* **Terrain Horizon**: 地形地平線の重ね表示の表示/非表示を切り替えます。CLI で `--terrain-horizon-opacity 0` を指定して起動した場合、その起動中はメニューから再有効化できません。
* **Urban Outline**: 都市アウトラインの重ね表示の表示/非表示を切り替えます。CLI で `--urban-outline-opacity 0` を指定して起動した場合、その起動中はメニューから再有効化できません。
* **Fullscreen**: フルスクリーン表示を切り替えます。
* **Exit**: アプリケーションを終了します。

ジャンプ/検索の確定後は約 3 秒間、マウスホバー時と同じ見た目（円マーカー + 名称ラベル）で対象星を強調表示します。

ウィンドウサイズ変更中も同じ簡易描画モードを使うため、表示の追随性を保ちます。

</details>

<details>
  <summary>都市アウトライン用データ</summary>

`zstarview` は都市アウトライン用の元データを Overture Maps から必要時に取得し、
アプリのキャッシュディレクトリに派生タイルとして保存します。新しい
地点・半径・建物高さ条件の組み合わせでは、初回だけ数秒待ってから都市
アウトラインが表示されます。

この機能には `overturemaps` CLI の別途インストールが必要です。次で確認できます。

```bash
overturemaps --help
```

起動例:

```bash
zstarview "Tokyo Tower" -r 2.5 -b 0
zstarview -p "Matsue Station" -r 2.0 -b 20
```

- `-r`, `--urban-outline-radius-km`: 取得半径（km）
- `-b`, `--urban-outline-min-building-height-m`: 建物の最小高さ（m）
- `--urban-outline-feature-type`: Overture キャッシュモード。既定値は `both`

キャッシュキーには地点・半径・最小建物高さが含まれるため、これらを変えると
別のキャッシュデータセットが作られます。

</details>

<details>
  <summary>トラブルシューティングとプラットフォーム補足</summary>


## トラブルシューティング

### X11（Ubuntu/Debian）

Qt の xcb プラグインが `libxcb-cursor0` を必要とする場合があります。
X11/Wayland を意識していないと分かりづらいですが、ターミナルから実行すると次のようなエラーが表示されます:

```sh
$ zstarview
qt.qpa.plugin: From 6.5.0, xcb-cursor0 or libxcb-cursor0 is needed to load the Qt xcb platform plugin.
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.

Available platform plugins are: eglfs, offscreen, wayland-egl, linuxfb, wayland, minimal, xcb, vkkhrdisplay, minimalegl, vnc.
```

この場合は以下で `libxcb-cursor0` をインストールしてください:

`sudo apt install libxcb-cursor0`

### Wayland 環境でウィンドウの影が表示されない

Wayland のデスクトップ環境では、通常のフレーム付きウィンドウとして
`zstarview --window-frame window` を起動しても、外側の影が表示されないことがあります。
これは zstarview 自身のウィンドウ設定というより、Wayland の window decoration /
compositor 側の挙動によることが多いです。

X11 系の見た目でウィンドウ影を出したい場合は、実用的な回避策として
XWayland 経由で起動できます:

```sh
QT_QPA_PLATFORM=xcb zstarview --window-frame window
```

`QT_QPA_PLATFORM=xcb` を付けると影が表示される場合は、その差分は
Wayland と X11 のウィンドウ装飾経路の違いによるものと考えてよいです。

### ネットワークが遅い / オフラインで使いたい
1. 惑星暦データ

   初回起動時は惑星暦データ（`de442s.bsp`）を自動ダウンロードします。  
   この一度だけはネットワーク接続が必要です。ダウンロード後はキャッシュを利用してオフラインでも動作します。

2. 雲衛星画像

   雲の描画は公開 S3 バケットから衛星画像を取得し、依存も比較的重めです。
   回線が細い、またはオフラインの場合は `-c 0` で雲描画を無効化してください。
   雲を無効化しても、恒星・惑星・空の色の表示は利用できます。

3. 地形地平線

   地形地平線表示は初回に Copernicus DEM タイルをダウンロードし、その後はローカルキャッシュを再利用します。
   回線が細い、またはオフラインの場合は `--terrain-horizon-opacity 0` で地形地平線表示を無効化してください。
   地形地平線を無効化しても、恒星・惑星・空の色の表示は利用できます。

4. 人工衛星データ

   人工衛星オーバーレイは実行時に ISS の軌道要素データを取得し、取得元は `wheretheiss.at` を優先し、失敗時だけ CelesTrak を使います。fresh な current cache は最大 24 時間まで再利用します。
   このレイヤーはリアルタイム表示でのみ利用でき、タイムシフト表示では人工衛星の取得も描画も行いません。
   回線が細い、またはオフラインの場合は `--satellite-opacity 0` で人工衛星レイヤーを無効化してください。
   新しいキャッシュがすでにあれば、ネットワークがなくても人工衛星オーバーレイを表示し続けられます。

5. 航空機データ

   航空機オーバーレイは実行時に OpenSky Network の state data を取得します。
   既定では 5 分ごとに再取得します。この間隔は、無料枠での利用や一時的な取得失敗・再試行に対して余裕を持たせるため、過度に短くせず保守的に設定しています。
   OpenSky への問い合わせ自体を避けたい場合は、`-a 0` で航空機レイヤーを無効化してください。

雲関連のステータス表示は `idle` / `downloading` / `partial` を使います:
- `downloading`: S3 から衛星ソースを取得中
- `partial`: 入手済みデータのみで描画（欠損領域は薄い黄色で表示）

GOES-East の参照先を GOES-19 に更新したことで、以前は単に「衛星のカバー外」となっていた場所でも、場所によっては部分カバーとして表示されるようになりました。ヨーロッパでも一部の地点でこの挙動になります。その場合、カバーされている部分には雲画像が描かれ、カバー外の部分には薄い黄色の欠損色が表示されます。77〜87% 程度のカバー率の例は [screenshot5](docs/images/screenshot5.png) を参照してください。

### 星空の更新間隔と CPU 負荷

CPU 性能によっては星空の自動更新が負荷になる場合があります。更新間隔を長くして負荷を下げてください（例: `-i 300` で 5 分ごと）。余裕があれば短くして構いません。

### ログの確認

ターミナルから `$ zstarview` で起動すると、起動メッセージやエラーを確認できます。
併せてログファイルにも出力されます（OS 依存）。例:
- Linux: `~/.cache/zstarview/logs/app.log`
- macOS: `~/Library/Logs/zstarview/app.log`
- Windows: `%LOCALAPPDATA%/tos-kamiya/zstarview/Logs/app.log`

通常起動だとすぐ閉じて確認しづらい場合は、ターミナルから `zstarview-debug` を試してください。
これは主に Windows でのトラブルシュート用です。Linux では `zstarview-debug` は `zstarview` と実質同じ動作です。

Windows では、Windows セキュリティにより Python 拡張モジュールの読み込みがブロックされ、起動時に止まることがあります。
その場合は、Windows セキュリティの `App & browser control` にある `Smart App Control` の設定を変更すると回避できることがあるようです。
ただし、セキュリティを弱くすることになるので、安全な環境以外での実行は推奨しません。
[Smart App Control の画面例](docs/images/windows-smart-app-control_ja.png) も参照してください。

</details>

<details>
  <summary>コード・データのライセンスとクレジット</summary>

## コード・データのライセンスとクレジット

このソフトウェアは [MIT](LICENSE.txt) の下で提供されています。

ただし、**同梱されているデータ** はそれぞれのライセンスに従って再配布されます。

以下のパスは `src/zstarview/data/` 配下を基準としています。

| ファイル | 内容 | 出典 | ライセンス |
| :--- | :--- | :--- | :--- |
| `cities1000.txt`, `admin1CodesASCII.txt` | 人口 1000 人以上の都市一覧 | [GeoNames](https://download.geonames.org/export/dump/) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| `viewpoints/tower_viewpoints.json` | タワー名起動用に同梱している展望塔/タワーデータ（Wikidata 由来の整形データ） | [Wikidata](https://www.wikidata.org/) をローカル整形したもの（手順は `dev-samples/` に記録） | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/)（Wikidata データ） |
| `viewpoints/mountain_viewpoints.json` | 山名起動用に同梱している山頂ビューポイントデータ（Wikipedia で収集した候補を Wikidata メタデータで正規化したデータ） | [Wikipedia](https://www.wikipedia.org/) での候補収集と [Wikidata](https://www.wikidata.org/) による正規化手順（`dev-samples/` に記録） | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/)（Wikidata データ） |
| 実行時に OpenStreetMap Nominatim へ送る `--place` ジオコーディング要求 | `--place` 指定時だけ使うオンライン地名検索 | [OpenStreetMap Nominatim](https://nominatim.openstreetmap.org/) | [Nominatim Usage Policy](https://operations.osmfoundation.org/policies/nominatim/) |
| アプリのキャッシュディレクトリ配下にオンデマンドで保存される都市アウトラインキャッシュ | ダウンロードした Overture 建物データから生成した派生建物タイルと `tile_index.json` | `overturemaps` CLI を通じて実行時に取得する [Overture Maps Buildings](https://docs.overturemaps.org/guides/buildings/) | [ODbL 1.0](https://opendatacommons.org/licenses/odbl/) |
| 実行時に `wheretheiss.at` から取得し、失敗時は CelesTrak を使う人工衛星オーバーレイ用データ | ISS 表示に使う軌道要素データ | [wheretheiss.at](https://wheretheiss.at/w/developer), [CelesTrak](https://celestrak.org/) | 利用条件やライセンスは各出典サイトを参照 |
| `dso.csv` | DSO（銀河/散開星団/球状星団）カタログ（OpenNGC 由来の生成データ） | [OpenNGC](https://github.com/mattiaverga/OpenNGC)（[PyOngc](https://github.com/mattiaverga/PyOngc) 経由で生成） | [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)（OpenNGC データベース） |
| アプリのキャッシュディレクトリ配下にオンデマンドで保存される地形 DEM キャッシュ | 地形地平線用の地形データ（Copernicus DEM GLO-90） | [Copernicus DEM / Copernicus Data Space Ecosystem](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM)（アプリは公開 AWS 配布を利用） | Copernicus Data Space Ecosystem の案内する Copernicus DEM GLO-90 の利用条件（"Licence for COP-DEM-GLO-90-F Global 90m Full, Free & Open" / "Licence for the use of the Copernicus WorldDEM™-90"） |
| `stars/IAU-Catalog of Star Names (always up to date).csv` | IAU 恒星名作業部会 (WGSN) による恒星固有名カタログ | [exopla.net](https://exopla.net/star-names/modern-iau-star-names/) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| `Noto_Sans/*` | テキスト表示フォント | [Google Fonts](https://fonts.google.com/) | [SIL Open Font License 1.1](https://openfontlicense.org) |

## クレジット

* 天文データを提供していただいている CDS Strasbourg および ESA Hipparcos Mission に感謝します。
* 都市データは GeoNames に基づいています。
* タワー/展望塔の起動データは Wikidata に基づく整形データであり、Wikidata の CC0 条件に従って再配布しています。
* 山頂ビューポイントの起動データは Wikipedia で収集した候補を Wikidata メタデータで正規化したものであり、ここでは Wikidata の CC0 条件に従って再配布しています。
* 都市アウトライン用の元データは **Overture Maps Buildings** から必要時に取得し、実行時利用向けに派生タイルへ変換したものです。
* 恒星の固有名は IAU 恒星名作業部会 (WGSN) による承認済みリスト（[exopla.net](https://exopla.net/star-names/modern-iau-star-names/) 経由）を使用しています。
* 雲データは気象衛星 **Himawari**（提供: JMA）および **NOAA GOES** シリーズ（提供: NOAA/NESDIS）による赤外線観測データを、それぞれの公開 S3 バケットから取得して利用しています。
* 人工衛星オーバーレイで使う軌道要素データ（TLE/OMM）は **wheretheiss.at** を優先し、失敗時は **CelesTrak** を fallback として利用します。
* `--place` による地名・駅名検索は公開の OpenStreetMap Nominatim サービスを使っており、[Nominatim Usage Policy](https://operations.osmfoundation.org/policies/nominatim/) の対象です。
* 地形地平線データは **Copernicus DEM GLO-90** に基づいており、欧州委員会のために ESA が管理するデータを、アプリでは公開 AWS 配布とローカルキャッシュを通じて利用しています。
* 大規模建物データを公開している Overture Maps とそのソースデータ提供者に感謝します。
* 雲画像や地形 DEM の取得に利用している公開 S3 配布/ミラーを提供している AWS および各データ提供者に感謝します。
* フォントは Google Noto Project を利用しています。
* ウィンドウタイトル「Zenith Star View」は ChatGPT の提案に由来します。
* Gemini および ChatGPT に、仕様の相談、コード生成、デバッグなど、多くの助力をいただきました。

</details>

## 付録

→ [開発者向けドキュメント](docs/developer/README.md)

→ [仕様書](docs/specification.md), [設計書](docs/design.md)

→ [2026〜2028年の月食（皆既・部分）, 2026〜2028年の皆既日食](docs/appendix-eclipses-ja_JP.md)
