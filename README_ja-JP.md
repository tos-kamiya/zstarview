# zstarview 🌌

雲があっても、太陽が出ていても、満天の星空を。

**Zenith Star View** は、地球上の任意の都市から見える星空を表示するアプリケーションです。
名前の *zenith*（天頂）は、観測者の真上の一点を意味し、その場に立って夜空を見上げる感覚を表しています。

**特徴:**

- 明るい恒星、惑星、天の赤道、黄道をリアルタイムで描画
- 名前付きのDSO（銀河/散開星団/球状星団）を薄い青系の領域として表示。DSOホバーは恒星ホバーと独立して動作
- 星座（アステリウム）表示に対応。暗い線で常時表示し、構成星にホバーすると該当アステリズムを明るく強調してラベルを表示。同一恒星に複数ある場合は3秒ごとに切り替え
- 太陽・月・主要惑星に対応。小惑星（アステロイド）は未対応です。
- 都市名で場所を指定可能（GeoNames に基づく）
- 表示中心を `-A`（高度）/`-Z`（方位）で調整可能
- リアルタイムの衛星雲画像（Himawari/GOES）を縞模様（ハッチ）のオーバーレイとして重ねて表示します。
- 雲データ取得と雲レイヤー生成は分離されており、カメラ（方位/仰角）変更時はキャッシュ済みソースから即座に再描画します。
- 衛星データが部分的な場合、欠損領域は薄い黄色ティントで表示します。
- Copernicus DEM を使った地形地平線オーバーレイを利用でき、観測地の地形に沿った控えめな黄土色がかった地形線を表示できます。
- 幾何学的地平線または地形地平線より下は、向きの把握を助けるため地面色で塗り分けます。
- 視線変更中やウィンドウリサイズ中は、一時的にビューポート操作用の簡易描画モードへ切り替わり、重いレイヤーの更新を待たずに明るい星・ガイドライン・地形地平線を素早く追随表示します。
- 観測地の緯度に対して昇らない天球領域（never-rises）を赤いティントで表示します。
- CPython 3.10, 3.11, 3.12, 3.13 で継続的にテストしています。

  ![](docs/images/screenshot1.png)


## インストール方法（推奨：`pipx`）

[`pipx`](https://pypa.github.io/pipx/) を使ってインストールする想定です。

```bash
pipx install zstarview
```

アップグレード:

```bash
pipx upgrade zstarview
```

> 注記: X11 ライブラリやネットワークが細い場合の回避策などは「トラブルシューティング」を参照してください。

## 使い方

```bash
zstarview [options] [location]
```

> 注記（Ubuntu/Wayland, GNOME）: ターミナル起動時にタスクバーのアイコンが表示されない場合は、後述の [`.desktop` ランチャーの生成（GNOME専用）](#desktop-ランチャーの生成gnome専用) を実行してください。

### GUI

#### キー操作

* **← / →**: 視線の方位を ±5° 回転
* **↑ / ↓**: 視線の高度を ±5° 変更（0°..90°にクランプ）
  方向キー入力が続く間と最後の入力から約0.7秒の間は、ビューポート操作用の簡易描画モードになります。この間は `Vmag <= 4.0` の恒星、天の赤道、黄道、地平線、地形地平線、方位ラベル、天頂マーカーのみを表示し、惑星、全星等の星空、空ディスク、雲、DSO、アステリウムは一時的に非表示になります。
* **M**: 月の5倍表示をトグル
* **D**: DSOオーバーレイの表示/非表示を切り替え
* **A**: アステリウムオーバーレイの表示/非表示を切り替え
* **S**: 空ディスク表示をグラデーションとフラットディスクで切り替え
* **C**: 雲オーバーレイの表示/非表示を切り替え
* **T**: 地形地平線オーバーレイの表示/非表示を切り替え
* **Ctrl+J**: Jump to Named Star を開く
* **Ctrl+F**: Search Stars and Asterisms を開く
* **F11**: フルスクリーン表示の切り替え
* **ESC**: フルスクリーンから復帰
* **Q**: 終了

#### メニュー操作

ハンバーガーメニュー（`☰`）から次を利用できます。

* **Jump to Named Star...**: 代表的な固有名星（`Vmag <= 2.0`）を北天 / 赤道付近 / 南天で選んで、視点中心をその星へ移動します。
* **Search Stars and Asterisms...**: 固有名付き恒星と対応アステリウムを横断検索し、選択した対象へ移動します。
* **Enlarge Moon**: 月の5倍表示を切り替えます。
* **DSO**: DSOオーバーレイの表示/非表示を切り替えます。
* **Asterisms**: 星座（アステリウム）オーバーレイの表示/非表示を切り替えます（有効時は暗い線を常時表示し、構成星にホバーすると該当アステリズムを明るく強調してラベルを表示します）。
* **Sky Color Disc**: 空ディスク表示を、空色グラデーション表示とフラットな暗色ディスク表示で切り替えます。
* **Clouds**: リアルタイム雲オーバーレイの表示/非表示を切り替えます。
* **Terrain Horizon**: 地形地平線オーバーレイの表示/非表示を切り替えます。CLI で `--terrain-horizon-opacity 0` を指定して起動した場合、その起動中はメニューから再有効化できません。
* **Fullscreen**: フルスクリーン表示を切り替えます。
* **Exit**: アプリケーションを終了します。

ジャンプ/検索の確定後は約3秒間、マウスホバー時と同じ見た目（円マーカー + 名称ラベル）で対象星を強調表示します。

ウィンドウリサイズ中も同じ簡易描画モードを使うため、サイズ変更中の追随性を保ちます。

### CLI

#### 引数

| 引数     | 説明                                                          | デフォルト              |
| :----- | :---------------------------------------------------------- | :----------------- |
| `location` | 表示する都市名、タワー名、山名、明示指定の `t/NAME`・`m/NAME`、または `"<lat>;<lon>"` 形式の緯度経度（十進度）を指定できます（例: `Tokyo`, `Tokyo Skytree`, `t/Tokyo Skytree`, `Mount Fuji`, `m/Mount Fuji`, `35.68;139.76`, `N35.68;E139.76`, `-35.68;139.76`）。省略時は前回起動時の location を使用します。初回起動時に省略すると `Tokyo` になります。 | 前回の location（初回は `Tokyo`） |

#### オプション

| オプション                                       | 説明                                                            | デフォルト |
| :------------------------------------------ | :-------------------------------------------------------------- | :------- |
| `-h`, `--help`                              | ヘルプメッセージを表示して終了します。                                 |          |
| `-Z`, `--view-center-az VIEW_CENTER_AZ`     | 表示中心の方位角を指定します。                                     | `180`    |
| `-A`, `--view-center-alt VIEW_CENTER_ALT`   | 表示中心の高度角を指定します（90=天頂、0=地平線）。                       | `90`     |
| `--observer-height-m METERS` | 観測者の地面からの高さをメートルで指定します。デフォルトは都市/緯度経度/山名入力で `1.7`、タワー名入力ではタワー高さです。 | location依存 |
| `-c`, `--cloud-opacity CLOUD_OPACITY`                 | 雲の不透明度を指定します（0.0〜1.0）。0.0で描画を無効化します。※2 | `0.2`   |
| `--cloud-missing-tint-opacity OPACITY` | 雲欠損領域の黄色ティント不透明度を指定します（0.0〜1.0）。 | `0.176` |
| `--sky-opacity SKY_OPACITY`                 | 空の色ディスクの不透明度を指定します（0.0〜1.0）。0.0で描画を無効化します。 | `0.2`   |
| `--terrain-horizon-opacity OPACITY` | 地形地平線ポリラインの不透明度を指定します（0.0〜1.0）。0.0で DEM ダウンロード・地形地平線計算・描画を無効化します。※4 | `0.05` |
| `--ground-tint-opacity OPACITY` | 幾何学的地平線または地形地平線より下の地面色塗りの強さを指定します（0.0〜1.0）。 | `0.1` |
| `-m`, `--enlarge-moon`                      | 月を5倍に拡大して表示します。                                      |          |
| `-s`, `--star-base-radius STAR_BASE_RADIUS` | 2等星の基本サイズを指定します。                                   | `4.0`    |
| `-w`, `--expected-render-width EXPECTED_RENDER_WIDTH` | 恒星をフル解像度で描画する想定ウィンドウ幅を指定します。天球幅がこの値を超える場合、恒星レイヤーは平方根スケーリングで描画します。 | `600` |
| `--window-geometry restore\|X,Y,W,H` | 初期ウィンドウ位置と大きさを指定します。`restore` で前回終了時の位置/サイズを復元し、`X,Y,W,H` で整数値を直接指定できます。注: Wayland ではウィンドウ位置の復元は利用できません（サイズ復元は有効です）。 |          |
| `-V`, `--vmag-limit V_MAG_LIMIT`            | 表示する恒星の等級（明るさ）の上限を指定します。                          | `6.0`    |
| `--vmag-brightness-multiplier MULTIPLIER`   | 等級1段階あたりの光量変化倍率（`1.58`〜`2.512`、デフォルト：`2.5`。Pogson の定義は `2.512` です）。 ※3 | `2.5`   |
| `-i`, `--sky-update-interval SKY_UPDATE_INTERVAL` | 星空を更新する時間間隔（秒） を指定します。 | `60` |
| `--show-dso-initial true\|false` | 起動時にDSOを表示するかを指定します。 | 自動（カタログがあれば表示） |
| `--show-asterisms-initial true\|false` | 起動時に星座（アステリウム）を表示するかを指定します。 | 表示 |
| `-t`, `--theme {night,day,white,black}` | 背景と星の見え方のテーマを指定します。 | `night` |
| `-H`, `--hours HOURS`                       | 現在時刻に加算する時間数を指定します。※1                                 | `0`      |
| `-D`, `--days DAYS`                         | 現在時刻に加算する日数を指定します。※1                                  | `0`      |
| `--datetime "YYYY-MM-DD HH[:MM[:SS]] [TZ]"` | 絶対的な日時を指定します。時刻は「時」「時:分」「時:分:秒」のいずれかで指定できます。タイムゾーン（TZ）を省略した場合はUTC。※1 |          |

※1 これらのオプションを指定してリアルタイムではない星空を表示した場合には、雲は描かれません。

※2 雲の描画は気象衛星（Himawari, GOES）の赤外線データを S3 バケットから取得して行います。ネットワーク関連の注意や回避策は「トラブルシューティング」を参照してください。

※3 Pogson の定義では 5 等級差 = 100 倍なので、等級1段階あたりの光量倍数は \(100^{1/5} \approx 2.512\) です。`--vmag-brightness-multiplier` の上限はこの値です。

※4 地形地平線表示は初回利用時に Copernicus DEM タイルをダウンロードし、以後はローカルキャッシュを再利用します。有効時は、ディスク内の地面/空の塗り分け境界にも地形プロファイルを使います。

#### 起動時のオーバーレイ表示設定

起動時の初期表示をメニュー操作なしで切り替えるには、次を使います。

```bash
# DSOは非表示、アステリウムは表示で起動
zstarview --show-dso-initial false --show-asterisms-initial true Tokyo
```

#### 表示中心オプションについて

`-Z`（方位角）と `-A`（高度角）のオプションで、画面の表示中心を指定できます。

デフォルトでは `-Z 180`（南向き）、`-A 90`（天頂）です。画面下が南、画面左が東で、天頂を見上げたような円形の表示になります。

例えば、`-Z 90`（東向き）および `-A 25`（高度25°、地平線から25°上）を指定すると、東の空を見上げる視野が生成されます。  
下のスクリーンショットは、この設定で表示した夏の大三角（ベガ、アルタイル、デネブ）を示しています。

[→ 例：高度25°の東の空に表示された夏の大三角](docs/images/screenshot2.png)

方位角は度数または16方位の方位記号（大小区別なし）で指定できます。
例: `-Z E`, `-Z ne`, `-Z SSW`（= 202.5°）。
（対応表: 0=北, 90=東, 180=南, 270=西）

#### 等級上限オプションについて

`-V 等級` で指定した等級までの明るさの星を描画します。
デフォルトは `-V 6.0` です。例えば 10.0 等級を指定すると、約32万4千個の星が描画されます。
この値を大きくすると処理が重くなる点に注意してください。

[→ 例： 10.0等級まで表示](docs/images/screenshot3.png)

#### 日時指定オプションについて

`--datetime "YYYY-MM-DD HH[:MM[:SS]] [TZ]"` で絶対的な日時を指定できます。
時刻部分は「時」だけ、「時:分」、「時:分:秒」のいずれも使用可能です。
タイムゾーン（TZ）を省略した場合は UTC として扱われます。

タイムゾーンは以下のいずれかの形式で指定できます：

- よく使われる略称（JST, UTC, GMT, KST, HKT, AWST, ACST, AEST, NZST, NZDT, MSK, EAT）
- IANA タイムゾーン名（例: `Asia/Tokyo`, `Europe/Moscow`）
- UTCオフセット（例: `UTC+9`, `UTC-07:30`）

例:

```bash
zstarview --datetime "2025-08-17 21:00:00 JST" Tokyo
zstarview --datetime "2025-09-12 9" Tokyo         # 9時ちょうど
zstarview --datetime "2025-09-12 09:00" Tokyo     # 9時0分
zstarview --datetime "2025-09-12 9:0:0 JST" Tokyo # JSTの9時
```

#### 緯度経度の直接指定

都市名の代わりに、`"<lat>;<lon>"` 形式で緯度経度（十進度）を指定できます。

- フォーマット: `緯度;経度`（セミコロン区切り）
- 許容例:
  - `35.68;139.76`
  - `N35.68;E139.76`
  - `-35.68;139.76`
  - `S35.68;W139.76`
- 緯度は `-90〜90`、経度は `-180〜180` の範囲でバリデーションされます。
- 方向記号 `N/S/E/W` を付けた場合は符号として解釈します（数値にマイナスが明示されている場合はその符号を優先）。
- 座標で起動した場合、**タイムゾーンはデフォルトで UTC** として扱われます（`--datetime` にタイムゾーンを明示すれば上書き可能）。

例:

```bash
zstarview "35.68;139.76"
zstarview "N35.68;E139.76" --datetime "2025-09-12 21 JST"
```

#### タワー名入力

Wikidata 由来の同梱タワー/展望地点 dataset から起動できます。

- 例:
  - `Tokyo Skytree`
  - `t/Tokyo Skytree`（タワー明示指定）
  - `Tsutenkaku`（`Tsūtenkaku` の ASCII フォールバック）
  - `Tokyo Tower`
  - `wikidata:Q57965`
- タワー名を使った場合、観測者高さのデフォルトはそのタワーの登録高です。
- `--observer-height-m` を使えば上書きできます。
- ダイアクリティカルマーク付き名称については、ASCII フォールバック表記でも解決できます。

例:

```bash
zstarview "Tokyo Skytree"
zstarview "Tokyo Tower" --observer-height-m 150
```

#### 山名入力

同梱山頂ビューポイント dataset から起動できます。

- 例:
  - `Mount Fuji`
  - `m/Mount Fuji`（山明示指定）
  - `Aconcagua`
  - `Snezka`（`Sněžka` の ASCII フォールバック）
  - `wikidata:Q39231`
- 山名を使った場合、観測者高さのデフォルトは `1.7m` です。
- ダイアクリティカルマーク付き名称については、ASCII フォールバック表記でも解決できます。

例:

```bash
zstarview "Mount Fuji"
zstarview "Snezka"
```

### ビューポイント dataset CLI 参照オプション

GUI を起動せずに、同梱タワー/展望地点 dataset と山頂ビューポイント dataset を参照できます。

| オプション                                       | 説明                                                            | デフォルト |
| :------------------------------------------ | :-------------------------------------------------------------- | :------- |
| `-h`, `--help`                              | ヘルプメッセージを表示して終了します。                                 |          |
| `--list-towers` | 同梱タワー/展望地点 dataset の一覧表示名を出力して終了します。利用可能な場合は ASCII フォールバック名を優先します。 | |
| `--list-tower-names` | 同梱タワー/展望地点 dataset の名前を、多言語名と ASCII フォールバック名込みで一覧出力して終了します。 | |
| `--show-tower-json NAME` | 指定名で同梱タワー/展望地点を解決し、利用可能な場合は `ascii_name` を含む JSON メタデータを出力して終了します。 | |
| `--list-mountains` | 同梱山頂ビューポイント dataset の一覧表示名を出力して終了します。利用可能な場合は ASCII フォールバック名を優先します。 | |
| `--list-mountain-names` | 同梱山頂ビューポイント dataset の名前を、多言語名と ASCII フォールバック名込みで一覧出力して終了します。 | |
| `--show-mountain-json NAME` | 指定名で同梱山頂ビューポイントを解決し、利用可能な場合は `ascii_name` を含む JSON メタデータを出力して終了します。 | |

```bash
zstarview --list-towers
zstarview --list-tower-names
zstarview --show-tower-json "Tokyo Skytree"
zstarview --list-mountains
zstarview --show-mountain-json "Mount Fuji"
```

これらのオプションは相互排他で、`location` 引数や時刻・描画オプションとは併用できません。
`--list-towers` では、利用可能な場合は ASCII フォールバック名を優先表示します。
`--list-tower-names` では、元の綴りと ASCII フォールバック綴りの両方を含みます。
`--list-mountains` では、利用可能な場合は ASCII フォールバック名を優先表示します。
`--list-mountain-names` では、元の綴りと ASCII フォールバック綴りの両方を含みます。

### 対応している星座（アステリウム）

ここでの「星座（アステリウム）」は、見かけ上の星の並びを結んだ通称パターンです。  
**IAU（国際天文学連合）が定義する正式な88星座の境界とは別概念**です。

- 冬: `Winter Triangle`（冬の大三角）, `Orion's Belt`（オリオンの三ツ星）, `Winter Hexagon`（冬のダイヤモンド）, `Southern Cross`（南十字）, `Southern Pointers`, `Diamond Cross`, `False Cross`
- 春: `Big Dipper`（北斗七星）, `Little Dipper`, `Spring Triangle`, `Arc to Arcturus`, `Leo Sickle`, `Southern Triangle`
- 夏: `Summer Triangle`（夏の大三角）, `Northern Cross`（北十字）, `Teapot`, `Keystone`
- 秋: `Great Square of Pegasus`（ペガススの四辺形）, `Circlet of Pisces`, `Water Jar of Aquarius`, `Cassiopeia W`, `House of Cepheus`, `Job's Coffin`

## `.desktop` ランチャーの生成（GNOME専用）

GNOME系デスクトップ環境（Ubuntu Dock や DockToPanel を含む）では、
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

## トラブルシューティング

### X11（Ubuntu/Debian）

Qt の xcb プラグインが `libxcb-cursor0` を必要とする場合があります。
X11/Wayland を意識していないと分かりづらいですが、ターミナル（コンソール）から実行すると次のようなエラーが表示されます：

```sh
$ zstarview
qt.qpa.plugin: From 6.5.0, xcb-cursor0 or libxcb-cursor0 is needed to load the Qt xcb platform plugin.
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.

Available platform plugins are: eglfs, offscreen, wayland-egl, linuxfb, wayland, minimal, xcb, vkkhrdisplay, minimalegl, vnc.
```

この場合は以下で `libxcb-cursor0` をインストールしてください：

`sudo apt install libxcb-cursor0`

### ネットワークが遅い/オフラインで使いたい
雲の描画は S3 から衛星画像を取得します。回線が細い/オフラインのときは `-c 0` で雲描画を無効化してください。
地形地平線表示は初回に Copernicus DEM タイルをダウンロードし、その後はローカルキャッシュを再利用します。DEM ダウンロードを避けたい場合は `--terrain-horizon-opacity 0` を指定してください。
雲や地形地平線を無効化しても恒星・惑星・空の色の表示は利用できます。

雲ステータス表示は `idle` / `downloading` / `partial` を使います。
- `downloading`: S3 から衛星ソースを取得中
- `partial`: 入手済みデータのみで描画（欠損領域は薄い黄色ティント）

> 注記: 初回起動時は惑星暦データ（`de440s.bsp`）を自動ダウンロードします。
> この一度だけはネットワーク接続が必要です。ダウンロード後はキャッシュを利用してオフラインでも動作します（特に雲を無効化した場合）。

### 星空の更新間隔とCPU負荷
CPU性能によっては星空の自動更新が負荷になる場合があります。様子を見ながら更新間隔を長く（数値を大きく）して負荷を下げてください（例: `-i 300` で5分ごと）。余裕があれば短くして構いません。

### ログの確認
ターミナル（端末）から `$ zstarview` により起動すると、起動メッセージやエラーを確認できます。
併せてログファイルにも出力されます（OSに依存）。例:
- Linux: `~/.cache/zstarview/logs/app.log`
- macOS: `~/Library/Logs/zstarview/app.log`
- Windows: `%LOCALAPPDATA%/tos-kamiya/zstarview/Logs/app.log`

## 星カタログ再生成（開発者向け）

カタログ生成スクリプトは以下で実行できます。

```bash
uv run -p .venv/bin/python src/zstarview/data/stars/generate_star_catalog.py
```

Tycho-2 入力や分割出力を含む詳細オプションは次を参照してください。

- `docs/developer/star-catalog-generation.md`

## ライセンス

このソフトウェアは [MIT](LICENSE.txt) の下で提供されています。

ただし、 **同梱されているデータ** はそれぞれのライセンスに従って再配布されます。

以下のパスは `src/zstarview/data/` 配下を基準としています。

| ファイル                                         | 内容                                               | 出典                                                                       | ライセンス                                                                                                                      |
| ------------------------------------------------ | -------------------------------------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `cities1000.txt`, `admin1CodesASCII.txt` | 人口1000人以上の都市一覧                           | [GeoNames](https://download.geonames.org/export/dump/)                     | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                       |
| `viewpoints/tower_viewpoints.json` | タワー名起動用に同梱している展望塔/タワー dataset（Wikidata 由来の整形データ） | [Wikidata](https://www.wikidata.org/) をローカル整形したもの（手順は `dev-samples/` に記録） | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/)（Wikidata データ） |
| `viewpoints/mountain_viewpoints.json` | 山名起動用に同梱している山頂ビューポイント dataset（Wikipedia で収集した候補を Wikidata メタデータで正規化したデータ） | [Wikipedia](https://www.wikipedia.org/) での候補収集と [Wikidata](https://www.wikidata.org/) による正規化手順（`dev-samples/` に記録） | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/)（Wikidata データ） |
| `dso.csv`                                                   | DSO（銀河/散開星団/球状星団）カタログ（OpenNGC 由来の生成データ） | [OpenNGC](https://github.com/mattiaverga/OpenNGC)（[PyOngc](https://github.com/mattiaverga/PyOngc) 経由で生成） | [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)（OpenNGC データベース） |
| アプリのキャッシュディレクトリ配下にオンデマンドで保存される地形 DEM キャッシュ | 地形地平線用の地形データ（Copernicus DEM GLO-90） | [Copernicus DEM / Copernicus Data Space Ecosystem](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM)（アプリは公開 AWS 配布を利用） | Copernicus DEM 向け ESA User Licence（Copernicus Contributing Mission data access terms） |
| `stars/IAU-Catalog of Star Names (always up to date).csv`             | IAU 恒星名作業部会 (WGSN) による恒星固有名カタログ | [exopla.net](https://exopla.net/star-names/modern-iau-star-names/)         | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                       |
| `Noto_Sans/*`     | テキスト表示フォント                             | [Google Fonts](https://fonts.google.com/)   | [SIL Open Font License 1.1](https://openfontlicense.org)                                                                        |

## クレジット

* 天文データを提供していただいている CDS Strasbourg および ESA Hipparcos Mission に感謝します。
* 都市データは GeoNames に基づいています。
* タワー/展望塔の起動データは Wikidata に基づく整形データであり、Wikidata の CC0 条件に従って再配布しています。
* 山頂ビューポイントの起動データは Wikipedia で収集した候補を Wikidata メタデータで正規化したものであり、ここでは Wikidata の CC0 条件に従って再配布しています。
* 恒星の固有名は IAU 恒星名作業部会 (WGSN) による承認済みリスト（[exopla.net](https://exopla.net/star-names/modern-iau-star-names/) 経由）を使用しています。
* 雲データは気象衛星 **Himawari**（提供: JMA）および **NOAA GOES** シリーズ（提供: NOAA/NESDIS）による赤外線観測データを、それぞれの S3 公開バケットから取得して利用しています。
* 地形地平線データは **Copernicus DEM GLO-90** に基づいており、欧州委員会のために ESA が管理するデータを、アプリでは公開 AWS 配布とローカルキャッシュを通じて利用しています。
* 雲画像や地形 DEM の取得に利用している公開 S3 配布/ミラーを提供している AWS および各データ提供者に感謝します。
* フォントは Google Noto Project を利用しています。
* ウィンドウタイトル「Zenith Star View」は ChatGPT の提案に由来します。
* Gemini および ChatGPT に、仕様の相談、コード生成、デバッグなど、多くの助力をいただきました。

## 付録

→ [開発者向けドキュメント](docs/developer/README.md)

→ [仕様書](docs/specification.md), [設計書](docs/design.md)

→ [2026〜2028年の月食（皆既・部分）, 2026〜2028年の皆既日食](docs/appendix-eclipses-ja_JP.md)
