# zstarview 🌌

雲があっても、太陽が出ていても、満天の星空を。

**Zenith Star View** は、地球上の任意の都市から見える星空を表示するアプリケーションです。
名前の *zenith*（天頂）は、観測者の真上の一点を意味し、その場に立って夜空を見上げる感覚を表しています。

**特徴:**

- 明るい恒星、惑星、天の赤道、黄道をリアルタイムで描画
- 太陽・月・主要惑星に対応。小惑星（アステロイド）は未対応です。
- 都市名で場所を指定可能（GeoNames に基づく）

  ![](docs/images/screenshot1.png)

- 表示中心を `-A`（高度）/`-Z`（方位）で調整可能
- リアルタイムの衛星雲画像（Himawari/GOES）を縞模様（ハッチ）のオーバーレイとして重ねて表示します。

  ![](docs/images/screenshot4.png)


## インストール方法（推奨：`pipx`）

[`pipx`](https://pypa.github.io/pipx/) を使ってインストールする想定です。

```bash
pipx install git+https://github.com/tos-kamiya/zstarview.git
```

> 注記: X11 ライブラリやネットワークが細い場合の回避策などは「トラブルシューティング」を参照してください。

## 使い方

```bash
zstarview [options] [city]
```

### 引数

| 引数     | 説明                                                          | デフォルト              |
| :----- | :---------------------------------------------------------- | :----------------- |
| `city` | 表示する都市名を指定します。**または** `"<lat>;<lon>"` 形式で緯度経度（十進度）を直接指定できます（例: `35.68;139.76`, `N35.68;E139.76`, `-35.68;139.76`）。省略時は前回起動時の都市 **または直前の座標** を使用します。初回起動時に省略すると `Tokyo` になります。 | 前回の都市/座標（初回は `Tokyo`） |

### オプション

| オプション                                       | 説明                                                            | デフォルト |
| :------------------------------------------ | :-------------------------------------------------------------- | :------- |
| `-h`, `--help`                              | ヘルプメッセージを表示して終了します。                                 |          |
| `-Z`, `--view-center-az VIEW_CENTER_AZ`     | 表示中心の方位角を指定します。                                     | `180`    |
| `-A`, `--view-center-alt VIEW_CENTER_ALT`   | 表示中心の高度角を指定します（90=天頂、0=地平線）。                       | `90`     |
| `-c`, `--cloud-opacity CLOUD_OPACITY`                 | 雲の不透明度を指定します（0.0〜1.0）。0.0で描画を無効化します。※2 | `0.2`   |
| `--sky-opacity SKY_OPACITY`                 | 空の色ディスクの不透明度を指定します（0.0〜1.0）。0.0で描画を無効化します。 | `0.2`   |
| `-m`, `--enlarge-moon`                      | 月を5倍に拡大して表示します。                                      |          |
| `-s`, `--star-base-radius STAR_BASE_RADIUS` | 星の基本サイズを指定します。                                       | `8.0`    |
| `-V`, `--vmag-limit V_MAG_LIMIT`            | 表示する恒星の等級（明るさ）の上限を指定します。                          | `6.0`    |
| `-i`, `--sky-update-interval SKY_UPDATE_INTERVAL` | 星空を更新する時間間隔（秒） を指定します。 | `180` |
| `-H`, `--hours HOURS`                       | 現在時刻に加算する時間数を指定します。※1                                 | `0`      |
| `-D`, `--days DAYS`                         | 現在時刻に加算する日数を指定します。※1                                  | `0`      |
| `--datetime "YYYY-MM-DD HH[:MM[:SS]] [TZ]"` | 絶対的な日時を指定します。時刻は「時」「時:分」「時:分:秒」のいずれかで指定できます。タイムゾーン（TZ）を省略した場合はUTC。※1 |          |

※1 これらのオプションを指定してリアルタイムではない星空を表示した場合には、雲は描かれません。

※2 雲の描画は気象衛星（Himawari, GOES）の赤外線データを S3 バケットから取得して行います。ネットワーク関連の注意や回避策は「トラブルシューティング」を参照してください。

**表示中心オプションについて**

`-Z`（方位角）と `-A`（高度角）のオプションで、画面の表示中心を指定できます。

デフォルトでは `-Z 180`（南向き）、`-A 90`（天頂）です。画面下が南、画面左が東で、天頂を見上げたような円形の表示になります。

例えば、`-Z 90`（東向き）、`-A 5`（高度5度＝地平線から5度見上げる）にすると、おおよそ半円型で星空が表示されます。

→ 東の空に夏の大三角（ベガ、アルタイル、デネブ）を捉えた表示 [<img width="40px" src="docs/images/screenshot2t.png" />](docs/images/screenshot2.png)

方位角は度数または16方位の方位記号（大小区別なし）で指定できます。
例: `-Z E`, `-Z ne`, `-Z SSW`（= 202.5°）。
（対応表: 0=北, 90=東, 180=南, 270=西）

**等級上限オプションについて**

`-V 等級` で指定した等級までの明るさの星を描画します。
デフォルトは `-V 6.0` です。例えば 9.0 等級を指定すると、約8万3千個の星が描画されます。
この値を大きくすると処理が重くなる点に注意してください。

→ 9.0等級まで表示 [<img width="40px" src="docs/images/screenshot3t.png" />](docs/images/screenshot3.png)

**日時指定オプションについて**

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

**緯度経度の直接指定**

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

### キー操作

* **← / →**: 視線の方位を ±5° 回転
* **M**: 月の5倍表示をトグル
* **F11**: フルスクリーン表示の切り替え
* **ESC**: フルスクリーンから復帰
* **Q**: 終了

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
$ zsterview
qt.qpa.plugin: From 6.5.0, xcb-cursor0 or libxcb-cursor0 is needed to load the Qt xcb platform plugin.
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.

Available platform plugins are: eglfs, offscreen, wayland-egl, linuxfb, wayland, minimal, xcb, vkkhrdisplay, minimalegl, vnc.
```

この場合は以下で `libxcb-cursor0` をインストールしてください：

`sudo apt install libxcb-cursor0`

### ネットワークが遅い/オフラインで使いたい
雲の描画は S3 から衛星画像を取得します。回線が細い/オフラインのときは `-c 0` で雲描画を無効化してください。
雲を無効化しても恒星・惑星・空の色の表示は利用できます。

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

## ライセンス

このソフトウェアは [MIT](LICENSE.txt) の下で提供されています。

ただし、 **同梱されているデータ** はそれぞれのライセンスに従って再配布されます。

以下のパスは `src/zstarview/data/` 配下を基準としています。

| ファイル                                         | 内容                                               | 出典                                                                       | ライセンス                                                                                                                      |
| ------------------------------------------------ | -------------------------------------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `cities1000.txt`, `admin1CodesASCII.txt` | 人口1000人以上の都市一覧                           | [GeoNames](https://download.geonames.org/export/dump/)                     | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                       |
| `stars/hip_main.dat.zip`                                   | Hipparcos および Tycho カタログ（ESA 1997）        | [CDS Strasbourg](https://cdsarc.cds.unistra.fr/ftp/I/239/)                 | [ODbL](https://www.data.gouv.fr/licences) または [CC BY-NC 3.0 IGO](https://creativecommons.org/licenses/by-nc/3.0/igo/)（非商用） |
| `stars/IAU-Catalog of Star Names (always up to date).csv`             | IAU 恒星名作業部会 (WGSN) による恒星固有名カタログ | [exopla.net](https://exopla.net/star-names/modern-iau-star-names/)         | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                       |
| `Noto_Sans/*`, `Noto_Sans_Symbols/*`     | テキスト / 惑星記号表示フォント                             | [Google Fonts](https://fonts.google.com/)   | [SIL Open Font License 1.1](https://openfontlicense.org)                                                                        |

## クレジット

* 天文データを提供していただいている CDS Strasbourg および ESA Hipparcos Mission に感謝します。
* 都市データは GeoNames に基づいています。
* 恒星の固有名は IAU 恒星名作業部会 (WGSN) による承認済みリスト（[exopla.net](https://exopla.net/star-names/modern-iau-star-names/) 経由）を使用しています。
* 雲データは気象衛星 **Himawari**（提供: JMA）および **NOAA GOES** シリーズ（提供: NOAA/NESDIS）による赤外線観測データを、それぞれの S3 公開バケットから取得して利用しています。
* フォントは Google Noto Project を利用しています。
* ウィンドウタイトル「Zenith Star View」は ChatGPT の提案に由来します。
* Gemini および ChatGPT に、仕様の相談、コード生成、デバッグなど、多くの助力をいただきました。

## 付録

→ [2025年の月食, 2026〜2028年の皆既日食](docs/appendix-eclipses-ja_JP.md)
