# zstarview 🌌

雲があっても、太陽が出ていても、満天の星空を。

Zenith Star View は、地球上の任意の都市を指定して、頭上の星空を描画するアプリケーションです。

- 明るい恒星、惑星、天の赤道、黄道をリアルタイムで描画
- 都市名で場所を指定可能（GeoNames 収録）

![](docs/images/screenshot1.png)

## インストール方法（推奨：`pipx`）

[`pipx`](https://pypa.github.io/pipx/) を使ってインストールする想定です。

```bash
pipx install git+https://github.com/tos-kamiya/zstarview.git
```

> 注記（X11 環境・Ubuntu/Debian）: X11 セッションでは Qt の xcb プラグインが実行時に `libxcb-cursor0` を必要とする場合があります。`sudo apt install libxcb-cursor0` によりインストールしてください。

## 使い方

```bash
zstarview [options] [city]
```

### 引数

| 引数     | 説明                                                          | デフォルト              |
| :----- | :---------------------------------------------------------- | :----------------- |
| `city` | 表示する都市名を指定します。省略時は前回起動時の都市を使用します。初回起動時に省略すると `Tokyo` になります。 | 前回の都市（初回は `Tokyo`） |

### オプション

| オプション                                       | 説明                                      | デフォルト |
| :------------------------------------------ | :---------------------------------------- | :------- |
| `-h`, `--help`                              | ヘルプメッセージを表示して終了します。             |          |
| `-H`, `--hours HOURS`                       | 現在時刻に加算する時間数を指定します。             | `0`      |
| `-D`, `--days DAYS`                         | 現在時刻に加算する日数を指定します。              | `0`      |
| `--datetime "YYYY-MM-DD HH:MM:SS [TZ]"`     | 絶対的な日時を指定します（TZを省略した場合はUTC）。 |          |
| `-m`, `--enlarge-moon`                      | 月を5倍に拡大して表示します。                  |          |
| `-s`, `--star-base-radius STAR_BASE_RADIUS` | 星の基本サイズを指定します。                   | `8.0`    |
| `-Z`, `--view-center-az VIEW_CENTER_AZ`     | 表示中心の方位角を指定します。                  | `180`    |
| `-A`, `--view-center-alt VIEW_CENTER_ALT`   | 表示中心の高度角を指定します。                  | `90`     |
| `-V`, `--vmag-threshold V_MAG_THRESHOLD`    | 表示する恒星の等級（明るさ）の上限を指定します。      | `6.0`    |

**日時指定オプションについて**

`--datetime "YYYY-MM-DD HH:MM:SS [TZ]"` で絶対的な日時を指定できます。  
タイムゾーン（TZ）を省略した場合は UTC として扱われます。

タイムゾーンは以下のいずれかの形式で指定できます：

- よく使われる略称（JST, UTC, GMT, KST, HKT, AWST, ACST, AEST, NZST, NZDT, MSK, EAT）
- IANA タイムゾーン名（例: `Asia/Tokyo`, `Europe/Moscow`）
- UTCオフセット（例: `UTC+9`, `UTC-07:30`）

例:

```bash
zstarview --datetime "2025-08-17 21:00:00 JST" Tokyo
````

**表示中心オプションについて**

`-Z`（方位角）と `-A`（高度角）のオプションで、画面の表示中心を指定できます。

デフォルトでは `-Z 180`（南向き）、`-A 90`（天頂）です。画面下が南、画面左が東で、天頂を見上げたような円形の表示になります。

例えば、`-Z 90`（東向き）、`-A 10`（高度10度＝地平線から10度見上げる）にすると、おおよそ半円型で星空が表示されます。
→ 東の空に [夏の大三角（ベガ、アルタイル、デネブ）](docs/images/screenshot2.png) を捉えた表示

方位角は度数または16方位の方位記号（大小区別なし）で指定できます。
例: `-Z E`, `-Z ne`, `-Z SSW`（= 202.5°）。
（対応表: 0=北, 90=東, 180=南, 270=西）

**等級上限オプションについて**

`-V 等級` で指定した等級までの明るさの星を描画します。
デフォルトは `-V 6.0` です。例えば 9.0 等級を指定すると、約8万3千個の星が描画されます。
この値を大きくすると処理が重くなる点に注意してください。

→ [9.0等級](docs/images/screenshot3.png) まで表示

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

## ライセンス

このソフトウェアは [MIT](LICENSE.txt) の下で提供されています。

ただし、 **同梱されているデータ** はそれぞれのライセンスに従って再配布されます。

| ファイル                                         | 内容                                               | 出典                                                                       | ライセンス                                                                                                                      |
| ------------------------------------------------ | -------------------------------------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `data/cities1000.txt`                            | 人口1000人以上の都市一覧                           | [GeoNames](https://download.geonames.org/export/dump/)                     | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                       |
| `data/stars/hip_main.dat`                                   | Hipparcos および Tycho カタログ（ESA 1997）        | [CDS Strasbourg](https://cdsarc.cds.unistra.fr/ftp/I/239/)                 | [ODbL](https://www.data.gouv.fr/licences) または [CC BY-NC 3.0 IGO](https://creativecommons.org/licenses/by-nc/3.0/igo/)（非商用） |
| `data/stars/IAU-Catalog-of-Star-Names.csv`             | IAU 恒星名作業部会 (WGSN) による恒星固有名カタログ | [exopla.net](https://exopla.net/star-names/modern-iau-star-names/)         | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                       |
| `data/Noto_Sans/*`, `data/Noto_Sans_Symbols/*`     | テキスト / 惑星記号表示フォント                             | [Google Fonts](https://fonts.google.com/)   | [SIL Open Font License 1.1](https://openfontlicense.org)                                                                        |

## クレジット

* 天文データを提供していただいている CDS Strasbourg および ESA Hipparcos Mission に感謝します。
* 都市データは GeoNames に基づいています。
* 恒星の固有名は IAU 恒星名作業部会 (WGSN) による承認済みリスト（[exopla.net](https://exopla.net/star-names/modern-iau-star-names/) 経由）を使用しています。
* フォントは Google Noto Project を利用しています。
* ウィンドウタイトル「Zenith Star View」は ChatGPT の提案に由来します。
* Gemini および ChatGPT に、仕様の相談、コード生成、デバッグなど、多くの助力をいただきました。

## 付録：2025年の月食のプレビュー

2025年に観測可能な**2回の皆既月食**を、日時と都市を指定してプレビューできます。

### 2025年3月13〜14日（アメリカ大陸で見られる皆既月食）

```bash
zstarview --datetime "2025-03-14 02:58:43 America/New_York" "New York City"
zstarview --datetime "2025-03-13 23:58:43 America/Los_Angeles" "US/Los Angeles"
zstarview --datetime "2025-03-13 20:58:43 HST" "Honolulu"
zstarview --datetime "2025-03-14 00:58:43 America/Mexico_City" "Mexico City"
zstarview --datetime "2025-03-14 03:58:43 America/Sao_Paulo" "BR/São Paulo"
```

### 2025年9月7〜8日（アジア・ヨーロッパ・アフリカ・オセアニアで見られる皆既月食）

```bash
zstarview --datetime "2025-09-08 03:12:00 JST" "Tokyo"
zstarview --datetime "2025-09-08 03:11:47 KST" "Seoul"
zstarview --datetime "2025-09-08 02:11:47 Asia/Shanghai" "Beijing"
zstarview --datetime "2025-09-08 01:11:47 Asia/Bangkok" "Bangkok"
zstarview --datetime "2025-09-07 23:41:47 Asia/Kolkata" "New Delhi"
zstarview --datetime "2025-09-07 21:11:47 Europe/Istanbul" "Istanbul"
zstarview --datetime "2025-09-07 21:11:47 EAT" "Nairobi"
zstarview --datetime "2025-09-08 04:11:47 Australia/Sydney" "Sydney"
zstarview --datetime "2025-09-08 06:11:47 Pacific/Auckland" "Auckland"
```
