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

## 使い方

```bash
zstarview [options] [city]
```

### 引数

| 引数     | 説明                                                          | デフォルト              |
| :----- | :---------------------------------------------------------- | :----------------- |
| `city` | 表示する都市名を指定します。省略時は前回起動時の都市を使用します。初回起動時に省略すると `Tokyo` になります。 | 前回の都市（初回は `Tokyo`） |

### オプション

| オプション                                       | 説明                                            | デフォルト  |
| :------------------------------------------ | :-------------------------------------------- | :----- |
| `-h`, `--help`                              | ヘルプメッセージを表示して終了します。                           |        |
| `-H`, `--hours HOURS`                       | 現在時刻に加算する時間数を指定します。                           | `0`    |
| `-D`, `--days DAYS`                         | 現在時刻に加算する日数を指定します。                            | `0`    |
| `-m`, `--enlarge-moon`                      | 月を3倍の大きさで表示します。                               |        |
| `-s`, `--star-base-radius STAR_BASE_RADIUS` | 星の基本サイズを指定します。                                | `15.0` |
| `-Z`, `--view-center-az VIEW_CENTER_AZ`     | 表示中心の方位角を度単位で指定します (0=北, 90=東, 180=南, 270=西)。 | `180`  |
| `-A`, `--view-center-alt VIEW_CENTER_ALT`   | 表示中心の高度を度単位で指定します (90=天頂, 0=地平線)。             | `90`   |
| `-V`, `--vmag-threshold V_MAG_THRESHOLD`    | 表示する恒星の等級（明るさ）の上限を指定します。                      | `7.0`  |

**表示中心に関するオプションについて**

`-Z`（方位角）と `-A`（高度）のオプションで、画面の表示中心を指定できます。

デフォルトでは、`-Z 180`（南向き）、`-A 90`（天頂）です。画面下が南、画面左が東で、天頂を見上げたような円形の表示になります。

例えば、`-Z 90` （東向き）、`-A 10` （高度10度＝地面から10度見上げる）にすると、おおよそ半円型で星空が表示されます。  
→ 東の空に [夏の大三角（ベガ、アルタイル、デネブ）](docs/images/screenshot2.png) を捉えた表示

**等級の上限に関するオプションについて**

`-V 等級` で指定した等級までの明るさの星を描画します。デフォルトは `-V 7.0` です。例えば、8.0等級を指定すると、約4万1千個の星が描画されます。この値が大きいと処理が重くなります。

→ [8.0等級](docs/images/screenshot3.png) まで表示

### キー操作

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

## 依存パッケージ

* Python 3.10+
* [`appdirs`](https://pypi.org/project/appdirs/)
* [`PyQt5`](https://pypi.org/project/PyQt5/)
* [`astropy`](https://pypi.org/project/astropy/)
* [`skyfield`](https://pypi.org/project/skyfield/)
* [`numpy`](https://pypi.org/project/numpy/)
* [`Pillow`](https://pypi.org/project/Pillow/)

## ライセンス

このソフトウェアは [MIT](LICENSE.txt) の下で提供されています。

ただし、 **同梱されているデータ** はそれぞれのライセンスに従って再配布されます。

| ファイル                                         | 内容                                               | 出典                                                                       | ライセンス                                                                                                                      |
| ------------------------------------------------ | -------------------------------------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `data/cities1000.txt`                            | 人口1000人以上の都市一覧                           | [GeoNames](https://download.geonames.org/export/dump/)                     | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                       |
| `data/stars/hip_main.dat`                                   | Hipparcos および Tycho カタログ（ESA 1997）        | [CDS Strasbourg](https://cdsarc.cds.unistra.fr/ftp/I/239/)                 | [ODbL](https://www.data.gouv.fr/licences) または [CC BY-NC 3.0 IGO](https://creativecommons.org/licenses/by-nc/3.0/igo/)（非商用） |
| `data/stars/IAU-Catalog-of-Star-Names.csv`             | IAU 恒星名作業部会 (WGSN) による恒星固有名カタログ | [exopla.net](https://exopla.net/star-names/modern-iau-star-names/)         | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                       |
| `data/NotoSansSymbols-VariableFont_wght.ttf`     | 惑星記号表示用フォント                             | [Google Fonts](https://fonts.google.com/noto/specimen/Noto+Sans+Symbols)   | [SIL Open Font License 1.1](https://openfontlicense.org)                                                                        |

## クレジット

* 天文データを提供していただいている CDS Strasbourg および ESA Hipparcos Mission に感謝します。
* 都市データは GeoNames に基づいています。
* 恒星の固有名は IAU 恒星名作業部会 (WGSN) による承認済みリスト（[exopla.net](https://exopla.net/star-names/modern-iau-star-names/) 経由）を使用しています。
* フォントは Google Noto Project を利用しています。
* ウィンドウタイトル「Zenith Star View」は ChatGPT の提案に由来します。
* Gemini および ChatGPT に、仕様の相談、コード生成、デバッグなど、多くの助力をいただきました。
