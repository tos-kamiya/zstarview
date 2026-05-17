# 星空と天体

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `--sky-opacity SKY_OPACITY` | 空の色ディスクの不透明度を指定します（0.0〜1.0）。0.0 で描画を無効化します。 | `0.15` |
| `--sky-disc-altaz-rings {off,dimalt,altaz}` | 常時表示の空ディスク方位/高度オーバーレイです。`dimalt` は控えめな高度リング、`altaz` はフルグリッドを表示します。 | `dimalt` |
| `--sky-disc-altaz-rings-hover {off,dimalt,altaz}` | ホバー時の空ディスク方位/高度オーバーレイです。意味は上記と同じです。 | `altaz` |
| `--bright-bodies {outline,fill}` | 明るい天体の描画モードを指定します。`outline` では明るい恒星をひし形輪郭、惑星を輪郭のみ、月を通常表示では輪郭のみで描画し、`--enlarge-moon` や月ホバー時は通常の月描画を使います。`fill` では従来どおり塗りつぶし表示にします。 | `outline` |
| `-m`, `--enlarge-moon` | 月を 5 倍に拡大して表示します。 | |
| `-s`, `--star-base-radius STAR_BASE_RADIUS` | 2 等星の基本サイズを指定します。 | `4.0` |
| `-w`, `--expected-render-width EXPECTED_RENDER_WIDTH` | 恒星をフル解像度で描画する想定ウィンドウ幅を指定します。天球幅がこの値を超える場合、恒星レイヤーは平方根スケーリングで描画します。 | `600` |
| `-V`, `--vmag-limit V_MAG_LIMIT` | 表示する恒星の等級上限を指定します。 | `7.0` |
| `--vmag-brightness-multiplier MULTIPLIER` | 等級 1 段階あたりの光量変化倍率（`1.58`〜`2.512`、デフォルト `2.5`。Pogson の定義は `2.512`）を指定します。※3 | `2.5` |
| `-i`, `--sky-update-interval SECONDS` | 星空を更新する時間間隔（秒）を指定します。 | `60` |
| `--show-dso-initial true\|false` | 起動時に DSO（deep-sky objects）を表示するかを指定します。 | 自動（カタログがあれば表示） |
| `--show-asterisms-initial true\|false` | 起動時にアステリウムを表示するかを指定します。 | `show` |

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
デフォルトは `-V 7.0` です。現在の同梱星表では最大 `-V 10.5` までサポートしており、その場合の描画候補となる恒星は約 536,000 個です。
この値を大きくすると描画時間も増えます。

#### テーマプリセットについて

`--theme` を使うと、背景とコントラストの見え方を切り替えられます。

* `night`: 標準の暗色テーマ
* `black`: より黒く不透明な背景
* `transparent`: 黒寄りの半透明背景。デスクトップ側も暗い場合に向く
* `day`: 明るい昼空寄りの背景表現
* `white`: 最も明るい淡色テーマ

### 対応しているアステリウム

ここでの「アステリウム」は、見かけ上の星の並びを結んだ通称パターンです。  
**IAU（国際天文学連合）が定義する正式な 88 星座の境界とは別概念**です。

- 冬: `Winter Triangle`（冬の大三角）, `Orion's Belt`（オリオンの三ツ星）, `Winter Hexagon`（冬のダイヤモンド）, `Southern Cross`（南十字）, `Southern Pointers`, `Diamond Cross`, `False Cross`
- 春: `Big Dipper`（北斗七星）, `Little Dipper`, `Spring Triangle`, `Arc to Arcturus`, `Leo Sickle`, `Southern Triangle`
- 夏: `Summer Triangle`（夏の大三角）, `Northern Cross`, `Teapot`, `Keystone`
- 秋: `Great Square of Pegasus`（ペガススの四辺形）, `Circlet of Pisces`, `Water Jar of Aquarius`, `Cassiopeia W`, `House of Cepheus`, `Job's Coffin`

#### 脚注

※3 最も明るい等級差の倍率は、古典的な Pogson 値 \(100^{1/5}\approx2.512\) を超えられません。
