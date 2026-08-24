# 天体

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `-V`, `--vmag-limit V_MAG_LIMIT` | 表示する恒星の等級上限を指定します。 | `7.0` |
| `--vmag-brightness-multiplier MULTIPLIER` | 等級 1 段階あたりの光量変化倍率（`1.58`〜`2.512`、デフォルト `2.5`。Pogson の定義は `2.512`）を指定します。※3 | `2.5` |
| `-m`, `--enlarge-moon` | `--moon-style sphere --moon-scale 5` の互換ショートカットです。`--moon-style` または `--moon-scale` と同時には指定できません。 | |
| `--moon-style {marker,sphere,image}` | 月の描画方式を指定します。`marker` は月相を示す簡潔な輪郭、`sphere` はLambert陰影による自前の球体月相、`image` はNASA Dial-A-Moon画像を表示します。画像を利用できない間は平面的な自前月相へフォールバックします。 | `marker` |
| `--moon-scale {1,2,3,4,5,6,7,8}` | 月の表示倍率を整数で指定します。倍率は `marker`、`sphere`、`image` のすべてに適用されます。 | `1` |
| `--bright-bodies {outline,fill}` | 明るい恒星、太陽、惑星の描画モードを指定します。`outline` は輪郭を強調し、`fill` は塗りつぶし表示にします。このオプションは月には影響しません。 | `outline` |
| `-s`, `--star-base-radius STAR_BASE_RADIUS` | 2 等星の基本サイズを指定します。 | `4.0` |
| `-w`, `--expected-render-width EXPECTED_RENDER_WIDTH` | 恒星をフル解像度で描画する想定ウィンドウ幅を指定します。天球幅がこの値を超える場合、恒星レイヤーは平方根スケーリングで描画します。 | `600` |
| `--show-dso-initial true\|false` | 起動時に DSO（deep-sky objects）を表示するかを指定します。 | 自動（カタログがあれば表示） |
| `--show-asterisms-initial true\|false` | 起動時にアステリウムを表示するかを指定します。 | `show` |
| `--asterism-opacity OPACITY` | 通常のアステリズム線の絶対opacityを指定します（0.0〜0.5）。0.0で通常線を非表示にします。ホバー強調は表示され、通常線より暗くなりません。 | `0.1`（Atlasはテーマ既定値を維持） |
| `--diffuse-sky-source {akari,gaia}` | 拡散全天レイヤーのソースを選択します。Gaiaは同梱のEDR3総合明るさ・色テクスチャ、AKARIは準備済み遠赤外線キャッシュを使用します。 | `gaia` |
| `--diffuse-sky-opacity OPACITY` | 選択した拡散全天レイヤーの不透明度を指定します（0.0〜1.0）。0.0で無効にします。 | Gaia・AKARIともに`0.30` |
| `--akari-ir-bands-opacity OPACITY` | `--diffuse-sky-source akari --diffuse-sky-opacity OPACITY` の互換ショートカットです。`--diffuse-sky-*` との併用はエラーになります。 | — |
| `--twinkle-count N` | 2秒ごとの表示更新で選ぶ、星の瞬き（シンチレーション）候補数を指定します。`0` で瞬きを無効にします。通常GUIでのみ利用できます。 | `30` |
| `--show-guidelines-initial true\|false` | 起動時に幾何学・天体ガイドラインを表示するかを指定します。 | `show` |
| `-i`, `--sky-update-interval SECONDS` | 恒星・DSOなどのsky snapshotを更新する時間間隔（秒）を指定します。空色ディスクはこの値から独立し、太陽高度 `+15`〜`-15` 度では15秒、それ以外では60秒ごとに更新します。空色ディスクの計算画像は基本的に縦横1/4で、1920px幅を超える大きなディスクでは平方根スケールでさらに縮小します。 | `60` |

#### 起動時のオーバーレイ表示設定

起動時の初期表示をメニュー操作なしで切り替えるには、次を使います。

```bash
# DSO は非表示、アステリウムは表示で起動
zstarview --show-dso-initial false --show-asterisms-initial true Tokyo
```

#### 月表示オプション

月の描画方式と倍率は独立して指定できます。

```bash
# Lambert陰影による自前の球体月相を3倍で表示
zstarview --moon-style sphere --moon-scale 3 Tokyo

# Dial-A-Moon画像を5倍で表示（取得できない間は自前月相）
zstarview --moon-style image --moon-scale 5 Tokyo
```

GUIでは、Mキーと **Moon Option** メニュー項目が、設定した月表示を一時的に切り替えます。
既定の `marker`、1倍では、自前の `sphere`、5倍との間を切り替えます。それ以外の方式または
倍率を設定した場合は、その設定を一時的に `marker`、1倍へ戻し、次の操作で元の設定を復元します。

<a id="about-magnitude-limit"></a>

#### 等級上限オプションについて

`-V 等級` で、指定した等級までの明るさの星を描画します。
デフォルトは `-V 7.0` です。現在の同梱星表では最大 `-V 10.5` までサポートしており、その場合の描画候補となる恒星は約 536,000 個です。
この値を大きくすると描画時間も増えます。

### 対応しているアステリウム

ここでの「アステリウム」は、見かけ上の星の並びを結んだ通称パターンです。  
**IAU（国際天文学連合）が定義する正式な 88 星座の境界とは別概念**です。

- 冬: `Winter Triangle`（冬の大三角）, `Orion's Belt`（オリオンの三ツ星）, `Winter Hexagon`（冬のダイヤモンド）, `Southern Cross`（南十字）, `Southern Pointers`, `Diamond Cross`, `False Cross`
- 春: `Big Dipper`（北斗七星）, `Little Dipper`, `Spring Triangle`, `Arc to Arcturus`, `Leo Sickle`, `Southern Triangle`
- 夏: `Summer Triangle`（夏の大三角）, `Northern Cross`, `Teapot`, `Keystone`
- 秋: `Great Square of Pegasus`（ペガススの四辺形）, `Circlet of Pisces`, `Water Jar of Aquarius`, `Cassiopeia W`, `House of Cepheus`, `Job's Coffin`

#### 脚注

※3 最も明るい等級差の倍率は、古典的な Pogson 値 \(100^{1/5}\approx2.512\) を超えられません。
