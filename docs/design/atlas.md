# zstarview Atlas 設計

最終更新: 2026-07-10

この文書は、Atlas の内部設計案をまとめる。
利用者向け仕様は `docs/atlas-specification.md` を参照する。

## 1. アプリケーション境界

Atlas は、`zstarview` パッケージ内の別アプリ入口として実装する。
新リポジトリは作らず、地点解決、天体計算、DEM、水面、都市アウトライン、航空機、人工衛星などのデータ取得・処理は共有する。
一方で描画は、通常の `zstarview` と同じ演出的パイプラインを既定値だけ変えて流用するのではなく、位置確認に向いた instrument presentation から共有 scene data を描画する。

想定する構成は次の通り。

- `zstarview.gui.viewer`
  - 既存の星空ビューア入口。
- `zstarview.gui.atlas`
  - Atlas の白ベース星図入口。
- `zstarview.gui.window_inputs`
  - 共有の `SkyWindowUserOptions` / `SkyWindowRuntimeOptions` を引き続き使う。
- `zstarview.render`
  - 共有 scene data を使う。
  - 通常ビューアは scenic presentation、Atlas は instrument presentation を使う。

### 1.1 Presentation 境界

通常の `zstarview` は scenic presentation として扱う。
この presentation は、空の雰囲気、夜間光、Earth guide、雲、にじむような重ね塗り、viewport interaction 中の fast rendering などを使う。

Atlas は instrument presentation として扱う。
この presentation は、存在する対象の位置を読むための表示に寄せる。
雲、夜間光、Earth guide glow、地形 ridge glow などの雰囲気用合成は base scene では通さず、地形主線、水面、都市アウトライン、恒星、太陽・月・惑星、航空機、人工衛星を安定した順序で描く。
viewport interaction 中も表示内容を fast-mode 用に減らさず、通常時と同じ presentation を使う。

恒星データの保持方針も presentation と独立した内部設定として分ける。
通常ビューアは `star_data_policy="scenic_view_scoped"` を使い、現在の視野に必要な恒星だけを sky worker 側で保持する。
Atlas は `star_data_policy="positional_static"` を使い、等級上限などで選ばれた恒星を視線方向ではマスクせずに保持し、描画時の FOV 判定に任せる。
これにより、矢印キーなどで視線方向を変更している間も、新しく画面に入る領域の恒星が一時的に欠けにくくなる。

## 2. 入口プロファイル

Atlas は、CLI オプションを個別に増やすより、まず入口専用の既定プロファイルで挙動を決める。

初期プロファイルでは次の値を上書きする。

- theme: Atlas 専用の白ベーステーマ。通常の `zstarview` からは選択できない。
- sky opacity: `0.0`。
- labels: 黒文字を基本にする。
- observation info: 画面下部に常時表示し、`ThemeStyle.text` から Atlas のラベル色を解決する。
- vmag limit: 既定 `4.0`、`-V` の上限 `6.0`。
- ground tint: `0.08`。
- night light: `0.0`。
- ridge glow: `0.0`。
- terrain horizon: 有効。
- urban outline: 有効。
- water surface: 有効。
- Earth guide: 有効。Atlas 専用の単一細線パスで描き、通常版や fast-mode の塗り線・多重下敷き線は使わない。
- DSO: 有効。通常版と同系色の暗い青による単一塗りとし、下敷き線は使わない。
- asterisms: 有効。通常版と同系色の暗い青緑による単一細線とし、ホバーしていない通常時は低アルファ・細線にする。
- guidelines: 有効。水平線、方位ラベル、天の赤道、黄道、天頂・天底マーカーを黒または濃いグレーのくっきりした細線で描く。
- direction grid: guidelines 有効時は常時表示する。
- Space キーの簡易表示では地面 tint と都市アウトラインを省略し、ラベルあり状態では共有の名前付き恒星ラベル描画を有効にする。
- aircraft: 有効。
- satellites: 有効。
- cloud: 有効。ただし白背景で見える青灰色にする。
- solar-system bodies: filled, opaque markers; dark opaque labels.

既存 CLI オプションで明示された値は、入口プロファイルの既定値より優先してよい。

## 3. テーマとラベル

Atlas 向けには、通常のテーマとは分離した白ベーステーマを first-class な内部スタイルとして持つ。
`object-white` は廃止し、通常の `zstarview` のテーマ選択肢にも Atlas の白ベーステーマを追加しない。

Atlas の白ベーステーマは次の性質を持つ。

- window background は完全な白とし、背景テクスチャや紙色の着色は行わない。
- sky disc は既定で無効。
- ラベルと HUD は黒または濃いグレー。
- ラベルアウトラインは原則不要。
- 補助レイヤーは黒、濃いグレー、青灰色、水色を使い、線幅は細く調整する。アルファを十分に確保し、中間色に見えないくっきりした線にする。

### 3.1 テーマとレイヤースタイルの責務

既存の `ThemeStyle` は、画面全体のテーマを表す責務に留める。
具体的には、背景、基本文字色、ステータス文字色、ウィンドウ chrome、sky disc、splash、明背景向けの大枠の振る舞いを扱う。

航空機、人工衛星、地形地平線、都市アウトライン、水面、Earth guide など、描画レイヤー固有の色、線幅倍率、不透明度倍率、図形アウトラインは、別のレイヤースタイルとして扱う。
このレイヤースタイルは `ThemeStyle` から参照できる構造にするのが望ましい。
描画パイプラインや GUI 側が航空機色、衛星色、水面色などを個別引数で持つ形にはしない。

水平線、方位ラベル、天の赤道、黄道、天頂・天底マーカーも `GuideStyle` として `ThemeStyle` から参照する。
Atlas は黒または濃いグレーの細い線を使い、黄道だけを点線、その他を実線とする。通常版の色付き・破線を使うガイドとは別の視認性設定にする。
ガイド線幅、グリッド線幅、マーカー線幅はスタイルから注入し、共有描画関数内の通常版向け固定値を直接変更しない。Earth guide も通常版の多重線パスとは別に、共有投影処理を使う単一細線モードで描く。
投影、クリッピング、描画順のルーチンは共有し、ガイド線のスタイルだけをテーマから注入する。

想定する構造は次の通り。

```python
@dataclass(frozen=True, slots=True)
class OverlayLayerStyle:
    rgb: tuple[int, int, int]
    alpha_scale: float = 1.0
    width_scale: float = 1.0
    outline_rgba: tuple[int, int, int, int] | None = None
    label_rgb: tuple[int, int, int] | None = None
    label_outline_rgba: tuple[int, int, int, int] | None = None
    line_alpha: int | None = None
    fill: bool = True


@dataclass(frozen=True, slots=True)
class OverlayStyles:
    aircraft: OverlayLayerStyle
    satellite: OverlayLayerStyle
    terrain_horizon: OverlayLayerStyle
    urban_outline: OverlayLayerStyle
    water: OverlayLayerStyle
    earth_guide: OverlayLayerStyle
    asterism: OverlayLayerStyle
    dso: OverlayLayerStyle
```

`ThemeStyle` へは、個別の `aircraft_rgb` や `satellite_rgb` を直接増やすのではなく、`overlays: OverlayStyles` のようにグループ化して追加する。
これにより、テーマ全体の責務と、各描画レイヤーの視認性調整を分離できる。

`OverlayStyles.asterism` と `OverlayStyles.dso` も同じグループで管理する。DSO は `fill=True`、アステリズムは `line_alpha` と `width_scale` を使って、通常版と同系色のまま暗く細くする。

### 3.2 文字色とアウトラインの分担

文字色だけは `ThemeStyle` とレイヤースタイルの責務が重なりやすい。
基本方針は次の通りにする。

- HUD、汎用ラベル、状態表示は `ThemeStyle.text` / `ThemeStyle.status_text` を使う。
- 航空機、人工衛星、惑星など対象色と紐づくラベルは、対応する `OverlayLayerStyle.label_rgb` が指定されていればそれを使う。
- `OverlayLayerStyle.label_rgb` が未指定なら、`ThemeStyle.text.foreground_rgb` を使う。
- 文字の縁取りは原則 `ThemeStyle.text.outline_rgba` / `outline_width` を使う。
- 対象ラベルだけ縁取りを変えたい場合に限り、`OverlayLayerStyle.label_outline_rgba` で上書きできる。

Atlas の専用テーマでは `ThemeStyle.text` が黒または濃いグレーを返すため、レイヤー側が `label_rgb` を指定しないラベルは自然に暗色になる。
描画ルーチンへ Atlas 用の色や線幅を直書きしない。

レイヤースタイルの `outline_rgba` は文字用ではなく、線、点、十字マーカー、面、リボンなど図形そのものの縁取りまたは下敷き線を表す。
例えば、白背景の航空機軌跡では、必要な場合だけ濃い半透明の下敷き線を描き、その上に本体色の細線を描くことで視認性を補う。

実行時のレイヤー opacity は CLI や UI の値を優先し、レイヤースタイル側の `alpha_scale` はテーマ由来の基礎倍率として扱う。
最終 alpha は概ね `object_alpha * runtime_layer_opacity * alpha_scale` のように合成する。

## 4. 恒星描画パス

Atlas では、恒星のアウトラインと本体を別パスで描く。

描画処理は概ね次の手順にする。

1. 等級上限で候補恒星を絞る。
2. 投影と見かけサイズを計算する。
3. 見かけ直径が 1.0 px 未満の恒星を除外する。
4. 残った恒星のアウトラインを黒または濃グレーで描く。
5. 同じ恒星の本体を恒星色で描く。

アウトラインは恒星本体より前に描く。
これにより、近接した恒星がある場合でも、後から描くアウトラインで既に描いた恒星本体の色を潰すことを避ける。

この描画パスは通常の `zstarview` へは既定適用しない。
通常の暗背景ビューアでは、暗い星のサブピクセル表現や淡い発光表現の価値が高いためである。

## 5. レイヤー構成

Atlas では、次の順序を基本にする。

1. 白背景
2. 稜線下の ground tint
3. 常時表示の方位・高度グリッド
4. 水平線、天の赤道、黄道、方位ラベル、天頂・天底マーカー
5. DSO の単一塗り
6. アステリズムの単一細線
7. 地形地平線、都市アウトライン、水面、Earth guide
8. 雲
9. 恒星アウトライン
10. 恒星本体
11. 太陽、月、惑星
12. 航空機、人工衛星
13. 検索マーカー、ラベル、HUD

地理補助レイヤーは、場所の特定に役立つため既定で含める。
ただし、主対象を妨げないように通常の `zstarview` より淡い配色にしてよい。

## 5.1 雲の width スタイル

Atlas の雲は、通常の `zstarview` で使っている `width` スタイルを基礎にする。
衛星データの取得、`CloudAltAzGrid` の生成・キャッシュ、部分カバーの欠損マスク、
視線変更後の世代管理は既存の雲パイプラインを共有する。Atlas 専用に別の雲データ
形式や別の衛星投影を追加しない。

### 描画責務

通常の `zstarview` の `SkyCompositorCache` は、sky disc を下地にして雲を加算合成し、
雲のある箇所を灰色化する責務を持つ。Atlas は sky disc を描かないため、この合成を
そのまま呼び出さない。代わりに instrument presentation の描画列へ、次の独立した
雲ストライプパスを挿入する。

1. `scene.cloud_altaz_grid` を screen-space の雲量マップへサンプリングする。
2. 既存 `width` 方式と同じ周期、中心対称のストライプ、雲量に応じた最大線幅、
   線の外側へ向かう ease-out 的な alpha 減衰を適用する。
3. 白いストライプ画像を作る既存実装を、色を注入できる RGBA レンダラーへ拡張する。
4. Atlas の雲色（淡い青灰色）を使い、白背景へ通常の alpha 合成を行う。
5. `cloud_missing_mask` がある場合は、通常ビューアと同じマスクでストライプを抑え、
   欠損色を別レイヤーとして重ねる。

Atlas の描画順は、既存の instrument context layers（地形、都市、水面、Earth guide）
の直後、恒星アウトラインの直前とする。雲は位置確認の補助情報なので、恒星本体や
太陽・月・惑星の輪郭を隠さない。雲の下に地面 tint を置き、雲の上に恒星を置くことで、
白背景上でも雲の面と天体の点を区別できる。

### 設定と既定値

- Atlas の雲表示は既定で有効にする。
- `cloud_opacity`、`cloud_stripe`、衛星データの時刻可否、部分カバーの欠損表示は、
  既存 CLI / 起動プロファイルの値を優先する。
- Atlas の入口プロファイルには、通常ビューアと同じ既定の `width` 設定を明示的に
  設定する。これにより、共有プロファイルの既定値変更で Atlas の見え方が予期せず
  変わらないようにする。
- `cloud_opacity == 0`、ストライプ本数が `0`、または最大線幅が `0` の場合は雲を
  描かない。衛星データ未取得・時刻が非リアルタイム・全域欠損の場合も、空の雲面を
  描かず、利用可能なら欠損マスクだけを表示する。
- Atlas の `cloud` スタイルは `OverlayLayerStyle` と同じく theme から参照できる
  専用スタイルにする。少なくとも `rgb`、`alpha_scale`、`width_scale`、欠損 tint の
  色と alpha を持たせ、描画関数へ Atlas 固有の色を直書きしない。

想定する最小のスタイル拡張は次の通りである。

```python
@dataclass(frozen=True, slots=True)
class CloudLayerStyle:
    rgb: tuple[int, int, int]
    alpha_scale: float = 1.0
    width_scale: float = 1.0
    missing_rgba: tuple[int, int, int, int] = (255, 220, 80, 45)
```

`ThemeStyle` に `cloud: CloudLayerStyle` または `OverlayStyles.cloud` として保持し、
通常テーマは現在の白 / 夜間の見え方を維持する。Atlas では、白背景上で識別できる
青灰色を設定する。雲の RGB と alpha は、既存の `compose_cloud_over_sky()` の
`cloud_tint_rgb` と同じ入力境界に集約し、通常ビューアと Atlas の色分岐を一箇所に
する。

### 実装境界

- `clouddisc`、cloud worker、`CloudAltAzGrid`、サンプリング、キャッシュは変更しない。
- `gui/composite.py` の width レンダラーは、出力 RGB を固定値にせず、色を引数で受ける
  か、アルファマスクと色塗りを分離する。既存呼び出しの既定色は白にして互換性を保つ。
- Atlas の instrument presentation は、sky disc を要求しない雲描画ヘルパーを呼ぶ。
  既存 `SkyCompositorCache.draw()` に `instrument=True` の分岐を追加して責務を混ぜない。
- 雲画像生成は UI スレッドで行わず、既存の cloud worker 完了結果を使い、再描画時には
  `(cloud grid identity, geometry, projection, stripe style, cloud style)` をキーにした
  画像キャッシュを使う。リサイズ・視線変更・新しい衛星データ到着時は旧画像を破棄する。

### 検証項目

- 同じ `CloudAltAzGrid` と同じ geometry で、通常ビューアと Atlas の width ストライプの
  周期・中心線・雲量による線幅順位が一致する。
- 雲量の低い領域が消え、高い領域ほど線幅が広くなる。alpha を濃くするだけの挙動に
  ならないことを確認する。
- Atlas の白背景で雲が白飛びせず、恒星・惑星・ラベルを隠さない。
- 部分カバーではカバー済み領域だけにストライプが出て、欠損色が雲色と混同されない。
- cloud opacity 0、空グリッド、全域欠損、未来時刻、リサイズ、視線変更中の世代不一致を
  回帰テストする。
- `tests/clouddisc/test_cloud_hatch.py` に色注入・Atlas 幅マッピングの純粋関数テストを
  追加し、instrument presentation の描画テストでは cloud helper が呼ばれる順序と
  Atlas の雲色を検証する。

## 6. 共有と分離

共有するもの。

- 地点解決。
- 時刻解釈。
- 天体計算。
- 雲、航空機、人工衛星、地形、都市アウトライン、水面の取得とキャッシュ。
- 検索とジャンプ先解決。
- GUI ウィンドウの基本操作。

分離するもの。

- アプリ入口。
- 既定表示プロファイル。
- Atlas 専用の白ベーステーマ。
- 白背景向けレイヤースタイル。
- 白背景向け雲色。
- Atlas 向け恒星アウトライン描画。
- 利用者向けドキュメント。

## 7. 実装順序案

1. 別入口とプロファイルだけを追加し、既存描画でAtlasアプリとして起動できるようにする。完了。
2. Atlas 専用テーマ、暗色のくっきりした細線、白ベース背景を追加する。完了。
3. sky disc 既定無効、地理補助レイヤー既定有効の挙動を固定する。完了。
4. Atlas 向けの恒星アウトライン別パスと 1.0 px 未満非表示を追加する。完了。
5. Atlas 向け雲色を調整し、instrument 描画へ接続する。完了。
6. `ThemeStyle` から参照するレイヤースタイルを追加し、航空機、人工衛星、地理補助レイヤー、ガイドのAtlas向け配色と図形アウトラインを整理する。完了。
7. 起動前ダイアログを共有するか、専用入口を追加するか判断する。未実施。
8. Atlas の雲 width パスを instrument presentation に接続する。完了。
