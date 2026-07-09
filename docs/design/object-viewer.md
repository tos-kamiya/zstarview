# 白背景オブジェクトビューア設計

最終更新: 2026-07-10

この文書は、白背景オブジェクトビューアの内部設計案をまとめる。
利用者向け仕様は `docs/object-viewer-specification.md` を参照する。

## 1. アプリケーション境界

白背景オブジェクトビューアは、`zstarview` パッケージ内の別アプリ入口として実装する。
新リポジトリは作らず、地点解決、天体計算、DEM、水面、都市アウトライン、航空機、人工衛星などのデータ取得・処理は共有する。
一方で描画は、通常の `zstarview` と同じ演出的パイプラインを既定値だけ変えて流用するのではなく、位置確認に向いた instrument presentation から共有 scene data を描画する。

想定する構成は次の通り。

- `zstarview.gui.viewer`
  - 既存の星空ビューア入口。
- `zstarview.gui.object_viewer`
  - 白背景オブジェクトビューア入口。
- `zstarview.gui.window_inputs`
  - 共有の `SkyWindowUserOptions` / `SkyWindowRuntimeOptions` を引き続き使う。
- `zstarview.render`
  - 共有 scene data を使う。
  - 通常ビューアは scenic presentation、白背景オブジェクトビューアは instrument presentation を使う。

### 1.1 Presentation 境界

通常の `zstarview` は scenic presentation として扱う。
この presentation は、空の雰囲気、夜間光、Earth guide、雲、にじむような重ね塗り、viewport interaction 中の fast rendering などを使う。

白背景オブジェクトビューアは instrument presentation として扱う。
この presentation は、存在する対象の位置を読むための表示に寄せる。
雲、夜間光、Earth guide glow、地形 ridge glow などの雰囲気用合成は base scene では通さず、地形主線、水面、都市アウトライン、恒星、太陽・月・惑星、航空機、人工衛星を安定した順序で描く。
viewport interaction 中も表示内容を fast-mode 用に減らさず、通常時と同じ presentation を使う。

恒星データの保持方針も presentation と独立した内部設定として分ける。
通常ビューアは `star_data_policy="scenic_view_scoped"` を使い、現在の視野に必要な恒星だけを sky worker 側で保持する。
白背景オブジェクトビューアは `star_data_policy="positional_static"` を使い、等級上限などで選ばれた恒星を視線方向ではマスクせずに保持し、描画時の FOV 判定に任せる。
これにより、矢印キーなどで視線方向を変更している間も、新しく画面に入る領域の恒星が一時的に欠けにくくなる。

## 2. 入口プロファイル

白背景オブジェクトビューアは、CLI オプションを個別に増やすより、まず入口専用の既定プロファイルで挙動を決める。

初期プロファイルでは次の値を上書きする。

- theme: `object-white`。
- sky opacity: `0.0`。
- labels: 黒文字を基本にする。
- vmag limit: 既定 `4.0`、`-V` の上限 `6.0`。
- ground tint: `0.08`。
- night light: `0.0`。
- ridge glow: `0.0`。
- terrain horizon: 有効。
- urban outline: 有効。
- water surface: 有効。
- Earth guide: 有効。
- guidelines: 有効。水平線、方位ラベル、天の赤道、黄道、天頂・天底マーカーを黒系の細線で描く。
- direction grid: guidelines 有効時は常時表示する。
- aircraft: 有効。
- satellites: 有効。
- cloud: 有効。ただし白背景で見える青灰色にする。
- solar-system bodies: filled, opaque markers; dark opaque labels.

既存 CLI オプションで明示された値は、入口プロファイルの既定値より優先してよい。

## 3. テーマとラベル

白背景向けには、完全な白背景を first-class なテーマとして持つのが望ましい。
既存の `white` テーマは白系ではあるが、星空ビューアの背景・sky disc・文字色の文脈を含むため、白背景オブジェクトビューアの恒久仕様とは分けてよい。

白背景テーマは次の性質を持つ。

- window background は完全な白。
- sky disc は既定で無効。
- ラベルと HUD は黒または濃いグレー。
- ラベルアウトラインは原則不要。
- 補助レイヤーは低アルファのグレー、青灰色、水色を使う。

### 3.1 テーマとレイヤースタイルの責務

既存の `ThemeStyle` は、画面全体のテーマを表す責務に留める。
具体的には、背景、基本文字色、ステータス文字色、ウィンドウ chrome、sky disc、splash、明背景向けの大枠の振る舞いを扱う。

航空機、人工衛星、地形地平線、都市アウトライン、水面、Earth guide など、描画レイヤー固有の色、線幅倍率、不透明度倍率、図形アウトラインは、別のレイヤースタイルとして扱う。
このレイヤースタイルは `ThemeStyle` から参照できる構造にするのが望ましい。
描画パイプラインや GUI 側が航空機色、衛星色、水面色などを個別引数で持つ形にはしない。

水平線、方位ラベル、天の赤道、黄道、天頂・天底マーカーも `GuideStyle` として `ThemeStyle` から参照する。
object viewer は黒の細い線を使い、黄道だけを点線、その他を実線とする。通常版の色付き・破線を使うガイドとは別の視認性設定にする。
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


@dataclass(frozen=True, slots=True)
class OverlayStyles:
    aircraft: OverlayLayerStyle
    satellite: OverlayLayerStyle
    terrain_horizon: OverlayLayerStyle
    urban_outline: OverlayLayerStyle
    water: OverlayLayerStyle
    earth_guide: OverlayLayerStyle
```

`ThemeStyle` へは、個別の `aircraft_rgb` や `satellite_rgb` を直接増やすのではなく、`overlays: OverlayStyles` のようにグループ化して追加する。
これにより、テーマ全体の責務と、各描画レイヤーの視認性調整を分離できる。

### 3.2 文字色とアウトラインの分担

文字色だけは `ThemeStyle` とレイヤースタイルの責務が重なりやすい。
基本方針は次の通りにする。

- HUD、汎用ラベル、状態表示は `ThemeStyle.text` / `ThemeStyle.status_text` を使う。
- 航空機、人工衛星、惑星など対象色と紐づくラベルは、対応する `OverlayLayerStyle.label_rgb` が指定されていればそれを使う。
- `OverlayLayerStyle.label_rgb` が未指定なら、`ThemeStyle.text.foreground_rgb` を使う。
- 文字の縁取りは原則 `ThemeStyle.text.outline_rgba` / `outline_width` を使う。
- 対象ラベルだけ縁取りを変えたい場合に限り、`OverlayLayerStyle.label_outline_rgba` で上書きできる。

`object-white` では `ThemeStyle.text` が黒または濃いグレーを返すため、レイヤー側が `label_rgb` を指定しないラベルは自然に黒文字になる。
描画ルーチンへ `object-viewer` 用の黒文字を直書きしない。

レイヤースタイルの `outline_rgba` は文字用ではなく、線、点、十字マーカー、面、リボンなど図形そのものの縁取りまたは下敷き線を表す。
例えば、白背景の航空機軌跡では、先に濃い半透明の太線を描き、その上に本体色の細線を描くことで視認性を補う。

実行時のレイヤー opacity は CLI や UI の値を優先し、レイヤースタイル側の `alpha_scale` はテーマ由来の基礎倍率として扱う。
最終 alpha は概ね `object_alpha * runtime_layer_opacity * alpha_scale` のように合成する。

## 4. 恒星描画パス

白背景では、恒星のアウトラインと本体を別パスで描く。

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

白背景オブジェクトビューアでは、次の順序を基本にする。

1. 完全白背景
2. 稜線下の ground tint
3. 常時表示の方位・高度グリッド
4. 水平線、天の赤道、黄道、方位ラベル、天頂・天底マーカー
5. 地形地平線、都市アウトライン、水面、Earth guide
6. 雲
7. 恒星アウトライン
8. 恒星本体
9. 太陽、月、惑星
10. 航空機、人工衛星
11. 検索マーカー、ラベル、HUD

地理補助レイヤーは、場所の特定に役立つため既定で含める。
ただし、主対象を妨げないように通常の `zstarview` より淡い配色にしてよい。

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
- 白背景専用テーマ。
- 白背景向けレイヤースタイル。
- 白背景向け雲色。
- 白背景向け恒星アウトライン描画。
- 利用者向けドキュメント。

## 7. 実装順序案

1. 別入口とプロファイルだけを追加し、既存描画で白背景アプリとして起動できるようにする。完了。
2. 完全白背景テーマと黒ラベルを追加する。完了。
3. sky disc 既定無効、地理補助レイヤー既定有効の挙動を固定する。完了。
4. 白背景向けの恒星アウトライン別パスと 1.0 px 未満非表示を追加する。完了。
5. 白背景向け雲色を調整する。完了。
6. `ThemeStyle` から参照するレイヤースタイルを追加し、航空機、人工衛星、地理補助レイヤー、ガイドの白背景向け配色と図形アウトラインを整理する。完了。
7. 起動前ダイアログを共有するか、専用入口を追加するか判断する。未実施。
