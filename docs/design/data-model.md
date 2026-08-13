# データモデル

この文書は、`zstarview`、`zstarview-gui`、`zstarview-export-image` で共有する主要データ構造と状態モデルをまとめる。
描画そのものは [rendering-pipeline.md](rendering-pipeline.md) を参照する。

## 1. モデルの考え方

- アプリ 3 種は同じドメインモデルを共有する。
- 違うのは、どの入口から状態を初期化し、どの永続設定を読むかである。
- そのため、データモデルは「共有 scene data」と「入口固有の launch / UI state」を分けて考える。

## 2. 主要な共有データ

### 2.1 地点・視点関連

- `ViewerData`
  - 観測地点
  - タイムゾーン
  - 表示中心の方位・高度
  - 画面投影用の edge FOV
  - 描画対象保持用の content FOV
  - 基準観測点の上に乗る追加高さ
  - 地盤標高と構造物高の表示用値
  - 追加高さを UI 表示するかどうかのフラグ
  - 画面描画に必要な視点情報
- `ScreenGeometry`
  - 描画キャンバスの中心と半径のみを表す
  - 視点中心や FOV は持たない
  - 半径は `edge_fov_deg` に対応するスクリーン半径として扱う
  - `content_fov_deg` に対応するスクリーン半径は `radius * content_fov_deg / edge_fov_deg` で求める

`ScreenGeometry` は画面サイズと view center altitude から draw-time に作る。
横長または正方形のクライアント領域では、従来通り sky disc の上端をクライアント領域上端に接する配置を使う。
縦長のクライアント領域では中心をクライアント領域中央に置き、幅ではなく高さを基準に半径を決める。

縦横比 `R = width / height` が `R < 1` のとき、半径は次の 2 つの状態を線形補間してよい。

- `R = 1` 側: `edge_fov_deg` のスクリーン半径を `height / 2` にする。
- `R = 1 / 2` 側: `content_fov_deg` のスクリーン半径を `height / 2` にする。つまり `edge_fov_deg` のスクリーン半径は `(height / 2) * edge_fov_deg / content_fov_deg` になる。

`R <= 1 / 2` では後者を固定してよい。
これにより、縦長ウィンドウでは sky disc が左右へはみ出しつつ、1:2 付近からより細長い領域では四隅に `content_fov_deg` 外側の領域が見える。

地点 dataset が持つ高さ情報と追加高さは別概念として扱う。

- mountain viewpoint
  - dataset 側の高さは山頂ビューポイントの海抜標高 `Elevation`
- tower viewpoint
  - dataset 側の高さは地表からのタワー高またはビューポイント高 `Tower height`
- `ViewerData.height_add_m`
  - どの入力種別でも、基準観測点から観測者の追加高さ
  - CLI `--height-add-m` はこの値を設定する
  - CLI `--observer-height-m` は互換エイリアスとして残る
  - 既定値は `1.7m`
- `ViewerData.ground_elevation_m`
  - DEM から求めた観測地点の地盤標高
  - DEM を取れない場合は内部計算用に `0.0m` へ正規化してよい。ただし、これは地形地平線レイヤーを表示可能とする意味ではなく、地形データなしの状態は別途保持する。
- `ViewerData.location_height_m`
  - 建物頂部、タワー高、または解決済み地点の構造物高を表す
  - 該当しない場合は `0.0m` へ正規化してよい

### 2.2 天体計算結果

- `CelestialData`
  - 恒星の描画用データ
  - 太陽、月、惑星の位置
  - 地平線、黄道、赤道などの線データ
  - DSO、ホバー対象、補助表示に必要な情報

### 2.3 航空機関連

- `AircraftSnapshot`
  - `icao24`
  - `callsign`
  - `latitude`
  - `longitude`
  - `baro_altitude_m`
  - `velocity_mps`
  - `heading_deg`
  - `vertical_rate_mps`
  - `on_ground`
  - `last_contact_unix`
  - OpenSky の 1 回の取得結果から正規化した機体状態
- `AircraftOverlayPoint`
  - `icao24`
  - `callsign`
  - `alt_deg`
  - `az_deg`
  - `trail_alt_az_points`
  - `distance_km`
  - `age_seconds`
  - `alpha_scale`
  - 描画直前まで絞り込んだ軽量軌跡モデル
  - 折れ線描画や帯状ポリゴン描画のための中間表現として使ってよい
- `AircraftState`
  - 最新スナップショットの機体一覧
  - 表示用折れ線列
  - 最終成功タイムスタンプ
  - 読込中状態
  - エラーバナー

### 2.3.1 GMNメテオ軌跡関連

- `MeteorObservation`
  - GMNの軌跡識別子と開始UTC時刻
  - WGS84楕円体高を含む発光区間の始点・終点
  - 任意のDuration、最大光度、初速度、流星群コード
  - 配布ファイルを解析した直後の地理座標モデル
- `CelestialMeteorTrail`
  - 軌跡識別子と開始UTC時刻
  - 発生時刻の観測地点から可視だった区間の始点・終点RA/Dec
  - 描画時刻や画面投影に依存しない天球固定モデル
- `GmnLoadResult`
  - 時刻窓に含まれる`MeteorObservation`列
  - 読み込めた配布ファイルと利用不能だった配布ファイル
  - staleな索引または配布ファイルを使ったかを示すフラグ
- `MeteorWindowResult`
  - 表示時刻以前の最新観測を終端とする24時間窓と`CelestialMeteorTrail`列
  - 取得元、部分欠損、stale利用状態

### 2.4 雲関連

- `CloudMeta`
  - 雲画像のデータ元メタ情報
- `CloudState`
  - 最新画像
  - 欠損マスク
  - 進行中状態
  - エラーバナー
- `CloudAltAzGrid`
  - 観測者中心の中間表現。衛星データを取得した直後に構築し、カメラ視点に依存しない。
  - `amount`: 形状 `(alt_bins, az_bins)` の `float32` 雲量配列（0.0〜1.0）。
  - `alt_min_deg` / `alt_max_deg` / `az_min_deg` / `az_max_deg`: グリッドの角度範囲。
  - `observer_lat` / `observer_lon`: グリッドを構築した観測地点。
  - `satellite` / `product` / `time_utc`: データ元メタ情報。
  - `shells_km`: 使用した雲シェルの半径リスト。
  - `source_key` / `coverage_ratio` / `source_completeness_ratio`: キャッシュ・状態管理用メタ情報。
  - `missing_mask`: 欠損セルを表す同形状の `uint8` 配列。
- `CloudImageState`
  - `altaz_grid`: 最新の `CloudAltAzGrid`（推奨）。カメラが変わってもこのグリッドを再利用する。
  - `image`: 必要に応じて `altaz_grid` から生成した雲レイヤーの `numpy RGBA`。カメラ依存のキャッシュとして使ってよい。
  - `missing_mask`: 欠損領域を表す 2D `uint8` alpha 配列。
  - `cloud_amount_field`: ストライプ線幅生成用の雲量場（後続で `CloudAltAzGrid` 由来に移行してよい）。
  - `meta`、`coverage_ratio`、`source_key`、`render_key`、`request_id`、`source_refreshed_at_utc` などの更新メタデータ

### 2.5 人工衛星関連

- `SatelliteOverlayPoint`
  - `group_key`
  - `satellite_name`
  - `alt_deg`
  - `az_deg`
  - `marker_scale`
  - 描画直前まで絞り込んだ人工衛星マーカーモデル
- `SatelliteState`
  - group ごとの最新軌道要素 records
  - 表示用マーカー列
  - 最終成功タイムスタンプ
  - 読込中状態
  - エラーバナー

### 2.6 地形・地平線関連

- `TerrainHorizonState`
  - 地形地平線の点列
  - 読込中状態
  - エラー表示状態
  - データソース状態（最新取得、stale cache、unavailable）
- `TerrainSecondaryRidgeState`
  - 主地平線とは別の補助稜線の点列群
  - 読込中状態は主地平線と共有してよい
  - 表示強度や on/off は独立してよい

### 2.6.1 水面レイヤー関連

- `WaterOverlayState`
  - 海面・内陸水域の表示点列と境界折れ線
  - 観測地点・要求範囲に対応する cache scope
  - 読込中状態
  - データソース状態（最新取得、stale cache、unavailable）

### 2.7 地球裏面ガイド関連

- `SubsurfaceEarthGuideState`
  - 粗い大陸ポリゴンまたはそれと等価な低解像度 land-mask の参照
  - 現在の観測地点と視線条件に対して投影済みの screen polyline / polygon 群
  - 最終合成用の bitmap cache または `QPainterPath`
  - 読込中状態やネットワーク取得状態は持たず、同梱 static data のみを前提としてよい
- `EarthGuideState`
  - Earth guide の表示強度
  - GUI トグルの on/off 状態
  - CLI で `--earth-guide-opacity` が `0` のときのロックアウト状態

### 2.8 ウィンドウ状態

- `SkyWindowState`
  - 現在の視点
  - 直近の描画用視点
  - `rotation_step` による 5° 刻みの粗調整
  - `viewport_interaction_mode`
  - `viewport_interaction_release_pending`
  - `viewport_interaction_stars`
  - ホバー対象
  - ハイライト対象
  - `inverted_city_enabled`: `Inverted City`のGUI一時表示状態。設定やexport-image入力には保存しない。
  - 各更新パイプラインの UI 反映状態
- `SkyWindow._frame_cache_image`
  - `paintEvent` のベース描画部分を保持する `QImage`
- `SkyWindow._frame_cache_key`
  - base frame cache の無効化条件をまとめたキー
  - ウィンドウサイズ、`render_view_center`、描画トグル、`inverted_city_enabled`、`CelestialData`、空ディスク画像、雲画像、地形/都市アウトラインなどを含む
- `SkyWindow._present_frame_cache_image`
  - base frame に通常オーバーレイを重ねた `QImage`
  - HUD、status line、mouse hover、search marker などの volatile overlay は含めず、paint 時に別途重ねる
- `SkyWindow._present_frame_cache_key`
  - present frame cache の無効化条件をまとめたキー
  - base frame key、短周期 overlay の時刻 bucket、航空機・人工衛星・台風などの入力識別子を含む

### 2.9 都市アウトライン関連

- `UrbanOutlineState`
  - 都市アウトライン点列
  - 必要に応じて base レイヤーと skyscraper レイヤーの件数を別々に保持する
  - 読込中 / 取得中バナー
  - 失敗表示状態
  - `cache` または `overture` などの現在ソース表示

### 2.10 台風関連

- `TropicalCycloneState`
  - source snapshot 群
  - refresh / cache metadata
  - projected geometry を状態として持たなくてよい

### 2.11 設定とキャッシュメタデータ

- `ConfigData`
  - `config.json` に保存する共有設定のまとまり
  - 前回地点、window geometry、共有する表示初期値などを保持する
- `GuiLaunchProfile`
  - `zstarview-gui` 用の起動前ダイアログ設定
  - GUI 固有の初期値や選択状態を保持する
- `CacheMetadata`
  - `fetched_at_utc`
  - `source`
  - `source_version` や `release`
  - `ttl_class`
  - `last_checked_at_utc`
  - `stale` / `invalid` の判定に使う属性
- 設定ファイルや cache sidecar の具体的な読み書きタイミングは runtime の責務であり、この文書では保存対象の形だけを扱う

## 3. 入口ごとの状態の違い

### 3.1 `zstarview`

- CLI で受け取った初期値を `ViewerData` と各 overlay state の初期値に変換する。
- 共有設定は `config.json` を使う。

### 3.2 `zstarview-gui`

- 起動前ダイアログの入力値と launch profile を保持する。
- 前回起動値は `gui-launch-profile.json` に保存してよい。

### 3.3 `zstarview-export-image`

- 永続的な UI 状態を最小限にし、1 回の描画に必要な scene data を集めて使う。
- 途中の preview state は持たず、最終画像と stderr 要約を優先する。

## 4. 設計上の境界

- `ViewerData` は scene の中核であり、表示や描画の共通入力に使う。
- `SkyWindowState` は UI 固有の transient state であり、renderer へはできるだけ直接漏らさない。
- `CloudState`、`TerrainHorizonState`、`AircraftState`、`SatelliteState`、`UrbanOutlineState` は、各パイプラインの責務に閉じ込める。
- 入口ごとの違いは「状態の初期化」と「保存先」に寄せ、データモデル自体は揃える。
- task queue、worker、timer、signal handoff、latest-request-wins の実行順は runtime の責務であり、この文書では保持する値だけを定義する。
- `request_id` や `refresh_at` のような更新制御用の値は、必要なら state の一部として持たせてよいが、投入順や破棄条件の説明は runtime に置く。
