# zstarview 設計書

最終更新: 2026-03-21

## 1. この文書の位置づけ

この文書は、`zstarview` の内部設計をまとめたものである。  
アーキテクチャ、モジュール責務、主要データ構造、処理フロー、スレッド分離、外部依存の境界を扱う。  
利用者向けの機能説明は `docs/specification.md`、時系列の作業記録や TODO は `docs/implementation-history.md` を参照する。

## 2. 設計方針

- UI と計算処理を分離する。
- 長時間処理は UI スレッドから切り離す。
- 外部 I/O は局所化し、失敗を UI 側で扱いやすい状態へ正規化する。
- 描画入力はできるだけ前処理済みデータとして渡し、描画側を単純化する。
- CLI 起動時の選択肢と GUI 実行中の状態を同じモデルで扱う。
- 雲、地形地平線、星図更新は相互に独立した補助パイプラインとして扱う。
- 更新頻度の低い補助レイヤーは、連続アニメーションよりも責務分離と UI 安全性を優先する。

## 3. 全体アーキテクチャ

システムは大きく次の層に分かれる。

- 起動・設定層
  - CLI 解析
  - 設定読込
  - 地点解決
  - 初期データロード
- ドメイン計算層
  - 天体位置計算
  - 星カタログ前処理
  - アステリズム、DSO、地形地平線の補助計算
- 描画層
  - 星空ディスク
  - 恒星、惑星、ラベル、補助線描画
  - 雲と地面ティントの合成
- UI 層
  - ウィンドウ管理
  - 入力イベント
  - バックグラウンド更新の起動と反映
- 外部データ連携層
- 衛星クラウドデータ取得
- DEM 取得
- Overture 建物データ取得
- キャッシュ管理

## 4. モジュール構成

### 4.1 起動・設定

- `src/zstarview/zstarview.py`
  - アプリケーションの主エントリポイント
  - 起動シーケンスの組み立て
- `src/zstarview/cli_args.py`
  - CLI オプション定義と値解釈
  - タワー一覧・タワー詳細 JSON 出力の即時終了オプションを扱う
  - `--place`、`--place-countrycode`、`--place-lang` の online 地点検索オプションを扱う
  - parser 構築は `add_location_arguments()`、`add_dataset_query_arguments()`、`add_time_arguments()`、`add_render_arguments()` の helper に分割し、将来の別 CLI からも再利用できるようにする
- `src/zstarview/startup.py`
  - 起動時の地点解決、設定復元、初期値決定
  - Nominatim による online 地点検索と結果正規化を扱う
- `src/zstarview/config.py`
  - ユーザー設定の保存と読込
  - 前回地点を legacy 文字列または構造化地点オブジェクトとして保存する
- `src/zstarview/paths.py`
  - 設定・キャッシュ・データのパス解決

### 4.2 ドメイン計算

- `src/zstarview/astro.py`
  - 恒星、太陽、月、惑星、補助線の計算
  - 可視判定と投影前データ生成
- `src/zstarview/catalog.py`
  - 星カタログの読込と描画用配列の前処理
- `src/zstarview/asterisms.py`
  - アステリズム定義
  - 恒星との対応付け
  - ローテーション選択
- `src/zstarview/tower_viewpoints.py`
  - タワー・展望地点データの解決
  - ASCII フォールバック名を含む名前解決
  - 一覧表示名、全名前一覧、単一タワー JSON 出力用の列挙ヘルパを持つ
- `src/zstarview/mountain_viewpoints.py`
  - 山頂ビューポイントデータの解決
  - ASCII フォールバック名を含む名前解決
  - 一覧表示名、全名前一覧、単一 mountain viewpoint JSON 出力用の列挙ヘルパを持つ
  - 同梱 `mountain_viewpoints.json` は、Wikipedia 起点の候補収集を経て Wikidata メタデータへ正規化した生成物を読む
- `src/zstarview/startup.py`
  - 都市、タワー、山、緯度経度の通常起動時地点解決
  - `--place` 指定時の Nominatim 検索による地点解決
  - タワー/山については最近傍都市からタイムゾーンを補完
  - Nominatim 検索結果についても最近傍都市からタイムゾーンを補完
  - `t/NAME` / `m/NAME` の明示プレフィックスを解釈し、都市名解決より優先する
  - 観測点の基準高さと観測者の目線高さを合成して、実効観測高さを初期化する
- `src/zstarview/types.py`
  - ドメインデータの共有型

### 4.2.1 ビューポイント dataset 参照専用 CLI

同梱 `tower_viewpoints.json` と `mountain_viewpoints.json` を直接参照する軽量 CLI 経路を持つ。

- 対象オプション
  - `--list-viewpoints KIND`
  - `--list-viewpoint-names KIND`
  - `--show-viewpoint-json NAME`
- この経路では `startup.py` の都市解決や GUI 初期化へ進まず、`zstarview.py` が `tower_viewpoints.py` / `mountain_viewpoints.py` を直接呼び出して標準出力を書き、即時終了する。
- 一覧出力はローカル JSON のみを参照し、GeoNames、設定保存、ネットワークアクセスには触れない。
- `KIND=t` は tower dataset、`KIND=m` は mountain dataset を選ぶ。
- 一覧出力は `t/NAME` / `m/NAME` 形式で prefix 付き表示を返す。
- `--show-viewpoint-json NAME` は、prefix 付き入力なら該当 kind の resolver だけを使う。
- prefix なし入力では tower / mountain の両方を試し、両方に exact match がある場合は曖昧一致エラーとして候補名を列挙する。
- `tower_viewpoints.py` は dataset の `name` を保持したまま、必要に応じて ASCII フォールバック名 `ascii_name` を算出する。
- `mountain_viewpoints.py` も同様に dataset の `name` を保持したまま、必要に応じて ASCII フォールバック名 `ascii_name` を算出する。
- `--list-viewpoints` は `ascii_name` がある場合それを表示名として優先する。
- `--list-viewpoint-names` は元の綴りと ASCII フォールバック綴りの両方を含む。
- mountain viewpoint dataset の生成は、Wikipedia 候補抽出、Wikidata raw query、review JSON、curated seed、最終 JSON という段階を持つ。詳細は `docs/developer/viewpoint-dataset-generation.md` を参照する。
- viewpoint dataset の別名クリーニング基準
  - ASCII フォールバックは、ダイアクリティカルマーク除去後に英字を含む場合だけ採用する。
  - 記号だけ、数字だけ、絵文字だけになる fallback は捨てる。
  - mountain の `names` では、長文フレーズや明らかなノイズ値を捨てる。
  - mountain の `names` では、`Mt.` / `Mtn.` を `Mount` / `Mountain` に展開できる場合は展開側へ寄せ、重複側を捨てる。
  - mountain の `names` では、`/` の前後に空白を含む表記は候補一覧から外し、必要なら空白なし表記だけを残す。
  - tower の短縮 alias は、固有名詞として成立するものだけを人手で追加する。
- tower / mountain とも、数字だけの部分一致では解決しない。

### 4.2.2 `--place` による Nominatim 検索起動

`--place QUERY` は、既存の `location` 引数とは別の明示的な online resolver 経路を持つ。

- 対象オプション
  - `--place QUERY`
  - `--place-countrycode CODE`
  - `--place-lang LANG`
- `--place` は位置引数 `location` と排他とする。
- `--place` が指定されない通常起動では、既存の offline-first 解決規則を維持する。
- `--place` が指定された場合、`startup.py` は Nominatim Search API を 1 回呼び出し、返却された候補を正規化する。
- 候補の正規化では少なくとも `name`、`lat`、`lon`、`category`、`type`、`importance` を保持する。
- 候補は `importance` 降順で扱い、先頭候補を起動地点として採用する。
- 複数候補が得られた場合でも起動時の対話選択は行わず、候補一覧を標準エラーまたは logger 出力へ流しつつ先頭候補を採用する。
- `name` には Nominatim の `display_name` を保持し、GUI の地点表示にはその全体をそのまま使う。
- タイムゾーンは、採用した座標から最近傍の GeoNames 都市を引いて補完する。補完できない場合は UTC を使う。
- HTTP エラー、レート制限、通信失敗、JSON 解析失敗、0 件結果は `StartupAbortError` 相当で起動中断とし、logger 経由でターミナルとスプラッシュへ表示する。
- Nominatim 利用は起動時の単発検索に限定し、候補列挙だけの反復照会経路は持たない。

### 4.2.3 前回地点の保存形式

前回地点保存は後方互換のため複数形式を許容する。

- legacy 形式
  - 既存どおり文字列 1 個を保存する。
  - 都市名、`CC/City`、`t/NAME`、`m/NAME`、`lat;lon` などを表す。
- Nominatim 形式
  - 構造化オブジェクトとして保存する。
  - 少なくとも resolver 種別、元の query、採用した正規化結果 JSON を含める。
  - 次回起動時は再検索せず、この保存済み結果をそのまま `ResolvedLocation` へ復元する。
- `config.py` は `str | object` を読めるようにし、未知形式は安全側で無視または起動失敗へ正規化する。
- 保存形式が異なっても、起動後に UI へ渡す地点モデルは共通の `ResolvedLocation` に揃える。

### 4.2.4 単発画像書き出し CLI

GUI 常駐とは別に、1 枚の画像を書き出して終了する headless CLI 経路を持つ。

- エントリポイントは `zstarview-export-image` とする。
- この経路は `zstarview.py` の GUI `main()` と別の `main` を持つ。
- parser は `cli_args.py` の helper 群を組み合わせて構築し、地点、時刻、視線、描画オプションは通常 CLI と共有する。
- 画像書き出し CLI 固有オプションは少なくとも次を想定する。
  - `--output`
  - `--image-size`
  - `--layer-timeout-seconds`
  - `--allow-partial-data`
- `--window-geometry` と `--sky-update-interval` は GUI 専用とし、この CLI では parser に載せない。
- dataset 参照専用オプション群も、この CLI では parser に載せない。
- 地点解決と天体計算は GUI と同じ下位ロジックを共有するが、レイヤー取得順序は GUI と共有しない。
- 出力画像には既定で hover/HUD を含めず、`RenderSceneData` と `RenderStyle` を中心にベース描画だけを行う想定とする。
- `guide` はベース描画に含めてよい。
- 外部依存レイヤーの取得は逐次でも並列でもよいが、CLI 側では「いつまで待つか」と「部分データを許容するか」を引数で決められるようにする。
- 既定は安全側として「部分データは保存しない」とし、明示的に `--allow-partial-data` を指定したときだけ部分出力を許可する。
- `opacity == 0` で無効化されたレイヤーは、取得キュー自体に積まず、layer timeout の待機対象からも外す。
- 実装では `SkyWindow` と GUI controller 群には依存せず、sky/cloud/terrain/urban/aircraft を同期的に順番に取得してから、shared pipeline で `QImage` へ 1 回だけ描画して保存する。
- Qt はフォント読込と `QImage` / `QPainter` 利用のためだけに初期化し、CLI 側ではバックグラウンド worker や signal ベースの寿命管理を持たない。
- `ui/sky_worker.py` の celestial / sky-disc 計算は pure helper `compute_sky_snapshot()` として切り出し、GUI worker と export CLI の両方から共有する。

### 4.3 描画

- `src/zstarview/render/draw.py`
  - 恒星、惑星、ラベル、補助線、地平線関連の描画
  - アステリズム線は大円弧をサンプルし、アステリズム専用の広い FOV 境界で円形クリップして描画
  - 天の赤道、黄道、地平線は `(alt, az)` サンプル列から描画時に `render_view_center` 基準で投影する
- `src/zstarview/render/draw_sky_disc.py`
  - sky color disc の生成
- `src/zstarview/ui/composite.py`
  - 星空、雲、欠損ティント、地面色の合成

#### 4.3.1 描画リファクタリング方針

- 将来の単発画像書き出し CLI に備え、描画処理は `SkyWindow` / Mixin 直結から、関数呼び出し中心の再利用可能なパイプラインへ寄せる方針とする。
- 共通化の中心は「最終的に何をどう描くか」であり、「どの順序でデータを取得するか」ではない。
- GUI 経路と将来の CLI 経路は、別々のオーケストレーションを持ってよい。
  - GUI 経路は、従来どおり星空の初期表示を優先し、雲・地形地平線・都市アウトライン・航空機を後段で非同期に反映してよい。
  - CLI 経路は、必要なレイヤーを逐次取得して 1 回だけ描画する構成を許容する。
- このため、描画層では少なくとも次の境界を意識する。
  - 描画入力: 観測地点、時刻、視線方向、画面サイズ、表示オプション
  - レイヤー入力: celestial、sky disc、cloud、terrain horizon、urban outline、aircraft overlay などの描画用データ
  - 描画関数群: `QPainter` または `QImage` に対して各レイヤーを順に描く純粋寄りの関数
- `ui/window_render.py` は最終的に、hover 判定、jump highlight、frame cache、interaction mode など GUI 固有の責務に寄せ、描画本体は薄いラッパに縮小するのが望ましい。

#### 4.3.2 現在の描画パイプライン到達点

- 共有描画本体は `src/zstarview/render/pipeline.py` に置く。
- 共有描画入力は次の 3 つに分ける。
  - `RenderSceneData`
  - `RenderStyle`
  - `RenderHudState`
- `render_scene_into_painter()` と下位の `draw_*` 関数群は、`geometry`、`viewport_rect`、`scene`、`style`、`hud` を明示的に受ける。
- `RenderPipelineState` のような中間ラッパ型は廃止し、shared pipeline 側では直接引数で依存関係を表す。
- `ui/window_render.py` は、`paintEvent()` 本線、scene/style/hud の組み立て、frame cache、jump highlight、hover 解決など GUI 固有処理に絞る。
- 現在の通常描画順は概ね次のとおり。
  - `background`
  - `sky-cloud`
  - `guide`
  - `terrain`
  - `stars`
  - `aircraft`
  - `planets`
  - `overlay`
  - `labels`
  - `hover-overlay`
  - `status`
- `guide` は方位ラベルと天頂マーカーを含む独立レイヤーであり、空色・雲合成の上、通常の hover/HUD オーバーレイより手前に置く。

#### 4.3.3 次段のリファクタリング方針

- 次段の主目的は、hover/HUD とベース描画の分離である。
- `guide` レイヤーは HUD ではなくベース描画側に残す。
- ベース描画には、少なくとも `background`、`sky-cloud`、`guide`、`terrain`、`stars`、`aircraft`、`planets`、`labels` を含める想定とする。
- HUD 側には、少なくとも次を寄せる方向で整理する。
  - 恒星ホバー
  - DSO ホバー
  - jump highlight
  - status line
- `paintEvent()` は最終的に「ベースフレームをキャッシュし、その上に hover/HUD を都度重ねる」形を目指す。
- この変更により、ベースフレーム cache key から `mouse_pos` や hover 対象名などの高頻度変化要素を外し、キャッシュ効率を改善する。

#### 4.3.4 hover/HUD 分離の現在位置

- shared pipeline は、`render_base_scene_into_painter()` と `render_hud_overlay_into_painter()` に分かれている。
- ベース描画は、`background`、`sky-cloud`、`guide`、`terrain`、`stars`、`aircraft`、`planets`、静的 overlay 情報、`labels` を担当する。
- HUD 描画は、少なくとも次を担当する。
  - アステリズムの hover 強調
  - 月 hover 時の拡大上書き
  - DSO hover
  - 恒星・惑星 hover 情報
  - jump highlight
  - status line
- 月の `5x` 拡大は、角半径の生値ではなく「通常時の見た目半径」を基準に適用する。
- `guide` レイヤーはベース側に残し、マウス位置によるラベル回避には依存しない安定描画として扱う。
- `ui/window_render.py` の frame cache はベース描画だけを保持し、hover/jump/status はキャッシュ後に都度上書きする。
- これにより、frame cache key から `mouse_pos`、hover 対象名、jump highlight 名、status message を外している。

### 4.4 UI

- `src/zstarview/ui/window.py`
  - メインウィンドウ
  - UI 状態と更新制御の集約点
- `src/zstarview/ui/window_state.py`
  - 画面状態の保持
- `src/zstarview/ui/window_inputs.py`
  - ユーザー指定値と実行時オプションの正規化
- `src/zstarview/ui/window_render.py`
  - 再描画とレンダリング関連の UI ロジック
  - 恒星レイヤは描画時に現在のウィンドウサイズから内部レンダリング面サイズを再計算する
  - 天球ディスク幅が `expected-render-width` 以下なら等倍描画し、それを超える場合は `expected-render-width * sqrt(disc_width / expected-render-width)` に従って内部描画面を縮小する
  - 縮小時は低解像度 `QImage` に恒星を描いてからウィンドウ全体へ拡大転写し、大型ウィンドウでの負荷を抑える
  - 最終フレームの `QImage` キャッシュを持ち、geometry、描画入力、hover/ハイライト状態、ステータスメッセージ、interaction mode が不変なら前回フレームをそのまま再利用する
- `src/zstarview/ui/window_updates.py`
  - バックグラウンド更新結果の反映
- `src/zstarview/ui/sky_worker.py`
  - 星空計算のバックグラウンド実行
- `src/zstarview/ui/famous_star_dialog.py`
  - 代表恒星ジャンプ UI
- `src/zstarview/ui/famous_star_search_dialog.py`
  - 星・アステリズム検索 UI
- `src/zstarview/ui/famous_star_shortcuts.py`
  - ジャンプ・検索用データの整形

### 4.5 雲データ処理

- `src/zstarview/ui/cloud_controller.py`
  - 雲更新の実行制御
  - キューイング
  - latest-request-wins の適用
- `src/zstarview/ui/cloud_state.py`
  - 雲画像、バナー、欠損状態の保持
- `src/zstarview/clouddisc/core.py`
  - クラウドディスク生成のオーケストレーション
  - `CloudDiscConfig.alt_min_deg` による可視高度下限の適用
- `src/zstarview/clouddisc/providers/*.py`
  - 衛星データ取得
- `src/zstarview/clouddisc/projectors/az.py`
  - 空ディスク向け投影
  - 既定 `alt_min_deg = 3.0°` により、地平線近傍の遠距離雲を可視マスクから外す
- `src/zstarview/clouddisc/render/grayscale.py`
  - 雲画像生成
- `src/zstarview/clouddisc/cache/*.py`
  - キャッシュ保存と清掃

### 4.6 地形地平線処理

- `src/zstarview/ui/terrain_controller.py`
  - 地形地平線更新の実行制御
- `src/zstarview/ui/terrain_state.py`
  - 地形地平線の状態保持
- `src/zstarview/terrain/dem.py`
  - DEM 取得、読込、サンプリング
- `src/zstarview/terrain/horizon.py`
  - 方位ごとの見かけ地平線計算

### 4.7 航空機オーバーレイ処理

- `src/zstarview/ui/aircraft_controller.py`
  - 航空機更新の実行制御
  - `5分` タイマー、明示更新要求、latest-request-wins の適用
  - OpenSky 取得結果の正規化と UI 反映
- `src/zstarview/aircraft_constants.py`
  - 航空機更新間隔と `bbox` 既定値の共有定数
  - `5分` 取得、`2秒` 予想再投影、bbox / fade / trail span 定数を定義
- `src/zstarview/ui/aircraft_state.py`
  - 最新スナップショットの航空機一覧
  - 最新の描画用折れ線データ
  - 最終成功時刻による cache age 判定
  - 読込中 / 取得中バナー
  - 失敗表示状態
  - 最終成功時刻
- `src/zstarview/aircraft/opensky.py`
  - OpenSky `states/all` へのアクセス
  - 観測地点由来 `bbox` の組み立て
  - 緯度依存の east-west coverage 調整と 1-credit 帯面積クリップ
  - 生レスポンス配列を名前付き内部モデルへ正規化
- `src/zstarview/aircraft/project.py`
  - 航空機の `lat/lon/alt` を観測地点基準の `alt/az` へ変換
  - `velocity` / `heading` / `vertical_rate` による短時間前進予測
  - `4秒前 -> 現在 -> 4秒後` を `2秒` 刻みでサンプリングした折れ線点列と age-based alpha の算出
- `src/zstarview/aircraft/types.py`
  - OpenSky state vector の正規化モデル
  - UI が保持する描画用航空機折れ線モデル

### 4.8 ユーティリティ

- `src/zstarview/utils/resolve_city.py`
  - 都市解決補助
- `src/zstarview/utils/timezone_parser.py`
  - `--datetime` 文字列のタイムゾーン解釈
- `src/zstarview/utils/image.py`
  - 画像変換補助
- `src/zstarview/utils/qt.py`
  - Qt 補助
- `src/zstarview/data/import_overture_buildings.py`
  - `lat/lon + radius` に対応する bbox を計算する
  - `overturemaps download` を呼び、`building` および必要に応じて `building_part` を取得する
  - ダウンロード結果を既存の derived tile JSON と `tile_index.json` 形式へ変換する
  - 出力先既定値は `CACHE_PATH/overture_buildings`

## 5. 主要データ構造

### 5.1 地点・視点関連

- `ViewerData`
  - 観測地点
  - タイムゾーン
  - 表示中心の方位・高度
  - 観測者の目線高さ
  - 地点表示用の補足高さラベルと値
  - 観測者高さを UI 表示するかどうかのフラグ
  - 画面描画に必要な視点情報

地点 dataset が持つ高さ情報と `ViewerData.observer_height_m` は別概念として扱う。

- mountain viewpoint
  - dataset 側の高さは山頂ビューポイントの海抜標高 `Elevation`
- tower viewpoint
  - dataset 側の高さは地表からのタワー高またはビューポイント高 `Tower height`
- `ViewerData.observer_height_m`
  - どの入力種別でも、基準観測点から観測者の目線までの高さ
  - CLI `--observer-height-m` はこの値だけを置き換える
  - 既定値は `1.7m`

### 5.2 天体計算結果

- `CelestialData`
  - 恒星の描画用データ
  - 太陽、月、惑星の位置
  - 地平線、黄道、赤道などの線データ
  - DSO、ホバー対象、補助表示に必要な情報

### 5.3 航空機関連

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
  - 描画直前まで絞り込んだ軽量折れ線モデル
  - 折れ線は既定で `-4, -2, 0, +2, +4秒` の `alt/az` サンプルを保持する
- `AircraftState`
  - 最新スナップショットの機体一覧
  - 表示用折れ線列
  - 最終成功タイムスタンプ
  - 読込中状態
  - エラーバナー
  - 最終成功時刻
  - snapshot は `5分` 単位で更新し、折れ線データは `2秒` ごとに再投影してよい

### 5.4 雲関連

- `CloudMeta`
  - 雲画像のデータ元メタ情報
- `CloudState`
  - 最新画像
  - 欠損マスク
  - 進行中状態
  - エラーバナー

雲パイプラインでは 2 種類のキーを分離して扱う。

- `SourceKey`
  - 衛星種別や時間スロットなど、取得元を識別するキー
- `RenderKey`
  - `SourceKey` に視点条件を加えた描画条件キー

この分離により、ソース取得をやり直さずに視点変更のみ再描画できる。

### 5.5 地形地平線関連

- `TerrainHorizonState`
  - 地形地平線の点列
  - 読込中状態
  - エラー表示状態
- 地形地平線プロファイル
  - `(altitude_deg, azimuth_deg)` の系列として保持する
  - 地点依存、時刻非依存のデータとして扱う

### 5.6 ウィンドウ状態

- `SkyWindowState`
  - 現在の視点
  - 直近の描画用視点
  - `viewport_interaction_mode` による簡易描画状態
  - `viewport_interaction_stars` による簡易描画用の明るい星テーブル
  - ホバー対象
  - ハイライト対象
  - 各更新パイプラインの UI 反映状態
- `SkyWindow._frame_cache_image`
  - `paintEvent` の最終出力を保持する `QImage`
  - 同一フレームの再利用に使う
- `SkyWindow._frame_cache_key`
  - 最終フレームキャッシュの無効化条件をまとめたキー
  - ウィンドウサイズ、`render_view_center`、描画トグル、`CelestialData`、空ディスク画像、雲画像、地形/都市アウトライン、hover、jump highlight、ステータス行文言などを含む
- `UrbanOutlineState`
  - 都市アウトライン点列
  - 読込中 / 取得中バナー
  - 失敗表示状態
  - `cache` または `overture` などの現在ソース表示
- `AircraftState`
  - 航空機描画折れ線列
  - 読込中 / 取得中バナー
  - 失敗表示状態
  - 最終成功時刻

## 6. 処理フロー

### 6.1 起動フロー

1. `zstarview.py` が CLI とログを初期化する。
2. 設定ファイルから前回の地点やウィンドウ状態を復元する。
3. 入力を、`--place` による online 検索地点または通常の都市・タワー・山・座標として解決する。
4. 星カタログや補助データを読み込む。
5. `SkyWindow` を生成し、初回描画を行う。
6. 必要に応じて雲更新、地形地平線更新、都市アウトライン更新、航空機更新をバックグラウンドで開始する。

### 6.2 星空更新フロー

1. UI 操作またはタイマーで再計算要求が発生する。
2. `SkyDataWorker` がバックグラウンドで天体計算と sky disc 生成を行う。
3. 計算結果を `CelestialData` と描画補助データとして UI へ返す。
4. `SkyWindow` が内部状態を更新し、再描画する。

視線変更とリサイズの連続入力時は例外的に `viewport_interaction_mode` を使う。

1. `render_view_center` を即時更新する。
2. 明るい星 (`vmag <= 4.0`) のみ同期的に再計算し、簡易描画に使う。
3. 補助線と地形地平線は、保持済みデータを `render_view_center` で再投影して追随させる。
4. 最後の入力から 0.7 秒経過後に通常更新を 1 回だけ開始する。

### 6.2.1 最終フレームキャッシュ

Qt はメニュー操作やボタン状態変化でも `paintEvent` を再発行するため、空の内容が変わらない repaint では重い描画をやり直さないようにする。

1. `window_render.py` は、geometry、描画入力、hover、jump highlight、interaction mode、ステータス行文言などからフレームキーを作る。
2. キーが前回と同じなら、保持済みの最終フレーム `QImage` を `drawImage()` するだけで終了する。
3. キーが変わったときだけ、背景、空ディスク、雲、地形、都市アウトライン、星、惑星、オーバーレイ、ラベル、ステータス行を一時 `QImage` に順に描いて新しいフレームとして保存する。

この最終フレームキャッシュは `SkyCompositorCache` の上位にある。`SkyCompositorCache` は sky/cloud 合成だけを再利用し、その出力を含むウィンドウ全体の描画結果をさらに `window_render.py` 側でキャッシュする二段構えにしている。

### 6.3 雲更新フロー

1. `CloudController` が新しい要求を受け取る。
2. ソース取得が必要か、既存ソースで再描画可能かを判定する。
3. 必要なら衛星データを取得し、CloudDisc 用ソースを生成する。
4. 視点条件に応じて雲画像を描画する。
5. `request_id` により古い結果を破棄し、最新結果のみ UI に反映する。
6. 欠損領域がある場合は欠損マスクも渡す。

### 6.4 地形地平線更新フロー

1. `TerrainHorizonController` が地点に応じた DEM 更新要求を受ける。
2. 必要な DEM タイルを取得またはキャッシュから読込する。
3. 方位ごとに見かけ地平線を計算する。
4. 結果を地形地平線プロファイルとして UI に返す。
5. 画面再投影時は既存プロファイルを使い回し、再取得はしない。

### 6.5 描画フロー

1. sky disc を生成する。
2. 恒星、惑星、月、補助線を重ねる。
3. 地形地平線があれば地平線関連描画へ反映する。
4. 都市アウトラインがあれば白線オーバーレイとして描画する。
5. 航空機オーバーレイがあればオレンジの折れ線オーバーレイとして描画する。
6. 雲画像と欠損ティントを合成する。
7. ラベル、オーバーレイ、ステータス行を描画する。

### 6.6 都市アウトライン更新フロー

1. `SkyWindow` が起動時またはトグル再有効化時に `UrbanOutlineController` へ更新要求を出す。
2. `UrbanOutlineController` は `lat/lon + radius + min_height + mode` から実行キーを作る。既定 mode は `both`、既定半径は `2.5km` である。
3. `mode=both` の場合、`UrbanOutlineController` は `building` 用と `building_part` 用の derived dataset をそれぞれ確認し、欠けている方だけを取得する。
4. `mode=building` の場合は `building` 用 derived dataset のみを確認・取得する。
5. 同じ `update()` サイクルの中で、`UrbanOutlineController` は `src/zstarview/data/skyscraper_tiles_z14.json` を読み、視点中心 `10km` 円に交差しつつ内側半径 `radius_km` だけには収まらない seed tile を選ぶ。既定 `radius_km` は `2.5km` である。
6. 選ばれた skyscraper seed tile については、専用 cache root 配下の derived dataset を確認し、未取得 tile だけを `overturemaps download --bbox=... --type building` で取得する。
7. skyscraper tile は import 時に `height_m >= 150` で前処理され、runtime では `resolve_urban_outline_layer_for_viewer(..., radius_km=10.0, min_distance_km=radius_km)` として読む。
8. runtime 側は通常 derived dataset 群を追加の高さフィルタなしで読む。skyscraper derived dataset 群は cache 自体は常に `150m` 下限で共有しつつ、`min_building_height_m > 150` の場合のみ runtime 側で追加高さフィルタをかけてよい。
9. runtime マージ時には、通常レイヤー側で `building_part` が持つ `parent_building_id` を参照し、対応する親 `building` 外形を除外する。
10. `compute_urban_outlines()` は建物ごとの `height_m` を保持した輪郭列を返し、`resolve_urban_outline_layer_for_viewer()` はそれを `UrbanOutlinePolyline` の列に変換する。
11. `UrbanOutlineController` は通常レイヤーと skyscraper レイヤーをマージして 1 回の `urban_ready` として反映する。skyscraper 取得が失敗した場合は、通常レイヤーだけで `urban_ready` してよい。
12. `--urban-outline-skyscraper-only` 指定時は、通常近距離 derived dataset の確認・取得・解決をスキップし、skyscraper レイヤーだけを解決する。
13. 描画時は `50m` 以上を CLI 指定 opacity の基準とし、`0m` ではその `25%` になるよう高さ比例で alpha を下げる。
14. 結果の outline 列は `UrbanOutlineState` と `SkyWindowState.urban_outlines` に反映し、再描画する。
15. 取得中や失敗時はバナー文字列を UI 状態へ反映する。

補足:
- 旧 `list[list[(alt, az)]]` 形式の runtime 互換コードは削除し、都市アウトライン描画は `UrbanOutlinePolyline` のみを受け付ける。

### 6.7 航空機オーバーレイ更新フロー

1. `SkyWindow` が起動時に `AircraftController` を生成し、初回更新要求を出す。
2. `SkyWindowUserOptions` は `aircraft_opacity` と `aircraft_gui_allowed` を持ち、`--aircraft-opacity 0.0` のときは起動直後から航空機問い合わせを行わない。
3. `AircraftController` は観測地点の `lat/lon` から OpenSky 用 `bbox` を作る。南北は既定 `±1.0°`、東西は緯度に応じて最低 `90km` を確保しつつ、面積は `25 square degrees` 以下へ抑える。
4. `AircraftController` は明示更新要求または次回 fetch timer に従い、OpenSky 取得を 1 回だけ開始する。
5. `aircraft/opensky.py` は `states/all?lamin=...&lamax=...&lomin=...&lomax=...` を呼び、state vector 配列を `AircraftSnapshot` 列へ正規化する。
6. 正規化時には `lat/lon` 欠損、`on_ground=true`、極端に低速な機体を落としてよい。
7. `aircraft/project.py` は各 `AircraftSnapshot` の `velocity`、`heading`、`vertical_rate`、`last_contact` を使って短時間前進予測し、`2秒前 -> 現在 -> 2秒後` の折れ線端点を含む `AircraftOverlayPoint` を作る。
8. `aircraft/project.py` は age に応じた alpha scale も計算し、`90秒` を超えた機体が次回取得まで徐々に薄くなるようにする。
9. `window.py` は API 取得とは別に `AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS` の UI タイマーを持ち、既定では `2秒` ごとに保持済み snapshot から折れ線データだけを再投影してよい。
10. `window.py` は fetch timer を single-shot で扱い、レイヤー再表示時には `AircraftState.last_success_utc` から cache age を計算して次回 API 呼び出し時刻を再調整してよい。
11. 保持済み snapshot が新しければ、GUI で航空機レイヤーを再表示しても API を再問い合わせせず、再投影だけで即時復帰してよい。
12. 保持済み snapshot が更新間隔を超えて古いときだけ、レイヤー再表示時に即時再取得してよい。
13. 描画時は観測地点から `50km` を超える機体を落とし、`10km` 以内かつ `90秒` 以内の機体だけに `callsign` を付けてよい。
14. 新しい取得が成功したときだけ `AircraftState.snapshots` と `AircraftState.overlay_points` をまとめて置き換える。
15. 取得失敗時は直前成功スナップショットを保持したまま、`AircraftState` の失敗表示だけを更新する。
16. 古い非同期結果で UI を巻き戻さないため、`AircraftController` も request id ベースの latest-request-wins を使ってよい。
17. `SkyWindow` は `aircraft_ready` または予想再投影完了を受けたら `SkyWindowState` を更新し、再描画する。
18. `viewport_interaction_mode` 中は航空機オーバーレイ描画を抑止してよい。モード終了後に保持済み `overlay_points` で通常描画へ戻る。

## 7. スレッドモデル

- UI スレッド
  - Qt イベントループ
  - 描画
  - メニュー、入力、状態反映
- バックグラウンドワーカー
  - 星空計算
  - 雲データ取得と描画
  - 地形地平線計算
  - 都市アウトライン取得と outline 生成
  - 航空機データ取得と可視折れ線生成
  - キャッシュ清掃の補助処理

設計上の原則は次の通り。

- UI オブジェクト更新は必ず UI スレッドで行う。
- バックグラウンド結果はシグナルで UI スレッドへ戻す。
- 終了開始後は新規要求を止め、破棄済み UI への通知を避ける。

## 8. 状態管理

### 8.1 状態の分離

- 長寿命の UI 状態は `SkyWindow` と `SkyWindowState` が持つ。
- 雲専用状態は `CloudState` に分離する。
- 地形地平線専用状態は `TerrainHorizonState` に分離する。
- 航空機専用状態は `AircraftState` に分離する。

### 8.2 latest-request-wins

視点変更やタイマー更新が短時間に連続した場合、古い非同期結果で UI を巻き戻さないことを優先する。  
そのため、雲更新や航空機更新では最新要求のみが採用される。

### 8.3 CLI と GUI の整合

- CLI で与えられた初期オプションは、起動後の GUI 状態の初期値になる。
- ただし CLI で明示的に無効化した機能の一部は、セッション中の GUI 再有効化を禁止する。
- 代表例が `--terrain-horizon-opacity 0` による地形地平線ロックアウトである。
- `--sky-opacity 0`、`--cloud-opacity 0`、`--terrain-horizon-opacity 0`、`--urban-outline-opacity 0` は、そのセッションで各 GUI トグルをロックアウトする。

### 8.4 航空機オーバーレイの更新粒度

- 航空機観測値の取得は `5分` 間隔とし、表示上の予想再投影は `2秒` 間隔で行ってよい。
- `aircraft_opacity <= 0.0` の間は、fetch timer と予想再投影 timer を止めてよい。
- 折れ線は `2秒前 -> 現在 -> 2秒後` の 3 点から構成し、機体の短時間進行方向を示す。
- 線幅は観測地点に近い機体ほど太く、遠い機体ほど細くしてよい。
- 描画は観測地点から `50km` 以内に限定し、近距離機体だけにコールサインを出してよい。
- 取得時刻から `90秒` を超えたら alpha を段階的に下げ、次回 API 更新まで視覚的に古さを示してよい。
- GUI から再表示したときは `last_success_utc` を見て fresh cache を優先し、不要な OpenSky 再問い合わせを避けてよい。

### 8.5 都市アウトライン描画の簡略化

- 都市アウトラインは derived tile から得た建物上端輪郭を描く。
- derived tile の各建物レコードは `building_id` を持ち、`building_part` 由来レコードは必要に応じて `parent_building_id` を持つ。
- `building` と `building_part` を併用する場合、`parent_building_id` を持つ part が存在する親 `building` は描画対象から外す。
- 建物高さのしきい値は derived tile 生成時の前処理パラメータとし、runtime 読込時の既定値では再適用しない。
- 遠距離スカイスクレーパー補助レイヤーは `building_part` を使わず `building` のみを扱い、runtime では `min_distance_km` を使って通常レイヤー半径より内側の建物を落とす。
- `UrbanOutlinePolyline` は `source` を持ち、通常レイヤーと skyscraper レイヤーの由来を区別できる。ただし現行描画色は共通である。
- ただし見かけの方位幅が `0.5°` 未満の輪郭は、細い polyline ではなく太い水平線に簡略化する。
- `viewport_interaction_mode` 中は都市アウトライン描画を抑止し、方向キー操作の負荷を下げる。
- ハンバーガーメニュー表示そのものでは `viewport_interaction_mode` へ入らない。メニュー起動で `paintEvent` が走っても、最終フレームキャッシュの再利用で余分な重い描画を避ける。

## 9. エラーモデル

### 9.1 起動エラー

- 必須データの欠落や起動時解決不能は起動中断対象とする。
- 起動中断は例外または明示的な abort 制御で扱う。

### 9.2 補助機能エラー

- 雲取得失敗は雲機能内で閉じる。
- 地形地平線取得失敗は terrain 機能内で閉じる。
- 航空機取得失敗は aircraft 機能内で閉じる。
- 本体表示を継続し、UI には状態表示のみ反映する。

### 9.3 再試行方針

- 雲は更新タイミングごとに再取得の機会がある。
- 地形地平線は同一セッション中の自動再試行を抑制する。
- 航空機は `5分` タイマーごとに再取得の機会がある。
- 航空機の予想再投影はセッション内 `2秒` タイマーで繰り返してよい。

## 10. データとキャッシュ

### 10.1 同梱データ

- 星カタログ
- DSO データ
- フォント
- タワー・展望地点データ
- スカイスクレーパー seed tile リスト
- 画像アセット

### 10.2 外部取得データ

- 衛星クラウドデータ
- DEM タイル
- Overture 建物データ
- OpenSky 航空機 state vector

### 10.3 キャッシュ方針

- キャッシュ対象は再利用価値の高い外部取得データとする。
- DEM データは永続キャッシュする。
- Overture 建物由来の derived tile は地点条件ごとに永続キャッシュする。
- 遠距離スカイスクレーパー補助レイヤー用 derived tile は `overture_skyscrapers/<tile-cache-key>/bldg` 配下に tile 単位で永続キャッシュする。
- 地形地平線の計算済みポリラインは永続化しない。
- 雲は取得ソースと中間成果物をキャッシュし、視点変更時の再利用を優先する。
- 初期実装の航空機スナップショットはセッション内メモリだけに保持し、永続キャッシュしない。

## 11. テスト観点と設計上の分離

テストしやすさのため、次を重視する。

- 数学・投影・判定ロジックは純粋関数に寄せる。
- UI 依存を持つ処理は state/controller に閉じ込める。
- 外部 I/O は provider/controller 経由に限定し、モック可能にする。

代表的なテスト対象は次の通り。

- 地点解決
- タイムゾーン解釈
- 星カタログ前処理
- 惑星・恒星描画ルール
- 雲更新キューイング
- 地形地平線計算
- 航空機 state vector 正規化
- 観測地点由来 `bbox` の生成
- 航空機の短時間前進予測
- 航空機の地平線上 / 視野内判定
- age-based alpha と折れ線端点の計算
- GUI トグルの状態遷移

## 12. 外部依存

- GUI: `PySide6`
- 天文計算: `astropy`, `skyfield`
- 数値計算: `numpy`, `polars`
- 雲データ処理: `xarray`, `satpy`, `boto3`, `botocore`
- 地形データ処理: `rasterio`, `pyproj`
- 画像処理: `Pillow`

## 13. 文書間の責務分担

- `docs/specification.md`
  - 利用者から見た機能と動作
- `docs/design.md`
  - 内部構造、責務分割、データ構造、処理フロー
- `docs/implementation-history.md`
  - 実装の時系列メモ、TODO、INPROGRESS、変更の履歴
