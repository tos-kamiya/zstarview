# zstarview 設計書

最終更新: 2026-03-14

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

### 4.3 描画

- `src/zstarview/render/draw.py`
  - 恒星、惑星、ラベル、補助線、地平線関連の描画
  - アステリズム線は大円弧をサンプルし、アステリズム専用の広い FOV 境界で円形クリップして描画
  - 天の赤道、黄道、地平線は `(alt, az)` サンプル列から描画時に `render_view_center` 基準で投影する
- `src/zstarview/render/draw_sky_disc.py`
  - sky color disc の生成
- `src/zstarview/ui/composite.py`
  - 星空、雲、欠損ティント、地面色の合成

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
- `src/zstarview/clouddisc/providers/*.py`
  - 衛星データ取得
- `src/zstarview/clouddisc/projectors/az.py`
  - 空ディスク向け投影
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

### 4.7 ユーティリティ

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
  - `overturemaps download` を呼び、`building` または `building_part` を取得する
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

### 5.3 雲関連

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

### 5.4 地形地平線関連

- `TerrainHorizonState`
  - 地形地平線の点列
  - 読込中状態
  - エラー表示状態
- 地形地平線プロファイル
  - `(altitude_deg, azimuth_deg)` の系列として保持する
  - 地点依存、時刻非依存のデータとして扱う

### 5.5 ウィンドウ状態

- `SkyWindowState`
  - 現在の視点
  - 直近の描画用視点
  - `viewport_interaction_mode` による簡易描画状態
  - `viewport_interaction_stars` による簡易描画用の明るい星テーブル
  - ホバー対象
  - ハイライト対象
  - 各更新パイプラインの UI 反映状態
- `UrbanOutlineState`
  - 都市アウトライン点列
  - 読込中 / 取得中バナー
  - 失敗表示状態
  - `cache` または `overture` などの現在ソース表示

## 6. 処理フロー

### 6.1 起動フロー

1. `zstarview.py` が CLI とログを初期化する。
2. 設定ファイルから前回の地点やウィンドウ状態を復元する。
3. 入力を、`--place` による online 検索地点または通常の都市・タワー・山・座標として解決する。
4. 星カタログや補助データを読み込む。
5. `SkyWindow` を生成し、初回描画を行う。
6. 必要に応じて雲更新、地形地平線更新、都市アウトライン更新をバックグラウンドで開始する。

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
5. 雲画像と欠損ティントを合成する。
6. ラベル、オーバーレイ、ステータス行を描画する。

### 6.6 都市アウトライン更新フロー

1. `SkyWindow` が起動時またはトグル再有効化時に `UrbanOutlineController` へ更新要求を出す。
2. `UrbanOutlineController` は `lat/lon + radius + min_height + feature_type` からキャッシュキーを作る。既定半径は `2.5km` である。
3. 対応する derived dataset がキャッシュ内にあれば、それを読込対象にする。
4. キャッシュが無ければ `import_overture_buildings.py` を通じて `overturemaps download` を実行し、GeoJSON 系の中間結果を derived tile と `tile_index.json` に変換する。
5. runtime 側は `resolve_urban_outline_layer_for_viewer()` を使って、その derived dataset を追加の高さフィルタなしで読む。
6. 結果の outline 点列は `UrbanOutlineState` と `SkyWindowState.urban_outlines` に反映し、再描画する。
7. 取得中や失敗時はバナー文字列を UI 状態へ反映する。

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

### 8.2 latest-request-wins

視点変更やタイマー更新が短時間に連続した場合、古い非同期結果で UI を巻き戻さないことを優先する。  
そのため、雲更新では最新要求のみが採用される。

### 8.3 CLI と GUI の整合

- CLI で与えられた初期オプションは、起動後の GUI 状態の初期値になる。
- ただし CLI で明示的に無効化した機能の一部は、セッション中の GUI 再有効化を禁止する。
- 代表例が `--terrain-horizon-opacity 0` による地形地平線ロックアウトである。
- `--sky-opacity 0`、`--cloud-opacity 0`、`--terrain-horizon-opacity 0`、`--urban-outline-opacity 0` は、そのセッションで各 GUI トグルをロックアウトする。

### 8.4 都市アウトライン描画の簡略化

- 都市アウトラインは derived tile から得た建物上端輪郭を描く。
- 建物高さのしきい値は derived tile 生成時の前処理パラメータとし、runtime 読込時の既定値では再適用しない。
- ただし見かけの方位幅が `0.5°` 未満の輪郭は、細い polyline ではなく太い水平線に簡略化する。
- `viewport_interaction_mode` 中は都市アウトライン描画を抑止し、方向キー操作の負荷を下げる。

## 9. エラーモデル

### 9.1 起動エラー

- 必須データの欠落や起動時解決不能は起動中断対象とする。
- 起動中断は例外または明示的な abort 制御で扱う。

### 9.2 補助機能エラー

- 雲取得失敗は雲機能内で閉じる。
- 地形地平線取得失敗は terrain 機能内で閉じる。
- 本体表示を継続し、UI には状態表示のみ反映する。

### 9.3 再試行方針

- 雲は更新タイミングごとに再取得の機会がある。
- 地形地平線は同一セッション中の自動再試行を抑制する。

## 10. データとキャッシュ

### 10.1 同梱データ

- 星カタログ
- DSO データ
- フォント
- タワー・展望地点データ
- 画像アセット

### 10.2 外部取得データ

- 衛星クラウドデータ
- DEM タイル
- Overture 建物データ

### 10.3 キャッシュ方針

- キャッシュ対象は再利用価値の高い外部取得データとする。
- DEM データは永続キャッシュする。
- Overture 建物由来の derived tile は地点条件ごとに永続キャッシュする。
- 地形地平線の計算済みポリラインは永続化しない。
- 雲は取得ソースと中間成果物をキャッシュし、視点変更時の再利用を優先する。

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
