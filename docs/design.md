# zstarview 設計書

最終更新: 2026-04-05

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
- ターミナル、console、log、CLI help、例外文言、subprocess の stdout/stderr に出る可能性がある文字列は、Windows 文字コード事故を避けるため ASCII-only を原則とする。
- 非 ASCII 文字を判定ロジックで扱う必要がある場合は、ソースコード中に文字を直書きせず `"\u2019"` のような Unicode escape を優先する。
- UI 専用文字列は、サポート環境で表示確認済みなら非 ASCII を許容する。

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
  - `--theme` は `night`、`day`、`white`、`black` の 4 preset のみを受け付ける
  - 同梱星表の実上限に合わせ、`-V` / `--vmag-limit` は `10.5` を超える指定を parse 時点で `10.5` へ丸める
  - parser 構築は `add_location_arguments()`、`add_dataset_query_arguments()`、`add_time_arguments()`、`add_render_arguments()` の helper に分割し、将来の別 CLI からも再利用できるようにする
  - ガイドライン表示は `--show-guidelines-initial` として扱い、DSO / アステリウムと同じ起動時 boolean 指定に揃える
- `src/zstarview/startup.py`
  - 起動時の地点解決、設定復元、初期値決定
  - Nominatim による online 地点検索と結果正規化を扱う
- `src/zstarview/config.py`
  - ユーザー設定の保存と読込
  - 前回地点を legacy 文字列または構造化地点オブジェクトとして保存する
- `src/zstarview/paths.py`
  - 設定・キャッシュ・データのパス解決
  - テーマ preset ごとの共有表示定義を持つ

### 4.2 ドメイン計算

- `src/zstarview/astro.py`
  - 恒星、太陽、月、惑星、補助線の計算
  - 可視判定と投影前データ生成
- `src/zstarview/catalog.py`
  - 星カタログの読込と描画用配列の前処理
  - 同梱の分割星カタログは `stars_base` (`vmag <= 6`)、`stars_extra7` (`6 < vmag <= 7`)、`stars_extra8` (`7 < vmag <= 8`)、`stars_extra9` (`8 < vmag <= 9`)、`stars_extra10` (`9 < vmag <= 10`)、`stars_extra_faint` (`10 < vmag <= 10.5`) を前提とする
  - loader 上の `extra_faint` バケット名は `vmag > 10` 用のままだが、同梱データでは実質的な上限は `10.5` として扱う
  - GUI / export 向けの前処理済み星カタログは、単一の数値配列セットと、`vmag <= 6.0` 用の index 配列を保持する
  - 恒星名とアステリズム用 `source_id` は、全件 object 配列ではなく、full catalog index に紐づく疎なメタデータとして別保持する
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
  - `--use-building-top` 指定時は、地点近傍の building / building_part を使って建物頂部基準の観測高さへ置き換える
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

### 4.2.2a `--use-building-top` による起動地点高さ補正

`--use-building-top` は、都市名、`--place`、緯度経度、Google Maps URL で解決した地点について、建物頂部を観測基準に使う補助経路である。

- この補助経路は tower / mountain viewpoint には適用しない。
- location resolver は、解決済みの緯度経度を中心に小半径の Overture building / building_part dataset を同期取得してよい。
- 初期実装では、起動前にこの建物取得を同期的に行う。そのため `--use-building-top` 指定時は、建物取得完了までウィンドウ表示を待ってよい。
- 取得した building footprint 群に対して、指定地点を厳密に含む建物だけでなく、指定地点から `5m` 以内の建物も候補とする。
- `building_part` がある場合は、親建物とその part 群を同一グループとして扱い、そのグループの最大 `height_m` を建物頂部高として採用してよい。
- `min_height_m` は浮いた底面の属性として保持するが、この経路で観測基準に使う値は上端であるため、観測基準高の決定には最大 `height_m` を使う。
- 候補建物が見つからない場合は、従来どおり地表基準の `observer_height_m` を使う。
- 利用者向けの地点情報には `Building height` を表示してよい。

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
- この経路は `zstarview.gui.viewer` の GUI `main()` と別の `main` を持つ。
- console script の実体は `zstarview.cli.export_image` とする。
- parser は `zstarview.cli.args` の helper 群を組み合わせて構築し、地点、時刻、視線、描画オプションは通常 CLI と共有する。
- 画像書き出し CLI 固有オプションは少なくとも次を想定する。
  - `--output`
  - `--image-size`
  - `--layer-timeout-seconds`
  - `--allow-partial-data`
  - `--sixel`
- `--window-geometry`、`--window-frame`、`--sky-update-interval` は GUI 専用とし、この CLI では parser に載せない。
- dataset 参照専用オプション群も、この CLI では parser に載せない。
- 地点解決と天体計算は GUI と同じ下位ロジックを共有するが、レイヤー取得順序は GUI と共有しない。
- 出力画像には既定で hover/HUD を含めず、`RenderSceneData` と `RenderStyle` を中心にベース描画だけを行う想定とする。
- `guide` はベース描画に含めてよい。
- export 画像では static overlay info も既定で描かない。
- export 画像では GUI 向けの外周背景グラデーションも描かず、`content_fov_deg` の外側は透明のままにする。
- 外部依存レイヤーの取得は逐次でも並列でもよいが、CLI 側では「いつまで待つか」と「部分データを許容するか」を引数で決められるようにする。
- 既定は安全側として「部分データは保存しない」とし、明示的に `--allow-partial-data` を指定したときだけ部分出力を許可する。
- `--sixel` は `--output` と併用可能とし、実装順序は「まずファイル保存、その後に端末出力」とする。
- `--sixel` 出力が失敗しても、ファイル保存済みなら警告扱いで成功終了としてよい。
- SIXEL 変換は、一時 PNG ファイルを経由せず、`QImage` を `QBuffer` 等で PNG bytes 化して `img2sixel -` の stdin へ流すパイプ方式を前提とする。
- `--sixel` 指定時は、重い初期化やレイヤー取得へ進む前に `shutil.which("img2sixel")` で存在確認を行い、見つからない場合は明示エラーで終了する。
- `--output -` は PNG bytes を stdout へ直接流す用途とし、stdout を SIXEL 出力でも使う `--sixel` とは併用不可とする。
- `opacity == 0` で無効化されたレイヤーは、取得キュー自体に積まず、layer timeout の待機対象からも外す。
- 実装では `SkyWindow` と GUI controller 群には依存せず、sky/cloud/terrain/urban/aircraft を同期的に順番に取得してから、shared pipeline で `QImage` へ 1 回だけ描画して保存する。
- export-image の雲経路は、GUI と同じく `numpy RGBA` と 2D missing-mask alpha を保持し、最終合成段まで `QImage` へ早期変換しない。
- Qt はフォント読込と `QImage` / `QPainter` 利用のためだけに初期化し、CLI 側ではバックグラウンド worker や signal ベースの寿命管理を持たない。
- `gui/sky_worker.py` の celestial / sky-disc 計算は pure helper `compute_sky_snapshot()` として切り出し、GUI worker と export CLI の両方から共有する。

### 4.3 描画

- `src/zstarview/render/stars.py`
  - 恒星描画と hover object 選択
- `src/zstarview/render/sky_disc.py`
  - sky color disc の生成
- `src/zstarview/render/text.py`
  - preset ごとの文字色、アウトライン色、アウトライン幅の解決
- `src/zstarview/render/asterisms.py`
  - アステリズム線の描画
  - 大円弧をサンプルし、アステリズム専用の広い FOV 境界で円形クリップして描画
- `src/zstarview/render/guides.py`
  - 天の赤道、黄道、地平線などの補助線を `(alt, az)` サンプル列から描画時に `render_view_center` 基準で投影する
- `src/zstarview/render/background.py`
  - テーマ定義に基づいてウィンドウ背景 radial gradient とウィンドウ枠を描画する
- `src/zstarview/splash.py`
  - テーマ定義に基づいてスプラッシュ背景、枠線、情報文字色を構成する
- `src/zstarview/gui/composite.py`
  - 星空、雲、欠損ティント、地面色の合成
  - 雲ハッチ、縞密度生成、欠損マスク適用は NumPy ベースで進め、合成結果の出力段で `QImage` に戻す

#### 4.3.0 テーマ表示定義

- `paths.py` には `ThemeStyle` を持つ。
- `ThemeStyle` は少なくとも次を含む。
  - 通常テキスト用 `TextStyle`
  - ステータス行用 `TextStyle`
  - ウィンドウ背景用 `WindowBackgroundStyle`
  - スプラッシュ用 `SplashStyle`
- `render/text.py` は `ThemeStyle` から文字色、アウトライン色、アウトライン幅を解決する。
- `render/background.py` は `ThemeStyle.window_background` を唯一の参照元として、背景 gradient と border 色を決める。
- `splash.py` は `ThemeStyle.splash` を色定義の参照元とし、背景 alpha は `ThemeStyle.window_background` から導出した平均 alpha を使う。
- 明るい preset (`day`, `white`) では、暗い preset (`night`, `black`) より広い文字アウトライン幅を持てる。

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
- `gui/window_render.py` は最終的に、hover 判定、jump highlight、frame cache、interaction mode など GUI 固有の責務に寄せ、描画本体は薄いラッパに縮小するのが望ましい。

#### 4.3.2 現在の描画パイプライン到達点

- 共有描画本体は `src/zstarview/render/pipeline.py` に置く。
- 共有描画入力は次の 3 つに分ける。
- `RenderSceneData`
- `RenderStyle`
- `RenderHudState`
- `render_scene_into_painter()` と下位の `draw_*` 関数群は、`geometry`、`viewport_rect`、`scene`、`style`、`hud` を明示的に受ける。
- `RenderStyle` は `show_guidelines` を持ち、guide レイヤーと viewport interaction 中の reference line 描画を同じ boolean で制御する。
- `RenderPipelineState` のような中間ラッパ型は廃止し、shared pipeline 側では直接引数で依存関係を表す。
- `RenderSceneData` の cloud image / cloud missing mask は `QImage` ではなく NumPy 配列を持ち、cloud path の変換回数を抑える。
- `gui/window_render.py` は、`paintEvent()` 本線、scene/style/hud の組み立て、frame cache、jump highlight、hover 解決など GUI 固有処理に絞る。
- 現在の通常描画順は概ね次のとおり。
  - `background`
  - `sky-cloud`
  - `guide`
  - `terrain`
  - `stars`
  - `planets`
  - `satellites`
  - `aircraft`
  - `labels`
  - `hover-overlay`
  - `overlay`
  - `status`
- ここでいう `overlay` は、地点名、時刻、Alt/Az、地点高さ、観測者高さなどの static observation overlay を指す。`Vmag limit` は GUI メニューの disabled 項目で示す。
- static overlay は hover 系 HUD と同じ更新タイミングに揃えるため、ベース描画ではなく HUD 描画側で重ねる。
- `guide` は方位ラベルと天頂マーカーを含む独立レイヤーであり、空色・雲合成の上、通常の hover/HUD オーバーレイより手前に置く。
- 幾何学的な地平線、天の赤道、黄道も `show_guidelines` に従う guide 系表示として扱う。
- `show_guidelines == False` のときは、guide レイヤー本体だけでなく、viewport interaction 中の sky reference line 描画もまとめて省略してよい。
- `show_overlay_info` は GUI 側の表示トグルとして保持し、既定では `True` でよい。
- CLI の `--show-observation-info-initial` は `show_overlay_info` の初期値上書きとして扱ってよい。

#### 4.3.3 次段のリファクタリング方針

- 次段の主目的は、hover/HUD とベース描画の分離である。
- `guide` レイヤーは HUD ではなくベース描画側に残す。
- ベース描画には、少なくとも `background`、`sky-cloud`、`guide`、`terrain`、`stars`、`planets`、`satellites`、`aircraft`、`labels` を含める想定とする。
- HUD 側には、少なくとも次を寄せる方向で整理する。
  - 恒星ホバー
  - DSO ホバー
  - static observation overlay
  - jump highlight
  - status line
- `paintEvent()` は最終的に「ベースフレームをキャッシュし、その上に hover/HUD を都度重ねる」形を目指す。
- この変更により、ベースフレーム cache key から `mouse_pos` や hover 対象名などの高頻度変化要素を外し、キャッシュ効率を改善する。

#### 4.3.4 hover/HUD 分離の現在位置

- shared pipeline は、`render_base_scene_into_painter()` と `render_hud_overlay_into_painter()` に分かれている。
- ベース描画は、`background`、`sky-cloud`、`guide`、`terrain`、`stars`、`planets`、`satellites`、`aircraft`、`labels` を担当する。
- HUD 描画は、少なくとも次を担当する。
  - アステリズムの hover 強調
  - 月 hover 時の拡大上書き
  - DSO hover
  - 恒星・惑星 hover 情報
  - static observation overlay
  - jump highlight
  - status line
- static observation overlay は `mouse_pos` に応じて左上と左下を切り替えるため、HUD 側で毎フレーム再評価する。
- 位置切替にはヒステリシスを入れ、上 1/3 で左下、下 1/3 で左上へ切り替え、中央 1/3 では直前位置を保持してよい。
- static observation overlay は、GUI の `show_overlay_info` が有効なときだけ HUD 側で描画する。
- overlay 行順は、地点名、時刻、Alt/Az を先頭に固定してよい。
- overlay anchor の保持状態は window state に持たせてよい。
- 月の `5x` 拡大は、角半径の生値ではなく「通常時の見た目半径」を基準に適用する。
- `guide` レイヤーはベース側に残し、マウス位置によるラベル回避には依存しない安定描画として扱う。
- `gui/window_render.py` の frame cache はベース描画だけを保持し、hover/jump/status はキャッシュ後に都度上書きする。
- これにより、frame cache key から `mouse_pos`、hover 対象名、jump highlight 名、status message を外している。

### 4.4 UI

- `src/zstarview/gui/window.py`
  - メインウィンドウ
  - UI 状態と更新制御の集約点
- `src/zstarview/gui/window_state.py`
  - 画面状態の保持
- `src/zstarview/gui/window_inputs.py`
  - ユーザー指定値と実行時オプションの正規化
- `src/zstarview/gui/window_render.py`
  - 再描画とレンダリング関連の UI ロジック
  - 恒星レイヤは描画時に現在のウィンドウサイズから内部レンダリング面サイズを再計算する
  - 天球ディスク幅が `expected-render-width` 以下なら等倍描画し、それを超える場合は `expected-render-width * sqrt(disc_width / expected-render-width)` に従って内部描画面を縮小する
  - 縮小時は低解像度 `QImage` に恒星を描いてからウィンドウ全体へ拡大転写し、大型ウィンドウでの負荷を抑える
  - ベースフレームの `QImage` キャッシュを持ち、geometry、描画入力、interaction mode などが不変なら前回ベースフレームをそのまま再利用する
  - hover 対象、jump highlight、status line、static observation overlay はキャッシュ後に HUD として重ねる
- `src/zstarview/gui/window_updates.py`
  - バックグラウンド更新結果の反映
- `src/zstarview/gui/sky_worker.py`
  - 星空計算のバックグラウンド実行
- `src/zstarview/gui/famous_star_dialog.py`
  - 代表恒星ジャンプ UI
- `src/zstarview/gui/famous_star_search_dialog.py`
  - 星・アステリズム・place 検索 UI
- `src/zstarview/gui/famous_star_shortcuts.py`
  - ジャンプ・検索用データの整形

### 4.4.1 GUI place 検索ターゲット

- GUI の検索機能は、既存の恒星・アステリズムに加えて `place` ターゲットを扱ってよい。
- `place` ターゲットは Nominatim 検索結果を正規化した固定地表点であり、少なくとも `name`、`display_name`、`latitude_deg`、`longitude_deg`、`kind=\"place\"` を持つ。
- `place` ターゲットは天体や人工衛星と異なり、時間変化で位置更新される対象としては扱わない。
- Nominatim への問い合わせ自体は既存の location resolver 系の責務に寄せてよい。
- GUI 検索用には、起動地点解決用の `--place` とは別に、検索候補を `place` ターゲットへ正規化する薄い変換層を持ってよい。
- この変換層は検索候補の `display_name`、緯度経度、および候補種別の整形だけを担当し、GUI の jump 制御や描画状態は持たない。
- GUI では、恒星・アステリズム検索ダイアログと place 検索ダイアログを分けてよい。
- place の Nominatim 問い合わせは、place 専用ダイアログ内の明示的な検索操作で開始してよく、入力変更のたびに自動問い合わせしなくてよい。
- 候補一覧では結果種別を識別できるようにし、`Place` を他の検索結果と区別してよい。
- place 検索で複数候補が返った場合は、利用者が候補を選択できる一覧を表示してよい。
- 候補選択後は、種別ごとの処理へ分岐するが、最終的には既存の検索 jump / highlight 経路へ統合してよい。
- `place` ターゲット選択時は、現在の観測地点と対象緯度経度から見かけの方位角・仰角を計算する。
- この変換は GUI イベントコードに直接書かず、純粋関数として分離してよい。
- 初期仕様では、地形遮蔽、建物遮蔽、地球曲率による可視・不可視の UI 分岐は導入しなくてよい。
- 仰角が `0°` 付近または `0°未満` であっても、その値をそのまま検索ターゲットの jump / highlight に使ってよい。
- place 選択時は、恒星検索と同様に検索ターゲット方向へ視界中心を移す。
- ハイライト表示には恒星検索と同じ通常の白丸マーカーを使ってよく、`place` 用に別の描画記号や別の jump モードを導入する必要はない。
- ラベル文字列は `place` の `name` または `display_name` を使ってよい。
- 検索ジャンプ状態は既存の jump highlight state を再利用してよく、必要なら `kind` を追加してラベル生成だけ分岐してよい。
- Nominatim 検索は GUI スレッドを塞がないよう非同期で行ってよい。
- 検索失敗、候補なし、タイムアウト、レート制限時は、検索 UI またはステータス表示へ簡潔に返してよい。
- 失敗時に描画 state や視界中心を変更してはならない。
- 検索ダイアログを開くメニュー項目は、debug-oriented な補助機能としてキーボードショートカットを持たなくてよい。
- 初期仕様では、place 検索結果の永続保存、検索履歴、CLI からの place ターゲット検索、place ターゲットの継続追跡は扱わない。

### 4.5 雲データ処理

- `src/zstarview/gui/cloud_controller.py`
  - 雲更新の実行制御
  - キューイング
  - latest-request-wins の適用
- `src/zstarview/gui/cloud_state.py`
  - 雲画像、バナー、欠損状態の保持
  - 雲画像は `numpy RGBA`、欠損マスクは 2D `uint8` alpha 配列として保持する
- `src/zstarview/clouddisc/core.py`
  - クラウドディスク生成のオーケストレーション
  - `CloudDiscConfig.alt_min_deg` による可視高度下限の適用
  - 同一 `CloudSourceData` に対する repeated render では brightness-temperature sampler を再利用してよい
- `src/zstarview/clouddisc/providers/*.py`
  - 衛星データ取得
  - GOES は `CMIPF C13` NetCDF を `xarray` で直接読み、`goes_imager_projection` から内部の geostationary `area` 定義を再構築する
  - Himawari は `ISatSS M1C13` タイル群を直接 stitch し、`fixedgrid_projection` から `area` を再構築する
  - Himawari の observer 指定時は、描画用の局所タイルに加え、warm-threshold 推定用の赤道帯タイルを少数だけ追加取得する
- `src/zstarview/clouddisc/projectors/az.py`
  - 空ディスク向け投影
  - 既定 `alt_min_deg = 1.0°` により、地平線より少し上まで雲投影を残しつつ、極端な低仰角ノイズは可視マスクから外す
  - cloud fetch/render では `content_fov_deg + 2.0°` の overscan を許してよい
- `src/zstarview/clouddisc/render/grayscale.py`
  - 雲画像生成
  - ランタイム本線では `numpy RGBA` を返す
- `src/zstarview/clouddisc/cache/*.py`
  - キャッシュ保存と清掃

### 4.6 地形地平線処理

- `src/zstarview/gui/terrain_controller.py`
  - 地形地平線更新の実行制御
- `src/zstarview/gui/terrain_state.py`
  - 地形地平線の状態保持
- `src/zstarview/terrain/dem.py`
  - DEM 取得、読込、サンプリング
- `src/zstarview/terrain/horizon.py`
  - 方位ごとの見かけ地平線計算

### 4.7 航空機オーバーレイ処理

- `src/zstarview/gui/aircraft_controller.py`
  - 航空機更新の実行制御
  - `5分` タイマー、明示更新要求、latest-request-wins の適用
  - OpenSky 取得結果の正規化と UI 反映
- `src/zstarview/aircraft_constants.py`
  - 航空機更新間隔と `bbox` 既定値の共有定数
  - `5分` 取得、`2秒` 予想再投影、bbox / fade / trail span 定数を定義
- `src/zstarview/gui/aircraft_state.py`
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
- `src/zstarview/aircraft/cache.py`
  - `bbox` 単位の航空機キャッシュ key 生成
  - OpenSky 正規化結果の永続保存と読込
  - fresh TTL / stale fallback 判定
  - 古い cache file の簡易清掃
- `src/zstarview/aircraft/project.py`
  - 航空機の `lat/lon/alt` を観測地点基準の `alt/az` へ変換
  - `velocity` / `heading` / `vertical_rate` による短時間前進予測
  - `4秒前 -> 現在 -> 4秒後` を `2秒` 刻みでサンプリングした折れ線点列と age-based alpha の算出
- `src/zstarview/aircraft/types.py`
  - OpenSky state vector の正規化モデル
  - UI が保持する描画用航空機折れ線モデル

### 4.8 人工衛星データ取得とキャッシュ

- `src/zstarview/satellite_constants.py`
  - 人工衛星の軌道要素取得間隔、fresh TTL、表示位置更新間隔、fetch timeout、失敗 backoff の共有定数
- `src/zstarview/satellites/fetch.py`
  - `wheretheiss.at` ISS TLE endpoint の組み立て
  - fallback 用 CelesTrak `stations` endpoint の組み立て
  - ISS TLE payload と fallback OMM/JSON payload の正規化
  - Skyfield へ渡す軌道要素表現への変換
- `src/zstarview/satellites/cache.py`
  - ISS 用人工衛星キャッシュ file path 管理
  - 軌道要素 JSON の永続保存と読込
  - `ISS` 用 fresh TTL 判定
  - 失敗メタデータと `failure_backoff_until_utc` の永続化
  - stale cache の backoff 再利用
- `src/zstarview/satellites/types.py`
  - cache 読込結果をまとめる内部モデル

### 4.9 人工衛星位置計算と描画

- `src/zstarview/satellites/project.py`
  - Skyfield `EarthSatellite` を観測地点基準の `alt/az` へ変換
  - 視野内判定前の人工衛星描画用軽量モデルを生成
  - `ISS` marker の scale とラベル情報を決定
- `src/zstarview/gui/satellite_state.py`
  - `ISS` の最新軌道要素
  - 最新の描画用人工衛星マーカー列
  - 最終成功時刻
  - 読込中 / 取得中バナー
  - 失敗表示状態
- `src/zstarview/gui/satellite_controller.py`
  - 人工衛星軌道要素更新の実行制御
  - 明示更新要求、latest-request-wins の適用
  - `wheretheiss.at` primary / CelesTrak fallback と失敗 backoff を前提にした fetch orchestration
  - 軌道要素取得結果の UI 反映
- `src/zstarview/render/satellites.py`
  - 人工衛星マーカーの描画
  - 航空機と同系統の紫色、小型クロス、`ISS` の marker を担当
- `src/zstarview/render/pipeline.py`
  - `satellites` レイヤーを `planets` の後、`aircraft` の前へ挿入

### 4.10 ユーティリティ

- `src/zstarview/utils/resolve_city.py`
  - 都市解決補助
- `src/zstarview/utils/timezone_parser.py`
  - `--datetime` 文字列のタイムゾーン解釈
- `src/zstarview/utils/image.py`
  - 月相 RGBA 画像などの NumPy ベース画像生成補助
- `src/zstarview/render/qt_image.py`
  - ランタイム本線で使う `NumPy <-> QImage` 変換補助
- `src/zstarview/utils/qt_pil.py`
  - 補助用途の `Pillow <-> Qt` 変換
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
- `BuildingFootprint`
  - `height_m` は Overture 建物属性から得た地表基準の建物高を表す
  - `ground_elevation_m` は DEM から求めた建物 footprint 代表点の地盤標高を表す
  - 都市アウトライン計算では `top_elevation_m = ground_elevation_m + height_m` を用いる

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
- `CloudImageState`
  - `image`: 雲レイヤーの `numpy RGBA`
  - `missing_mask`: 欠損領域を表す 2D `uint8` alpha 配列
  - `cloud_amount_field`: ストライプ線幅生成用の雲量場
  - `meta`、`coverage_ratio`、`source_key`、`render_key`、`request_id` などの更新メタデータ

雲パイプラインでは 2 種類のキーを分離して扱う。

- `SourceKey`
  - 衛星種別や時間スロットなど、取得元を識別するキー
- `RenderKey`
  - `SourceKey` に視点条件を加えた描画条件キー

この分離により、ソース取得をやり直さずに視点変更のみ再描画できる。

### 5.5 人工衛星関連

- `SatelliteOverlayPoint`
  - `group_key`
  - `satellite_name`
  - `alt_deg`
  - `az_deg`
  - `marker_scale`
  - `show_label`
  - 描画直前まで絞り込んだ人工衛星マーカーモデル
  - 初期実装では軌跡線を持たず、現在時刻の位置だけを保持する
- `SatelliteState`
  - group ごとの最新軌道要素 records
  - 表示用マーカー列
  - 最終成功タイムスタンプ
  - 読込中状態
  - エラーバナー
  - 軌道要素は `24時間` 単位で更新し、マーカー列は `5秒` ごとに再計算してよい

### 5.6 地形地平線関連

- `TerrainHorizonState`
  - 地形地平線の点列
  - 読込中状態
  - エラー表示状態
- 地形地平線プロファイル
  - `(altitude_deg, azimuth_deg)` の系列として保持する
  - 地点依存、時刻非依存のデータとして扱う

### 5.7 ウィンドウ状態

- `SkyWindowState`
  - 現在の視点
  - 直近の描画用視点
  - `viewport_interaction_mode` による簡易描画状態
  - `viewport_interaction_stars` による簡易描画用の明るい星テーブル
  - ホバー対象
  - ハイライト対象
  - 各更新パイプラインの UI 反映状態
- `SkyWindow._frame_cache_image`
  - `paintEvent` のベース描画部分を保持する `QImage`
  - 同一ベースフレームの再利用に使う
- `SkyWindow._frame_cache_key`
  - 最終フレームキャッシュの無効化条件をまとめたキー
  - ウィンドウサイズ、`render_view_center`、描画トグル、`CelestialData`、空ディスク画像、雲画像、地形/都市アウトラインなどを含む
  - `mouse_pos`、hover 対象、jump highlight、ステータス行文言、static observation overlay の anchor 状態は含めない
- jump highlight は恒星・ISS に加えて GUI place 検索で選ばれた固定地表点も含めてよい
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
- `SatelliteState`
  - 人工衛星描画マーカー列
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
6. 必要に応じて雲更新、地形地平線更新、都市アウトライン更新、人工衛星更新、航空機更新をバックグラウンドで開始する。

### 6.2 星空更新フロー

1. UI 操作またはタイマーで再計算要求が発生する。
2. `SkyDataWorker` がバックグラウンドで天体計算と sky disc 生成を行う。
3. 計算結果を `CelestialData` と描画補助データとして UI へ返す。
4. `SkyWindow` が内部状態を更新し、再描画する。

視線変更とリサイズの連続入力時は例外的に `viewport_interaction_mode` を使う。

1. `render_view_center` を即時更新する。
2. 明るい星 (`vmag <= 4.0`) のみ同期的に再計算し、簡易描画に使う。
3. 補助線と地形地平線は、保持済みデータを `render_view_center` で再投影して追随させる。
4. 雲はこのモード中に旧視線の bitmap を描かず、一時的に非表示としてよい。
5. 最後の入力から 0.7 秒経過後に通常更新を 1 回だけ開始する。

### 6.2.0 雲レイヤー内部表現

- `clouddisc` のランタイム出力は `numpy RGBA` と 2D missing-mask 配列を基本形とする。
- `CloudController` は `QImage` を先に作らず、そのまま `CloudImageState` へ渡してよい。
- `SkyCompositorCache` は cloud image / missing mask / cloud amount field を NumPy ベースで扱い、ストライプ描画、masking、missing tint 適用をそのまま進めてよい。
- 雲投影モデルは、単一球殻だけに固定せず、同じ source から複数の代表高度球殻へ再投影してから混合してよい。
- 初期候補として、`5km` 単層の代わりに `3km` と `8km` の 2 層球殻を `0.5 / 0.5` で混合してよい。
- この複数高度モデルは物理的な雲頂高度推定ではなく、単一高度投影で生じる視差由来の不自然な穴を緩和するための視覚補正として位置付ける。
- `cloud_amount_field` は、雲 RGBA の alpha から `(u, v)` 正規化座標上へ集約した 2D 雲量場として扱ってよい。
- ストライプ描画では、`width` モードと `alpha` モードの 2 方式を持ってよい。
- `width` モードでは、基準線から片側へ 1px 単位で白線を積み上げ、整数本数に加えて次の 1 本だけ小数部相当の alpha を与えてよい。
- `alpha` モードでは、白線幅は固定とし、雲量に応じて白線 alpha を変えてよい。
- 外周境界の見た目を和らげるため、cloud fetch/render 側に小さな overscan を持たせてよい。
- 雲量場の再正規化では、非ゼロ値の下側 `8 percentile` と上側 `92 percentile` を使ってよい。
- `--cloud-opacity` は追加の内部係数を掛けず、そのまま最終 cloud 合成 opacity として使ってよい。
- `QImage` 化は、合成済み画像または最終描画に必要になった段階でのみ行う。
- これにより、cloud path の `QImage <-> NumPy` 往復と、missing mask の不要な 4ch 展開を避ける。

### 6.2.0a ビューポート幾何

- `get_screen_geometry()` は sky/cloud disc 用に固定 `10px` の内側余白を持たず、ウィンドウ全体を基準に geometry を計算してよい。
- border は別レイヤーとして後段で重ねるため、sky/cloud disc は border の下まで描かれてよい。

### 6.2.1 最終フレームキャッシュ

Qt はメニュー操作やボタン状態変化でも `paintEvent` を再発行するため、空の内容が変わらない repaint では重い描画をやり直さないようにする。

1. `window_render.py` は、geometry、描画入力、hover、jump highlight、interaction mode、ステータス行文言などからフレームキーを作る。
2. キーが前回と同じなら、保持済みの最終フレーム `QImage` を `drawImage()` するだけで終了する。
3. キーが変わったときだけ、背景、空ディスク、雲、地形、都市アウトライン、星、惑星、オーバーレイ、ラベル、ステータス行を一時 `QImage` に順に描いて新しいフレームとして保存する。

この最終フレームキャッシュは `SkyCompositorCache` の上位にある。`SkyCompositorCache` は sky/cloud 合成だけを再利用し、その出力を含むウィンドウ全体の描画結果をさらに `window_render.py` 側でキャッシュする二段構えにしている。

### 6.3 人工衛星更新フロー

1. `SkyWindow` が起動時に `SatelliteController` を生成し、初回更新要求を出す。
2. `SatelliteController` は、まず `satellites/cache.py` に `ISS` の fresh cache があるかを確認してよい。
3. fresh cache がある場合、その records をそのまま採用してよい。
4. fresh cache が無い場合、`satellites/fetch.py` は最初に `wheretheiss.at` の ISS TLE endpoint を取得してよい。
5. `wheretheiss.at` 取得に成功した場合は、その TLE、`element_epoch_utc`、取得試行時刻、成功状態、source=`wheretheiss` を `satellites/cache.py` で永続保存してよい。
6. `wheretheiss.at` 取得に失敗した場合だけ、fallback として CelesTrak `stations` endpoint を試し、その中から `ISS (ZARYA)` だけを抽出してよい。
7. fallback 取得に成功した場合は、その records、`element_epoch_utc`、取得試行時刻、成功状態、source=`celestrak` を `satellites/cache.py` で永続保存してよい。
8. primary / fallback のどちらも失敗した場合は、失敗試行時刻、失敗理由、`failure_backoff_until_utc` を cache file に保存してよい。
9. stale cache が残っていて `failure_backoff_until_utc` が未来の場合は、再起動後であっても network fetch を行わず `cache-backoff` として再利用してよい。
10. `satellites/project.py` は `ISS` の正規化済み軌道要素を Skyfield `EarthSatellite` へ変換し、観測地点と現在時刻から `alt/az` を計算する。
11. `satellites/project.py` は視野内にある人工衛星を `SatelliteOverlayPoint` へ落とし込み、地平線下も保持してよい。
12. `satellites/project.py` は `ISS` に `marker_scale` とラベルを与えてよい。
13. `window.py` は API 取得とは別に人工衛星位置再計算を行い、既定では `ISS` を有効対象として扱ってよい。
14. 人工衛星と航空機の位置再計算は、`2秒` 間隔の共通 overlay projection timer で駆動してよい。
15. 人工衛星レイヤーを GUI で再表示したとき、records が fresh なら外部再取得を行わず marker 再計算だけで即時復帰してよい。
16. records が stale でも失敗 backoff 中なら、cache fallback により marker 再計算だけで復帰してよい。
17. `SkyWindow` は `satellite_ready` または位置再計算完了を受けたら `SkyWindowState` を更新し、再描画する。
18. 通常描画では人工衛星を `planets` の後、`aircraft` の前に描く。
19. 初期実装では人工衛星と月・惑星の接近時に特別な隠蔽処理は行わなくてよい。

### 6.4 時刻モードによる補助レイヤー可否判定

1. 補助レイヤーの可否判定は、対象時刻を `now` と比較して `past`、`present`、`future` に分類してよい。
2. `present` では、雲、航空機、人工衛星をすべて有効候補としてよい。
3. `past` では、雲、航空機、人工衛星をすべて無効化してよい。
4. `future` では、雲、航空機、人工衛星をすべて無効化してよい。
5. 地形地平線と都市アウトラインは地点依存レイヤーであり、この時刻分類では無効化しない。
6. 恒星、太陽、月、惑星、DSO、アステリズム、sky disc、ガイド線は `delta_t` または絶対日時から計算した対象時刻で通常どおり描画してよい。
7. GUI と export-image は同じ可否判定を使うのが望ましい。

### 6.5 雲更新フロー

1. `CloudController` が新しい要求を受け取る。
2. ソース取得が必要か、既存ソースで再描画可能かを判定する。
3. 必要なら衛星データを取得し、CloudDisc 用ソースを生成する。
   - GOES は単一 `CMIPF C13` ファイルから直接 BT field を作る。
   - Himawari は `ISatSS M1C13` タイルを stitch して BT field を作る。
   - Himawari の部分タイル取得では、observer 周辺の描画用タイルに加え、赤道帯の warm-threshold 推定用タイルを和集合で取得する。
4. 視点条件に応じて雲画像を描画する。
5. `request_id` により古い結果を破棄し、最新結果のみ UI に反映する。
6. 欠損領域がある場合は欠損マスクも渡す。

### 6.6 地形地平線更新フロー

1. `TerrainHorizonController` が地点に応じた DEM 更新要求を受ける。
2. 必要な DEM タイルを取得またはキャッシュから読込する。
3. 方位ごとに見かけ地平線を計算する。
4. 結果を地形地平線プロファイルとして UI に返す。
5. 画面再投影時は既存プロファイルを使い回し、再取得はしない。

### 6.7 描画フロー

1. sky disc を生成する。
2. 恒星、惑星、月、補助線を重ねる。
3. 地形地平線があれば地平線関連描画へ反映する。
4. 都市アウトラインがあれば白線オーバーレイとして描画する。
5. 人工衛星があれば紫色の小型クロスマーカーとして描画する。
6. 航空機オーバーレイがあれば紫色の折れ線オーバーレイとして描画する。
7. 雲画像と欠損ティントを合成する。
8. ラベル、オーバーレイ、ステータス行を描画する。

### 6.8 都市アウトライン更新フロー

1. `SkyWindow` が起動時またはトグル再有効化時に `UrbanOutlineController` へ更新要求を出す。
2. `UrbanOutlineController` は `lat/lon + radius + min_height + mode` から実行キーを作る。既定 mode は `both`、既定半径は `2.5km` である。
3. `mode=both` の場合、`UrbanOutlineController` は `building` 用と `building_part` 用の derived dataset をそれぞれ確認し、欠けている方だけを取得する。
4. `mode=building` の場合は `building` 用 derived dataset のみを確認・取得する。
5. 同じ `update()` サイクルの中で、`UrbanOutlineController` は `src/zstarview/data/skyscraper_tiles_z14.json` を読み、視点中心 `60km` 円に交差しつつ内側半径 `radius_km` だけには収まらない seed tile を選ぶ。既定 `radius_km` は `2.5km` である。
6. 選ばれた skyscraper seed tile については、専用 cache root 配下の derived dataset を確認し、未取得 tile だけを `overturemaps download --bbox=... --type building` で取得する。
7. skyscraper tile は import 時に `height_m >= 150` で前処理され、runtime では `resolve_urban_outline_layer_for_viewer(..., radius_km=60.0, min_distance_km=radius_km)` として読む。
8. runtime 側は通常 derived dataset 群を追加の高さフィルタなしで読む。skyscraper derived dataset 群は cache 自体は常に `150m` 下限で共有しつつ、`min_building_height_m > 150` の場合のみ runtime 側で追加高さフィルタをかけてよい。
9. runtime マージ時には、通常レイヤー側で `building_part` が持つ `parent_building_id` を参照し、対応する親 `building` 外形を除外する。
10. runtime は都市アウトライン生成前に、観測地点の DEM 地盤標高 `observer_ground_elevation_m` を解決する。海上など `nodata` は terrain horizon と同様に `0.0m` として扱ってよい。
11. runtime は各建物について DEM から footprint 代表点の `ground_elevation_m` を解決し、derived tile 読込結果へ付与する。代表点は polygon centroid を第一候補とし、centroid が polygon 外に出る場合は bbox center または外周サンプルの中央値へフォールバックしてよい。
12. `compute_urban_outlines()` は `observer_ground_elevation_m + observer_height_m` を観測者標高、`ground_elevation_m + height_m` を建物頂部標高として見かけ仰角を計算する。
13. `building_part` の `min_height` は derived tile では `min_height_m` として保持してよい。これは底面の持ち上がり量であり、頂部標高計算では `height_m` を ground-to-top として優先する。
14. 遠距離スカイスクレーパー補助レイヤーで DEM が必要な場合も、専用 cache ではなく `copernicus-dem` の長寿命 cache を terrain horizon と共用してよい。
15. `compute_urban_outlines()` は建物ごとの `height_m` を保持した輪郭列を返し、`resolve_urban_outline_layer_for_viewer()` はそれを `UrbanOutlinePolyline` の列に変換する。
16. `UrbanOutlineController` は通常レイヤーと skyscraper レイヤーをマージして 1 回の `urban_ready` として反映する。skyscraper 取得が失敗した場合は、通常レイヤーだけで `urban_ready` してよい。
17. `--urban-outline-skyscraper-only` 指定時は、通常近距離 derived dataset の確認・取得・解決をスキップし、skyscraper レイヤーだけを解決する。
18. 描画時は `50m` 以上を CLI 指定 opacity の基準とし、`0m` ではその `25%` になるよう高さ比例で alpha を下げる。
19. 結果の outline 列は `UrbanOutlineState` と `SkyWindowState.urban_outlines` に反映し、再描画する。
20. 取得中や失敗時はバナー文字列を UI 状態へ反映する。

補足:
- 旧 `list[list[(alt, az)]]` 形式の runtime 互換コードは削除し、都市アウトライン描画は `UrbanOutlinePolyline` のみを受け付ける。

### 6.8 航空機オーバーレイ更新フロー

1. `SkyWindow` が起動時に `AircraftController` を生成し、初回更新要求を出す。
2. `SkyWindowUserOptions` は `aircraft_opacity` と `aircraft_gui_allowed` を持ち、`--aircraft-opacity 0.0` のときは起動直後から航空機問い合わせを行わない。
3. `AircraftController` は観測地点の `lat/lon` から OpenSky 用 `bbox` を作る。南北は既定 `±1.0°`、東西は緯度に応じて最低 `90km` を確保しつつ、面積は `25 square degrees` 以下へ抑える。
4. `AircraftController` は明示更新要求または次回 fetch timer に従い、OpenSky 取得を 1 回だけ開始する。
5. `AircraftController` は `bbox` を確定したら、まず `aircraft/cache.py` に同じ `bbox` の fresh cache があるかを確認してよい。
6. fresh cache が存在する場合、`AircraftController` は OpenSky API を呼ばずにその `AircraftSnapshot` 列を使って `aircraft_ready` してよい。
7. fresh cache が無い場合だけ、`aircraft/opensky.py` は `states/all?lamin=...&lamax=...&lomin=...&lomax=...` を呼び、state vector 配列を `AircraftSnapshot` 列へ正規化する。
8. 正規化時には `lat/lon` 欠損、`on_ground=true`、極端に低速な機体を落としてよい。
9. OpenSky 取得が成功したときは、正規化済み `AircraftSnapshot` 列、取得時刻、`bbox`、source 名を `aircraft/cache.py` で永続保存してよい。
10. OpenSky 取得が失敗した場合でも、同じ `bbox` の stale cache が fallback 上限以内なら、それを `aircraft_ready` 相当として返しつつ UI には stale 利用中であることを示してよい。
11. `aircraft/project.py` は各 `AircraftSnapshot` の `velocity`、`heading`、`vertical_rate`、`last_contact` を使って短時間前進予測し、`2秒前 -> 現在 -> 2秒後` の折れ線端点を含む `AircraftOverlayPoint` を作る。
12. `aircraft/project.py` は age に応じた alpha scale も計算し、`90秒` を超えた機体が次回取得まで徐々に薄くなるようにする。
13. `window.py` は API 取得とは別に `AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS` の UI タイマーを持ち、既定では `2秒` ごとに保持済み snapshot から折れ線データだけを再投影してよい。
14. `window.py` は fetch timer を single-shot で扱い、レイヤー再表示時には `AircraftState.last_success_utc` から cache age を計算して次回 API 呼び出し時刻を再調整してよい。
15. 保持済み snapshot が新しければ、GUI で航空機レイヤーを再表示しても API を再問い合わせせず、再投影だけで即時復帰してよい。
16. 保持済み snapshot が更新間隔を超えて古いときだけ、レイヤー再表示時に即時再取得してよい。
17. 描画時は観測地点から `50km` を超える機体を落とし、`10km` 以内かつ `90秒` 以内の機体だけに `callsign` を付けてよい。
18. 新しい取得または fresh/stale cache 復元が成功したときだけ `AircraftState.snapshots` と `AircraftState.overlay_points` をまとめて置き換える。
19. 取得失敗かつ有効な stale cache も無い場合は、直前成功スナップショットを保持したまま、`AircraftState` の失敗表示だけを更新する。
20. 古い非同期結果で UI を巻き戻さないため、`AircraftController` も request id ベースの latest-request-wins を使ってよい。
21. `SkyWindow` は `aircraft_ready` または予想再投影完了を受けたら `SkyWindowState` を更新し、再描画する。
22. `viewport_interaction_mode` 中は航空機オーバーレイ描画を抑止してよい。モード終了後に保持済み `overlay_points` で通常描画へ戻る。

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
  - 人工衛星データ取得と可視マーカー生成
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
- 人工衛星専用状態は `SatelliteState` に分離する。
- 航空機専用状態は `AircraftState` に分離する。

### 8.2 latest-request-wins

視点変更やタイマー更新が短時間に連続した場合、古い非同期結果で UI を巻き戻さないことを優先する。  
そのため、雲更新、人工衛星更新、航空機更新では最新要求のみが採用される。

### 8.3 CLI と GUI の整合

- CLI で与えられた初期オプションは、起動後の GUI 状態の初期値になる。
- ただし CLI で明示的に無効化した機能の一部は、セッション中の GUI 再有効化を禁止する。
- 代表例が `--terrain-horizon-opacity 0` による地形地平線ロックアウトである。
- `--sky-opacity 0`、`--cloud-opacity 0`、`--terrain-horizon-opacity 0`、`--urban-outline-opacity 0` は、そのセッションで各 GUI トグルをロックアウトする。

### 8.4 人工衛星レイヤーの更新粒度

- 人工衛星の current 用軌道要素 cache の fresh 判定は `element_epoch_utc` 基準とし、現在の実装では `ISS=24時間` を用いてよい。表示上の位置再計算は `2秒` 間隔で行ってよい。
- 人工衛星描画は realtime view の現在時刻だけを描き、タイムシフト表示では描かない。
- 初期実装では軌跡線を持たない。
- `satellite_opacity <= 0.0` または `ISS` 表示が無効の間は、人工衛星 fetch timer と位置再計算 timer を止めてよい。
- 描画は視野内に限定し、地平線下も表示してよい。`ISS` は少し大きい marker を使ってよい。
- GUI 既定の有効対象は `ISS` としてよい。
- 航空機と人工衛星の位置再計算は、共通 overlay projection timer で同期させてよい。
- GUI から再表示したときは `last_success_utc` を見て fresh cache を優先し、不要な `wheretheiss.at` / CelesTrak 再取得を避けてよい。
- stale cache は通常は再取得優先でよいが、失敗 backoff 中は表示継続に使ってよい。
- タイムシフト表示では人工衛星レイヤー自体を無効化してよく、current cache の過去探索や archive snapshot は持たなくてよい。

### 8.5 航空機オーバーレイの更新粒度

- 航空機観測値の取得は `5分` 間隔とし、表示上の予想再投影は `2秒` 間隔で行ってよい。
- 航空機オーバーレイは現在時刻近傍だけを対象とし、過去時刻と未来時刻では取得も描画も行わない。
- `aircraft_opacity <= 0.0` の間は、fetch timer と予想再投影 timer を止めてよい。
- 折れ線は `2秒前 -> 現在 -> 2秒後` の 3 点から構成し、機体の短時間進行方向を示す。
- 線幅は観測地点に近い機体ほど太く、遠い機体ほど細くしてよい。
- 描画は観測地点から `50km` 以内に限定し、近距離機体だけにコールサインを出してよい。
- 取得時刻から `90秒` を超えたら alpha を段階的に下げ、次回 API 更新まで視覚的に古さを示してよい。
- GUI から再表示したときは `last_success_utc` を見て fresh cache を優先し、不要な OpenSky 再問い合わせを避けてよい。

### 8.6 都市アウトライン描画の簡略化

- 都市アウトラインは derived tile から得た建物上端輪郭を描く。
- derived tile の各建物レコードは `building_id` を持ち、`building_part` 由来レコードは必要に応じて `parent_building_id` を持つ。
- `building` と `building_part` を併用する場合、`parent_building_id` を持つ part が存在する親 `building` は描画対象から外す。
- 建物高さのしきい値は derived tile 生成時の前処理パラメータとし、runtime 読込時の既定値では再適用しない。
- 遠距離スカイスクレーパー補助レイヤーは `building_part` を使わず `building` のみを扱い、runtime では `min_distance_km` を使って通常レイヤー半径より内側の建物を落とす。
- 建物の高さ属性 `height_m` は Overture 建物属性の地表基準高さとして保持し、見かけ仰角計算では DEM から求めた `ground_elevation_m` を必ず加える。
- `building_part` の `min_height_m` は、底面が地表より上から始まる場合のオフセットとして保持する。現在の上端輪郭描画では頂部計算を変えないが、将来の底面表現や厚み表現に利用してよい。
- 観測者側も `observer_height_m` 単独ではなく、観測地点の DEM 地盤標高を加えた絶対標高で扱う。
- 仰角計算式は `atan2((ground_elevation_m + height_m) - (observer_ground_elevation_m + observer_height_m), distance_m)` を正本とする。
- 現行 derived tile に `ground_elevation_m` を永続化しない場合でも、runtime 解決時に DEM サンプリングして同等の結果を得られることを優先する。
- 遠距離スカイスクレーパー補助レイヤー用の DEM も、観測地点中心の `copernicus-dem` cache に保存して共用する。スカイスクレーパー専用 DEM sidecar や専用 DEM root は現時点では持たない。
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
- 航空機は OpenSky 再取得に失敗しても、同じ `bbox` の stale cache が fallback 上限以内ならそれを使って継続してよい。

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
- `wheretheiss.at` ISS TLE
- CelesTrak `stations` fallback OMM JSON

### 10.3 キャッシュ方針

- キャッシュ対象は再利用価値の高い外部取得データとする。
- DEM データは永続キャッシュするが、保存期限の既定値は `90日` とする。
- Overture 建物由来の derived tile は地点条件ごとに永続キャッシュするが、保存期限の既定値は `30日` とする。
- 遠距離スカイスクレーパー補助レイヤー用 derived tile は `overture_skyscrapers/<tile-cache-key>/bldg` 配下に tile 単位で永続キャッシュし、保存期限の既定値は通常建物キャッシュと同じ `30日` とする。
- 地形地平線の計算済みポリラインは永続化しない。
- 雲は取得ソースと中間成果物をキャッシュし、視点変更時の再利用を優先する。
- 雲の source から作る sampler は source 単位で再利用してよい。Alt/Az 変更で source が同じ場合、sampler の再構築は避ける。
- 航空機スナップショットは `bbox` 単位の少数 JSON file として短寿命永続キャッシュしてよい。
- 航空機 cache file には少なくとも `bbox`、`fetched_at_utc`、`source`、`snapshots` を保持する。
- cache key は観測地点そのものではなく、実問い合わせに使う OpenSky `bbox` から導出する。
- clean up は GOES/Himawari のような時刻ディレクトリ走査ではなく、航空機 cache root 配下の古い file を `fetched_at_utc` 基準で削除する簡易方式でよい。
- 人工衛星の軌道要素 cache は current 1 層だけでよい。
- `current` は `ISS` 用の少数 JSON file としてよく、cache key は `iss` 固定でよい。
- 人工衛星の current cache file には少なくとも `element_epoch_utc`、`fetched_at_utc`、`source`、`records`、`last_fetch_attempt_utc`、`last_fetch_failed`、`last_fetch_error`、`last_fetch_failure_utc`、`failure_backoff_until_utc` を保持する。
- 人工衛星の current cache の fresh 判定は `element_epoch_utc` 基準とし、現在の実装では `ISS=24時間` を使う。
- fetch 失敗後の `2時間` backoff は cache file 側にも保存し、アプリ再起動後も継続してよい。
- DEM / Overture 建物キャッシュは、各取得単位ごとに `fetched_at_utc` をメタデータとして持たせ、利用時に TTL 超過かどうかを判定できるようにする。
- TTL 超過時は「即削除」ではなく「stale として再取得対象」に落とし、再取得成功までは既存キャッシュをフォールバック利用できるようにする。
- 別系統の clean up は任意とし、長期間使われない stale キャッシュだけを後段で物理削除してよい。初期方針としては `TTL x 3` 超過を clean up 候補としてよい。
- `--clear-long-lived-cache` は別系統の明示的削除手段として扱い、TTL 判定とは独立に `copernicus-dem`、`overture_buildings`、`overture_skyscrapers` を削除してよい。
- ただし常用防止のため、cache root 直下に `clear_long_lived_cache_meta.json` を置き、`last_cleared_at_utc` を記録して `3日` のクールダウンを設ける。
- クールダウン中に `--clear-long-lived-cache` が再度指定された場合は、削除を実行せず splash と通常ログの両方に拒否理由を表示して終了してよい。
- CLI からの強制再実行オプションは持たず、どうしても必要な場合は cache directory の手動削除を README で案内する。

### 10.4 DEM / 建物キャッシュ再取得設計

- 目的は、毎回の再ダウンロードを避けつつ、長期間放置された DEM / 建物キャッシュが無期限に固定化されることを防ぐことである。
- DEM は更新頻度が低いため、`fresh=90日`、`stale>90日` とする。
- Overture 建物由来 cache は DEM より更新頻度が高いため、通常 derived dataset / skyscraper tile ともに `fresh=30日`、`stale>30日` とする。
- stale 判定に使う基準時刻は、すべての取得単位で `fetched_at_utc` に統一する。
- TTL 判定はファイルの `mtime` に依存しない。`mtime` は補助的な clean up 用ヒントに留め、fresh/stale の一次判定は明示メタデータで行う。
- stale entry は利用時に非同期で再取得を試みる。UI スレッドは既存 cache を読める限りそれを使い、取得完了後に差し替える。
- 再取得に成功した場合だけ cache payload / tile directory を原子的に入れ替え、成功時刻メタデータを更新する。
- 再取得に失敗した場合、既存 stale cache が読めるなら表示継続を優先し、状態表示だけを warning 系へ更新する。
- stale cache も読めない場合に限り、当該レイヤーを unavailable 扱いにしてよい。
- `fetched_at_utc` を付ける粒度は、TTL 判定と再取得の実行単位に一致させる。
- DEM では tile ごとに `fetched_at_utc` を持たせる。
- 通常の Overture derived dataset では dataset directory ごとに `fetched_at_utc` を持たせる。
- 遠距離スカイスクレーパー補助レイヤーでは tile directory ごとに `fetched_at_utc` を持たせる。
- この方針により、fresh/stale 判定と再取得単位を一致させ、ディレクトリ構造の大規模変更を避ける。
- 後方互換のため、既存 cache に `fetched_at_utc` が無い場合は初回読込時の現在時刻を暫定 `fetched_at_utc` として補完し、その metadata を書き戻してよい。
- この補完時刻は元の実取得時刻ではなく移行時刻として扱う。
- そのため既存 cache は移行後の最初の TTL 周期だけ実際より新しく見えるが、大量再取得を避けるため許容してよい。

#### 10.4.1 DEM フロー

1. `TerrainHorizonController` が必要 tile 一覧を確定する。
2. 各 tile について sidecar metadata から `fetched_at_utc` を読む。無い場合は現在時刻で補完して書き戻す。
3. fresh なら既存 tile をそのまま使う。
4. stale なら既存 tile を入力に含めつつ、バックグラウンドで同じ key を再取得する。
5. 再取得成功時は `*.tmp` + `replace()` で tile 本体と metadata を更新する。
6. 再取得失敗時は stale tile を使い続け、UI には stale 利用中であることを示してよい。

#### 10.4.2 Overture 建物フロー

1. `UrbanOutlineController` は通常 derived dataset と skyscraper tile cache をそれぞれ独立に fresh/stale 判定する。
2. `mode=both` の場合、`building` と `building_part` は別々に TTL 判定し、片方だけ stale ならその片方だけ再取得する。
3. 通常 derived dataset は dataset directory 単位で metadata を持ち、`cache_meta.json` などの sidecar file に `fetched_at_utc`、半径、mode、最小高さなどの cache key 情報を保持する。`fetched_at_utc` が無い旧 cache は現在時刻で補完して書き戻してよい。
4. skyscraper 補助レイヤーは tile directory 単位で metadata を持ち、seed tile ごとに stale 判定する。
5. fresh cache があればそれを即時読込する。
6. stale cache があればそれを即時読込しつつ、欠けている dataset / stale dataset だけをバックグラウンドで `overturemaps download` し直す。
7. 再取得成功時は新 directory を一時パスに組み立て、整合性確認後に `replace()` 相当で切り替える。
8. 再取得失敗時は stale dataset を使い続け、次回要求時に再試行してよい。

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
- 雲データ処理: `xarray`, `boto3`, `botocore`
- 地形データ処理: `rasterio`, `pyproj`
- 画像補助: `Pillow`
  - 現在の主要ランタイム経路 (`zstarview`, `zstarview-export-image`) は、月相と雲を含めて `NumPy -> QImage` を基本とする。
  - `Pillow` は主に補助変換と `src/zstarview/data/` 配下の画像生成スクリプトのために残している。

## 13. 文書間の責務分担

- `docs/specification.md`
  - 利用者から見た機能と動作
- `docs/design.md`
  - 内部構造、責務分割、データ構造、処理フロー
- `docs/implementation-history.md`
  - 実装の時系列メモ、TODO、INPROGRESS、変更の履歴
