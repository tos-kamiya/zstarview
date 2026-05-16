# zstarview 設計書

最終更新: 2026-05-16

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
  - 地平線下の地球裏面ガイド合成
- UI 層
  - ウィンドウ管理
  - 入力イベント
  - バックグラウンド更新の起動と反映
- 外部データ連携層
  - 衛星クラウドデータ取得
  - DEM 取得
  - 水面レイヤー取得
  - Overture 建物データ取得
  - 夜間光 GeoTIFF 取得
  - キャッシュ管理

## 4. モジュール構成

### 4.1 起動・設定

- `src/zstarview/gui/viewer.py`
  - GUI アプリケーションの主エントリポイント
  - 起動シーケンスの組み立て
  - startup overlay は初期 sky データと terrain DEM の両方が揃うまで維持し、ready/fail の両経路で初回表示の切り替えを行う
- `src/zstarview/cli/args.py`
  - CLI オプション定義と値解釈
  - タワー一覧・タワー詳細 JSON 出力の即時終了オプションを扱う
  - `--place`、`--place-countrycode`、`--place-lang` の online 地点検索オプションを扱う
  - `--search`、`--list`、`--search-keep-marker` の検索オプションを扱い、GUI 起動時検索と export-image で同じ引数定義を共有する
- `--theme` は `night`、`day`、`white`、`black`、`transparent` の 5 preset を受け付ける
  - `--edge-fov-deg` は起動時の画面投影スケールを、`--content-fov-deg` は描画対象の保持範囲を制御する
  - `--night-light-opacity` は夜間光オーバーレイの基礎強度を制御し、GUI の Layers メニュー上の Night Lights トグルと同じ意味を持つ
  - `--water-surface-opacity` は水面レイヤーの初期表示強度を制御し、GUI の Layers メニュー上の `Water Surface` トグルと同じ意味を持つ
  - 同梱星表の実上限に合わせ、`-V` / `--vmag-limit` は `10.5` を超える指定を parse 時点で `10.5` へ丸める
  - parser 構築は `add_location_arguments()`、`add_dataset_query_arguments()`、`add_time_arguments()`、`add_render_arguments()` の helper に分割し、将来の別 CLI からも再利用できるようにする
  - ガイドライン表示は `--show-guidelines-initial` として扱い、DSO / アステリウムと同じ起動時 boolean 指定に揃える
- `src/zstarview/startup.py`
  - 起動時の地点解決、設定復元、初期値決定
  - Nominatim による online 地点検索と結果正規化を扱う
  - `auto` 指定時の IP-API による現在地自動取得を扱う
- `src/zstarview/location_resolver/ip_api.py`
  - `ip-api.com` へのリクエストと結果のパースを担当する
- `src/zstarview/config.py`
  - ユーザー設定の保存と読込
  - 前回地点を legacy 文字列または構造化地点オブジェクトとして保存する
- `src/zstarview/paths.py`
  - 設定・キャッシュ・データのパス解決
  - テーマ preset ごとの共有表示定義を持つ
  - 夜間光 cache root の解決を持つ

### 4.2 ドメイン計算

- `src/zstarview/astro.py`
  - 恒星、太陽、月、惑星、補助線の計算
  - 可視判定と投影前データ生成
  - `altaz_to_normalized_xy()` と `is_in_fov()` は画面投影用の edge FOV と保持対象の content FOV を分け、前者で正規化半径を決め、後者で描画可否を決める
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

#### 4.2.1 ビューポイント dataset 参照専用 CLI

同梱 `tower_viewpoints.json` と `mountain_viewpoints.json` を直接参照する軽量 CLI 経路を持つ。

- 対象オプション
  - `--list-viewpoints KIND`
  - `--list-viewpoint-names KIND`
  - `--show-viewpoint-json NAME`
- この経路では `startup.py` の都市解決や GUI 初期化へ進まず、`zstarview` コマンドが `tower_viewpoints.py` / `mountain_viewpoints.py` を直接呼び出して標準出力を書き、即時終了する。
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

#### 4.2.2 `--place` による Nominatim 検索起動

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
- 緯度経度を直接入力した場合の GUI 表示は、小数点以下 5 桁で整形してよい。表示用の丸め規則は `utils/latlon_format.py` にまとめ、表示用定数と cache 再利用用定数を分けて管理してよい。
- タイムゾーンは、採用した座標から最近傍の GeoNames 都市を引いて補完する。補完できない場合は UTC を使う。
- HTTP エラー、レート制限、通信失敗、JSON 解析失敗、0 件結果は `StartupAbortError` 相当で起動中断とし、logger 経由でターミナルとスプラッシュへ表示する。
- Nominatim 利用は起動時の単発検索に限定し、候補列挙だけの反復照会経路は持たない。

#### 4.2.3 `--use-building-top` による起動地点高さ補正

`--use-building-top` は、都市名、`--place`、緯度経度、Google Maps URL で解決した地点について、建物頂部を観測基準に使う補助経路である。

- この補助経路は tower / mountain viewpoint には適用しない。
- location resolver は、解決済みの緯度経度を中心に小半径の Overture building / building_part dataset を同期取得してよい。
- 初期実装では、起動前にこの建物取得を同期的に行う。そのため `--use-building-top` 指定時は、建物取得完了までウィンドウ表示を待ってよい。
- 取得した building footprint 群に対して、指定地点を厳密に含む建物だけでなく、指定地点から `5m` 以内の建物も候補とする。
- `building_part` がある場合は、親建物とその part 群を同一グループとして扱い、そのグループの最大 `height_m` を建物頂部高として採用してよい。
- `min_height_m` は浮いた底面の属性として保持するが、この経路で観測基準に使う値は上端であるため、観測基準高の決定には最大 `height_m` を使う。
- 候補建物が見つからない場合は、従来どおり地表基準の `observer_height_m`を使う。
- 利用者向けの地点情報には、地点名がある場合は 1 行目に名前、2 行目に `Lat: ..., Lon: ... | Ground: ..., Building: ...` の compact summary を表示してよい。

#### 4.2.4 `auto` による現在地自動取得起動

地点引数として `auto`（大文字小文字不問）が指定された場合、外部サービスを利用して現在地を推定する。

- 対象引数: `location` に `auto` を指定
- 解決フロー:
  - `location_resolver/ip_api.py` が `http://ip-api.com/json` を呼び出す。
  - レスポンスから緯度、経度、タイムゾーン、都市名、国名を取得する。
  - 取得した座標を基に `ResolvedLocation` を構築する。
  - タイムゾーンは API から取得した値を優先するが、失敗時は最近傍都市からの補完を試みる。
  - この経路は **ip-api.com** の利用条件に従う。非商用利用の制限と 1 分あたり 45 リクエストの上限がある。
- 失敗時の挙動:
  - ネットワークエラーや API エラー（レート制限等）が発生した場合は `LocationResolveError` として起動を中断する。
  - エラー内容はターミナルとスプラッシュスクリーンに表示する。
- データの保存:
  - 解決された地点は、他の指定方法と同様に次回起動用設定として保存する。
  - 保存形式は緯度経度ベース（`lat;lon`）または構造化オブジェクトであってよい。

#### 4.2.5 前回地点の保存形式

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

#### 4.2.6 単発画像書き出し CLI

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
  - `--include-direction-grid`
  - `--sixel`
- `--window-geometry`、`--window-frame`、`--sky-update-interval` は GUI 専用とし、この CLI では parser に載せない。
- dataset 参照専用オプション群も、この CLI では parser に載せない。
- 地点解決と天体計算は GUI と同じ下位ロジックを共有するが、レイヤー取得順序は GUI と共有しない。
- 出力画像には既定で hover/HUD を含めず、`RenderSceneData` と `RenderStyle` を中心にベース描画だけを行う想定とする。
- `guide` はベース描画に含めてよい。
- `--include-direction-grid` が指定された場合、export-image 専用の静的方向グリッドを guide の上に重ねてよい。
- この方向グリッドは GUI の hover guidance とは別実装で、単発の書き出し画像にだけ適用してよい。
- 30 度刻みの主線は通常の線で描き、10 度刻みの交点は小さな十字で示してよい。
- export 画像では static overlay info を画像に焼き込まず、代わりに stderr へ compact summary を出してよい。
  - summary には、地点名がある場合は名前行と `Lat: ..., Lon: ... | Ground: ..., Building: ...` の 1 行要約、時刻、Alt/Az を含めてよい。
- export 画像では GUI 向けの外周背景グラデーションも描かず、`content_fov_deg` の外側は透明のままにする。
- 外部依存レイヤーの取得は逐次でも並列でもよいが、CLI 側では「いつまで待つか」と「部分データを許容するか」を引数で決められるようにする。
- 既定は安全側として「部分データは保存しない」とし、明示的に `--allow-partial-data` を指定したときだけ部分出力を許可する。
- `--sixel` は `--output` と併用可能とし、実装順序は「まずファイル保存、その後に端末出力」とする。
- `--sixel` 出力が失敗しても、ファイル保存済みなら警告扱いで成功終了としてよい。
- SIXEL 変換は、一時 PNG ファイルを経由せず、`QImage` を `QBuffer` 等で PNG bytes 化して `img2sixel -` の stdin へ流すパイプ方式を前提とする。
- `--sixel` 指定時は、重い初期化やレイヤー取得へ進む前に `shutil.which("img2sixel")` で存在確認を行い、見つからない場合は明示エラーで終了する。
- `--sixel` 指定時は、`lsix` と同様に `/dev/tty` へ `ESC[c` を送り、応答の device attributes に `4` が含まれる場合だけ SIXEL 対応とみなす。応答がない、または `4` を含まない場合は明示エラーで終了する。
- `--output -` は PNG bytes を stdout へ直接流す用途とし、stdout を SIXEL 出力でも使う `--sixel` とは併用不可とする。
- 検索関連オプションは 1 つの `--search QUERY` と `--search-keep-marker` を GUI 起動時検索と `zstarview-export-image` で共有してよく、`--list` は `zstarview-export-image` 専用の探索モードとして扱ってよい。
- 単純文字列の `--search QUERY` は label と id の両方を候補にして解決してよい。
- `--list` は検索結果を 1 行 1 件で列挙して終了する探索モードとし、`zstarview` GUI では受け付けずエラーにしてよい。
- `--search-keep-marker` は、1 件に解決できた場合にその対象を描画へ持ち込むためのフラグとして扱ってよい。
- `-A` と `-Z` は、それぞれ仰角と方位の固定値として扱い、未指定の軸だけ検索結果の `alt/az` で補完してよい。
- 検索が 0 件または複数件のとき、`--list` が無ければ候補一覧を stderr へ出して非 0 終了し、`--list` があれば stdout へ 1 行 1 件で列挙して 0 終了してよい。
- `opacity == 0` で無効化されたレイヤーは、取得キュー自体に積まず、layer timeout の待機対象からも外す。
- 実装では `SkyWindow` と GUI controller 群には依存せず、sky/cloud/terrain/urban/aircraft を同期的に順番に取得してから、shared pipeline で `QImage` へ 1 回だけ描画して保存する。
- export-image の雲経路は、GUI と同じく `numpy RGBA` と 2D missing-mask alpha を保持し、最終合成段まで `QImage` へ早期変換しない。
- `zstarview-export-image` の検索解決は headless ルールを優先し、0 件または複数件では表示だけを行って終了し、1 件に解決できた場合だけ描画へ進む。
- 1 件に解決できた場合は、検索結果の `alt/az` を使って未指定の軸だけを補完し、固定軸はそのまま残す。
- `--search-keep-marker` が有効な場合、export-image はその 1 件の marker/label を描画状態へ持ち込んでよいが、プロセス外へ永続状態を残す必要はない。
- Qt はフォント読込と `QImage` / `QPainter` 利用のためだけに初期化し、CLI 側ではバックグラウンド worker や signal ベースの寿命管理を持たない。
- `gui/sky_worker.py` の celestial / sky-disc 計算は pure helper `compute_sky_snapshot()` として切り出し、GUI worker と export CLI の両方から共有する。

#### 4.2.7 検索リクエストの共有解決

- この節では、`search/models.py` / `search/query.py` / `search/jpl.py` / `search/resolver.py` で共有する検索リクエスト設計を定義する。
- `--search` は、GUI 起動時検索と `zstarview-export-image` の両方で共通の単一検索リクエストとして扱ってよい。
- 検索リクエストは少なくとも raw query、`list_only`、`keep_marker` を持ち、必要なら `label` / `id` の限定解釈を表す selector を持ってよい。
- `QUERY` に `=` が含まれない場合は label と id の両方を検索対象にし、`label=...` と `id=...` は片側だけを検索してよい。
- 共有解決層は候補列挙だけを担当し、GUI で dialog を開くか、export-image で stderr へ流して終了するかは上位オーケストレーションで決めてよい。
- 候補の正規化結果は `label`、`kind`、`id`、`command`、`source`、および必要なら当該時刻の `alt/az` を含む `SearchJumpTarget` に落としてよい。
- 候補が 1 件だけに解決した場合は、その 1 件だけを後段の jump / render に渡してよい。

### 4.3 描画

- `src/zstarview/render/stars.py`
  - 恒星描画と hover object 選択
- `src/zstarview/render/sky_disc.py`
  - sky color disc の生成
  - 連続的な sky disc を生成する
  - 実描画は render surface 寸法、view center、sun alt/az、FOV、opacity をキーにキャッシュしてよい
- `src/zstarview/gui/composite.py`
  - sky-disc 合成の直前または直後に Alt/Az オーバーレイの always モードを適用する
  - always モードの既定は `dimalt`、hover モードの既定は `altaz`
  - `dimalt` は控えめな altitude ring、`altaz` は guide レイヤー相当のフルグリッドとして扱う
  - hover モードは HUD レイヤーで同じモード語彙を使い、mouse hover 時に追加表示してよい
- `src/zstarview/render/text.py`
  - preset ごとの文字色、アウトライン色、アウトライン幅の解決
- `src/zstarview/render/asterisms.py`
  - アステリズム線の描画
  - 大円弧をサンプルし、アステリズム専用の広い FOV 境界で円形クリップして描画
- `src/zstarview/render/guides.py`
  - 天の赤道、黄道、地平線などの補助線を `(alt, az)` サンプル列から描画時に `render_view_center` 基準で投影する
- `src/zstarview/render/background.py`
  - テーマ定義に基づいてウィンドウ背景 radial gradient を描画する
  - 端へ行くほど急に抜けないよう、四隅の透明度は内部より少しだけ高く保つ
  - 独自ウィンドウ枠は frameless host のときだけ描画する
  - 右下の resize handle は `QSizeGrip` ではなく専用の子ウィジェットで実装し、その `paintEvent()` で斜線マーカーを描く
- `src/zstarview/render/night_lights.py`
  - 夜間光 GeoTIFF から方位ごとの glow band を生成する
  - 距離帯ごとの区間積分プロファイルを作り、主稜線と副稜線の各 band に対応づける
  - 地形地平線や副稜線の少し上に重ねるための、固定角度帯の強度マップを作る
  - 描画層は 4 本の band を `NIGHT_LIGHTS_BAND_SPECS = ((alpha_scale, width_scale), ...)` の形でまとめ、外側から順に重ねている
  - band の並びは外側から内側へ向かって `alpha` を強く、`width` を細くしていく
  - 帯の切れ目は観測者の真後ろを seam として扱い、固定の `az=0°` では切らない
  - 地形地平線が無い場合は、水平線を基準にする
- `src/zstarview/splash.py`
  - テーマ定義に基づいてスプラッシュ背景、枠線、情報文字色を構成する
- `src/zstarview/gui/composite.py`
  - 星空、雲、欠損ティント、地面色の合成
  - 雲ハッチ、縞密度生成、欠損マスク適用は NumPy ベースで進め、合成結果の出力段で `QImage` に戻す
- `src/zstarview/night_lights/*.py`
  - NASA Earth at Night/Black Marble 2016 Grayscale 500m tiled GeoTIFF の取得、cache、検証、簡易投影を担当する

#### 4.3.1 テーマ表示定義

- `paths.py` には `ThemeStyle` を持つ。
- `ThemeStyle` は少なくとも次を含む。
  - 通常テキスト用 `TextStyle`
  - ステータス行用 `TextStyle`
  - ウィンドウ背景用 `WindowBackgroundStyle`
  - ウィンドウ chrome 用 `WindowChromeStyle`
  - sky-disc 用 `SkyDiscStyle`
  - スプラッシュ用 `SplashStyle`
- `render/text.py` は `ThemeStyle` から文字色、アウトライン色、アウトライン幅を解決する。
- `render/background.py` は `ThemeStyle.window_background` と `ThemeStyle.window_chrome` を参照元として、背景 gradient と frameless 用 border 色を決める。
- `gui/sky_worker.py` は `ThemeStyle.sky_disc` を参照して sky-disc opacity を決める。
- `gui/sky_worker.py` は sky-disc の render surface 寸法を含めて生成し、描画結果を worker 側の入力で再利用してよい。
- `splash.py` は `ThemeStyle.splash` を色定義の参照元とし、背景 alpha は `ThemeStyle.window_background` から導出した平均 alpha を使う。
- 明るい preset (`day`, `white`) では、暗い preset (`night`, `black`) より広い文字アウトライン幅を持てる。

#### 4.3.2 テーマ文字色メモ

`ThemeStyle.text` は HUD、マウスホバー中の名前ラベル、惑星ラベルなどの通常文字列に使われる。  
`ThemeStyle.status_text` はステータス行などの補助文字列に使われる。  
白系テーマは、読みやすさを少し上げるため、暖色寄りのままやや白に近い値へ寄せてある。

- `night`
  - text: `(180, 180, 180)` / `#B4B4B4`
  - status_text: `(190, 190, 160)` / `#BEBEA0`
- `white`
  - text: `(229, 163, 100)` / `#E5A364`
  - status_text: `(220, 155, 94)` / `#DC9B5E`
- `day`
  - text: `(233, 148, 112)` / `#E99470`
  - status_text: `(224, 142, 106)` / `#E08E6A`
- `black`
  - text: `(246, 249, 255)` / `#F6F9FF`
  - status_text: `(255, 220, 220)` / `#FFDCDC`
- `transparent`
  - text: `(242, 245, 250)` / `#F2F5FA`
  - status_text: `(255, 224, 224)` / `#FFE0E0`
  - sky disc opacity: `0.4`
  - window chrome: menu square uses theme-contrast fill and hamburger lines;
    grip uses a dedicated resize-handle widget that draws a single diagonal stroke itself

#### 4.3.3 カラーパレット参照

`color-palette.png` の RGB 値は、描画色の調整時に参照する固定メモとして以下を残す。

- 1色目: `(240, 173, 122)` / `#F0AD7A`
  - never-rises tint
- 2色目: `(206, 122, 240)` / `#CE7AF0`
  - aircraft overlay
  - satellite marker
- 3色目: `(206, 240, 122)` / `#CEF07A`
  - horizon line
  - direction labels
- 4色目: `(122, 226, 240)` / `#7AE2F0`
  - asterisms
- 5色目: `(216, 206, 192)` / `#D8CEC0`
  - terrain horizon
- 6色目: `(112, 99, 89)` / `#706359`
  - earth guide

#### 4.3.4 描画境界と責務分離

- 描画処理は `SkyWindow` / Mixin 直結ではなく、関数呼び出し中心の再利用可能なパイプラインへ分離している。
- 共通化の中心は「最終的に何をどう描くか」であり、「どの順序でデータを取得するか」ではない。
- GUI 経路と CLI 経路は別々のオーケストレーションを持つ。
  - GUI 経路は、星空の初期表示を優先し、雲・地形地平線・都市アウトライン・航空機を後段で非同期に反映する。
  - GUI の frameless window drag は OS / compositor の移動処理へ委ね、viewport interaction とは分離する。
  - CLI 経路は、必要なレイヤーを逐次取得して 1 回だけ描画する。
- 描画層では少なくとも次の境界を持つ。
  - 描画入力: 観測地点、時刻、視線方向、画面サイズ、表示オプション
  - レイヤー入力: celestial、sky disc、cloud、terrain horizon、urban outline、aircraft overlay などの描画用データ
  - 描画関数群: `QPainter` または `QImage` に対して各レイヤーを順に描く純粋寄りの関数
- `gui/window_render.py` は、hover 判定、jump highlight、frame cache、interaction mode など GUI 固有の責務を担い、描画本体への依存を薄く保つ。
- `terrain` レイヤーは、地形地平線の主線と距離帯ごとの副稜線を担当する。
  - 主線は距離に応じた線幅と alpha を使い、近い側ほど強く、遠い側ほど弱く描いてよい。
  - 距離帯ごとの副稜線はまず従来色の連続ポリラインとして全体を描き、その後に見えている部分だけを別色の細い線で上塗りしてよい。
  - 可視/不可視の判定は、方位を粗い bin にまとめたうえで、近い稜線から見た最大仰角を保持する単純な occlusion heuristic でよい。
  - 上塗り線は、下地線より細く、alpha も弱くして空の光のような縁取りとして扱ってよい。
  - 上塗り線の alpha は距離帯で変えず、全距離帯で一定としてよい。

#### 4.3.5 描画パイプライン構成

- 共有描画本体は `src/zstarview/render/pipeline.py` に置く。
- 共有描画入力は次の 3 つに分ける。
  - `RenderSceneData`
  - `RenderStyle`
  - `RenderHudState`
- `render_scene_into_painter()` と下位の `draw_*` 関数群は、`geometry`、`viewport_rect`、`scene`、`style`、`hud` を明示的に受ける。
- `RenderStyle` は `show_guidelines` を持ち、guide レイヤーと viewport interaction 中の reference line 描画を同じ boolean で制御する。
- never-rises の円弧は guide 系表示の一部として扱い、`show_guidelines` が False のときは表示しない。
- `RenderStyle` は `bright_bodies_mode` を持ち、`vmag < 2.0` の恒星と太陽・月・惑星の outline/fill 表示を共有パイプラインへ伝える。
- 恒星レイヤーは `size_px` に応じて描画を切り替え、`1px` は単一ピクセル+微弱ぼかし、`2px` は 2x2、`3px` から `6px` は塗りつぶし矩形、`7px` 以上は 2px 枠線矩形として扱ってよい。
- `vmag < 2.0` の恒星は、通常モードではダイヤ形の上書き強調を重ねてよく、`bright_bodies_mode == "outline"` ではその輪郭のみを使ってよい。
- `RenderPipelineState` のような中間ラッパ型は使わず、shared pipeline 側では直接引数で依存関係を表す。
- `RenderSceneData` の cloud image / cloud missing mask は `QImage` ではなく NumPy 配列を持ち、cloud path の変換回数を抑える。
- shared pipeline は星レイヤーの縮小レンダリング面サイズを一度計算し、cloud stripe density の参照値としても再利用してよい。
- `gui/window_render.py` は、`paintEvent()` 本線、scene/style/hud の組み立て、frame cache、jump highlight、hover 解決など GUI 固有処理に絞る。
- `gui/window_render.py` の present frame cache key には Alt リングの有無を含めてよく、リング表示の切り替えで present frame を再生成してよい。
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
- ここでいう `overlay` は、地点名、Lat/Lon、`Ground ...; Building ...` の 1 行要約、時刻、Alt/Az などの static observation overlay を指す。`Vmag limit` は GUI メニューの disabled 項目で示す。
- static overlay は hover 系 HUD と同じ更新タイミングに揃えるため、ベース描画ではなく HUD 描画側で重ねる。
- `guide` は方位ラベルと天頂マーカーを含む独立レイヤーであり、空色・雲合成の上、通常の hover/HUD オーバーレイより手前に置く。
- `guide` レイヤーには天頂マーカーの「×」が含まれる。
- 幾何学的な地平線、天の赤道、黄道、never-rises の円弧も `show_guidelines` に従う guide 系表示として扱う。
- 幾何学的な地平線、天の赤道、黄道の線分は、固定刻みの点列をそのまま描くのではなく、画面空間の誤差に応じて再帰分割してから `drawPolyline()` へ渡してよい。
- 再帰分割はスクリーン座標系で誤差を評価し、ウィンドウが大きいほど細かく、縮小時は粗くなるようにしてよい。
- `show_guidelines == False` のときは、guide レイヤー本体だけでなく、viewport interaction 中の sky reference line 描画もまとめて省略してよい。
- `show_overlay_info` は GUI 側の表示トグルとして保持し、既定では `True` でよい。
- 文字ラベルは、星、DSO、アステリズム、衛星、航空機、検索結果などの候補をいったん集約してから、HUD 直前に一括描画してよい。
- その一括描画では、先に配置したラベルの予約矩形を共有し、後から描くラベルほど既存ラベルを避けるようにしてよい。

#### 4.3.6 ベース描画と HUD 分離の現在状態

- shared pipeline は `render_base_scene_into_painter()` と `render_hud_overlay_into_painter()` に分かれている。
- `guide` レイヤーはベース描画側に残す。
- ベース描画は `background`、`sky-cloud`、`guide`、`terrain`、`stars`、`planets`、`satellites`、`aircraft`、`labels` を担当する。
- HUD 側は少なくとも次を担当する。
  - 恒星ホバー
  - DSO ホバー
  - static observation overlay
  - jump highlight
  - status line
- `bright_bodies_mode == "outline"` の場合、`stars` は `vmag < 2.0` の恒星をダイヤの輪郭のみで描き、`planets` は太陽をクロスカーソルのみ、惑星を輪郭のみで描く。
- `bright_bodies_mode == "outline"` の場合、月は通常の占有面積を抑えるため輪郭のみで描くが、`enlarge_moon` や hover による拡大表示では通常の月レンダリングを優先してよい。
- 月の phase 塗りや惑星円盤の内部塗りは `bright_bodies_mode == "outline"` で抑止してよいが、拡大月については塗りつぶしを残してよい。
- `stars` の通常描画では、`7px` 以上の恒星のみ 2px 枠線矩形へ切り替え、それ未満は塗りつぶし矩形を維持してよい。`bright_bodies_mode == "outline"` の明るい星ダイヤはこの切り替えより優先してよい。
- `paintEvent()` はベースフレームをキャッシュし、その上に hover/HUD を都度重ねる構成になっている。
- `paintEvent()` は `viewport_interaction_mode` に応じて fast / normal の描画入口を切り替える。
- fast 側では、星・地平線・Earth guide の簡易入口を使い、夜間光グローや詳細 overlay を normal 側に寄せる。
- ベースフレーム cache key から `mouse_pos`、hover 対象名、jump highlight 名、status message を外し、キャッシュ効率を上げている。
- `labels` レイヤーは、各候補の既定配置矩形を見て局所グループを作り、グループ内では重心に近い候補から先に配置する軽量なクラスタ方式を使ってよい。
- 文字色は、星名や惑星名は既定のテーマ文字色を使い、DSO とアステリズムのホバーラベルはそれぞれ本体色に寄せてよい。
- そのグループ内の要素が 2 件だけの場合は、上側の要素を先に処理する。ラベルは対象の左上に出るため、下側の候補は少し下へ押し出しても元の対象からの距離変化が小さい。
- その際も、完全な全体最適化ではなく、既存の priority 順と予約矩形の仕組みを温存し、密集時だけ配置順を補正する方針としてよい。
- status line は短い記号接頭辞を使い、レイヤーが CLI/GUI の設定で無効化されている場合は `---` を付けて不在であることを明示してよい。
- status line の失敗文言は、完全な例外文字列ではなく、`failed`、`timed out`、`unavailable`、`partial` などの短い要約へ正規化してよい。

#### 4.3.7 hover/HUD 分離の現在位置

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
- overlay 行順は、地点名または `Lat/Lon`、`Ground ...; Building ...` の 1 行要約、時刻、Alt/Az を先頭に固定してよい。
- overlay anchor の保持状態は window state に持たせてよい。
- 月の `5x` 拡大は、角半径の生値ではなく「通常時の見た目半径」を基準に適用する。
- `guide` レイヤーはベース側に残し、マウス位置によるラベル回避には依存しない安定描画として扱う。
- `gui/window_render.py` の frame cache はベース描画だけを保持し、hover/jump/status はキャッシュ後に都度上書きする。
- これにより、frame cache key から `mouse_pos`、hover 対象名、jump highlight 名、status message を外している。

### 4.4 UI

- `src/zstarview/gui/window.py`
  - UI 状態と更新制御の集約点
  - 共通ロジックは `SkyWindowCoreMixin` に置く
  - frameless window のドラッグ移動は `DraggableWindow` に委ね、通常の再描画抑制とは分離する
  - 描画対象のクライアント領域は `SkyWindowClientWidget` として分離する
  - ホストウィンドウは `FramelessSkyWindow` と `StandardSkyWindow` に分ける
  - `FramelessWindowFrame` は frameless 専用の外装であり、独自外枠、ハンバーガーボタン、専用 resize handle を持つ
  - `StandardSkyWindow` は OS 標準の `QMainWindow` とメニューバーを使い、独自外枠は描かない
  - 共通 action は `File`、`Search`、`Layers`、`View Direction` の 4 系統に編成する
  - `Layers` は空側の `Sky Color`、`Clouds`、`Satellites`、`Aircraft` と、地面側の `Night Lights`、`Urban Outline`、`Terrain Horizon`、`Earth Guide` に分ける
  - `Sky Guides` は地平線・赤道・黄道・never-rises 領域・方位ラベル・天頂マーカーをまとめた表示群として扱う
  - frameless ではハンバーガーメニューの 1 階層目に同じ 4 系統の submenu を並べる
- `src/zstarview/gui/window_state.py`
  - 画面状態の保持
- `src/zstarview/gui/window_inputs.py`
  - ユーザー指定値と実行時オプションの正規化
- `src/zstarview/gui/window_render.py`
  - 再描画とレンダリング関連の UI ロジック
  - `viewport_interaction_mode` に応じて fast / normal の描画入口を切り替える
  - fast 側では、星・地平線・Earth guide・夜間光などの専用の簡易入口を使い、詳細 overlay を後回しにする
  - 恒星レイヤは描画時に現在のウィンドウサイズから内部レンダリング面サイズを再計算する
  - 天球ディスク幅が `expected-render-width` 以下なら等倍描画し、それを超える場合は `expected-render-width * sqrt(disc_width / expected-render-width)` に従って内部描画面を縮小する
  - 縮小時は低解像度 `QImage` に恒星を描いてからウィンドウ全体へ拡大転写し、大型ウィンドウでの負荷を抑える
  - 惑星と月のマーカーも同じ内部アップスケール係数を `marker_scale` として受け取り、円盤や月面の見た目半径を恒星レイヤと揃える
- `bright_bodies_mode` は、`vmag < 2.0` の恒星を塗りつぶしなしのダイヤ輪郭へ切り替え、太陽をクロスのみ、惑星を輪郭のみのマーカーへ切り替える表示モードとして扱う
  - `bright_bodies_mode == "outline"` でも、`enlarge_moon` か月 hover での拡大表示は通常の月レンダリングを使ってよい
  - アステリズムと都市アウトラインの線幅はこの係数で膨らませず、内部描画面上では固定幅として扱う
  - ベースフレームの `QImage` キャッシュを持ち、geometry、描画入力、interaction mode などが不変なら前回ベースフレームをそのまま再利用する
  - このフレームキャッシュは描画時の実行時キャッシュであり、永続キャッシュの設計は `10.x` に分離する
  - hover 対象、jump highlight、status line、static observation overlay はキャッシュ後に HUD として重ねる
- `src/zstarview/gui/window_updates.py`
  - バックグラウンド更新結果の反映
- `src/zstarview/gui/sky_worker.py`
  - 星空計算のバックグラウンド実行
- `src/zstarview/gui/famous_star_dialog.py`
  - 代表恒星ジャンプ UI
- `src/zstarview/gui/famous_star_search_dialog.py`
  - 星・アステリズム・place・JPL 小天体検索 UI
- `src/zstarview/gui/famous_star_shortcuts.py`
  - ジャンプ・検索用データの整形

#### 4.4.1 GUI place 検索ターゲット

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

#### 4.4.2 検索ダイアログの持続表示オプション

- `gui/famous_star_search_dialog.py` は、`対象検索...` ダイアログとして、検索語入力、候補一覧、確定/取消に加えて、結果を継続表示するためのチェックボックス領域を持ってよい。
- チェックボックスは `Keep marker` の 1 つだけを持ち、マーカーとラベルの両方を継続表示する意味にしてよい。
- `対象検索...` と place 検索ダイアログは、`Keep CLI-specified Alt/Az` のチェックボックスを別に持ってよい。
- そのチェックボックスは、起動時に `-A` または `-Z` が明示されている場合だけ有効とし、初期状態ではチェック済みにしてよい。
- `-A` / `-Z` のうち指定された軸だけを保持し、未指定の軸は検索対象の `alt/az` で補完してよい。
- `SearchJumpTarget` に `preserve_cli_view_center` を持たせ、ダイアログの選択結果から検索ジャンプへ保持意図を伝えてよい。
- これらの設定は検索 UI のフィルタ条件ではなく、ジャンプ後の表示保持フラグとして扱う。
- 継続表示フラグは `SkyWindowState` の一時 jump highlight とは別の状態として保持してよい。
- 一時 jump highlight は従来どおり 3 秒程度で消えるが、継続表示フラグで選ばれたマーカーとラベルは、利用者が明示的に解除するまで残してよい。
- 継続表示対象は積み増さず、検索のたびに現在の 1 件へ置き換える実装としてよい。
- 継続表示対象は、星、アステリズム、place、衛星、JPL 小天体のいずれでもよいが、描画側ではターゲット種別に応じてラベル文言だけを切り替える。
- `対象検索...` は local first で実行し、恒星、アステリズム、place を先に調べてよい。
- local 候補が見つからない場合だけ、`ISS` 専用の衛星検索経路を現在位置基準で試してよい。
- `ISS` として認識されたのに現在位置ベースの位置解決に失敗した場合は、その検索を失敗として扱ってよく、同じ検索語を JPL へ自動フォールバックしてはならない。
- local / `ISS` の両方で候補が見つからない場合だけ JPL 小天体検索へフォールバックしてよい。`JWST` / `Voyager 1` / `Voyager 2` / `Parker` はこの JPL 経路で扱い、以後の継続表示も同一の追跡状態で維持してよい。
- 検索・継続表示の分類は次の対応として扱ってよい。
  - `local` = 恒星 / アステリズム / place
  - `satellite` = `ISS`
  - `JPL-backed spacecraft` = `JWST` / `Voyager 1` / `Voyager 2` / `Parker`
  - `solar-system bodies` = `Sun` / `Moon` / 惑星
  - `JPL small bodies` = 小惑星 / 彗星など
- JPL フォールバックは major body と small body の両方を検索対象に含めてよい。
- 検索欄の直下には、JPL データベースを明示して探すボタンを置いてよい。
- `Sun` / `Moon` / 惑星 は JPL フォールバック対象から除外し、solar-system 側に委ねてよい。
- JPL 小天体候補は、ネットワーク I/O を含む可能性があるため、検索・候補解決は UI スレッドを塞がない経路に寄せてよい。
- JPL フォールバックの結果は major body と small body を含み、keep フラグが有効なら major body でも持続表示してよい。
- JPL 検索結果の `kind` は、既存の `star` / `asterism` / `place` / `satellite` に加えて、小天体を識別できる値を持たせてよい。
- 結果を継続表示する際のマーカーは、衛星マーカーと同程度の見かけサイズを再利用してよい。
- 継続ラベルは、通常の hover ラベルより高い優先度で扱い、重なり回避で消えすぎないようにしてよい。

#### 4.4.3 検索ジャンプと描画状態

- `SkyWindowCoreMixin._jump_to_search_target()` は、選択した対象の現在時刻における見かけ方向を計算し、`viewer_data.view_center` を更新する役割を持つ。
- `preserve_cli_view_center` が `False` のときだけ、`_search_view_center_alt_specified` / `_search_view_center_az_specified` による固定を外してよい。
- `preserve_cli_view_center` が `True` か未指定のときは、既存の CLI 固定フラグの意味を維持してよい。
- 既存の `jump_highlight_name` / `jump_highlight_altaz` / `jump_highlight_until_ms` は、一時的な視覚フィードバック専用として保持する。
- 持続マーカーと持続ラベルは、上記の一時 jump highlight とは別の状態にし、HUD か補助オーバーレイとして再描画できるようにする。
- 持続マーカーと持続ラベルは単一対象の状態として保持し、複数対象の集合を管理しなくてよい。
- 描画キャッシュのキーには、持続マーカー/ラベルの有無と対象の識別子を含めてよい。
- これにより、結果の永続表示を切り替えたときだけ frame cache を破棄し、通常のマウス hover や短い jump highlight では base frame を再利用できる。

#### 4.4.4 永続検索状態モデル

- 永続表示は、`SkyWindowState` に次の2系統を持たせると扱いやすい。
  - `persistent_search_target`
  - `persistent_search_keep_marker`
- `persistent_search_target` は、最後に選ばれた 1 件の検索対象だけを保持する。
- `persistent_search_target` には `SearchJumpTarget` をそのまま使ってよいが、将来の拡張のために `kind`、`object_key`、`command` を正として扱い、描画用の一時座標や短期追跡状態は別途導出してよい。
- `persistent_search_keep_marker` は、検索ダイアログ下部のチェックボックス状態をそのまま反映する。
- 永続表示は「対象の識別」と「見かけ位置の再計算」を分ける。
  - 識別は `persistent_search_target`
  - 位置は `celestial_data` 更新時に再計算し、Horizons-backed spacecraft は短期追跡状態からローカル再投影する
- 位置再計算は、対象種別ごとに既存ロジックへ寄せてよい。
  - `star` / `asterism` はカタログや既存投影を使う
  - `place` は緯度経度から固定変換する
  - `satellite` は `ISS` 専用の既存衛星投影を使う
  - `JPL-backed spacecraft` は Horizons 由来の短期追跡状態や trajectory cache を使う
  - `JPL 小天体` は SBDB/Horizons 由来の時刻依存位置を使う
- `JPL-backed spacecraft` の trajectory cache は、`persistent_search_target` とは別の描画補助として扱ってよい。
- この trajectory cache は検索状態そのものではなく、`SkyWindowState` の描画補助として扱う。
- 画面に見せるラベルは、`persistent_search_keep_marker` が有効なときに同時に描画する。
- 画面に見せるマーカーは、`persistent_search_keep_marker` が有効なときだけ描画する。
- 検索を再実行した場合は、積み増しではなく `persistent_search_target` を上書きする。
- 明示解除は、`persistent_search_target = None` に加えて、`persistent_search_keep_marker = False` に戻す経路を持ってよい。

#### 4.4.5 JPL 小天体検索サービスとキャッシュ

- JPL 小天体検索は、`famous_star_search_dialog.py` 直下でネットワーク処理を書かず、検索サービス関数か controller 経由で行ってよい。
- サービス層は少なくとも次の2段に分けると扱いやすい。
  - **候補検索**: SBDB 由来の名前・番号・別名の解決
  - **位置取得**: Horizons 由来の時刻依存位置取得、または再投影可能な追跡状態の取得
- 候補検索は、検索ボックスで Enter / Search が押されたときだけ実行してよい。文字入力のたびに再問い合わせしなくてよい。
- 候補検索結果は、`normalized_query` をキーにした短期キャッシュを持ってよい。
  - 同一セッション内の再検索を減らす目的で、`24h` 程度の TTL を使ってよい。
  - 更新追従を少し早くしたいなら、実装上は `6h` 程度へ短縮してもよい。
- 候補一覧には少なくとも次を含めてよい。
  - 表示名
  - 小天体番号または designation
  - 種別
  - 検索元
- 候補から 1 件を選んだ後は、その 1 件だけを `persistent_search_target` に反映する。
- 位置取得は、選択対象・観測地点・高度・時刻をキーにした別キャッシュを持ってよい。
  - 小天体の見かけ位置は時刻で変化するため、候補検索キャッシュと同じ寿命で扱わない。
  - Horizons-backed spacecraft は、位置取得結果を `alt/az` の単発値ではなく、短期 trajectory / state vector として保持してよい。
  - 短期 trajectory は 2 秒描画 tick ではローカル補間に使い、Horizons 再問い合わせは粗い更新周期に限ってよい。
  - 1 時間単位で更新するなら、`[現在時刻の前後 30 分]` のような窓をまとめて取ってよい。
- 位置取得キャッシュは、少なくとも次の情報でキー化してよい。
  - 対象の一意 ID または command
  - 観測者の緯度
  - 観測者の経度
  - 観測者の高度
  - 観測時刻の hour bucket
- 位置取得結果は、描画用の `alt/az` 列または短期 trajectory/state vector として保持してよい。JPL-backed spacecraft は後者を優先してよい。
- 既存の衛星キャッシュと同様に、検索中は最新リクエスト優先の更新規則を使ってよい。
- 失敗時は、候補検索失敗と位置取得失敗を別メッセージに分けてよい。
- 候補検索に成功しても位置取得に失敗した場合は、候補一覧は残し、選択後の jump だけ失敗させてよい。

#### 4.4.6 起動時検索オーケストレーション

- `zstarview` の GUI 起動時検索は、`startup.py` で通常の地点解決と設定復元が終わったあとに、共通の検索リクエスト解決層へ渡す。
- GUI 側は、候補が 1 件に解決した場合だけその対象へ視線を向け、必要なら `persistent_search_target` と `persistent_search_keep_marker` を初回描画前に設定する。
- `--search-keep-marker` が付いている場合、GUI 起動直後の初期フレームから marker と label を表示する。
- GUI の検索ダイアログで CLI 視線保持チェックが無効な場合は、検索結果に `preserve_cli_view_center=False` を持たせて `_jump_to_search_target()` へ渡してよい。
- `JWST` / `Voyager 1` / `Voyager 2` / `Parker` の場合は、初回結果の `alt/az` を固定座標として保存せず、後続の描画 tick で追跡状態から再投影してよい。
- 候補が 0 件または複数件の場合、GUI は終了せず、`対象検索...` ダイアログへ検索語と候補一覧を渡して再選択を促す。
- GUI 側の初期検索は、対話ダイアログからの検索と同じ候補モデルを使い、ポスト処理だけを分ける。
- `zstarview-export-image` 側は同じ共通解決層を使うが、0 件または複数件のときは描画へ進まず、`--list` の有無に応じて列挙終了かエラー終了を選ぶ。`--list` は export-image にだけ存在する。
- 起動時検索と手動検索で同じ `SearchJumpTarget` を使い、GUI は jump highlight と persistent overlay を更新し、export-image は 1 回の描画入力として消費する。

#### 4.4.7 検索共通化

- 検索処理は、GUI ダイアログ、GUI 起動時検索、`zstarview-export-image` で共通利用するため、UI から切り出した resolver/service 群に分離している。
- 検索仕様の詳細は `4.2.7` にまとめ、ここでは配置だけを述べる。
- 現在の構造は次の通りである。
  - `search/models.py` は共有データモデルを持つ。
  - `search/query.py` は検索文字列の正規化を担当する。
  - `search/jpl.py` は JPL 候補生成と小天体位置取得を担当する。
  - `search/resolver.py` は共通解決戦略を担当する。
  - `render/search_overlay.py` は GUI と export-image の両方で使う marker / label 重ね描きを担当する。
  - `window.py` は、共通 resolver から返った 1 件結果に対して、視点移動と永続 overlay の state 更新だけを行う薄い orchestrator とする。

#### 4.4.8 ウィンドウドラッグと更新抑制

- `src/zstarview/gui/draggable_window.py`
  - frameless window のドラッグ開始・移動・終了をまとめて扱う
  - 実際の移動は Qt の system move か、環境に応じた代替移動処理へ委ねてよい
- GUI のウィンドウドラッグは viewport interaction mode へ入らない。
- Wayland では system move 中の再描画が compositor によって抑えられる場合があり、Windows では drag 開始時の fast-mode 再描画が操作感を悪化させる場合があるため、drag と fast-mode は結合しない。
- 方位/仰角変更とリサイズは引き続き viewport interaction mode を使ってよい。

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
  - Himawari の slot 採用は、88 タイルの完全一致を前提にせず、observer 描画と warm-threshold 推定に必要なタイルがそろっているかで判定してよい
  - Himawari の観測者向け描画では、観測地点からおおむね `50 km` より遠い欠損タイルを clear-sky 扱いで補完してよい
  - Himawari の warm-threshold 推定は赤道帯タイルを独立に扱い、欠損時は clear-sky 補完ではなく直近の有効な warm-threshold を再利用してよい
  - `_s3_io.py` は共通 S3 ダウンロード helper として、cache hit と download 完了後の両方を provider 固有 validator で検証してよい
  - 検証に失敗した cached file は破棄し、再取得パスへ戻してよい
  - 検証に失敗した freshly downloaded file は atomic replace せず破棄してよい
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
  - GeoTIFF として読み込めるかを cache reuse 前に確認し、壊れた tile は破棄して再取得に戻す
- `src/zstarview/terrain/horizon.py`
  - 方位ごとの見かけ地平線計算
- `src/zstarview/water_overlay.py`
  - 観測地点近傍の sea-mask tile を点群化するための座標正規化を担う
  - 海は `OSM Water Polygons` 由来のローカル sea-mask tile 群を主入力にする
  - 川・湖・運河・水路などの inland water は Overpass API 経由の OSM フットプリントを別系統で取得する
  - `water_tiles_125m`、`water_tiles_250m`、`water_tiles_500m` を距離帯で切り替え、遠方 `500m` 帯は in-memory の代表点へ縮約する
  - 海面は観測地点中心の地理座標を `target_height_m = 0.0` に投影した点として扱い、輪郭ポリゴンの再構成は行わない
  - 水面ドットは terrain horizon と同じ sky-dome 合成段へ重ねるが、海の内外判定や ring reconstruction は行わない
- `src/zstarview/gui/water_overlay_state.py`
  - 水面レイヤーの取得状態、sea-mask 由来の point set、status line 用の banner を保持する
- `src/zstarview/gui/water_overlay_controller.py`
  - 海マスク更新の実行制御
  - latest-request-wins と TTL 判定を適用し、band cache と active point set の切り替えを管理する
  - terrain horizon が有効なときだけ水面ドットを表示し、地形地平線が非表示のときは水面ドットも抑止する
  - `Terrain Horizon` が OFF の間は `Water Surface` の QAction を無効化して、依存関係を GUI で明示する
  - sea-mask の取得結果と inland water の取得結果は別レイヤーとして保持し、入手経路も分けて扱う
  - 取得完了後の active point set は band stats と共に GUI / export-image へ渡す
- `src/zstarview/render/terrain.py`
  - `WaterOverlayPoint` を小さな青色点として描画する
  - 点の alpha は `water_overlay_opacity` を基準にし、距離減衰を反映して terrain horizon の描画と同じ sky-dome 合成段へ重ねる

### 4.6.1 水面レイヤー処理

- `src/zstarview/water_overlay.py`
  - `water_tiles_125m`、`water_tiles_250m`、`water_tiles_500m` を距離帯ごとに読み分ける
  - `water_tiles_500m` は遠距離になるほど in-memory の代表点化を行い、点数を抑える
  - ローカル sea mask の `1` を水、`0` を地面として扱い、水ピクセル中心を水面ドットへ変換する
  - 点の投影は観測地点からの距離と方位を使って地理座標を復元し、`target_height_m = 0.0` の地表貼り付けとして扱う
  - 海の輪郭生成や point-in-polygon 判定は行わず、海マスクの画素をそのまま点群化する
- `src/zstarview/gui/water_overlay_state.py`
  - sea-mask 由来の dot set と、現在の active dots、表示中の mode/source/banners を保持する
- `src/zstarview/gui/water_overlay_controller.py`
  - 海マスク更新の実行制御
  - latest-request-wins と TTL 判定を適用し、band cache と active point set の切り替えを管理する
  - terrain horizon の ON/OFF や DEM ready/fail に応じて、active dots を再選択してよい
  - 取得済みの band stats は保持してよく、terrain horizon の OFF は表示切り替えだけに使ってよい
- `src/zstarview/render/terrain.py`
  - `WaterOverlayPoint` を小さな青色点として描画する
  - 点の alpha は `water_overlay_opacity` と距離減衰を反映し、terrain horizon の描画と同じ sky-dome 合成段へ重ねる
- `src/zstarview/gui/window.py`
  - `Water Surface` メニュー項目、`W` ショートカット、初期 opacity の入力を GUI 操作へ接続する
- `src/zstarview/gui/window_updates.py`
  - 水面レイヤーの status line と ready/fail banner をまとめる
  - status line は `W ---` / `W <count>` の最小表現にして、mode/source を出さない

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
  - JPL Horizons lookup/observer API を使う overlay 用 spacecraft ephemeris fetch
  - Skyfield へ渡す軌道要素表現への変換
- `src/zstarview/satellites/cache.py`
  - ISS と Horizons spacecraft 用人工衛星キャッシュ file path 管理
  - 軌道要素 JSON の永続保存と読込
  - `ISS` と Horizons spacecraft の fresh TTL 判定
  - 失敗メタデータと `failure_backoff_until_utc` の永続化
  - stale cache の backoff 再利用
- `src/zstarview/satellites/types.py`
  - cache 読込結果をまとめる内部モデル

### 4.9 人工衛星位置計算と描画

- `src/zstarview/satellites/project.py`
  - Skyfield `EarthSatellite` を観測地点基準の `alt/az` へ変換
  - Horizons spacecraft の short-term trajectory / observer projection を生成
  - 視野内判定前の人工衛星描画用軽量モデルを生成
  - `ISS` と Horizons spacecraft の marker scale と hover 表示名を決定
- `src/zstarview/gui/satellite_state.py`
  - `ISS` と Horizons spacecraft の最新軌道要素
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
  - 航空機と同系統の紫色、小型クロス、hover 時の線幅増加を含む人工衛星 marker を担当
- `src/zstarview/render/overlay_info.py`
  - 衛星 hover 名の表示
  - 星 hover と独立した衛星 hover 表示の組み立て
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
- `src/zstarview/utils/location_display.py`
  - スプラッシュ、静的 overlay、export-image stderr 向けの地点要約整形
  - `Lat: ..., Lon: ... | Ground: ..., Building: ...` の表示を共通化し、座標指定時は名前行を省略して重複表示を避ける
- `src/zstarview/data/import_overture_buildings.py`
  - `lat/lon + radius` に対応する bbox を計算する
  - `overturemaps download` を呼び、`building` および必要に応じて `building_part` を取得する
  - ダウンロード結果を既存の derived tile JSON と `tile_index.json` 形式へ変換する
  - 出力先既定値は `CACHE_PATH/overture_buildings`
  - 生成した cache metadata に `fetched_at_utc` と `overture_release` を保持する
  - `overture_release` は `2026-03-18.0` のような release version をそのまま保存する
  - release 照合のための直近確認時刻は cache root 直下の別 metadata に保持する

## 5. 主要データ構造

### 5.1 地点・視点関連

- `ViewerData`
  - 観測地点
  - タイムゾーン
  - 表示中心の方位・高度
  - 画面投影用の edge FOV
  - 描画対象保持用の content FOV
  - 観測者の目線高さ
  - 地盤標高と構造物高の表示用値
  - 観測者高さを UI 表示するかどうかのフラグ
  - 画面描画に必要な視点情報
  - `edge_fov_deg` は起動時固定値として扱い、実行中のリサイズや hover 状態では変えない

地点 dataset が持つ高さ情報と `ViewerData.observer_height_m` は別概念として扱う。

- mountain viewpoint
  - dataset 側の高さは山頂ビューポイントの海抜標高 `Elevation`
- tower viewpoint
  - dataset 側の高さは地表からのタワー高またはビューポイント高 `Tower height`
- `ViewerData.observer_height_m`
  - どの入力種別でも、基準観測点から観測者の目線までの高さ
  - CLI `--observer-height-m` はこの値だけを置き換える
  - 既定値は `1.7m`
- `ViewerData.ground_elevation_m`
  - DEM から求めた観測地点の地盤標高
  - DEM を取れない場合は `0.0m` へ正規化してよい
  - terrain horizon の表示ON/OFFとは独立に保持してよく、起動時の DEM 解決結果を観測者基準の絶対高度計算へ流用してよい
- `ViewerData.location_height_m`
  - 建物頂部、タワー高、または解決済み地点の構造物高を表す
  - 該当しない場合は `0.0m` へ正規化してよい
  - public な地点要約では、名前行を別にした上で `Lat: ..., Lon: ... | Ground: ..., Building: ...` として示してよい
- `BuildingFootprint`
  - `height_m` は Overture 建物属性から得た地表基準の建物高を表す
  - `ground_elevation_m` は DEM から求めた建物 footprint 代表点の地盤標高を表す
  - 都市アウトライン計算では `top_elevation_m = ground_elevation_m + height_m` を用いる
  - raw footprint は `rings_lonlat` のまま保持し、観測者基準へ変換した中間表現は別の vectorized cache として持ってよい
  - 中間表現は観測地点、`observer_height_m`、`ground_elevation_m`、建物高さ条件に依存し、`view_center` は含めなくてよい

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

雲レイヤーの内部表現は、`CloudState` と `CloudImageState` が担う。  
`CloudController` は取得結果を `QImage` に先に変換せず、NumPy ベースの雲バッファとして保持してよい。  
`SkyCompositorCache` は cloud image / missing mask / cloud amount field を NumPy ベースで扱い、ストライプ描画、masking、missing tint 適用をそのまま進めてよい。  
雲投影モデルは、単一球殻だけに固定せず、同じ source から複数の代表高度球殻へ再投影してから混合してよい。

内部表現の詳細は次の通り。

- `clouddisc` のランタイム出力は `numpy RGBA` と 2D missing-mask 配列を基本形とする。
- 3 層混合を使う場合、scene cloud amount に応じて `5km` 単層側から `3km` / `5km` / `7km` の既定配分へ滑らかに遷移してよい。
- 低雲量側と高雲量側の境界は、概ね `0.25` と `0.65` の補間帯として扱ってよい。
- この複数高度モデルは物理的な雲頂高度推定ではなく、単一高度投影で生じる視差由来の不自然な穴を緩和するための視覚補正として位置付ける。
- `cloud_amount_field` は、雲 RGBA の alpha から `(u, v)` 正規化座標上へ集約した 2D 雲量場として扱ってよい。
- ストライプ描画では、`width` モードと `alpha` モードの 2 方式を持ってよい。
- `width` モードでは、基準線を中心に左右対称へ白線を積み上げ、整数本数に加えて次の 1 本だけ小数部相当の alpha を与えてよい。
- `width` モードでは、ストライプ中心間隔を基準密度の 1/2 に詰めてよい。これにより、中心対称のまま修正前と同程度の本数を維持してよい。
- `width` モードの alpha 減衰は線形に限らず、基準線付近を強く保って外側でゆるやかに落ちる ease-out カーブとしてよい。
- `alpha` モードでは、白線幅は固定とし、雲量に応じて白線 alpha を変えてよい。
- 外周境界の見た目を和らげるため、cloud fetch/render 側に小さな overscan を持たせてよい。
- 雲量場の再正規化では、非ゼロ値の下側 `8 percentile` と上側 `92 percentile` を使ってよい。
- `--cloud-opacity` は追加の内部係数を掛けず、そのまま最終 cloud 合成 opacity として使ってよい。ここでの合成には、雲ストライプの加算成分と、その下地として使うグレー化の混合率の両方を含めてよい。
- `QImage` 化は、合成済み画像または最終描画に必要になった段階でのみ行う。
- これにより、cloud path の `QImage <-> NumPy` 往復と、missing mask の不要な 4ch 展開を避ける。

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
  - 描画直前まで絞り込んだ人工衛星マーカーモデル
  - 初期実装では軌跡線を持たず、現在時刻の位置だけを保持する
- `SatelliteState`
  - group ごとの最新軌道要素 records
  - 表示用マーカー列
  - 最終成功タイムスタンプ
  - 読込中状態
  - エラーバナー
  - 軌道要素は `24時間` 単位で更新し、マーカー列は `5秒` ごとに再計算してよい

#### 5.6 地形地平線関連

- `TerrainHorizonState`
  - 地形地平線の点列
  - 読込中状態
  - エラー表示状態
- `TerrainSecondaryRidgeState`
  - 主地平線とは別の補助稜線の点列群
  - 読込中状態は主地平線と共有してよい
  - 表示強度や on/off は独立してよい
- 地形地平線プロファイル
  - `(altitude_deg, azimuth_deg)` の系列として保持する
  - 地点依存、時刻非依存のデータとして扱う
- 距離帯ごとの副稜線プロファイル
  - 観測者からの距離帯ごとに、各方位で最も高い点を抽出した polyline 群として保持してよい
  - 距離帯の境界は `0.5 / 1 / 2 / 4 / 8 / 16 / 32 / 64 / 128 km` を初期値として扱ってよい
  - 主地平線のような 1 点/方位 の包絡線ではなく、近景・中景・遠景の層を分けて見せる補助表現として扱う
  - 距離帯ごとの副稜線は、距離の対数に基づいて単調に減衰させ、遠景ほど少し細くしてよい
  - 遠景のぼけは、距離の対数に基づく薄く太い下地線を先に描いてから本線を重ねることで近似してよい
- 副稜線の可視補助色は、距離帯ごとに alpha を変えず、全距離帯で一定としてよい
- 描画方針
  - fast-mode では主地平線の単一シルエット線を維持し、通常描画では距離帯ごとの副稜線の帯線のみを出してよい
  - 距離帯ごとの副稜線は別レイヤーとして描画し、同色でよい
  - 帯線の色は主地平線と同じ色を使い、遠景ほど本線を細くしてよい
  - 遠い帯ほど下地線の太さを少し強めて、ぼやけた遠景に見せてよい
  - 補助稜線が 0 件でも主地平線の表示は変えない
  - 補助稜線は主線より意味を弱く扱い、地形の「見えやすい複数の山並み」を示すための視覚補助として位置付ける

#### 5.6.1 地平線下地球ガイド関連

- `SubsurfaceEarthGuideState`
  - 粗い大陸ポリゴンまたはそれと等価な低解像度 land-mask の参照
  - 現在の観測地点と視線条件に対して投影済みの screen polyline / polygon 群
  - 最終合成用の bitmap cache または `QPainterPath`
  - 読込中状態やネットワーク取得状態は持たず、同梱 static data のみを前提としてよい
- `EarthGuideState`
  - Earth guide の表示強度
  - GUI トグルの on/off 状態
  - CLI で `--earth-guide-opacity` が `0` のときのロックアウト状態
- 地球裏面ガイドの表示は、精密地図ではなく「利用者向け方位ガイド」であることを前提にしてよい。
- 初期実装データは、大陸単位まで簡略化した閉ポリゴン集合でよいが、現在の同梱データは `build_octilinear_earth_guide_svg.py` で生成した runtime JSON を使う。
- その生成は、GeoJSON の ring をまず Douglas-Peucker で簡略化し、その後 SVG 空間で octilinear 化する 2 段処理として扱ってよい。
- 全体の基準簡略化は `2.0 deg` 近傍でよく、島・小島のような小さい ring は面積バケットごとに少し細かく残してよい。
- 面積バケットは、概ね `2 / 4 / 8 deg^2` の3段とし、より大きな landmass はより粗くしてよい。
- 地図データは地球固定座標系で持ち、少なくとも各頂点の `lat/lon` を保持する。
- 実行時には各頂点を Earth-centered な 3D ベクトルへ変換して扱ってよい。
- 観測点 `O` は、観測地点の `lat/lon` と観測者高さ `h` から球面地球半径 `R + h` 上の 3D ベクトルとして構成してよい。
- 観測点の局所基底は `east / north / up` を採用してよい。
- 各地表頂点 `P` について、観測点からの視線 `v = normalize(P - O)` と観測点天頂 `up = normalize(O)` を計算し、`dot(v, up) < 0` を満たす頂点だけを地平線下候補として扱ってよい。
- 地表法線 `n = normalize(P)` は表示の可否の必須条件ではなく、地平線近傍の減衰や線色調整など style 用補助量として使ってよい。
- 粗ポリゴンの辺が地平線境界をまたぐ場合は、両端の `dot(v, up)` 符号差を用いた線形補間で近似クリップしてよい。
- 地図ガイドは厳密な地表可視判定ではなく、観測地点依存で破綻しない近似投影を優先してよい。
- 地形地平線を使う場合でも、まず球面地球基準で地平線下候補を求め、その後に terrain horizon を最終クリップとして適用する二段構えにしてよい。
- `observer_height_m` の変化は、観測点半径の増減として直接反映してよい。塔や展望台の数百メートル級高さ差も同じ式系で扱える。
- earth-guide の描画は terrain horizon の前景線と同じ RGB / alpha カーブを共有し、前景線のみを使う単一ストロークでよい。
- earth-guide の見かけ形状は、球面上の粗リングを観測者ローカルへ投影し、スクリーン空間で再帰分割して近傍だけ細かく、遠景は粗く描いてよい。
- 再帰で得た短い断片は、連続するものだけ再結合して 1 本の polyline として描画してよい。
- 小さい ring は、大きい ring より少し高い解像度を残してよい。特に台湾や UK などの島嶼が消えないよう、面積に応じて `2 / 4 / 8` のような適応的な簡略化と最小面積判定を使ってよい。
- Earth guide は terrain horizon と独立した表示レイヤーとして扱い、`opacity == 0` のときは取得・計算・描画をスキップしてよい。
- Earth guide を CLI で `0` にした場合は、そのセッションの GUI トグルをロックアウトしてよい。

#### 5.6.2 近傍除外と地平線クリップ

- 自己帰属しやすい足元周辺の land 折れ線は、閉ポリゴン化よりも近傍除外を優先してよい。
- 近傍除外は地表距離ベースで扱い、観測者高度に応じて伸縮させてよい。
- 目安として、地表上の除外半径 `dead_zone_km` は次で与えてよい。
  ```text
  h_km = observer_height_m / 1000
  horizon_km = sqrt(2 * R_km * h_km)
  dead_zone_km = clamp(20, 0.25 * horizon_km, 80)
  ```
- `R_km` は地球半径の km 表記とし、`clamp(min, value, max)` は下限と上限で挟む。
- 近傍除外を角度で扱う場合は、`dead_zone_deg = degrees(dead_zone_km / R_km)` へ変換してよい。
- 追加の altitude クリップは、`horizon_dip_deg = degrees(acos(R / (R + h)))` に小さな余裕角 `margin_deg` を加えた値を使ってよい。
- 高所観測では `dead_zone_km` をやや広げてよいが、地平線下ガイド全体が消えないよう上限を設ける。
- この段では「自分のいる陸地を正しく塗る」ことを狙わず、破綻しない方向案内を優先する。
- 実装上は、`dead_zone_km` で足元周辺を落とした後に、`clip_alt_deg` で地平線近傍をさらに絞る二段構えでよい。

#### 5.6.3 点群による塗り表現

- Earth guide の塗り表現は、厳密な polygon fill ではなく、陸地内部に置いた点群ハッチングとして扱ってよい。
- 塗り用点群と輪郭用 polyline は別経路にし、輪郭の見た目と塗りの安全性を分離してよい。
- 塗り用は、まず地図固定座標系でサンプル点を生成し、陸地ポリゴン内部判定を通した点だけを採用してよい。
- その後に、観測者からの地表距離、地平線下判定、`dead_zone_km` による近傍除外、`clip_alt_deg` による地平線近傍クリップを点単位で適用してよい。
- 点群表示は、最終的に小さな円形または矩形のスタンプとして描画してよい。
- 点の配置は見た目の印象に直結するため、sampler は差し替え可能にしてよい。
- ベースラインの sampler は、おおむね等面積寄りの分布を優先してよい。
- 具体例としては、`Fibonacci sphere` 由来の球面点群を使う方法、緯度帯ごとの等面積近似を使う方法、`lat/lon` グリッドにジッタを加える方法を比較対象としてよい。
- `Fibonacci sphere` は地球全体で密度が均一になりやすく、極域の偏りが少ないため、初期の基準 sampler として採用してよい。
- 緯度帯ベースの sampler は、地理座標での直感的な調整がしやすく、見た目のチューニング候補として保持してよい。
- `lat/lon` グリッドは実装が単純で比較用としては便利だが、高緯度で密度ムラが出やすいため、本命よりは検証用に向く。
- 塗り用点群は厳密な面積再現を目標にせず、陸地の存在感を示すハッチングとして十分な密度だけを持てばよい。
- 塗り用点群が細かすぎて描画負荷が増える場合は、fast mode で sampler を間引いてよい。
- 輪郭用 polyline は従来どおり地平線近傍の再帰分割と近似投影を使ってよいが、塗り用点群の簡略化に引きずられなくてよい。

#### 5.6.4 点群の線分化

- 塗り用点群は、そのままの点描として使う代わりに、同じ緯度帯内の隣接点を短い線分で結んでよい。
- ここでの「同じ緯度」は厳密一致ではなく、数十分の 1 度から数分の 1 度程度の緯度帯にまとめた近傍群として扱ってよい。
- 各緯度帯では経度順に並べ、隣接点同士だけを結ぶ。
- 線分が長くなりすぎる場合は、さらに短い断片へ分割してよい。
- 線分化した塗りは、点描よりも恒星の光点と区別しやすく、地表ハッチングとして自然に見せやすい。
- 線分は輪郭線と紛れないよう、短く、薄く、細く描いてよい。
- 緯度帯ごとに少し位相をずらすか、偶数帯と奇数帯でオフセットを変え、碁盤目感を弱めてよい。
- 既存の点群 sampler を流用し、線分生成はその後段の軽い整形として扱ってよい。
- 線分化後も、観測者からの地表距離、地平線下判定、`dead_zone_km`、`clip_alt_deg` の点単位または線分単位のクリップ規則は維持してよい。
- 南極のように極点をまたぐ ring は、線分化の対象から除外してよい。
- fast mode では、点群の間引きに加えて、線分数もさらに減らしてよい。

### 5.7 ウィンドウ状態

- `SkyWindowState`
  - 現在の視点
  - 直近の描画用視点
  - `viewport_interaction_mode` による簡易描画状態
  - `viewport_interaction_release_pending` による release 後の通常更新待ち状態
  - `viewport_interaction_stars` による簡易描画用の明るい星テーブル
  - `viewport_interaction_mode` は fast レンダリングの入口を選ぶ状態で、入力停止後に通常描画へ戻るまで保持する
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

1. `zstarview` が CLI とログを初期化する。
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

この節は sky disc と恒星・太陽・月・惑星・補助線の共通更新フローを扱う。  
雲、人工衛星、地形地平線、都市アウトライン、航空機は、それぞれ `6.3` から `6.8` の個別フローで扱う。

視線変更とリサイズの連続入力時は例外的に `viewport_interaction_mode` を使う。

1. `render_view_center` を即時更新する。
2. 明るい星 (`vmag <= 4.0`) のみ同期的に再計算し、簡易描画に使う。
3. 補助線と地形地平線は、保持済みデータを `render_view_center` で再投影して追随させる。
4. 雲はこのモード中に旧視線の bitmap を描かず、一時的に非表示としてよい。
5. 矢印キーの `release` で簡易モードを終え、sky 更新を要求する。
6. 新しい sky disc が返ってきた時点で、`view-change-release` 起点の cloud / terrain 再開を行ってよい。
7. そのため、簡易モードの終了は固定時間ではなく、最後の物理キー解放と sky 更新完了の 2 段階で扱ってよい。

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
20. `対象検索...` と CLI 検索の `satellite` 解決は `ISS` だけを対象にしてよく、`JWST` / `Voyager 1` / `Voyager 2` / `Parker` は JPL 検索結果として扱ってよい。継続表示は `command` と追跡状態を保持し、表示時に再投影してよい。
21. `build_search_jump_targets(..., include_satellites=False)` を使う経路では、ローカル検索候補へ衛星ショートカットを混ぜなくてよい。

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
   - Himawari の slot 選定では、全 88 タイルの完全性よりも、描画と warm-threshold 推定に必要なタイルの充足を優先してよい。
   - Himawari の描画用タイルに限り、観測地点から 50 km より遠い欠損は clear-sky 相当として扱ってよい。
   - Himawari の赤道帯タイルは別の warm-threshold 推定入力として扱い、欠損時は前回有効値へのフォールバックを優先してよい。
4. 視点条件に応じて雲画像を描画する。
5. `request_id` により古い結果を破棄し、最新結果のみ UI に反映する。
6. 欠損領域がある場合は欠損マスクも渡す。
7. 雲取得または再描画が進行中の間は、viewport interaction 開始時でも既存の雲バッファを保持してよい。確定結果が届いた時点でのみ、必要に応じて更新・破棄する。

### 6.6 地形地平線更新フロー

1. `TerrainHorizonController` が地点に応じた DEM 更新要求を受ける。
2. `TerrainHorizonController` と export-image の地形処理は、`max_distance_km=128.0` を基準に距離サンプルを作る。
3. DEM 取得側は `max_distance_km + 10.0 km` のマージンを足した bbox を使うため、現行実装では `138 km` 相当までのタイルを対象にする。
4. 必要な DEM タイルを取得またはキャッシュから読込する。
5. 方位ごとに見かけ地平線を計算する。
6. 結果を地形地平線プロファイルとして UI に返す。
7. 画面再投影時は既存プロファイルを使い回し、再取得はしない。

### 6.7 描画フロー

1. sky disc を生成する。
2. 恒星、惑星、月、補助線を重ねる。
3. 地形地平線があれば地平線関連描画へ反映する。
4. 地平線下地球ガイドがあれば、地形地平線とは別の補助レイヤーとして描画する。
5. 都市アウトラインがあれば白線オーバーレイとして描画する。
6. 人工衛星があれば紫色の小型クロスマーカーとして描画する。
7. 航空機オーバーレイがあれば紫色の折れ線オーバーレイとして描画する。
8. 雲画像と欠損ティントを合成する。
9. ラベル、オーバーレイ、ステータス行を描画する。

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
15. `compute_urban_outlines()` は建物ごとの `height_m` を保持した輪郭列を返し、必要なら observer-centric の `numpy` 配列や同等のベクトル化配列としてキャッシュしてよい。`resolve_urban_outline_layer_for_viewer()` はそれを `UrbanOutlinePolyline` の列に変換し、描画時の `view_center` 回転へ渡してよい。
16. `UrbanOutlineController` は通常レイヤーと skyscraper レイヤーをマージして 1 回の `urban_ready` として反映する。skyscraper 取得が失敗した場合は、通常レイヤーだけで `urban_ready` してよい。
17. `--urban-outline-skyscraper-only` 指定時は、通常近距離 derived dataset の確認・取得・解決をスキップし、skyscraper レイヤーだけを解決する。
18. 描画時は `50m` 以上を CLI 指定 opacity の基準とし、`0m` ではその `25%` になるよう高さ比例で alpha を下げる。
19. 高層建物の見やすさを上げるため、`100m` から `600m` の間で下地線の線幅だけを線形に太くしてよい。前景の濃い線は固定幅のまま維持してよい。
20. 結果の outline 列は `UrbanOutlineState` と `SkyWindowState.urban_outlines` に反映し、再描画する。`render_view_center` が変わった場合は、既存の observer-centric 配列を再投影するだけでよい。
21. 取得中や失敗時はバナー文字列を UI 状態へ反映する。

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
- `--sky-opacity 0`、`--cloud-opacity 0`、`--terrain-horizon-opacity 0`、`--earth-guide-opacity 0`、`--urban-outline-opacity 0` は、そのセッションで各 GUI トグルをロックアウトする。
- `--visibility-boost` のような視認性補正は、CLI 層で `SkyWindowUserOptions` の opacity 群へ変換し、下流の描画コードには最終値だけを渡してよい。
- 補助レイヤーと小さい図形レイヤーは同じ倍率で持ち上げ、主役レイヤーは据え置きという profile として扱ってよい。
- 変換の初期案としては、地形地平線、Earth guide、都市アウトライン、cloud missing tint のような補助レイヤーは `visibility_boost` をそのまま適用し、sky disc、cloud disc、航空機、人工衛星、アステリズム、月マーカー、短いラベル、ground tint のような小さい図形や薄い補助表示も同じ倍率で適用してよい。地形地平線、Earth guide、アステリズムのような線主体レイヤーは、必要に応じて線幅だけでなく alpha も同じ倍率で持ち上げてよい。ground tint の既定濃さは `0.04` 程度、never-rises tint は `0.06` 程度でよい。
- 星、sky disc、雲本体、背景グラデーションのような主役レイヤーは、原則として boost しなくてよい。
- ただし sky disc / cloud disc は視認性低下が気になる環境向けに、補正対象として扱ってよい。
- 具体的には、`prepare_window_user_options` で `visibility_boost` を各レイヤーの opacity 値へ同じ倍率で展開し、`prepare_window_runtime_options` 以降は通常の opacity 値として扱ってよい。

### 8.4 人工衛星レイヤーの更新粒度

- 人工衛星の current 用軌道要素 cache の fresh 判定は `element_epoch_utc` 基準とし、現在の実装では `ISS` と Horizons spacecraft の両方に `24時間` を用いてよい。表示上の位置再計算は `2秒` 間隔でローカル再投影してよい。
- 人工衛星描画は realtime view の現在時刻だけを描き、タイムシフト表示では描かない。
- 初期実装では軌跡線を持たない。
- `satellite_opacity <= 0.0` または `ISS` 表示が無効の間は、人工衛星 fetch timer と位置再計算 timer を止めてよい。
- 描画は視野内に限定し、地平線下も表示してよい。
- GUI 既定の有効対象は `ISS` としてよい。
- `対象検索...` の `satellite` 経路も `ISS` 専用としてよく、Horizons spacecraft は検索時には JPL small-body / major-body 経路へ寄せてよい。継続表示時は `command` と追跡状態を保持し、2 秒 tick ではローカル再投影してよい。
- 航空機と人工衛星の位置再計算は、共通 overlay projection timer で同期させてよい。
- GUI から再表示したときは `last_success_utc` を見て fresh cache を優先し、不要な `wheretheiss.at` / CelesTrak 再取得を避けてよい。
- stale cache は通常は再取得優先でよいが、失敗 backoff 中は表示継続に使ってよい。
- 人工衛星マーカーは常時ラベルを描かず、hover 名表示を前提としてよい。
- hover の表示名は `SatelliteOverlayPoint.satellite_name` から生成してよい。
- 人工衛星 hover は恒星 hover と独立して扱い、同時に成立してよい。
- 人工衛星 hover の強調は丸囲みではなく、クロス記号の線幅を 2 倍にして示してよい。
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
- 描画の見やすさのため、高層建物の下地線だけは `height_m` に応じて太くしてよい。`100m` までは既定幅、`600m` 以上で 2 倍相当まで線形に増やしてよいが、最後の前景線は固定幅にしてよい。
- 観測者側も `observer_height_m` 単独ではなく、観測地点の DEM 地盤標高を加えた絶対標高で扱う。
- 仰角計算式は `atan2((ground_elevation_m + height_m) - (observer_ground_elevation_m + observer_height_m), distance_m)` を正本とする。
- 現行 derived tile に `ground_elevation_m` を永続化しない場合でも、runtime 解決時に DEM サンプリングして同等の結果を得られることを優先する。
- 遠距離スカイスクレーパー補助レイヤー用の DEM も、観測地点中心の `copernicus-dem` cache に保存して共用する。スカイスクレーパー専用 DEM sidecar や専用 DEM root は現時点では持たない。
- `UrbanOutlinePolyline` は `source` を持ち、通常レイヤーと skyscraper レイヤーの由来を区別できる。ただし現行描画色は共通である。
- 都市アウトラインの内部表現は三層に分けてよい。
  - raw footprint: `lon/lat` の ring 群
  - observer-centric cache: 観測者基準の `numpy` 配列または同等のベクトル化配列
  - render polyline: `UrbanOutlinePolyline.points` に入る描画直前の点列
- `view_center` は render polyline を画面へ回転させるための条件であり、raw footprint や observer-centric cache の生成キーには必須ではない。
- 時刻は建物アウトラインの生成キーに含めなくてよい。建物は地表固定物として扱い、時刻で変わるのは天体や移動体だけでよい。
- `observer_height_m`、`ground_elevation_m`、建物高さ条件が不変なら、観測者基準の配列は再利用してよい。`view_center` が変わった場合は再投影だけを行えばよい。
- 観測者基準の配列は、`numpy` の角度配列や 2D 平面座標配列、あるいは 3D ローカルベクトル配列として持ってよい。描画時には `view_center` に応じた回転行列を一括適用してから最終投影してよい。

#### 8.6.1 都市アウトラインの中間表現と責務分離

都市アウトラインは、次の 3 層に分けて扱ってよい。

- raw footprint layer
  - `BuildingFootprint.rings_lonlat` をそのまま保持する
  - derived tile の正本であり、地理データの更新や再取得の起点になる
  - cache key は観測地点 `lat/lon`、建物条件、`radius_km`、`feature_type`、`min_height_m`、skyscraper 系の選択条件を含めてよい
- observer-centric layer
  - raw footprint を観測者原点のローカル座標へ変換した中間表現
  - 実装は `numpy` ベースの配列でよく、`x/y` 平面配列、`alt/az` 配列、または 3D ローカルベクトル配列のいずれでもよい
  - 生成条件は `lat/lon`、`observer_height_m`、`ground_elevation_m`、建物高さ条件、DEM 解決結果であり、`view_center` は含めなくてよい
  - この層は「観測者が原点」という意味で topocentric な正規形として扱ってよい
- render polyline layer
  - observer-centric layer を `view_center` に応じて再投影した描画直前の点列
  - 実装上は `UrbanOutlinePolyline.points` に `alt/az` サンプルを保持し、`render/terrain.py` 側で最終スクリーン投影してよい
  - `view_center` が変わったときは raw footprint から作り直すのではなく、この層だけを再投影してよい

実装時の責務は次のように分けてよい。

- `compute_urban_outlines()`
  - 建物の選別、距離判定、穴リングの省略、線分化判定を行う
  - observer-centric layer を作る責務を持ってよい
  - `view_center` を受け取る場合は、必要に応じて線分化判定にのみ使ってよい
- `resolve_urban_outline_layer_for_viewer()`
  - derived tile 読込結果と viewer state を受け取り、UI 向けの outline 列へ正規化する
  - observer-centric cache がある場合はそれを使い、なければ `compute_urban_outlines()` へフォールバックしてよい
- `UrbanOutlineController`
  - derived dataset の取得、cache 判定、dataset merge を担う
  - observer-centric cache の保存先または再利用単位を決めてよい
- `SkyWindowRenderMixin` / `render/terrain.py`
  - render polyline layer を受け取り、`view_center` に応じてスクリーンへ投影する
  - 画面回転以外の geospatial 再計算を持たなくてよい

中間表現の更新条件は次のとおりとする。

- raw footprint layer の更新条件
  - derived tile の更新、再ダウンロード、feature type の変更
  - 建物高さの前処理条件の変更
- observer-centric layer の更新条件
  - 観測地点 `lat/lon` の変更
  - `observer_height_m` の変更
  - `ground_elevation_m` の変更
  - 建物高さ条件や DEM 由来の補助値の変更
- render polyline layer の更新条件
  - `view_center` の変更
  - 表示 FOV の変更
  - 線分化や thin-run 判定の閾値変更

この分離により、`view_center` の頻繁な変更は最終投影にだけ作用し、建物データの地理変換や距離判定を繰り返さずに済む。

- ただし、投影後の screen-space 縦幅が十分小さい輪郭は、細い polyline ではなく 2 点線分に簡略化する。
- 遠距離で外周リングの投影幅が小さい場合は、穴リングを描画対象から外してよい。外周リングは残し、`rings_lonlat[0]` を主要輪郭として使い、`rings_lonlat[1:]` は距離と見かけサイズに応じて省略してよい。
- 穴リングの省略判定は、描画段階ではなく `compute_urban_outlines()` 側で行ってよい。これにより、投影・サンプル・フラグメント分割の前に不要な ring を落とせる。
- 判定条件は距離だけに固定せず、外周リングの見かけ方位幅や投影 bbox の大きさを併用してよい。遠距離・小面積ほど穴リングを省略し、近距離では従来どおり保持する。
- さらに、輪郭の screen-space 縦幅が十分小さい場合は、`compute_urban_outlines()` 側で polyline を 2 点線分へ落としてよい。これは `draw_urban_outlines()` ではなく、都市アウトライン生成結果の表現選択として扱う。
- 縦幅のしきい値は、画面上で見分けがつかない程度の薄さを基準にし、距離や FOV に応じて固定ピクセル閾値ではなく投影後の高さで判定してよい。初期実装では、正規化スクリーン座標系での縦幅がごく小さい run を対象にしてよい。
- `viewport_interaction_mode` 中は都市アウトライン描画を抑止し、方向キー操作の負荷を下げる。
- ハンバーガーメニュー表示そのものでは `viewport_interaction_mode` へ入らない。メニュー起動で `paintEvent` が走っても、最終フレームキャッシュの再利用で余分な重い描画を避ける。

## 9. エラーモデル

### 9.1 起動エラー

- 必須データの欠落や起動時解決不能は起動中断対象とする。
- 起動中断は例外または明示的な abort 制御で扱う。

### 9.2 補助機能エラー

- 雲取得失敗は雲機能内で閉じる。
- 雲取得または再描画が in-flight の間は、既存の雲バッファを維持してよい。
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
- NASA Earth at Night / Black Marble 2016 Grayscale 500m tiled GeoTIFF
- OpenSky 航空機 state vector
- `wheretheiss.at` ISS TLE
- CelesTrak `stations` fallback OMM JSON
- `ip-api.com` 現在地データ ( `auto` 指定時のみ。非商用利用制限と 1 分あたり 45 リクエストの上限あり)

### 10.3 キャッシュ方針

- キャッシュ対象は再利用価値の高い外部取得データとする。
- DEM データは永続キャッシュするが、保存期限の既定値は `90日` とする。
- 水面レイヤーの OSM 由来 cache も long-lived cache として扱い、保存期限の既定値は `90日` とする。
- Overture 建物由来の derived tile は地点条件ごとに永続キャッシュするが、保存期限の既定値は `30日` とする。
- 遠距離スカイスクレーパー補助レイヤー用 derived tile は `overture_skyscrapers/<tile-cache-key>/bldg` 配下に tile 単位で永続キャッシュし、保存期限の既定値は通常建物キャッシュと同じ `30日` とする。
- 地形地平線の計算済みポリラインは永続化しない。
- 雲は取得ソースと中間成果物をキャッシュし、視点変更時の再利用を優先する。
- 雲の source から作る sampler は source 単位で再利用してよい。Alt/Az 変更で source が同じ場合、sampler の再構築は避ける。
- 夜間光 GeoTIFF は versioned static cache として保持し、2016 Grayscale 500m tiled 版を初回利用時にのみ取得してよい。
- 夜間光 GeoTIFF は TTL による定期再取得を前提にせず、欠損、破損、明示的な cache 削除、または別版への切り替え時だけ再取得してよい。
- 夜間光 cache は `science.nasa.gov` の Earth at Night flat maps から direct GeoTIFF を落とすだけでよく、Earthdata login token を要求しなくてよい。
- 夜間光は観測地点ごとに重い基礎プロファイルをキャッシュし、sun_alt は連続係数として描画時に掛けてよい。
- 距離帯は `0.5 / 1 / 2 / 4 / 8 / 16 / 32 / 64 / 128 km` を初期値として扱ってよい。
- 距離帯ごとの band profile は、その band 区間だけの積分値を使ってよく、遠い band は近距離分を引き継がない。
- 地形地平線を描かない場合は、水平線を fallback にしてよい。
- 夜間光の描画は、観測者の真後ろ方向を seam として持ち、`az=0°` に固定して切らない。
- 描画レイヤーは `CORE / MID_NEAR / MID_FAR / OUTER` の 4 段として扱い、`MID` と `OUTER` を細かく分けて自然なにじみを作ってよい。
- 航空機スナップショットは `bbox` 単位の少数 JSON file として短寿命永続キャッシュしてよい。
- 航空機 cache file には少なくとも `bbox`、`fetched_at_utc`、`source`、`snapshots` を保持する。
- cache key は観測地点そのものではなく、実問い合わせに使う OpenSky `bbox` から導出する。
- clean up は GOES/Himawari のような時刻ディレクトリ走査ではなく、航空機 cache root 配下の古い file を `fetched_at_utc` 基準で削除する簡易方式でよい。
- 人工衛星の軌道要素 cache は current 1 層だけでよい。
- `current` は `ISS` 用と Horizons spacecraft 用の 2 系列 JSON file としてよく、cache key は `iss` と `horizons` を使ってよい。
- 人工衛星の current cache file には少なくとも `element_epoch_utc`、`fetched_at_utc`、`source`、`records`、`last_fetch_attempt_utc`、`last_fetch_failed`、`last_fetch_error`、`last_fetch_failure_utc`、`failure_backoff_until_utc` を保持する。
- 人工衛星の current cache の fresh 判定は `element_epoch_utc` 基準とし、現在の実装では `ISS` と Horizons spacecraft の両方に `24時間` を使う。
- `iss` と `horizons` はキー名だけでなく payload 形状も分けてよく、`iss` は TLE/OMM 系 records、`horizons` は target command / epoch / short-term trajectory を保持してよい。observer-specific の `alt/az` は描画時にローカル導出してよい。
- fetch 失敗後の `2時間` backoff は cache file 側にも保存し、アプリ再起動後も継続してよい。
- DEM / Overture 建物キャッシュは、各取得単位ごとに `fetched_at_utc` をメタデータとして持たせ、利用時に TTL 超過かどうかを判定できるようにする。
- TTL 超過時は「即削除」ではなく「stale として再取得対象」に落とし、再取得成功までは既存キャッシュをフォールバック利用できるようにする。
- 別系統の clean up は任意とし、長期間使われない stale キャッシュだけを後段で物理削除してよい。初期方針としては `TTL x 3` 超過を clean up 候補としてよい。
- cache hit のたびに利用可能性を検証し、壊れた DEM / 建物 / 雲キャッシュは stale 扱いではなく invalid として破棄してよい。
- cache hit のたびに利用可能性を検証し、壊れた夜間光 GeoTIFF も invalid として破棄してよい。
- freshly downloaded payload も、ファイル名が期待どおりでも実体が壊れていればキャッシュへ昇格させず再取得失敗として扱ってよい。
- `--clear-long-lived-cache` は別系統の明示的削除手段として扱い、TTL 判定とは独立に `copernicus-dem`、`overture_buildings`、`overture_skyscrapers` を削除してよい。
- ただし常用防止のため、cache root 直下に `clear_long_lived_cache_meta.json` を置き、`last_cleared_at_utc` を記録して `3日` のクールダウンを設ける。
- クールダウン中に `--clear-long-lived-cache` が再度指定された場合は、削除を実行せず splash と通常ログの両方に拒否理由を表示して終了してよい。
- CLI からの強制再実行オプションは持たず、どうしても必要な場合は cache directory の手動削除を README で案内する。

### 10.4 DEM / 建物キャッシュ再取得設計

- 目的は、毎回の再ダウンロードを避けつつ、長期間放置された DEM / 建物キャッシュが無期限に固定化されることを防ぐことである。
- DEM は更新頻度が低いため、`fresh=90日`、`stale>90日` とする。
- 地形地平線の距離帯は `0.5 / 1 / 2 / 4 / 8 / 16 / 32 / 64 / 128 km` を初期値として扱い、最遠帯までの計算は `128 km` を基準に行ってよい。
- Overture 建物由来 cache は DEM より更新頻度が高いため、通常 derived dataset / skyscraper tile ともに `fresh=30日`、`stale>30日` とする。
- 夜間光 cache は static versioned dataset として扱い、TTL ベースの fresh/stale 判定を持たず、未取得・破損・版切り替え時のみ再取得してよい。
- 上記の TTL とは別に、可能な場合は Overture release の差分確認を追加の freshness シグナルとして使ってよい。
- release 照合は起動時または都市アウトライン初回有効化時に行ってよいが、前回照合から `24時間` 以内なら省略してよい。
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
- Overture 建物 cache には各取得単位ごとに `overture_release` も持たせ、起動時の release 照合に使う。
- release 照合のネットワーク負荷を抑えるため、`last_release_check_utc` は `overture_buildings_release_check_meta.json` のような cache root 直下のローカル metadata に保存する。
- この方針により、fresh/stale 判定と再取得単位を一致させ、ディレクトリ構造の大規模変更を避ける。
- 後方互換のため、既存 cache に `fetched_at_utc` が無い場合は初回読込時の現在時刻を暫定 `fetched_at_utc` として補完し、その metadata を書き戻してよい。
- この補完時刻は元の実取得時刻ではなく移行時刻として扱う。
- そのため既存 cache は移行後の最初の TTL 周期だけ実際より新しく見えるが、大量再取得を避けるため許容してよい。

#### 10.4.1 DEM フロー

1. `TerrainHorizonController` が必要 tile 一覧を確定する。
2. 各 tile について sidecar metadata から `fetched_at_utc` を読む。無い場合は現在時刻で補完して書き戻す。
3. fresh なら既存 tile をそのまま使う。
4. stale なら既存 tile を入力に含めつつ、バックグラウンドで同じ key を再取得する。
5. 距離帯や描画上限を変えても、既存 tile のうち要求範囲に重なるものは再利用し、外周で新たに必要になった tile だけを追加取得してよい。
5. 再取得成功時は `*.tmp` + `replace()` で tile 本体と metadata を更新する。
6. 再取得失敗時は stale tile を使い続け、UI には stale 利用中であることを示してよい。

#### 10.4.1.1 水面フロー

1. 水面レイヤーは sea と inland water を別レイヤーとして扱ってよい。
2. sea は観測地点中心の sea-mask tile 群を読み、海域だけを点群化してよい。
3. inland water は OSM の河川、湖、運河、水路、riverbank 系輪郭から別途点群化してよい。
4. sea mask は `125m`、`250m`、`500m` の 3 段で読み分け、遠方 `500m` 帯は in-memory の代表点化で密度を下げてよい。
5. sea mask の `1` を水、`0` を地面として扱い、水ピクセル中心を観測点基準の `alt/az` に投影してよい。
6. `WaterOverlayState` は sea と inland の dot set を別保持し、terrain horizon の ON/OFF や DEM ready/fail に応じて active dots を切り替えてよい。
7. sea だけが先に計算できた場合は sea-only の中間更新を emit し、その後 inland 完了時に sea + inland の合成結果を emit してよい。
8. 水面ドットは `target_height_m = 0.0` の地表貼り付けとして扱ってよく、地形地平線が非表示のときは水面ドットも非表示にしてよい。
9. 点群化後の alpha は `water_overlay_opacity` を基準にし、距離に応じた指数減衰をかけてよい。
10. `WaterOverlayPoint` は `render/terrain.py` で小さな点として描画し、海と inland で色を分けてよい。
11. cache key は観測地点中心の `lat/lon` と距離帯、tile root から導出し、`bbox` はその派生値として扱ってよい。inland water の OSM footprint は `CACHE_PATH/water_overlay` 配下に JSON cache として保持してよい。
12. `Water Surface` の GUI トグルや `--water-surface-opacity 0` は初期表示の有無を変えてよく、tile cache を破棄する必要はない。

#### 10.4.2 Overture 建物フロー

1. `UrbanOutlineController` は通常 derived dataset と skyscraper tile cache をそれぞれ独立に fresh/stale 判定する。
2. `mode=both` の場合、`building` と `building_part` は別々に TTL 判定し、片方だけ stale ならその片方だけ再取得する。
3. 通常 derived dataset は dataset directory 単位で metadata を持ち、`cache_meta.json` などの sidecar file に `fetched_at_utc`、半径、mode、最小高さなどの cache key 情報を保持する。`fetched_at_utc` が無い旧 cache は、`cache_meta.json` の `mtime` を優先し、次に `bldg` directory の `mtime` を使って推定してよい。推定した時刻は書き戻してよい。
   `overture_release` も sidecar に保持し、release 照合時に同一 release かどうかを比較する。
   `overture_release` が無い、または sidecar の読込に失敗した場合は release 不明として TTL 判定に落とす。
4. skyscraper 補助レイヤーは tile directory 単位で metadata を持ち、seed tile ごとに stale 判定する。
5. fresh cache があればそれを即時読込する。
6. 起動時または都市アウトライン初回有効化時にネットワーク接続がある場合は、Overture の最新 release を確認してよい。
7. ただし前回照合から `24時間` 以内なら release 照合を省略してよい。
8. 取得済み cache の release と最新 release が異なる場合は、TTL 期限前でも stale 扱いにして再取得してよい。
9. stale cache があればそれを即時読込しつつ、欠けている dataset / stale dataset / release 差分がある dataset だけをバックグラウンドで `overturemaps download` し直す。
10. 再取得成功時は新 directory を一時パスに組み立て、整合性確認後に `replace()` 相当で切り替える。
11. 再取得失敗時は stale dataset を使い続け、次回要求時に再試行してよい。
12. release 照合の結果は `last_release_check_utc` と `last_seen_overture_release` として root metadata に保存し、次回起動時の 24 時間スキップ判定に使う。

#### 10.4.3 Overture release metadata 形式

- dataset / tile sidecar
  - `cache_meta.json` などの sidecar file に次を持たせてよい。
  - `dataset_name`: cache directory 名
  - `fetched_at_utc`: cache を最後に取得した UTC ISO 8601 string
  - `overture_release`: その cache を生成した Overture release version
  - `feature_type`: `building` または `building_part`
  - `min_building_height_m`: 最小高さフィルタ
  - `bbox`: `west` / `south` / `east` / `north`
  - `query_lat_deg`, `query_lon_deg`, `query_radius_km`: 生成元 query を復元するための補助情報
- root-level release check metadata
  - cache root 直下に `overture_buildings_release_check_meta.json` を置いてよい。
  - `last_release_check_utc`: 最新 release を最後に照会した UTC ISO 8601 string
  - `last_seen_overture_release`: 最後に確認した release version
  - `checked_source`: `stac` など照会元の識別子
  - `checked_success`: 最後の照会が成功したかどうか
  - `checked_error`: 任意の失敗理由
- 起動時の release 照合では、root-level metadata の `last_release_check_utc` を読み、`24時間` 未満なら照会を省略してよい。
- 照会成功時は root-level metadata の `last_release_check_utc` と `last_seen_overture_release` を更新し、必要なら対象 cache の `overture_release` も更新してよい。
- 照会失敗時は `checked_success=false` を記録してよいが、既存 cache の利用可否は TTL 側の判定に委ねてよい。
- root-level metadata や sidecar file が読めない場合でも cache 本体を直ちに削除せず、release 照合不能として TTL 判定へフォールバックしてよい。
- sidecar の `overture_release` が欠けている場合も cache を保持し、次回更新時に補完してよい。
- `fetched_at_utc` が欠けている legacy cache では、`cache_meta.json` の `mtime` を最優先の推定値として使ってよい。
- `cache_meta.json` が無い、または読めない場合に限り、`bldg` directory の `mtime` を次点の推定値として使ってよい。
- 推定した `fetched_at_utc` は書き戻してよいが、初回移行時の推定値であることを意識して扱ってよい。

#### 10.4.4 Overture release 照合フロー

1. `UrbanOutlineController` または `zstarview-export-image` 側が release 照合を開始する前に root-level metadata を読む。
2. `last_release_check_utc` が存在し、現在時刻との差が `24時間` 未満なら release catalog への問い合わせを省略する。
3. root-level metadata が無い、壊れている、または古い形式でも、照合対象 cache の削除は行わず TTL 判定に進む。
4. 照会が必要な場合は Overture STAC catalog か Python client の最新 release 情報を 1 回だけ問い合わせる。
5. 最新 release version が取得できたら、該当 cache の sidecar `overture_release` と比較する。
6. sidecar の `overture_release` が一致しない cache は TTL 期限前でも stale 扱いとし、その dataset / tile を再取得候補に入れる。
7. sidecar の `overture_release` が欠けている場合は release unknown とみなし、TTL 判定だけで fresh/stale を決める。
8. 照会成功時は root-level metadata を更新し、`last_release_check_utc` と `last_seen_overture_release` を書き戻す。
9. root-level metadata 更新は `.tmp` ファイルへ書いてから `replace()` してよい。
10. dataset / tile sidecar 更新は、再取得成功時に新しい payload と同時に `overture_release` を書き戻してよい。
11. 照会失敗時は root-level metadata の `checked_success=false` と `checked_error` を記録してよいが、既存 cache の利用は継続してよい。
12. 同一起動内で複数の urban-outline 更新要求が来た場合は、1 回の release 照合結果を共有してよい。

#### 10.4.5 実装モジュールの責務

- `src/zstarview/data/import_overture_buildings.py`
  - dataset / tile sidecar の読み書きを担当する。
  - `fetched_at_utc`、`overture_release`、cache key を同じ sidecar に保存する。
  - 再取得成功時に release version を sidecar へ反映する。
  - `lat/lon` を埋め込む dataset 名は、cache 再利用用の 4 桁丸め規則を使い、表示用の 5 桁丸めとは分離して扱う。
- `src/zstarview/cache_maintenance.py`
  - 長寿命 cache 削除とは独立した root-level release check metadata の読み書き補助を持ってよい。
  - ただし cache 削除ロジック自体は release 照合に依存させない。
- `src/zstarview/utils/latlon_format.py`
  - GUI 表示用の 5 桁フォーマッタと、建物 dataset 名 / cache 再利用判定用の 4 桁フォーマッタを持つ。
  - 緯度経度の小数桁数をそれぞれ 1 箇所で変更できるようにする。
- `src/zstarview/gui/urban_outline_controller.py`
  - 起動時 / 初回有効化時の release 照合をトリガーする。
  - release 照合の失敗は警告ログに留め、表示中の cache は TTL で扱う。
- `src/zstarview/cli/export_image.py`
  - headless 経路でも同じ release 照合補助を共有する。
  - GUI と同じ root-level metadata を読むことで、release チェックの 24 時間抑制を共通化する。

### 10.5 夜間光データ

- 夜間光データは NASA `Earth at Night/Black Marble: Flat Maps` の `2016 Grayscale` 500m tiled GeoTIFF を前提とする。
- 取得元は `science.nasa.gov` の flat maps ページに置かれた direct GeoTIFF リンクでよく、Earthdata 認証を要求しなくてよい。
- アプリは cache root 配下に versioned directory を持ち、初回利用時だけ対象 tile をダウンロードしてよい。
- tile 名は `A1` 〜 `D2` の 8 枚を前提にしてよい。
- 起動時に cache があれば再利用し、無ければバックグラウンドでダウンロードを開始してよい。
- ダウンロード失敗時は夜間光レイヤーだけ unavailable とし、他のレイヤーは継続してよい。
- 観測地点ごとの基礎 profile は重いので、`observer_lat/lon + terrain horizon` 単位でキャッシュしてよい。
- `sun_alt` は profile そのものには含めず、連続係数として描画時に掛けてよい。
- 画面上の効果は、地形地平線または副稜線の少し上に置く固定角度帯の glow で表現し、角度の中心は演出として固定してよい。
- 描画時の seam は観測者の真後ろを基準にし、`az=0°` で固定して切らない。
- 描画層は `CORE / MID_NEAR / MID_FAR / OUTER` の 4 段として重ね、`MID` と `OUTER` は段階的に薄くしてよい。
- glow の強度は、方位ごとのサンプルを距離帯ごとに区間積分した値とし、各 band にはその区間だけの brightness を割り当ててよい。
- 地形地平線や副稜線が無い場合は、水平線を fallback にしてよい。
- 500m グリッド由来の段差を避けるため、GeoTIFF 読込時に双一次補間や方位方向の平滑化を行ってよい。

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
