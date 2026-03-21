# zstarview 実装履歴

最終更新: 2026-03-21

## 1. この文書の位置づけ

この文書は、`zstarview` の時系列の実装ノートを記録する。  
設計の正本ではなく、何をいつ実装したか、どこが継続課題か、現在何を進めているかを残すための文書である。

記録対象は次の通り。

- 実装済み変更の要約
- TODO
- INPROGRESS
- 実装上の判断メモ

利用者向け仕様は `docs/specification.md`、内部設計の正本は `docs/design.md` を参照する。

## 2. 書き方のルール

- 時系列を保つ
- 変更の背景と結果を短く残す
- 将来課題は `TODO` にまとめる
- 着手中の項目は `INPROGRESS` に置く
- 仕様として確定した内容は必要に応じて `design.md` または `specification.md` へ昇格させる

## 3. INPROGRESS

- 山頂ビュー用 dataset の作成フローを設計する
  - 目的は「山そのものの完全地理データ」ではなく、「山名から山頂ビュー用の代表点へ解決する curated dataset」を作ることとする。
  - 候補抽出は Wikipedia を起点に行い、地域代表性と知名度を重視して少数の山を選ぶ。
  - 正規化には Wikidata を使い、`qid`、多言語名、別名、座標、標高を補う。
  - 座標は厳密測量値ではなく、星空表示用の山頂ビュー代表点として妥当なものを採用する。
  - 初版は `mountain_viewpoints.json` を別 dataset とし、`tower_viewpoints.json` とは分けて管理する。
  - 将来的には mountain も tower と同じ viewpoint CLI から参照できる形を想定する。

## 4. TODO

- GUI 上で時刻を前後できるダイナミックなタイムシフト操作を追加する
- CLI 仕様と内部データ構造の対応表を設計書へ追加する
- 必要に応じて `CloudController` と `SkyWindow` の責務境界を再評価する
- 航空機オーバーレイの OpenSky 問い合わせ節約策を検討する
  - 日本の深夜帯では観測地点周辺の航空機が少ないため、固定 `5分` 更新より問い合わせ回数を減らせる可能性がある。
  - 候補は、夜間の更新間隔延長、`0機` 連続時のバックオフ、あるいは時刻帯ベースの可変更新間隔。
  - ただし、更新間隔を空けた分だけ `bbox` を無条件に広げる案は、国境越えや遠距離機体の取り込み増加で credit 効率を悪化させる可能性がある。
  - 当面は「bbox は保守的に固定し、更新間隔だけ可変にする」案を優先候補として扱う。

## 5. 実装履歴

### 2026-03-21

- Windows 起因の文字コード安全性ルールを明文化
  - `AGENTS.md` と `docs/design.md` に、terminal / console / log / CLI help / exception / subprocess 出力へ出る可能性がある文字列は ASCII-only を原則とする方針を追加した。
  - 非 ASCII を判定ロジックで扱う場合は、ソース中に直書きせず Unicode escape を優先する方針を追加した。

- 惑星暦カーネルを `de442s.bsp` へ更新
  - 既定の ephemeris filename を `de440s.bsp` から `de442s.bsp` へ切り替えた。
  - `de442s.bsp` は Skyfield の既定 JPL download base では 404 になるため、runtime ではこのファイルだけ NAIF の明示 URL を使うようにした。
  - README の「初回起動時にダウンロードされる惑星暦データ」の説明も `de442s.bsp` に更新した。

- 雲データ取得の satpy 依存撤去
  - GOES は `CMIPF C13` NetCDF を `xarray` で直接読み、`goes_imager_projection` から `area` を再構築する経路へ置き換えた。
  - Himawari は Satpy reader をやめ、`ISatSS M1C13` タイル群を直接 stitch して `fixedgrid_projection` から `area` を再構築する経路へ置き換えた。
  - これに合わせて `satpy` の直接依存を削除し、雲データ処理の依存は `xarray` と `pyresample` ベースへ整理した。

- Himawari warm-threshold の安定化
  - Himawari の部分タイル最適化後、observer 周辺タイルだけでは赤道帯サンプルが消え、`eq_samples=0` により warm-threshold が固定 fallback へ落ちる問題が出た。
  - このため、observer 周辺の描画用タイルに加えて、観測経度付近の赤道帯を横切る少数タイルを追加取得するようにした。
  - これにより、全 `88` タイル取得へ戻さずに warm-threshold 推定を復旧し、雲が白っぽく出すぎる症状を抑えた。

- 描画パイプラインの shared pipeline 化
  - 将来の単発画像書き出し CLI を見据え、`ui/window_render.py` の Mixin に集まっていた描画知識を `src/zstarview/render/pipeline.py` へ切り出した。
  - 描画入力は `RenderSceneData`、`RenderStyle`、`RenderHudState` に分離し、旧 `RenderPipelineState` は廃止した。
  - shared pipeline の関数群は `geometry`、`viewport_rect`、`scene`、`style`、`hud` を直接受ける形に揃えた。
  - `SkyWindowRenderMixin` は、`paintEvent()` 本線、scene/style/hud 組み立て、frame cache、jump highlight、hover 解決など GUI 固有処理に絞った。

- 描画順とガイドレイヤーの整理
  - 方位ラベルと天頂マーカーは `guide` として独立レイヤー化し、`sky-cloud` 合成の直後に置くようにした。
  - DSO hover は terrain 段から外し、通常描画では後段オーバーレイ側へ寄せた。

- hover/HUD とベース描画の分離
  - shared pipeline を `render_base_scene_into_painter()` と `render_hud_overlay_into_painter()` に分け、`paintEvent()` はベースをキャッシュ描画した後に HUD を都度重ねる構成にした。
  - `guide` はベース描画側に残し、方位ラベルのマウス回避はやめて安定ガイド扱いにした。
  - ベースフレームの cache key から `mouse_pos`、hover 対象名、jump highlight 名、status message を外した。

- アステリズム強調と月拡大の後段上書き化
  - 通常のアステリズム線と通常サイズの月はベース描画に残し、hover 時の強調だけを HUD 側で上書きする構成へ切り替えた。
  - 月の `5x` 拡大は、角半径の生値ではなく通常時の見た目半径を基準に適用するよう修正した。
  - 月 hover についても `name == "moon"` を持つ hover object なら拡大上書きに入るようにした。

- CLI 引数 builder の分割
  - `cli_args.py` の parser 構築を、地点、dataset query、時刻、描画の helper 関数群へ分割した。
  - main app 側は `build_main_argument_parser()` を使う形に整理した。
  - dataset query の整合性検証は、その parser に存在するオプションだけを見るようにして、将来の画像書き出し CLI parser でも再利用できるようにした。

- 単発画像書き出し CLI の追加
  - `zstarview-export-image` を追加し、地点・時刻・視線・描画オプションの大半を既存 `zstarview` CLI と共有する形で 1 枚の PNG を書き出して終了できるようにした。
  - 画像書き出し固有オプションとして `--output`、`--image-size`、`--layer-timeout-seconds`、`--allow-partial-data` を追加した。
  - `--layer-timeout-seconds` の既定値は、Himawari 分割ダウンロードを考慮して `30` 秒から `90` 秒へ引き上げた。
  - 将来拡張として、`--sixel` を `--output` と併用可能な端末出力オプションとして扱い、保存を先に行ってから `img2sixel -` へ PNG bytes を流すパイプ方式で端末表示を試みる方針を整理した。
  - `--sixel` を実装し、`--output` と併用可能にした。`--sixel` 指定時は事前に `img2sixel` の存在を確認し、生成済み `QImage` を PNG bytes 化して `img2sixel -` の stdin へ流す。
  - `-o -` を実装し、PNG bytes を stdout へ直接流せるようにした。stdout 競合を避けるため、`--output -` と `--sixel` の併用は弾く。
  - GUI 専用の `--sky-update-interval` と `--window-geometry`、dataset 参照専用 CLI オプションは export parser に載せない。
  - 実装では `SkyWindow` と Qt signal ベースの controller 群を使わず、sky/cloud/terrain/urban/aircraft を同期的に順番に取得してから、hover/HUD なしのベース描画を `QImage` へ保存する。
  - `opacity == 0` のレイヤーは取得キューと timeout 待機対象から外す。
  - export 画像では、地点名・時刻などの静的 overlay 情報と、FOV 外の GUI 向け背景グラデーションを描かないようにした。

### 2026-03-19

- OpenSky 航空機オーバーレイ
  - OpenSky `states/all` を使った航空機オーバーレイを本体へ追加した。
  - 取得は `5分` 間隔、表示上の予想再投影は `2秒` 間隔とした。
  - 機体は当初 `2秒前 -> 現在 -> 2秒後` のオレンジ折れ線で描画し、近い機体ほど線を太くした。
  - 後に視認性を上げるため、折れ線は `4秒前 -> 現在 -> 4秒後` へ延長した。
  - さらに近距離機体の大きな折れ角を避けるため、折れ線は `-4, -2, 0, +2, +4秒` の 5 点サンプルを結ぶ方式へ変更した。
  - 観測地点から `50km` を超える機体は描かず、`10km` 以内かつ `90秒` 以内の機体だけ `callsign` を表示するようにした。
  - `90秒` を超えた機体は次回更新まで徐々に薄くし、予想表示であることを視覚的に示すようにした。
  - bbox は南北 `±1.0°` を基準にしつつ、東西方向は緯度依存で最低 `90km` を確保し、OpenSky の `1 credit` 帯から外れないよう `25 square degrees` 以下へクリップする構成にした。

- 航空機レイヤーの opacity / toggle 制御
  - CLI に `-a` / `--aircraft-opacity` を追加し、`0.0` のときは起動時の航空機問い合わせ自体を行わないようにした。
  - GUI には `Aircraft` トグルを追加し、`I` ショートカットで表示・非表示を切り替えられるようにした。
  - 航空機レイヤーの再表示時は `last_success_utc` を見て fresh cache を優先し、`5分` 以内の保持データがあれば API 再問い合わせなしで再投影だけ行うようにした。
  - 描画側には航空機レイヤー opacity を追加し、折れ線と `callsign` ラベルの alpha に乗算するようにした。

- 雲オーバーレイの地平線近傍カットオフ見直し
  - `CloudDiscConfig.alt_min_deg` を `0.0` から `2.5` へ引き上げた。
  - 地平線ぎりぎりの遠距離雲を抑え、観測補助として意味の薄い領域の取得・表示を減らす方向にした。
  - その後、既定値を `3.0` に再調整し、水平線直上の雲をさらに控えめに扱うようにした。

### 2026-03-17

- スカイスクレーパー遠距離都市アウトライン補助レイヤー
  - 通常の `0-2.5km` urban outline とは別に、`2.5-10km` の遠距離スカイスクレーパー補助レイヤーを追加した。
  - 遠距離レイヤーは `150m` 以上の `building` のみを対象とし、`building_part` は使わない。
  - seed source は Scraperbase を基にした 30 都市 curated list で、Web Mercator `zoom 14` の `z/x/y` tile と bbox を同梱 JSON として持つ構成にした。
  - `zoom 14` は Overture 固有タイルではなく、`overturemaps download --bbox=...` に戻すための便宜的な scan 単位として採用した。
  - `UrbanOutlineController` は通常レイヤー更新と同じ `update()` サイクルで skyscraper tile を選定・取得・解決し、結果を 1 回の `urban_ready` で UI に反映する。
  - skyscraper tile の取得はオンデマンドとし、`overture_skyscrapers/<tile-cache-key>/bldg` 配下へ永続キャッシュする。
  - skyscraper tile 取得が失敗した場合は通常レイヤーだけで描画を継続する safe fallback を入れた。
  - 確認用に `--urban-outline-skyscraper-only` を追加し、通常 `0-2.5km` レイヤーを描かずに遠距離 skyscraper レイヤーだけを描画できるようにした。

### 2026-03-14

- 都市アウトラインの Overture パイプライン化
  - 都市アウトラインの既定 runtime ソースを bundled PLATEAU derived data から `CACHE_PATH/overture_buildings` へ切り替えた。
  - `UrbanOutlineController` と `UrbanOutlineState` を追加し、起動時またはトグル再有効化時に Overture の取得またはキャッシュ読込をバックグラウンドで開始するようにした。
  - `lat/lon + radius + min_height + feature_type` を正規化したキャッシュキーを導入した。
  - `overturemaps` CLI を呼ぶ `import_overture_buildings.py` と `zstarview-import-overture-buildings` を追加し、ダウンロード結果を既存 derived tile 形式へ変換できるようにした。

- 都市アウトライン用建物高さしきい値の見直し
  - runtime 側の都市アウトライン読込では 40m 再フィルタをやめ、derived tile に残っている建物をそのまま使うようにした。
  - 建物高さしきい値は import / 生成時の `--min-building-height-m` に集約し、低層中心都市では `0` を指定できるようにした。

- PLATEAU 実装の撤去
  - 旧 `plateau_derived` 同梱データ、CityGML 取込スクリプト、PLATEAU 名の runtime helper を削除し、都市アウトライン経路を Overture に一本化した。
  - derived tile 読込や tile index 生成のユーティリティは、データソース非依存の汎用名へ整理した。

- 都市アウトライン描画の高さ連動 alpha
  - 都市アウトラインは建物ごとに `height_m` を保持する runtime 形式へ変更した。
  - `50m` 以上は CLI 指定 opacity をそのまま使い、`0m` ではその `25%` まで線形に下げる描画則を追加した。
  - 旧 `list[list[(alt, az)]]` 形式を受ける後方互換コードと古いテスト前提を削除した。

### 2026-03-09

- ドキュメント整理
  - `specification.md` を利用者視点の機能仕様へ整理した。
  - `design.md` をアーキテクチャ、責務、データ構造、処理フロー中心へ整理した。
  - `implementation-history.md` を時系列ノート専用に整理し、TODO と INPROGRESS を明示した。

- アステリズム線の FOV クリップ調整
  - アステリズム用のクリップ境界を独立定数化した。
  - 片端だけが視野内にある線分も表示できるようにした。
  - 線分は広い表示境界で円形クリップし、境界外へはみ出す部分を描かないようにした。

- キーボード視線変更時の簡易表示モード
  - 方向キー操作とウィンドウリサイズ中は `viewport_interaction_mode` で `render_view_center` を即時更新し、全レイヤー再計算ではなく再投影ベースで追随するようにした。
  - 最後の入力から 0.7 秒間は、`vmag <= 4.0` の明るい星とガイドライン、地形地平線、方位ラベル、天頂マーカーを表示する。
  - この簡易モード中は惑星、DSO、アステリズムを表示しない。
  - 簡易モード中の明るい星は新しい FOV に合わせて同期的に再計算する。
  - 補助線は worker 時点の正規化点列ではなく、`(alt, az)` サンプル列から描画時に `render_view_center` 基準で投影する。
  - 簡易モード終了時に初めて、全星等の星空、sky disc、雲、地形地平線の通常更新を開始する。

### 2026-03-20

- 共有 overscan 視野 (`--content-fov-deg`)
  - `--content-fov-deg` の既定値を `100°` とし、ウィンドウ端 `90°` 固定のまま全レイヤーの描画対象を広げられるようにした。
  - 背景グラデーション、sky disc、雲合成、雲ハッチ、地面ティント、ガイド線、都市アウトライン、航空機、明るい星の簡易表示経路も同じ overscan 視野へ揃えた。
  - 暗い星だけを旧来の `90°` 境界で追加的に落とす最適化は削除した。
  - `content_fov_deg` の外側には、ウィンドウ背景へ自然につなぐ短い背景フェード帯を持たせる方針にした。
  - その後、低仰角の横長ウィンドウで sky disc と雲レイヤーが `radius*2` の固定正方形に見えていた問題を修正し、実ウィンドウ geometry に沿って overscan を広げるようにした。

- タワー参照専用 CLI
  - `--list-towers` で一覧表示名を標準出力へ出して終了するようにした。
  - `--list-tower-names` で英語名以外を含む入力候補名一覧を標準出力へ出して終了するようにした。
  - `--show-tower-json NAME` で既存のタワー解決規則を使って 1 件解決し、JSON を標準出力へ出して終了するようにした。
  - これら 3 オプションは相互排他とし、位置引数や描画・時刻指定オプションとは併用不可にした。
  - この経路では GUI を起動せず、GeoNames 読込や設定保存も行わない。
  - ダイアクリティカルマーク付き主表示名には ASCII フォールバック名を別フィールド `ascii_name` として算出し、`--list-towers` とタワー名解決で利用するようにした。

- mountain viewpoint 参照専用 CLI
  - `--list-mountains` で一覧表示名を標準出力へ出して終了するようにした。
  - `--list-mountain-names` で英語名以外を含む入力候補名一覧を標準出力へ出して終了するようにした。
  - `--show-mountain-json NAME` で mountain viewpoint 解決規則を使って 1 件解決し、JSON を標準出力へ出して終了するようにした。
  - この経路では GUI を起動せず、GeoNames 読込や設定保存も行わない。
  - ダイアクリティカルマーク付き主表示名には ASCII フォールバック名を別フィールド `ascii_name` として算出し、`--list-mountains` と mountain viewpoint 名解決で利用するようにした。
  - 通常起動の `location` 引数でも mountain viewpoint 名とその `wikidata:Q...` 指定を受け付けるようにした。
  - 通常起動の `location` 引数に `t/NAME` と `m/NAME` を追加し、メインウィンドウの地点名表示と保存キーも `t/...` / `m/...` 形式に揃えた。

- viewpoint 参照 CLI の一本化
  - 旧 `--list-towers` / `--list-mountains` 系を廃止し、`--list-viewpoints KIND`、`--list-viewpoint-names KIND`、`--show-viewpoint-json NAME` に統一した。
  - 一覧出力は `t/NAME` / `m/NAME` の prefix を常に付け、通常起動の `t/NAME` / `m/NAME` ルールと揃えた。
  - `--show-viewpoint-json` は prefix 付き入力ならその kind だけを解決し、prefix なし入力で tower / mountain の両方に exact match がある場合は曖昧一致エラーにして候補名を列挙するようにした。
  - 文書上も mountain viewpoint dataset の出典を、Wikipedia 候補収集 + Wikidata 正規化の流れとして README と docs に明記した。

- viewpoint 名クリーニングの整理
  - ASCII フォールバックは、ダイアクリティカルマーク除去後に英字を含むものだけを採用するようにした。
  - その結果、`=` や `==` のような記号だけの fallback は一覧と解決候補から除外した。
  - mountain の別名では、長文フレーズ、数字だけ、絵文字だけのノイズ値を除外した。
  - mountain の別名では、`Mt.` / `Mtn.` を `Mount` / `Mountain` に寄せ、`/` 前後に空白を含む表記は候補一覧から外した。
  - `Aoraki / Mount Cook` は主表示名としては残しつつ、別名候補としては `Aoraki/Mount Cook` を残す方針にした。
  - tower の短縮 alias は個別判断とし、`i360` は残し、`138` のような数字だけの短縮は許可しない方針にした。

### 2026-03-04

- 高密度恒星表示の操作応答改善
  - インタラクション中の更新とアイドル後更新を分離した。
  - 星カタログ前処理データを保持し、操作中の負荷を下げた。
  - クラウド更新との干渉を抑え、描画割り込みを減らした。

- 星カタログ生成の正本化
  - Hipparcos + Tycho-2 を統合したカタログ生成フローを整理した。
  - `vmag <= 10` 運用に対応した。

- 全天球の円形ビュー化と地平線下ティント表示
  - ビューを円ディスク前提に統一した。
  - 横長ウィンドウ時の円配置を整理した。
  - 地平線下に地面色ティントを導入した。

- 恒星描画の高速化
  - 恒星を矩形ベースで描画する方式に整理した。
  - NumPy バッファへの集約で描画負荷を下げた。

- 雲の地平線下マージン整理
  - 雲表示は観望補助用途とし、地平線下マージンを設けない方針へ統一した。

- 地平線方位マーカーの見直し
  - 画面固定の縦線ではなく、地平線に対して垂直な固定長マーカーへ変更した。

- 惑星表示スタイルの更新
  - emoji ベース表示を廃止し、専用の円マーカーとブルーム表現へ変更した。
  - 惑星名を常時表示にした。

### 2026-03-01

- 地形地平線オーバーレイ統合
  - Copernicus DEM を用いた地形地平線計算を本体へ統合した。
  - `Terrain Horizon` トグルと opacity 制御を導入した。
  - DEM 取得失敗時の status line 表示とセッション内自動再試行抑止を導入した。

### 2026-02-28

- 雲更新パイプライン分離
  - ソース取得と描画を分離した。
  - `SourceKey` と `RenderKey` を分け、視点変更時の再利用性を高めた。
  - partial coverage の欠損ティント表示を追加した。

### 2026-02-27

- アステリズム対応
  - 補助線表示としてアステリズムを追加した。
  - ホバー時強調と複数所属時のローテーション表示を導入した。

### 2025-09-16

- 実装履歴文書の前身を作成
  - 設計書から実装メモを分離し始めた。
