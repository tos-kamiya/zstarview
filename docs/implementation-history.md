# zstarview 実装履歴

最終更新: 2026-03-22

## 1. この文書の位置づけ

この文書は、`zstarview` の時系列の実装ノートを記録する。  
設計の正本ではなく、何をいつ実装したか、どこが継続課題かを残すための文書である。

記録対象は次の通り。

- 実装済み変更の要約
- TODO
- 実装上の判断メモ

利用者向け仕様は `docs/specification.md`、内部設計の正本は `docs/design.md` を参照する。

## 2. 書き方のルール

- 時系列を保つ
- 変更の背景と結果を短く残す
- 将来課題は `TODO` にまとめる
- 仕様として確定した内容は必要に応じて `design.md` または `specification.md` へ昇格させる

## 3. TODO

- 山頂ビュー用 dataset の作成フローを設計する
  - 目的は「山そのものの完全地理データ」ではなく、「山名から山頂ビュー用の代表点へ解決する curated dataset」を作ることとする。
  - 候補抽出は Wikipedia を起点に行い、地域代表性と知名度を重視して少数の山を選ぶ。
  - 正規化には Wikidata を使い、`qid`、多言語名、別名、座標、標高を補う。
  - 座標は厳密測量値ではなく、星空表示用の山頂ビュー代表点として妥当なものを採用する。
  - 初版は `mountain_viewpoints.json` を別 dataset とし、`tower_viewpoints.json` とは分けて管理する。
  - 将来的には mountain も tower と同じ viewpoint CLI から参照できる形を想定する。
- GUI 上で時刻を前後できるダイナミックなタイムシフト操作を追加する
- CLI 仕様と内部データ構造の対応表を設計書へ追加する
- 必要に応じて `CloudController` と `SkyWindow` の責務境界を再評価する
- 航空機オーバーレイの OpenSky 問い合わせ節約策を検討する
  - 日本の深夜帯では観測地点周辺の航空機が少ないため、固定 `5分` 更新より問い合わせ回数を減らせる可能性がある。
  - 候補は、夜間の更新間隔延長、`0機` 連続時のバックオフ、あるいは時刻帯ベースの可変更新間隔。
  - ただし、更新間隔を空けた分だけ `bbox` を無条件に広げる案は、国境越えや遠距離機体の取り込み増加で credit 効率を悪化させる可能性がある。
  - 当面は「bbox は保守的に固定し、更新間隔だけ可変にする」案を優先候補として扱う。

## 4. 実装履歴

### 2025-09-16

- 実装履歴文書の前身を作成
  - 設計書から実装メモを分離し始めた。

### 2026-02-27

- アステリズム対応
  - 補助線表示としてアステリズムを追加した。
  - ホバー時強調と複数所属時のローテーション表示を導入した。

### 2026-02-28

- 雲更新パイプライン分離
  - ソース取得と描画を分離した。
  - `SourceKey` と `RenderKey` を分け、視点変更時の再利用性を高めた。
  - partial coverage の欠損ティント表示を追加した。

### 2026-03-01

- 地形地平線オーバーレイ統合
  - Copernicus DEM を用いた地形地平線計算を本体へ統合した。
  - `Terrain Horizon` トグルと opacity 制御を導入した。
  - DEM 取得失敗時の status line 表示とセッション内自動再試行抑止を導入した。

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

### 2026-03-20

- メニューボタン表示の見直し
  - 右上メニューボタンの glyph は、縦三点案を試したうえでハンバーガー (`☰`) に戻し、共有ボタン font size は `18px` から `16px` へ下げた。
  - 視認性を改善しつつ、既存 UI の認知コストを増やさない方針にした。

- 共有 overscan 視野 (`--content-fov-deg`)
  - `--content-fov-deg` の既定値を `100°` とし、ウィンドウ端 `90°` 固定のまま全レイヤーの描画対象を広げられるようにした。
  - 背景グラデーション、sky disc、雲合成、雲ハッチ、地面ティント、ガイド線、都市アウトライン、航空機、明るい星の簡易表示経路も同じ overscan 視野へ揃えた。
  - 暗い星だけを旧来の `90°` 境界で追加的に落とす最適化は削除した。
  - `content_fov_deg` の外側には、ウィンドウ背景へ自然につなぐ短い背景フェード帯を持たせる方針にした。
  - その後、低仰角の横長ウィンドウで sky disc と雲レイヤーが `radius*2` の固定正方形に見えていた問題を修正し、実ウィンドウ geometry に沿って overscan を広げるようにした。

- リサイズ直後の stale disc 描画を抑止
  - `resizeEvent` 直後に sky-disc / cloud-disc の古い画像キャッシュを破棄し、新しい window state で再生成するようにした。
  - sky/cloud worker の payload には resize generation を持たせ、現在の window generation と render size に合わない結果は破棄するようにした。
  - これにより、リサイズ後に古い楕円状ディスク画像が遅れて上書きされる症状を防いだ。

- 海上観測地点向けの DEM fallback
  - Copernicus DEM の `nodata` は海面高度 `0.0m` として扱い、海上座標でも terrain horizon 計算を継続できるようにした。
  - 取得 bbox に DEM tile が 1 枚も無い場合も ocean-only horizon とみなし、失敗ではなく空 profile を返すようにした。

- GOES-East の参照先更新
  - GOES-East は `G16` ではなく `G19` を使うようにし、既存設定中の `G16` も `G19` として扱う互換処理を追加した。
  - README には、部分カバー時の黄色い欠損表示と `screenshot5` への参照を追加した。

- README の offline-use 節整理
  - README / README_ja-JP の「Slow or Unstable Network / Offline Use」を、cloud satellite imagery と terrain horizon で別項目に分けた。
  - disable 方法を機能ごとに読めるようにし、トラブルシュート導線を整理した。

- パッチ版更新
  - `1.2.0`、`1.2.1`、`1.2.2` の version bump を行った。

### 2026-03-21

- 直接依存ライブラリのライセンス互換性確認
  - 2026-03-21 時点の直接依存ライブラリについて、配布上のライセンス互換性に大きな問題がないかを確認した。

- Windows 起因の文字コード安全性ルールを明文化
  - `AGENTS.md` と `docs/design.md` に、terminal / console / log / CLI help / exception / subprocess 出力へ出る可能性がある文字列は ASCII-only を原則とする方針を追加した。
  - 非 ASCII を判定ロジックで扱う場合は、ソース中に直書きせず Unicode escape を優先する方針を追加した。

- 惑星暦カーネルを `de442s.bsp` へ更新
  - 既定の ephemeris filename を `de440s.bsp` から `de442s.bsp` へ切り替えた。

### 2026-03-22

- 直接座標入力の Google Maps URL 対応
  - `lat;lon` に加えて、`@lat,lon` と、現在広く観測される Google Maps shared URL 形式の一部を直接地点入力として受け付けるようにした。
  - 対象 URL は `maps.google.com/` または `www.google.com/maps/` で始まるものに絞り、`!3dLAT!4dLON` があればそれを優先し、なければ `@LAT,LON` を使うようにした。
  - Google Maps URL に含まれる zoom、高度風パラメータ、heading などは観測地点や観測者高さには使わない方針にした。

- 直接座標入力のタイムゾーン解決見直し
  - `lat;lon`、`@lat,lon`、Google Maps URL のいずれも、地点座標からタイムゾーンを自動解決するようにした。
  - これにより、従来の `lat;lon` 入力で暗黙に `UTC` 相当になっていた挙動は変わる。
  - 互換用途や明示上書きのために `--timezone` を追加した。

- `zstarview-export-image` のメタデータ出力追加
  - GUI 左上に表示される場所、地点高さ、観測者高さ、時刻、Alt/Az、`Vmag limit` を、export 画像には焼き込まず `stderr` へ出力するようにした。
  - `--sixel` 指定時は、端末画像の直前に同情報を出すようにした。

- マイナーバージョン更新
  - 直接座標入力のタイムゾーン挙動が後方互換ではないため、バージョンを `1.4.1` から `1.5.0` へ更新した。
  - `de442s.bsp` は Skyfield の既定 JPL download base では 404 になるため、runtime ではこのファイルだけ NAIF の明示 URL を使うようにした。
  - README の「初回起動時にダウンロードされる惑星暦データ」の説明も `de442s.bsp` に更新した。

- 雲データ取得の satpy 依存撤去
  - GOES は `CMIPF C13` NetCDF を `xarray` で直接読み、`goes_imager_projection` から `area` を再構築する経路へ置き換えた。
  - Himawari は Satpy reader をやめ、`ISatSS M1C13` タイル群を直接 stitch して `fixedgrid_projection` から `area` を再構築する経路へ置き換えた。
  - これに合わせて `satpy` の直接依存を削除し、雲データ処理の依存は `xarray` と `pyresample` ベースへ整理した。

- 雲データ取得の `pyresample` 依存撤去
  - GOES/Himawari で使っていた `pyresample.AreaDefinition` は、実際には geostationary 投影定義と `lon/lat -> pixel` 変換の最小機能しか使っていなかった。
  - そのため、同等の最小 API を持つ内部 `GeoArea` へ置き換え、雲ディスクの sampler と warm-threshold 推定はそのまま動くようにした。
  - これにより、`pyresample` の配布事情に引きずられずに cloud-disc 機能を維持できる構成になった。
  - `dev-samples/compare_pyresample_cloud_disc.py` を使って Himawari / GOES の複数ケースを比較し、確認した範囲では `bt_warm`、`bt_cold`、BT 配列、可視 mask、最終画像の差分はすべて `0` だった。

- Windows Arm64 でのインストール再確認
  - `pyresample` とその連鎖依存を削除した後、Windows Arm64 環境で `pipx install` が通ることを確認した。
  - 一方で、実行時には Windows Security の Smart App Control により Python 拡張モジュールの読込がブロックされる場合があり、そのケースでは設定変更で回避できることを README のトラブルシュートに追記した。

- Himawari warm-threshold の安定化
  - Himawari の部分タイル最適化後、observer 周辺タイルだけでは赤道帯サンプルが消え、`eq_samples=0` により warm-threshold が固定 fallback へ落ちる問題が出た。
  - このため、observer 周辺の描画用タイルに加えて、観測経度付近の赤道帯を横切る少数タイルを追加取得するようにした。
  - これにより、全 `88` タイル取得へ戻さずに warm-threshold 推定を復旧し、雲が白っぽく出すぎる症状を抑えた。

- 描画パイプラインの shared pipeline 化
  - 将来の単発画像書き出し CLI を見据え、`gui/window_render.py` の Mixin に集まっていた描画知識を `src/zstarview/render/pipeline.py` へ切り出した。
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

- 航空機データの短寿命ディスクキャッシュ導入
  - `zstarview-export-image` を短時間に連続実行した場合でも OpenSky API 問い合わせ回数を抑えるため、航空機レイヤーに `bbox` 単位の短寿命永続キャッシュを導入した。
  - キャッシュ対象は OpenSky から正規化した `AircraftSnapshot` 列、取得時刻、観測地点由来 `bbox`、source 名とした。
  - fresh 判定は既存の `AIRCRAFT_REFRESH_INTERVAL_SECONDS` と同じ `5分` を使い、stale fallback は `10分` とした。
  - GUI と `zstarview-export-image` の両方が同じ cached fetch 経路を使うようにした。
  - fresh cache が存在する場合は API 問い合わせを省略し、再取得失敗時のみ stale cache を警告付きで再利用する。
  - `--aircraft-opacity 0` のセッションでは、現行仕様どおり問い合わせもキャッシュ読込も必須ではない。
  - 実装は `aircraft/opensky.py` に直接混ぜず、キャッシュ責務を `aircraft/cache.py` へ分離した。
  - テストでは、fresh hit、stale hit 後の再取得成功、再取得失敗時の stale fallback、古い cache file の cleanup、export-image 経路への影響有無を確認した。

- 人工衛星データ取得とキャッシュ基盤の追加
  - Skyfield `EarthSatellite.from_omm()` を前提に、CelesTrak GP JSON を取得して group 単位でキャッシュする `satellites/` パッケージを追加した。
  - 初期対象 group は `ISS` と `Starlink` とし、CelesTrak group 名 `stations`、`starlink` へ解決するようにした。
  - 軌道要素 cache は `group_key` 単位の少数 JSON file とし、fresh 判定は `24時間` とした。
  - fresh を外れた cache は再取得優先とし、初版では stale fallback を行わない方針にした。
  - 取得結果は raw OMM JSON のまま永続保存し、runtime で Skyfield `EarthSatellite` へ変換する構成にした。
  - テストでは、CelesTrak URL 組み立て、OMM payload 正規化、Skyfield `EarthSatellite` 生成、fresh cache hit、stale 時の再取得失敗伝播、cache overwrite、cleanup を確認した。

- 人工衛星レイヤーの位置計算と描画統合設計
  - 人工衛星の位置計算は `satellites/project.py` に分離し、group ごとの raw OMM records を Skyfield `EarthSatellite.from_omm()` で satellite object へ変換して、観測地点基準の `alt/az` を計算する方針にした。
  - 描画用内部モデルとして `SatelliteOverlayPoint` を導入し、初期実装では軌跡線なしの「現在位置 marker のみ」を保持する構成にした。
  - 描画順は `planets` の後、`aircraft` の前とし、人工衛星は惑星・月より近く、航空機よりは遠い層として扱う。
  - marker は航空機と同系統の紫色、小型クロス、`ISS` だけ少し大きい scale を持つ方針にした。
  - GUI 側は `SatelliteController` と `SatelliteState` を追加し、軌道要素は `24時間`、位置再計算は `5秒` で更新する前提にした。

- 人工衛星レイヤーの位置計算と描画統合
  - `satellites/project.py` を追加し、cached OMM records から Skyfield `EarthSatellite` を生成して観測地点基準の `alt/az` を計算し、描画用 `SatelliteOverlayPoint` 列へ射影するようにした。
  - 初版の表示対象は `ISS` を既定有効 group とし、marker は航空機と同系統の紫色小型クロス、`ISS` だけ少し大きい scale とラベル表示を持つ構成にした。
  - 描画順は `planets` の後、`aircraft` の前へ組み込み、render pipeline と export-image の両方で同じ人工衛星レイヤーを描けるようにした。
  - GUI 側に `SatelliteController` と `SatelliteState` を追加し、軌道要素 fetch は `24時間`、位置再計算は `5秒` の timer で分離した。
  - GUI メニューに `Satellites` トグルを追加し、静止時のみ有効化する既存航空機トグルに近い挙動へそろえた。
  - テストでは、group 順序と marker 属性を含む射影結果、render order、export-image 経路、および周辺 GUI 同期への影響を確認した。
