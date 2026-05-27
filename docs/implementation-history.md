# zstarview 実装履歴

最終更新: 2026-05-27

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
- 実装履歴は末尾に追記し、新しい項目ほど後ろに並ぶようにする
- 変更の背景と結果を短く残す
- 将来課題は `TODO` にまとめる
- 仕様として確定した内容は必要に応じて `design.md` または `specification.md` へ昇格させる

## 3. TODO

- GUI 上で時刻を前後できるダイナミックなタイムシフト操作を追加する
- CLI 仕様と内部データ構造の対応表を設計書へ追加する
- 必要に応じて `CloudController` と `SkyWindow` の責務境界を再評価する
- 航空機オーバーレイの OpenSky 問い合わせ節約策を検討する
  - 日本の深夜帯では観測地点周辺の航空機が少ないため、固定 `5分` 更新より問い合わせ回数を減らせる可能性がある。
  - 候補は、夜間の更新間隔延長、`0機` 連続時のバックオフ、あるいは時刻帯ベースの可変更新間隔。
  - ただし、更新間隔を空けた分だけ `bbox` を無条件に広げる案は、国境越えや遠距離機体の取り込み増加で credit 効率を悪化させる可能性がある。
  - 当面は「bbox は保守的に固定し、更新間隔だけ可変にする」案を優先候補として扱う。
- 都市アウトラインの仰角計算を DEM 地盤標高基準へ切り替える
  - 現状は Overture の `height_m` をそのまま観測者高さと比較しており、建物側と観測地点側の地盤標高差を無視している。
  - 今後は `building_top_masl = building_ground_dem_m + height_m`、`observer_masl = observer_ground_dem_m + observer_height_m` を正本とする。
  - `building_part` の `min_height` / `height` の扱い、および derived tile へ `ground_elevation_m` を永続化するか runtime 解決に留めるかは実装時に再検討する。

## 4. 2026-05-15

- Water surface classification and dependency cleanup
  - Added a shared water-category classifier to `src/zstarview/water_overlay.py` so sea, river-like, and lake-like water footprints use the same tag-based rule in core code and dev samples.
  - Kept river-like surfaces as the only diagnostic 300m grid split target.
  - Removed `shapely` and `mapbox-earcut` from the project dependency set and kept water geometry handling in pure Python.
  - Updated the dev-sample SVG previews and the water-run scan demo to avoid geometry-library dependencies outside SVG path generation.

## 4. 実装履歴

### 2026-05-22

- 雲の取得と投影の分離
  - 気象衛星の更新を source fetch と projection に分けた。
  - source fetch 中は `Downloading`、projection 中は `Projecting` を表示できるようにした。
  - `cloud_source_ready` を導入し、取得完了後の投影更新を他の更新処理と分離した。
  - source / render の latest-wins を別系統で扱うようにし、古い視点の再投影結果を捨てやすくした。

### 2026-05-23

- ステータス行の白系テーマ文字色調整
  - `white` のステータス文字を `day` と同じグレーへ揃え、alpha を 255 にしてより不透明にした。
  - ステータス行の可読性を保ちながら、`white` だけが `day` よりはっきり見えるようにした。
- `white` 通常文字の `day` 統一
  - `white` の通常文字を `day` と同じ RGB に戻し、`white` 側だけ alpha を 255 にした。
  - 色相や明るさの差をなくし、`white` は不透明度だけで区別する方針に揃えた。
- ステータス行と通常文字の白系テーマアウトライン透明化
  - `day`、`white`、`black` の `ThemeStyle.text.outline_rgba` と `ThemeStyle.status_text.outline_rgba` を `night` と同じ黒の outline に揃えた。
  - 文字色はそのまま維持し、アウトラインだけを `night` 基準へ統一して見た目の差を減らした。

- version bump
  - `__version__` を `1.28.3` に上げた。

### 2026-05-17

- version bump
  - `__version__` を `1.27.0` に上げた。

### 2026-05-17

- 追加高さオプションの命名整理
  - 既存の `--observer-height-m` を、観測基準の上に足す追加高さという意味で整理し、仕様上の主名を `--height-add-m` に寄せる方針をまとめた。
  - `--use-building-top` は観測基準面を建物上端へ切り替えるモードとして維持する方針にした。
  - GUI の地点表示は `Height: ground ..., [building ... , ]add ...` の compact summary で最終高度の内訳が分かるようにした。

- 地点要約の短縮表示
  - 地点要約を `Height: ground ..., [building ... , ]add ...` に短縮し、`building` は 0 の場合だけ省く形へ整理した。
  - 主ウィンドウ、スプラッシュ、export の要約表示を同じ compact 形式に揃えた。

### 2026-05-13

- Water surface naming and docs sync
  - Renamed the user-facing water layer label to `Water Surface` and made `--water-surface-opacity` the documented CLI option.
  - Updated the specification and design docs to describe the current point-cloud water layer, including coastline handling, sea-level / DEM variants, and the GUI menu / status line terminology.

### 2026-05-12

- 起動ログのメインウィンドウ移行
  - 起動時のログ表示を splash ではなくメインウィンドウ内の一時ログオーバーレイへ切り替えた。
  - ロケーション解決・長寿命キャッシュ処理・初期データロードのログを、ウィンドウ表示直後から同じオーバーレイへ流すようにした。
  - `SkyWindow` は初期ロードを遅延できるようにし、解決後に `ViewerData` と `delta_t` を差し替えてから初期ロードを開始する形にした。
  - ロケーション解決失敗時は、例外種別を含む起動エラーをログとしてオーバーレイへ出し、短時間表示して終了するようにした。

- 起動失敗時のログ保持と自動終了切替
  - 起動失敗時はメインウィンドウ内の起動ログをそのまま残し、利用者が内容を確認できるようにした。
  - 起動失敗時だけ従来の自動クローズへ戻したい場合のために、`--close-on-startup-error` を追加した。
  - 起動ログオーバーレイは warning 以上を色分けして表示するようにした。
  - 起動失敗のときは、成功時と異なり起動ログハンドラをすぐ外さず、アプリ終了時までメインウィンドウ内に表示を残すようにした。
  - 起動失敗時は右上のハンバーガーメニューを前面へ上げ、ログを見ながらメニュー操作できるようにした。

- version bump
  - `__version__` を `1.24.2` に上げた。

### 2026-05-03

- `transparent` テーマの背景フラット化
  - `transparent` の `window_background` を平坦化し、`delta_rgb` を `0` にした。
  - `outer_alpha`、`edge_alpha`、`inner_rgba` の alpha を揃えて、半径方向の alpha 変化をなくした。
  - 背景描画側ではテーマ名ではなく `WindowBackgroundStyle.flat_background` を見て、単色の半透明塗りへ切り替えるようにした。
  - 回帰テストとして、`transparent` の背景が中心と隅で同一 RGBA になることを追加した。

- テーマ依存パラメータの追加整理
  - `WindowBackgroundStyle.draw_outer_border` を追加し、白系テーマの外枠描画をテーマ設定へ移した。
  - `ThemeStyle.star_visibility_boost` を追加し、`viewer` と `export_image` の `white` / `day` 分岐をテーマ値参照へ置き換えた。
  - `WindowChromeStyle.menu_fill_rgba` は各テーマの明示値として定義する形に整理した。
  - `RenderStyle.theme` を追加し、`window_render` と `export_image` で解決したテーマ実体を `pipeline` へ渡すようにした。
  - `draw_radial_background()` と `draw_window_border()` は `ThemeStyle` を直接受け取るようにして、描画側からのテーマ引き直しを外した。
  - これで、レンダラー側のテーマ名による条件分岐をさらに減らした。

### 2026-05-02

- 検索ダイアログの CLI Alt/Az 保持チェック
  - GUI の `対象検索...` と place 検索ダイアログに `Keep CLI-specified Alt/Az` を追加した。
  - `-A` / `-Z` の指定がある場合だけ有効化し、初期状態はチェック済みとした。
  - 指定された軸だけを CLI 視線として保持するようにした。
  - `preserve_cli_view_center` を `SearchJumpTarget` に載せ、選択結果から `_jump_to_search_target()` へ保持意図を渡すようにした。
  - place 検索と対象検索の両方で、保持フラグを無効にしたときは次回検索で CLI 視線へ戻さないことを明示した。

### 2026-04-25

- 恒星・惑星・月の輪郭のみ描画オプション
  - `vmag < 2.0` の恒星はダイヤ型の輪郭のみ、太陽・月・惑星は輪郭のみで描く表示モードを先に定義した。
  - `star_base_radius` を大きくしたときや高解像度表示時に、明るい天体の塗りつぶしが肥大化して目立ちすぎるため、占有面積を抑えるのが目的。
  - 実装時は星レイヤと solar-system レイヤの両方へ同じフラグを通す前提で進める。
- 恒星・惑星・月の輪郭のみ描画オプションの実装
  - `--outline-bright-bodies` を追加し、`RenderStyle` と `SkyWindowUserOptions` を通して共有描画へ配線した。
  - 恒星レイヤでは `vmag < 2.0` の星を塗りつぶしから除外し、ダイヤ輪郭だけを別描画するようにした。
  - 輪郭が低解像度の内部星面で潰れないよう、`outline_bright_bodies` 有効時は星レイヤの縮小レンダリングを止めてフル解像度で描画するようにした。
  - フル解像度化で星が小さく見えすぎないよう、輪郭モードでも通常の星レイヤ拡大率をサイズ計算へ戻した。
  - 恒星の輪郭はさらに細く見えるよう、輪郭ダイヤの表示スケールと強度を通常の星強調より下げた。
  - 惑星・月レイヤでは、太陽をクロスのみ、惑星を輪郭のみへ切り替え、月は通常表示では輪郭のみとしつつ、`enlarge_moon` と hover による拡大表示は通常の月レンダリングを使うようにした。
  - 回帰テストとして、星の center pixel が空くこと、惑星・月の描画経路が outline-only へ切り替わること、window render の配線が通ることを確認した。

- version bump
  - `__version__` を `1.21.7` に上げた。

### 2026-05-01

- 都市アウトラインの穴リング簡略化計画
  - 遠距離で建物の見かけサイズが十分小さい場合、外周リングだけを残して穴リングを省略する方針にした。
  - 判定は描画段階ではなく `compute_urban_outlines()` 側へ寄せ、不要な ring の投影・サンプリング・分割を先に避ける。
  - まずは外周リングの投影 bbox と見かけ方位幅を使って coarse 判定し、近距離では従来どおり穴リングを維持する。
  - 実装順は、判定ヘルパー追加 -> 回帰テスト追加 -> urban outline 生成経路へ接続 -> しきい値調整、の順で進める。

- 都市アウトラインの穴リング簡略化実装
  - `compute_urban_outlines()` に外周リング優先の ring 処理を入れ、穴リングは距離と投影サイズが十分小さいときだけ省略するようにした。
  - 省略判定は `compute_urban_outlines()` 内で行い、描画層に不要な ring を渡す前に止めるようにした。
  - `distance_km` の保持は建物ごとの最短距離に戻し、表示メタデータの意味も保った。
  - 回帰テストとして、遠距離で小さい穴リングが落ちることと、近距離では穴リングが残ることを追加した。

- 次期都市アウトライン最適化案
  - 画面上で縦方向に十分薄い輪郭を、生成段階で 2 点の線分へ圧縮する案を検討する。
  - この圧縮は描画段階ではなく `compute_urban_outlines()` 側で行い、`draw_urban_outlines()` を単純な表示器として保つ。
  - 目標は、遠距離の密集地で polyline 点数と Qt 描画コストをさらに減らしつつ、見た目の変化を最小化すること。

- 次期都市アウトライン最適化の実装
  - `compute_urban_outlines()` で、視点中心に対する正規化縦幅が十分小さい run を 2 点の線分へ圧縮するようにした。
  - ただし中距離の小さな建物まで折りたたまないよう、距離下限を追加して遠距離 run に限定した。
  - 線分化の端点は、run の normalized-screen vertical extrema を使うようにして、塔のような細い輪郭が同一点へ潰れないようにした。
  - その後、線分化のしきい値を `0.02` から `0.01` に下げて、台形に見える輪郭が過剰に折りたたまれないようにした。
  - 線分化の判定は描画層に入れず、生成段階で完了させることで Qt 側の polyline 構築を減らした。
  - 回帰テストとして、縦方向に薄い run が 2 点へ落ちることと、縦幅の大きい run が従来どおり残ること、さらに thin run の端点が別座標として保持されることを追加した。

### 2026-04-14

- Himawari partial-slot acceptance policy
  - Himawari ISatSS の slot 採用を、全 `88` タイルの完全一致ではなく、observer 描画と warm-threshold 推定に必要なタイルの充足基準へ寄せる方針に整理した。
  - 既存実装は slot 全体の完全性を先に見ていたが、必要タイルがそろっている場合まで捨ててしまうので、仕様・設計を先に合わせることにした。
- Himawari far-miss clear-sky policy
  - 観測者向け描画については、観測地点からおおむね `50 km` より遠い欠損タイルを `clear sky` として扱う方針に整理した。
  - 赤道帯の warm-threshold 推定用タイルは別扱いとし、欠損時は前回の有効 warm-threshold を再利用するフォールバックを優先する。

### 2026-04-24

- Moon phase screen rotation
  - 天頂を見上げた視点で観測方位を回転したとき、月相画像が空の回転と逆向きに回る問題を修正した。
  - `calculate_moon_render_data()` の画面回転角は、投影上の月ローカル上方向に同じ向きで追随する符号に戻した。
  - 天頂視点で `view_center` の方位を変えたとき、月の画面回転差も同じ向きになることを単体テストで固定した。
- Direction grid cleanup
  - GUI の hover grid と CLI export の static grid が共有する `_draw_direction_grid()` 内で、同一スタイルの緯線描画ループを一本化した。
  - GUI/CLI の入口は維持し、描画結果を変えずに内部重複だけを減らした。

- version bump
  - `__version__` を `1.20.17` に上げた。

### 2026-05-26

- Equidistant Conic validation asset
  - `dev-samples/fit_equidistant_conic_image_mapping.py` を追加し、`latlonmap.txt` の RGB 制御点から Equidistant Conic 仮説を当てはめて `eqdc_lonlat.npz` を生成できるようにした。
  - `dev-samples/draw_equidistant_conic_latlon_grid.py` を追加し、生成済みの `eqdc_lonlat.npz` を使って同じ 10 度グリッドを画像へ重ねられるようにした。
  - 検証対象は MET Norway Geo-Satellite API の `https://api.met.no/weatherapi/geosatellite/1.4/?area=europe&type=infrared` から取得した Europe infrared 画像とした。
  - `raw-data/geosatellite/README.md` にも同じ取得元 URL を追記し、検証アセットの出所を明示した。

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

- 山頂ビューポイント dataset と解決経路の実装
  - `mountain_viewpoints.json` を同梱 dataset として扱い、山名から山頂ビュー用の代表点へ解決できるようにした。
  - mountain 側の dataset 生成は、Wikidata raw query をそのまま使わず、まず review 用候補 JSON に畳んでから curated seed と最終 viewpoint JSON へ進める方針にした。
  - mountain viewpoint dataset は CLI 参照対象とし、一覧表示、名前一覧、JSON 出力を tower 系と同じ流儀で扱えるようにした。
  - 通常起動の `location` 引数でも mountain viewpoint 名と Wikidata ID を受け付けるようにし、最近傍都市からタイムゾーンを補完する構成にした。
  - 都市名との衝突を避けるため、`m/NAME` による明示プレフィックス解決を追加した。
  - `tests/test_wikidata_mountain_candidates.py`、`tests/test_wikidata_mountain_viewpoints.py`、`tests/test_cli_mountain_options.py`、`tests/test_mountain_viewpoints.py`、startup resolution 系テストで回帰を確認した。

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

### 2026-03-27

- スカイスクレーパー遠距離レイヤーの打ち切り距離見直し
  - 同梱の tower / mountain viewpoint `227` 地点に対して、現行 `select_skyscraper_seed_tiles_for_viewer(...)` で外側半径 `10km` と `60km` を比較した。
  - 増分は平均 `+0.83` tile、中央値 `0`、90 パーセンタイル `0` で、`+20` tile 以上増える地点は `4` 件、最大増分は Macau Tower の `+32` tile だった。
  - 実測上は増分が一部の沿岸都市・高層集積地に偏っており、全体としては急増しなかったため、遠距離スカイスクレーパー補助レイヤーの既定上限を `10km` から `60km` へ引き上げた。
  - 続いて `--urban-outline-skyscraper-radius-km` を追加し、遠距離スカイスクレーパー補助レイヤーの外側半径を CLI から上書きできるようにした。
  - 新オプションは `--urban-outline-radius-km` 以上を必須とし、GUI 本体と `zstarview-export-image` の両方で同じ検証と runtime 伝搬を行うようにした。

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
  - `--content-fov-deg` の既定値を `115°` とし、ウィンドウ端 `90°` 固定のまま全レイヤーの描画対象を広げられるようにした。
  - 背景グラデーション、sky disc、雲合成、雲ハッチ、地面ティント、ガイド線、都市アウトライン、航空機、明るい星の簡易表示経路も同じ overscan 視野へ揃えた。
  - 暗い星だけを旧来の `90°` 境界で追加的に落とす最適化は削除した。
  - `content_fov_deg` の外側には、ウィンドウ背景へ自然につなぐ短い背景フェード帯を持たせる方針にした。
  - その後、低仰角の横長ウィンドウで sky disc と雲レイヤーが `radius*2` の固定正方形に見えていた問題を修正し、実ウィンドウ geometry に沿って overscan を広げるようにした。

- version bump
  - `__version__` を `1.20.8` に上げた。
  - 前回の hover 月描画の `marker_scale` 伝搬修正を含む一連の変更を、区切りのよいパッチとして記録するための更新。

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

- 2026-04-07
  - GUI の描画対象をホストウィンドウ本体からクライアント領域ウィジェットへ分離し、同じ描画ロジックを frameless host と通常 decorated host の両方で再利用できるようにした。
  - 起動オプション `--window-frame {frameless,window}` を追加し、既定は従来互換の `frameless` とした。
  - export-image parser では `--window-frame` を GUI 専用オプションとして引き続き非対応にした。

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

### 2026-03-22

- バージョン更新
  - 時刻シフト時の補助レイヤー制御と人工衛星 cache 見直しに合わせて、バージョンを `1.6.1` から `1.6.2` へ更新した。

- 時刻シフト時の補助レイヤー方針を整理
  - 仕様として、雲と航空機は現在時刻近傍のみ、人工衛星は現在と過去のみ、未来時刻では 3 つとも非表示とする方針を文書化した。
  - 地形地平線と都市アウトラインは地点依存レイヤーとして扱い、時刻シフトだけでは無効化しない整理にした。
  - `specification.md` と `design.md` に、過去・現在・未来でのレイヤー可否を明記した。

- 人工衛星キャッシュ方針の見直し
  - 人工衛星の current cache の fresh 判定を `24時間` から `6時間` へ見直す方針を文書化した。
  - 過去表示向けに `archive` 層を持ち、`3日` 保持、対象時刻との差 `6時間` 以内だけ採用する設計へ更新した。
  - current と past の有効判定を同じ `6時間` にそろえ、説明と実装を単純化する方針にした。
  - current cache 更新時は古い snapshot を archive へ移し、archive 名と探索キーには `element_epoch_utc`、保持期間の cleanup には `fetched_at_utc` を使う整理にした。

- 時刻モード別レイヤー可否と人工衛星 archive 実装
  - `overlay_time.py` を追加し、`past / present / future` 判定と雲・航空機・人工衛星の可否判定を GUI / export-image で共有するようにした。
  - GUI では過去表示と未来表示の両方で、雲・航空機・人工衛星を無効化するようにした。
  - `satellites/cache.py` を `current + archive` 方式へ更新し、current refresh 時に古い snapshot を archive へ移し、archive cleanup を `3日` 保持へ変更した。
  - 人工衛星 cache は `element_epoch_utc` と `fetched_at_utc` を分け、利用判定は前者、retention/cleanup は後者を使うようにした。
  - `zstarview-export-image` も同じ時刻モード判定を使うようにし、過去表示と未来表示の両方で 3 補助レイヤーを skip するようにした。

### 2026-03-23

- バージョン更新
  - ISS 専用化と `wheretheiss.at` 優先取得への切り替えに合わせて、バージョンを `1.6.2` から `1.6.3` へ更新した。

- ISS API source policy update
  - Artificial satellite support was narrowed to `ISS` only, with `wheretheiss.at` as the primary TLE source and CelesTrak `stations` as the fallback source.
  - This keeps the runtime dependency focused on a single visible target while preserving a secondary source when the primary API is unavailable.
  - The specification, design docs, and README files were updated to match the implemented behavior.

- タイムシフト時の人工衛星レイヤー無効化
  - 人工衛星レイヤーは realtime view 専用とし、`--hours`、`--days`、`--datetime` による time-shifted view では取得も表示も行わない方針へ変更した。
  - 人工衛星 cache は current 1 層へ戻し、archive snapshot の保存・探索・cleanup を撤去した。
  - これに合わせて README、仕様書、設計書、衛星 cache/controller テストを更新した。

- 人工衛星取得失敗時の backoff 見直し
  - CelesTrak 取得 timeout の既定値を `20秒` から `60秒` へ延長した。
  - GUI では人工衛星取得が失敗した場合、即時再試行ではなく `2時間` 後に再試行するように変更した。
  - これにより、短時間に連続 timeout が起きた場合でも CelesTrak への再試行頻度を抑える方針にした。

### 2026-03-25

- バージョン更新
  - time-shift 時の realtime overlay 可否を GUI / export-image / tests でそろえた変更に合わせて、バージョンを `1.6.4` から `1.6.5` へ更新した。

- 人工衛星レイヤーの対象と地平線下表示の整理
  - 既定の対象は `ISS` のみとした。
  - `ISS` は検索ジャンプ対象としても扱い、地平線下であっても cache された軌道要素から位置解決できるようにした。
  - `satellites/project.py` と renderer の両方で地平線下の marker を保持・描画できるようにして、検索ジャンプ結果と通常 overlay の不一致を解消した。

- 人工衛星 stale cache の表示継続と失敗 backoff の永続化
  - stale な `ISS` cache が残っている場合、取得失敗時でも描画側は cache fallback で表示継続できるようにした。
  - 軌道要素 cache JSON に `last_fetch_attempt_utc`、`last_fetch_failed`、`last_fetch_error`、`last_fetch_failure_utc`、`failure_backoff_until_utc` を保存するようにした。
  - これにより、取得失敗後の `2時間` backoff はアプリ再起動後も継続し、短時間の再起動で CelesTrak へ再度アクセスしないようにした。

### 2026-03-24

- 航空機スクリーンショット取得用のデバッグフック
  - `ZSTARVIEW_DEBUG_SAVE_AIRCRAFT_READY_FRAME` を追加した。
  - 航空機更新直後の描画フレームを PNG として保存する、スクリーンショット取得用のデバッグオプションである。
  - 通常運用では未設定のままにし、README などの利用者向け文書には載せない。

### 2026-03-26

- DEM と都市アウトラインの長寿命キャッシュ管理
  - `DEM` は `90日`、都市アウトラインは `30日` の保存期限を持つ方針にした。
  - 保存期限は削除期限ではなく再取得を試みる目安とし、期限切れでも既存 cache はオフライン利用を継続できるようにした。
  - `DEM` は tile ごとの sidecar metadata、都市アウトラインは dataset directory ごとの metadata で `fetched_at_utc` を管理する構成にした。
  - 既存 cache に `fetched_at_utc` が無い場合は、初回読込時の現在時刻を暫定値として書き戻す移行方針にした。
  - stale 判定で再取得に入ったことがターミナルから分かるよう、`DEM` と都市アウトラインの refresh / stale fallback / legacy metadata 移行ログを追加した。

- 長寿命キャッシュの手動クリア支援
  - `zstarview` と `zstarview-export-image` に `--clear-long-lived-cache` を追加し、`copernicus-dem`、`overture_buildings`、`overture_skyscrapers` を起動前に削除できるようにした。
  - 誤用防止のため、最後の実行日時を cache root に記録し、`3日` 以内の再実行は拒否するクールダウンを入れた。
  - GUI では拒否理由を splash に表示するようにした。
  - 手動削除先を OS 非依存で確認できるよう、`zstarview-export-image --print-cache-dir` を追加した。

- mountain viewpoint dataset の補強
  - Wikidata 由来の座標と標高を使って、著名な山を 12 件追加した。
  - 既存の source policy を崩さず、同じ dataset 形式へ正規化できる山だけを追加する方針にした。
  - mountain resolver と startup location resolution の両方で追加項目が解決できることをテストで確認した。

### 2026-03-27

- 地点名解決まわりのディレクトリ分離
  - 名前やクエリ文字列から緯度経度または viewpoint を解決する機能を、`location_resolver/` サブディレクトリへ段階的に移した。
  - 第1段階として `nominatim_search.py` の本体を `location_resolver/nominatim.py` へ移した。
  - 第2段階として `viewpoints.py`、`tower_viewpoints.py`、`mountain_viewpoints.py` の本体を `location_resolver/` 配下へ寄せた。
  - 第3段階として `launch_location_time.py` から地点解決ロジックを分離し、最終的に時刻解釈は `launch_time.py`、地点解決は `location_resolver/resolve.py` へ整理した。
  - 段階移行の完了に合わせて旧ラッパーは解消し、参照先を新しいモジュール構成へ切り替えた。

- テスト構成の整理
  - `aircraft`、`clouddisc`、`satellites`、`terrain` のテストをサブディレクトリへまとめた。
  - 地点解決に直接属するテストを `tests/location_resolver/` に移した。
  - Wikidata、catalog、derived tile、建物 import などのデータ寄りテストを `tests/data/` に移した。

### 2026-03-30

- GUI 全体 opacity の試行と撤回
  - Wayland 環境で GUI ウィンドウ全体へ opacity を掛ける手段として、Qt の top-level window opacity を使う案を試した。
  - まず `WA_TranslucentBackground` を維持したまま whole-window opacity を検討し、その後、`WA_TranslucentBackground` を外して `setWindowOpacity()` を使う変種も試した。
  - いずれも対象の Wayland 環境では見た目上の効果が確認できず、compositor / platform 実装依存の制約が強いと判断した。
  - render 結果に後段 alpha を掛ける代替案も検討したが、popup や child widget との一貫性、既存の透過前提 UI との整合を考えると、その場での導入は見送った。
  - 最終的に、この実験に関する未コミット実装は `git reset --hard` で破棄し、仕様書・設計書に一時的に追加した `--window-opacity` 記述も取り下げた。

### 2026-04-01

- 直接座標入力のコンマ区切り緩和
  - CLI の直接座標入力で `lat, lon` と `@lat, lon` を受け付けるようにした。
  - Google Maps の右クリックメニューからコピーした `lat, lon` 文字列を、そのまま起動引数に貼り付けても解決できるようにした。

### 2026-04-02

- 都市アウトライン高さ基準の設計更新
  - Overture 建物高は DEM 由来の絶対標高ではなく、地表基準の建物高として扱う方針を明文化した。
  - 都市アウトラインの見かけ仰角は、建物 footprint 側 DEM 地盤標高と観測地点 DEM 地盤標高を加えた絶対標高差から計算する設計へ更新した。
  - 現行実装は未追随であり、設計先行の変更として `docs/design.md` に反映した。

### 2026-04-06

- 地平線下 earth guide の初期導入
  - まずは線のみで地平線下の大陸輪郭を描く earth guide を追加し、地形地平線とは別の補助レイヤーとして使い始めた。
  - 当初は terrain horizon と色・opacity を共有する単純な構成にして、見た目の確認を優先した。
  - 後続の調整で、観測者近傍の除外や GUI / CLI の独立制御へ発展していく前段階の実装である。

- ### 2026-04-29

- 地形地平線の距離帯ごとの副稜線と遮蔽ルール
  - 地形地平線に、`500m`、`1km` を含む距離帯ごとの副稜線を補助レイヤーとして追加し、山並みの重なりや視界の遷移が読み取れるようにした。
  - 同一方位では、より近い距離帯ごとの副稜線がより高く見えている場合に、より遠くて低い距離帯ごとの副稜線を表示しないルールを採用した。
  - これは、前景の山塊が背後の稜線を物理的に隠すという見え方を優先し、納得度の低い線の重なりを抑えるためである。
  - 距離帯ごとの副稜線は主シルエットの可読性を損なわない範囲で、低 opacity の補助情報として扱う。

### 2026-04-11

- 地形地平線と earth-guide のスタイル整理
  - 地形地平線の前景線 alpha 計算を共通関数へ切り出し、earth-guide も同じ alpha カーブを使うようにした。
  - `--terrain-horizon-opacity` の既定値を `0.028` に合わせ、renderer 側の `0.7` 係数は削除した。
  - earth-guide は terrain horizon と同じ RGB を使いつつ、前景線のみの単一ストロークへ整理した。
  - terrain horizon の見た目を earth-guide より少し強くするため、terrain 側の線幅を 2 割ほど増やした。
  - earth-guide は引き続きスクリーン空間で再帰分割し、遠景は粗く・近景は細かく見えるようにしている。

- 地平線下地球ガイドの近傍除外設計
  - 観測者足元の陸地は、自己帰属の難しさと閉ポリゴン化の破綻を避けるため、地球裏面ガイドから省略してよい方針を明文化した。
  - 近傍除外は固定角ではなく、`observer_height_m` から求める地平線距離に比例して広がる地表距離ベースの dead zone として扱う方針にした。
  - 追加の地平線近傍クリップは、地平線 dip と小さな余裕角に基づく二段目の絞り込みとして設計し、`docs/specification.md` と `docs/design.md` に反映した。

### 2026-04-12

- Earth guide の terrain horizon からの分離と README 反映
  - `Earth Guide` を terrain horizon とは独立したレイヤーとして扱い、`--earth-guide-opacity` と GUI トグル `E` / `Earth Guide` メニューを追加した。
  - `terrain_horizon_opacity` は DEM と地形地平線にのみ関与させ、Earth guide は別の opacity と cache key で制御するようにした。
  - Earth guide の色は terrain horizon と同じ系統を保つため、別名定数 `EARTH_GUIDE_LINE_COLOR` を導入した。
  - 利用者向けには README 英語版・日本語版へ反映し、`--earth-guide-opacity 0` で Earth guide のみを無効化できることを明記した。

- `zstarview-export-image` の SIXEL 端末可否チェック強化
  - `--sixel` 指定時に `img2sixel` の存在確認に加えて、`lsix` と同様の device attributes 問い合わせ (`ESC[c`) を行うようにした。
  - 応答に SIXEL 示唆の `4` が含まれない端末では、重い初期化やレイヤー取得へ進む前に明示エラーで終了するようにした。
  - 仕様文書と内部設計に端末可否チェックの前提を追記し、export-image の単体テストも追加した。

- `auto` による現在地取得の導入
  - 位置引数を `auto` にしたとき、`ip-api.com` の API を使って IP アドレスから現在地を取得する機能を追加した。
  - 取得結果は緯度、経度、タイムゾーン、地名を含む `ResolvedLocation` として扱い、他の地点指定と同じように起動時の観測地点へ反映するようにした。
  - その後、`ip-api.com` への自動問い合わせには 3 秒のレート制限を入れて、連続起動時の過剰アクセスを避けるようにした。

- `--observation-info` の起動モード拡張
  - 観測情報オーバーレイの起動モードを `auto / top / bottom / off` で明示的に指定できる `--observation-info` を追加した。
  - GUI では `top / bottom` を pinned mode として扱い、マウス位置による HUD の自動移動を抑止するようにした。
  - `zstarview-export-image` ではこの新しい起動モードを露出せず、従来どおり `--show-guidelines-initial` だけを headless 向けに残した。

- 赤道・黄道の点線間隔調整
  - 天の赤道と黄道の dash pattern を少し粗くして、通常ズームで点線がくっついて見えにくくなるようにした。
  - 地平線の点線は従来どおり維持し、赤道・黄道だけを対象にした。

- バージョン更新
  - パッチレベルの更新として `__version__` を `1.9.7` に上げた。

- バージョン更新
  - パッチレベルの更新として `__version__` を `1.9.8` に上げた。

- バージョン更新
  - パッチレベルの更新として `__version__` を `1.9.9` に上げた。

- バージョン更新
  - パッチレベルの更新として `__version__` を `1.9.10` に上げた。

### 2026-04-13

- 雲ストライプ `width` モードの減衰カーブ調整
  - `width` モードの alpha 減衰を線形から ease-out 寄りの二乗カーブへ変更し、基準線付近の明るさを保ったまま遠側だけをゆるやかに薄くするようにした。
  - 実装に合わせて `docs/specification.md` と `docs/design.md` へ、`width` モードの alpha 減衰は線形に限らずゆるやかなカーブとしてよい、という説明を追記した。

### 2026-04-17

- 人工衛星レイヤーの Horizons spacecraft 拡張
  - ISS の既存 TLE 経路は維持しつつ、JPL Horizons の observer ephemeris を使って `JWST`, `Voyager 1`, `Voyager 2`, `Parker` を追加した。
  - Horizons 側は observer site 別に cache を分け、CSV 出力から直接 `alt/az` を読み取って既存の overlay model に流し込むようにした。
  - ISS と Horizons spacecraft の cache TTL はともに `24h` に揃えた。
  - `README.md`、`docs/specification.md`、`docs/design.md` を新しい表示対象とデータ経路に合わせて更新した。

### 2026-04-18

- 対象検索拡張の設計メモ
  - JPL 小天体検索を既存の恒星・アステリズム検索 UI に統合し、選択対象へ視界中心を合わせる方針を確認した。
  - 結果の永続表示は、検索ダイアログ下部のチェックボックスで制御する方針にした。
  - 永続マーカーとラベルは、一時 jump highlight と別の状態として扱う前提で整理した。

- JPL small-body 本体接続
  - `対象検索` ダイアログを local first 方式にして、恒星・アステリズムの検索結果が無い場合だけ JPL small body へフォールバックするようにした。
  - 検索欄の直下に JPL database 検索ボタンを追加し、Enter キー依存ではなく明示的に JPL へ問い合わせられるようにした。
  - JPL フォールバックは major body と small body の両方を問い合わせるようにした。
  - `Keep marker` は JPL 結果にも反映し、major body の場合も持続表示だけは残るようにした。
  - `Sun` と `Moon` は JPL フォールバックから除外し、solar-system 側の表示に委ねるようにした。
  - `SearchJumpTarget` に JPL 用の `command` / `alt_deg` / `az_deg` / 永続フラグを追加し、`SkyWindowState` に永続ターゲットを保持するようにした。
  - 永続マーカーは衛星クロスと同じ大きさで描画し、ラベルは単独のアウトラインテキストとして重ねるようにした。
  - JPL 検索は Horizons lookup と observer ephemeris を使い、現在時刻の `alt/az` を検索結果として返す実装にした。

- JPL small-body 永続更新の retry policy
  - 永続表示対象の JPL small body について、検索時刻を起点に 1h ごとの再問い合わせを行う one-shot timer を追加した。
  - 再問い合わせに失敗した場合は即時再試行せず、失敗時刻から 1h 後へ次回 retry を延期するようにした。
  - retry の成功/失敗は `SkyWindowState` の `persistent_search_last_error` / `persistent_search_next_refresh_utc` で追跡し、HUD の status line に反映するようにした。

- 対象検索 / export-image overlay follow-up
  - `zstarview-export-image` でも検索結果の marker / label を出力画像へ重ねるようにし、検索で得た alt/az は GUI と同じ `-45°` 下限へ揃えた。
  - 検索解決の共有設計は `search/` サブパッケージへ移したが、UI との切り分けや overlay 連携はまだ仮設計なので、`docs/design.md` に provisional note を追加した。

- Validation note

### 2026-04-22

- Topic: FOV default alignment
  - Decision: Change the default projected edge FOV to `95°` and the default shared content FOV to `110°`.
  - Rationale: The defaults used by the CLI, runtime viewer data, and GUI input layer should match so the app starts with a consistent projection scale and overscan budget.

- Voyager 1 検索経路の整理
  - `対象検索...` と CLI 検索で使う `satellite` 経路を `ISS` のみに限定し、`JWST`, `Voyager 1`, `Voyager 2`, `Parker` は JPL 検索結果として扱うようにした。
  - ローカル検索候補の組み立てから衛星ショートカットを外し、`kind="satellite"` の曖昧な検索結果が Voyager 系へ混入しないようにした。
  - これにより、Voyager 1 の検索結果は `JPL [Voyager 1 (spacecraft)]` として扱われ、Horizons observer ephemeris の `alt/az` と比較しやすくした。

- JPL / Horizons デバッグ情報の強化
  - JPL 解決ログに `target_time_utc` と observer の緯度・経度・標高を含め、手動ダウンロードした Horizons CSV と GUI の比較条件を明示できるようにした。
  - 永続検索ターゲット更新ログには `kind`, `group`, `command`, `alt`, `az` を含め、`satellite` 経路と `jpl_body` 経路の取り違えを追跡しやすくした。
  - JPL small-body status line では、エラー文字列が無い場合に `None` を表示せず、`retry pending` などの状態だけを見せるようにした。
  - Modified GUI, renderer, and test files passed `py_compile` and `ruff check`.
  - Full `pytest` execution was not possible in the active `.venv` because `pytest` is not installed there.

### 2026-04-21

- Cloud stripe pre-composite debug dump
  - Added an optional debug dump for the pre-composite cloud stripe image so the width/alpha stripe result can be inspected before sky composition.
  - The dump is gated by `ZSTARVIEW_DEBUG_SAVE_CLOUD_STRIPE_FRAME`; the default destination is `CACHE_PATH/debug/cloud-stripe`, and an explicit path can be supplied instead.
  - Added regression coverage to ensure the snapshot is written only when the debug flag is enabled.
- 2026-04-22
  - Patch-level release bump to `1.20.15`.

- Topic: enlarged moon in outline mode
  - Decision: let `enlarge_moon` and hovered moon rendering use the normal filled moon artwork even when `outline_bright_bodies` is enabled.
  - Rationale: the outline-only mode should reduce footprint for ordinary bright bodies, but the enlarged moon view is explicitly meant to show the rendered moon, not a hollow marker.
  - Result: moon enlargement now bypasses the outline-only fallback while the rest of the outline-only behavior remains unchanged.

### 2026-04-25

- Bright-bodies mode rename plan
  - `--outline-bright-bodies` is being replaced by `--bright-bodies {outline,fill}` so the default state can be `outline` without encoding the mode as a boolean flag.
  - The new form keeps the intent explicit, removes the legacy flag entirely, and leaves room for future non-boolean body rendering modes.
- Star size threshold refinement
  - The large-star outline treatment was narrowed to stars whose rendered size is `>= 7px`, and the outline stroke for that path is `2px`.
  - Stars below that threshold remain filled rectangles, which keeps the outlined-star footprint restricted to the largest markers.
  - The implementation now batches full-fit outline rectangles by size to reduce Python-loop overhead, with a scalar fallback only for clipped edge cases.
- Documentation sync
  - Updated `docs/specification.md` and `docs/design.md` to describe the current size-based star rendering behavior and the `bright_bodies_mode` outline interaction.

### 2026-04-27

- Celestial reference line smoothing
  - Moved celestial equator and ecliptic refinement into the renderer and subdivided by screen-space error instead of relying on fixed angular spacing.
  - Kept the astro layer at coarse control points and added regression coverage for both the control-point contract and adaptive drawing behavior.
  - This matches the Earth-guide drawing strategy more closely and prevents visibly angular lines at fullscreen sizes.
- Arrow-key viewport interaction handoff
  - Kept viewport-interaction mode active while arrow keys are physically held, then ended the simplified mode on key release instead of the previous idle-timer handoff.
  - Deferred the normal cloud and terrain restart until the refreshed sky disc arrives, so release-time repaints do not show an intermediate black partial sphere.

### 2026-04-29

- Distance-band terrain ridge preview
- Replaced the previous chained secondary-ridge experiment with concentric distance-band secondary ridges at 3.75, 7.5, 15, 30, 60, and 120 km.
  - Kept the primary terrain horizon only for fast-mode interaction rendering and rendered the banded overlay as the normal-mode terrain visualization.
  - Aligned the nearest normal-mode band with the fast-mode horizon width and moved its alpha closer to the fast-mode look, then tapered farther bands down from there.
  - Updated the horizon extraction, renderer, state plumbing, and regression tests to carry banded ridge polylines instead of visibility-pruned secondary peaks.

### 2026-05-03

- Object-colored hover labels
  - Kept star and planet labels on the theme text color, but changed DSO and asterism hover labels to use the corresponding object colors.
  - This makes the auxiliary labels visually match the item being described without changing the conventional white labels for stars and planets.
  - Validation: `pytest tests/test_asterism_render.py tests/test_planet_marker_style.py tests/test_deep_sky_objects_render.py tests/test_window_render_sync.py -q`

- Version bump
  - Bumped `__version__` from `1.21.14` to `1.21.15`.
  - This follows the sky-opacity default restoration and keeps the patch-level version in sync with the current code state.

- Version bump
  - Bumped `__version__` from `1.21.16` to `1.21.17`.
  - This records the final location-summary layout update and keeps the patch-level version aligned with the source tree.

### 2026-05-09

- Night light band tuning and seam handling
  - Split the night-light glow into four draw layers, keeping the core band and dividing the mid-range glow into `MID_NEAR` and `MID_FAR`.
  - Moved the night-light seam to the azimuth behind the observer instead of hard-coding a seam at `az=0°`.
  - Kept the underlying distance-band integration and cached base profile behavior unchanged.

### 2026-05-11

- Cloud opacity compositing clarification
  - Adjusted the cloud composite path so `--cloud-opacity` also scales the gray underlay mixed behind cloud strokes, not only the additive cloud layer and final alpha.
  - Updated the specification and design docs to describe `--cloud-opacity` as the final cloud-compositing strength, including the gray desaturation blend.

### 2026-05-15

- Cloud download cancellation on shutdown
  - Added a cooperative abort event to cloud fetches so GUI shutdown can cancel in-flight GOES and Himawari downloads instead of waiting for every S3 transfer to finish.
  - Wired the abort event through the cloud controller, cloud-disc fetch path, and S3 download helper; progress callbacks now raise `KeyboardInterrupt` when shutdown is requested so boto3 can cancel the managed transfer.
  - Added regression coverage for canceled S3 downloads, Himawari fetch aborts, and cloud-controller shutdown signaling.

- Terrain and urban-outline cancellation on shutdown
  - Added the same cooperative abort event to Copernicus DEM fetches and Overture building imports so terrain-horizon and urban-outline shutdown can stop multi-file downloads early.
  - Copernicus DEM downloads now raise a cancellation error from the S3 progress callback, and Overture downloads terminate the `overturemaps` subprocess when shutdown is requested.
  - Added regression coverage for terrain-controller and urban-outline-controller shutdown cancellation plus direct function-level cancel tests.

- Water overlay startup responsiveness
  - Added cooperative cancellation to water overlay fetch and sampling paths so the controller can stop cleanly when shutdown is requested.
  - Made the water overlay geometry loops yield periodically while they scan Overpass elements and project sample points, which keeps the GUI thread responsive during long startup updates.
  - Added regression coverage for water-overlay cancellation and controller shutdown behavior.

- Water overlay startup order
  - Deferred the first water-overlay update until the terrain horizon has either resolved or been skipped, so the water layer no longer starts while terrain is still pending.
  - When terrain is disabled, the initial water update still starts immediately; this preserves the previous behavior for that configuration.
  - Added regression coverage for the startup ordering in both the terrain-enabled and terrain-disabled cases.

- Version bump
  - Bumped `__version__` from `1.25.7` to `1.25.8`.
  - This keeps the patch-level version aligned with the current source tree after the water-mask documentation updates and sea-mask rendering adjustments.

## 4. 2026-05-16

- Water overlay sea/inland split and two-stage redraw
  - Split the water overlay runtime state into separate sea and inland layers, so the sea-mask result and OSM inland-water result can be tracked independently.
  - Made the controller emit a sea-only intermediate update first, then a combined sea + inland update after inland sampling completes.
  - Updated the specification and design docs to match the current behavior, including the separate inland-water cache and distinct sea/inland colors.

- Water data source split in docs
  - Clarified in the specification and design docs that sea-mask tiles come from OSM Water Polygons, while rivers, lakes, ponds, and similar inland water still come from Overpass-fetched OSM footprint data.
  - Rationale: the user-facing and internal docs should reflect that the sea and inland acquisition paths are separate, even though they are presented together as one optional water overlay.

- Water overlay CLI docs alignment
  - Updated the CLI overlay docs in both English and Japanese to explain that sea points come from OSM Water Polygons sea-mask tiles, while inland water points come from Overpass-fetched OpenStreetMap features.
  - Rationale: the user-facing overlay help should match the README and the main specification so the different acquisition paths are not conflated.

## 2026-05-17

- Shared GUI worker pool
  - Routed long-running GUI background work through a single shared `ThreadPoolExecutor(max_workers=1)` and replaced per-task `threading.Thread` launches with `submit_gui_work(...)`.
  - The migration covers sky, cloud, terrain, satellite, aircraft, water, urban-outline, startup bootstrap, and search dialog work so native-heavy tasks do not overlap.
  - Added a lazy re-creation path after shutdown so the pool can be torn down during window close and still be used again inside the same process if needed.

- Water overlay crash mitigation
  - Removed the GUI water-overlay path's Copernicus DEM refetch and kept the update path on the terrain controller's existing ground information instead.
  - Added progress logging around sea-mask sampling, Overpass footprint fetch, and inland-water sampling, plus `faulthandler` enablement at startup for native crash diagnostics.
  - Added regression coverage to ensure the GUI water overlay no longer re-enters the DEM-fetch path during normal updates.

- Water overlay cache simplification
  - Simplified water-overlay footprints before saving or reloading them from the per-scope cache so distant, dense polygon geometry does not balloon the cache file.
  - The cache now writes simplified snapshots to `*_simplified.json` and keeps the legacy `*.json` as an input-only migration source.
  - When both files exist, the loader compares mtimes and regenerates the simplified cache from the newer legacy file if needed, which keeps the simplified cache current after a user updates from an older version.
  - Added info-level start/finish logging around the simplification path so slow legacy-cache migrations can be observed in the terminal or log file.
  - Switched the simplification design toward distance-band grid snapping with `2^k m` fixed cells, a `0.5°` apparent-angle band threshold, and ring rejection when the snapped geometry collapses to a single point or near-zero area.
  - Added regression coverage for the simplified cache path so the controller keeps storing the reduced geometry instead of the raw Overpass polygon set.

- Validation
  - Confirmed the worker pool serializes tasks, the startup overlay tests still pass, the controller shutdown waits still behave, and the water overlay regression coverage passes after the GUI simplification.

## 2026-05-19

- Update scheduler design consolidation
  - Decision: Consolidate the GUI update model around `fast-mode`, `normal-mode`, and deadline-based periodic refreshes. Fast-mode keeps only the latest view-change replay, while periodic cloud/satellite/aircraft refreshes are scheduled from the previous completion time instead of the previous timer tick.
  - Rationale: This matches the desired user-visible behavior more closely than the previous timer-driven enqueue model. It prevents worker busy periods from accumulating stale ticks, keeps the splash phase deterministic, and makes the priority of menu, search, and view-change updates explicit.
  - Status: The design document was updated to reflect the new model. The first implementation pass now routes timer-driven refreshes through an idle-gated scheduler tick and adds regression coverage for one-task-per-idle-tick behavior.

## 2026-05-24

- Viewport interaction mode restores solar-system bodies
  - Updated the fast viewport-interaction rendering path so it still draws the Sun, Moon, and planets instead of dropping the entire solar-system layer.
  - Kept labels disabled in that path so the interaction frame stays lightweight while the user is panning, rotating, or resizing.
  - Added regression coverage around `_draw_viewport_interaction_layers()` and synced the README, specification, and design docs with the current behavior.

- version bump
  - `__version__` を `1.28.6` に上げた。

## 2026-05-27

- Fixed-view export-image culling方針
  - `zstarview-export-image` は単発の固定視線書き出しなので、urban outline と water は cache を GUI 互換の完全版のまま保持しつつ、cache 読込後の in-memory 計算では背面半球を省略してよい方針にした。
  - 省略判定は az の単純な裏表だけでなく `view_center` と表示 FOV を基準にし、境界をまたぐ ring や footprint は必要に応じて見える側だけを残す。
  - GUI は起動後に az が変化しうるため、同じ省略を前提にした cache にはしない。
  - この方針は、urban outline resolve が支配的になっている export-image の固定視線経路を軽くするための設計メモとして残す。
