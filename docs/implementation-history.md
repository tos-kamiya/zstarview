# zstarview 実装履歴

最終更新: 2026-03-18

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

- 都市アウトラインを独立レイヤーとして扱う
  - 現在は都市アウトラインを描画経路の一部として常時表示しているが、他のオーバーレイと同じく独立したレイヤーとして扱えるように整理する。
  - 最小実装では、内部状態上で都市アウトラインの on/off を独立させ、既定値は on のままにする。
  - まずは state と描画経路の分離を優先し、UI トグルや設定保存の扱いは後段で追加する。
  - これにより、地形地平線や他レイヤーとの責務境界を明確にし、将来の表示切替追加を安全に進められるようにする。
  - 現在は `Urban Outline` メニュー項目と `U` ショートカット、`--urban-outline-opacity` を持つ。
  - 細い輪郭は太い水平線に簡略化し、方向キー操作中の簡易描画モードでは描かない。

- OpenSky ベースの航空機オーバーレイを設計する
  - 目的は「正確な航空管制表示」ではなく、「観測地点の空に動きがあることを見せる補助レイヤー」を追加することとする。
  - データソースは OpenSky Network `states/all` を想定し、観測地点由来の `bbox` 付き取得を前提にする。
  - 無料枠前提のため、取得更新間隔は `5分` を既定とする。
  - 初期実装では API 応答の間の補間・外挿は行わず、`5分` ごとのスナップショット表示として扱う。
  - 各機体の状態は少なくとも `lat`、`lon`、`alt`、`velocity`、`heading`、`vertical_rate`、`last_update_time` を持つ方針とする。
  - `bbox` の初期値は観測地点中心 `±1.0°` とし、`4 square degrees` で `1 credit` 帯に収める方針とする。
  - OpenSky 側の現行認証条件、匿名利用可否、credit 消費条件は実装着手時に再確認する。
  - 最小段階では「取得」「地平線上かつ視野内の点表示」「取得失敗時の前回状態維持」までを範囲とし、ラベルや軌跡は後段とする。

## 4. TODO

- GUI 上で時刻を前後できるダイナミックなタイムシフト操作を追加する
- CLI 仕様と内部データ構造の対応表を設計書へ追加する
- 必要に応じて `CloudController` と `SkyWindow` の責務境界を再評価する

## 5. 実装履歴

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
  - この簡易モード中は惑星、全星等の星空、sky disc、雲、DSO、アステリズムを表示しない。
  - 簡易モード中の明るい星は新しい FOV に合わせて同期的に再計算する。
  - 補助線は worker 時点の正規化点列ではなく、`(alt, az)` サンプル列から描画時に `render_view_center` 基準で投影する。
  - 簡易モード終了時に初めて、全星等の星空、sky disc、雲、地形地平線の通常更新を開始する。

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
