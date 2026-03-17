# zstarview 実装履歴

最終更新: 2026-03-17

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

- INPROGRESS: スカイスクレーパー専用広域タイルの追加レイヤー案
  - 背景として、視点中心半径 `2.5km` を `7.5km` へ単純に広げると、Skytree 周辺の Overture `building` ダウンロード量は約 `6.5x` に増えた。
  - 既定の通常レイヤーは引き続き視点中心 `2.5km` の derived tile を使い、全地点での取得コストは増やさない方針を優先する。
  - 代替案として、`150m` 超級の building だけを収録した広域 Overture derived tile 群を別ディレクトリへ事前生成し、runtime で補助レイヤーとして合成する構成を検討する。
  - ここでいう `zoom 14` は Overture 固有タイルではなく、Web Mercator の `z/x/y` 分割を便宜的な scan 単位として使う定義とする。各 tile は必要に応じて `west,south,east,north` の bbox へ戻して `overturemaps download --bbox=...` に渡す。
  - runtime では、通常 `2.5km` レイヤーを先に読み、その後スカイスクレーパー専用ディレクトリから視点周辺に関係するタイルだけを追加選択する。
  - 重複回避はまずタイル責務と距離リングで扱う。近距離 `0-2.5km` は通常レイヤー、遠距離 `2.5-10km` はスカイスクレーパーレイヤーを採用し、必要なら最終段で建物距離により内側成分を落とす。
  - 初期検討では `building_part` は広域高層専用レイヤーに含めず、`building` のみを対象としてよい。
  - Scraperbase `completed` かつ `150m` 以上の building を 30 都市 seed に突き合わせた試算では、scan 粒度を `zoom 12` から `zoom 14` へ上げるとユニーク tile 数は `199 -> 662` に増えた。一方で `zoom 12` は東京付近で 1 tile が現行 `2.5km` bbox よりかなり大きく、不要領域が広い。
  - そのため現時点の第一候補は `zoom 14` とし、scan tile 数の増加と引き換えに、広すぎる bbox を避ける方を優先する。
  - 外側半径は初期案の `7.5km` ではなく `10km` を第一候補とする。理由は、日本の代表例として東京タワーと東京スカイツリーの直線距離が約 `8.2km` あり、`7.5km` では主要ランドマーク間の遠景関係を取りこぼすためである。
  - 実装の第一段階では、skyscraper scan tile の seed リストをアプリ同梱データとして持つ。seed は少なくとも `z/x/y`、bbox、都市名、最大高さ帯、代表 building 名を保持し、Scraperbase から生成した JSON を人手レビュー後に取り込む前提とする。
  - ダウンロード方針は「初回起動で全件取得」ではなく「オンデマンド取得 + 永続キャッシュ」とする。seed リストだけは同梱し、実データは視点更新時に必要 tile が見つかったときだけ取得する。
  - キャッシュ配置は通常 urban outline 用 derived root とは別に、例えば `CACHE_PATH/overture_skyscrapers/z14/{z}_{x}_{y}/bldg/*.json` のような専用 root を持たせる。tile ごとに既存 derived tile 互換 JSON と `tile_index.json` を置ける形を優先する。
  - 1 tile の取得ジョブは `overturemaps download --bbox=... --type building` を使い、結果を import 時点で `height_m >= 150` に再フィルタして保存する。これにより runtime 側は高さ条件を持たず、通常 derived tile と同様に扱える。
  - runtime の tile 選択条件は、「tile bbox が視点中心 `10km` 円と交差し、かつ `2.5km` 内側だけには完全に収まらないこと」を一次条件とする。建物単位ではさらに `2.5km < distance <= 10km` を満たすものだけを採用する。
  - 都市 seed に含まれない場所では skyscraper tile の探索自体を行わないのではなく、視点中心 `10km` 円に交差する seed tile が存在しない限り何もしない構成とする。これにより世界全体の追加負荷を抑えつつ、都市名ベースの特別扱いを runtime へ持ち込まない。
  - 取得のトリガは `UrbanOutlineController` に寄せ、通常 `2.5km` レイヤー更新と同じ worker 経路で扱う。ただし skyscraper tile ダウンロードは低優先度の補助フェーズとして分離し、通常レイヤーの ready を待たずに UI を描画可能にする。
  - UI 状態としては、通常 urban outline と別トグルを増やさず、既存 `Urban Outline` の一部として扱う。必要なら status/banner で `Urban outline: downloading skyscraper cache...` のように補助取得を明示する。
  - 失敗時の扱いは safe-by-default とする。skyscraper tile 取得が失敗しても通常 `2.5km` urban outline は維持し、追加レイヤーだけを silently skip または banner 通知に留める。
  - データ更新は Overture upstream の変動を考慮し、通常 derived cache と同様に「手動削除すれば再取得される」運用を初期案とする。自動 refresh や世代管理は第一段階の対象外とする。
  - 最低限の検証観点として、(1) `10km` リングでの tile 選択が期待通りか、(2) `2.5km` 内側建物が混入しないか、(3) cache 未作成時でも通常レイヤー描画を阻害しないか、(4) skyscraper cache がない都市で余分なジョブを起こさないか、をテスト対象にする。
  - この案は未確定であり、実際に見え方差分が十分ある都市だけに適用する curated 運用を前提とする。

## 4. TODO

- GUI 上で時刻を前後できるダイナミックなタイムシフト操作を追加する
- CLI 仕様と内部データ構造の対応表を設計書へ追加する
- 必要に応じて `CloudController` と `SkyWindow` の責務境界を再評価する

## 5. 実装履歴

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
