# zstarview 設計書

最終更新: 2026-03-09

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
  - キャッシュ管理

## 4. モジュール構成

### 4.1 起動・設定

- `src/zstarview/zstarview.py`
  - アプリケーションの主エントリポイント
  - 起動シーケンスの組み立て
- `src/zstarview/cli_args.py`
  - CLI オプション定義と値解釈
  - タワー一覧・タワー詳細 JSON 出力の即時終了オプションを扱う
- `src/zstarview/startup.py`
  - 起動時の地点解決、設定復元、初期値決定
- `src/zstarview/config.py`
  - ユーザー設定の保存と読込
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
- `src/zstarview/types.py`
  - ドメインデータの共有型

### 4.2.1 ビューポイント dataset 参照専用 CLI

同梱 `tower_viewpoints.json` と `mountain_viewpoints.json` を直接参照する軽量 CLI 経路を持つ。

- 対象オプション
  - `--list-towers`
  - `--list-tower-names`
  - `--show-tower-json NAME`
  - `--list-mountains`
  - `--list-mountain-names`
  - `--show-mountain-json NAME`
- この経路では `startup.py` の都市解決や GUI 初期化へ進まず、`zstarview.py` が `tower_viewpoints.py` / `mountain_viewpoints.py` を直接呼び出して標準出力を書き、即時終了する。
- 一覧出力はローカル JSON のみを参照し、GeoNames、設定保存、ネットワークアクセスには触れない。
- `--show-tower-json NAME` の名前解決規則は、通常起動時のタワー解決と同じ `resolve_tower_viewpoint()` を再利用する。
- `--show-mountain-json NAME` は `resolve_mountain_viewpoint()` を再利用し、mountain viewpoint dataset を 1 件に解決して JSON 出力する。
- `tower_viewpoints.py` は dataset の `name` を保持したまま、必要に応じて ASCII フォールバック名 `ascii_name` を算出する。
- `mountain_viewpoints.py` も同様に dataset の `name` を保持したまま、必要に応じて ASCII フォールバック名 `ascii_name` を算出する。
- `--list-towers` は `ascii_name` がある場合それを表示名として優先し、通常のタワー名解決でも `ascii_name` を一致候補へ含める。
- `--list-mountains` も同様に `ascii_name` がある場合それを表示名として優先し、mountain viewpoint 名解決でも `ascii_name` を一致候補へ含める。

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

## 5. 主要データ構造

### 5.1 地点・視点関連

- `ViewerData`
  - 観測地点
  - タイムゾーン
  - 表示中心の方位・高度
  - 観測者高さ
  - 画面描画に必要な視点情報

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

## 6. 処理フロー

### 6.1 起動フロー

1. `zstarview.py` が CLI とログを初期化する。
2. 設定ファイルから前回の地点やウィンドウ状態を復元する。
3. 入力文字列を都市、タワー、座標のいずれかとして解決する。
4. 星カタログや補助データを読み込む。
5. `SkyWindow` を生成し、初回描画を行う。
6. 必要に応じて雲更新と地形地平線更新をバックグラウンドで開始する。

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
4. 雲画像と欠損ティントを合成する。
5. ラベル、オーバーレイ、ステータス行を描画する。

## 7. スレッドモデル

- UI スレッド
  - Qt イベントループ
  - 描画
  - メニュー、入力、状態反映
- バックグラウンドワーカー
  - 星空計算
  - 雲データ取得と描画
  - 地形地平線計算
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

### 10.3 キャッシュ方針

- キャッシュ対象は再利用価値の高い外部取得データとする。
- DEM データは永続キャッシュする。
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
