# zstarview 設計書（初版）

最終更新: 2026-03-05

## 1. 設計方針

- UI責務と計算責務を分離する。
- 長時間処理はバックグラウンドスレッドで実行する。
- 失敗しやすい外部I/O（S3）を局所化し、例外をドメイン例外に正規化する。
- 描画は純粋関数寄りに保ち、状態管理は `SkyWindow` に集約する。

## 2. モジュール構成

### 2.1 エントリ/起動

- `src/zstarview/zstarview.py`
  - CLI引数解析
  - ログ初期化
  - 都市解決・時刻解釈・星カタログ読み込み
  - `SkyWindow` 起動

### 2.2 ドメイン計算

- `src/zstarview/astro.py`
  - 天体位置計算
  - 可視判定
  - 補助線（地平線/赤道/黄道）点列生成
  - 食情報計算
  - 恒星カタログの前処理で RA/Dec のラジアン・三角関数・ICRS 単位ベクトルをキャッシュし、`SkyCoord.transform_to` を毎回呼ばず Alt/Az 行列変換で高速化

### 2.3 描画

- `src/zstarview/render/draw.py`
  - 背景、星、惑星、ラベル、オーバーレイ描画
- `src/zstarview/asterisms.py`
  - アステリウム定義（季節/表示名/恒星 `SourceId` の経路）
  - ホバー恒星に対するアステリウム逆引き
  - 時刻スロットによる回転表示対象選択
- `src/zstarview/render/draw_sky_disc.py`
  - 太陽高度・食係数を用いた空色ディスク生成

### 2.4 UI

- `src/zstarview/ui/window.py`
  - メインウィンドウ
  - 入力イベント
  - タイマー管理
  - ワーカー/コントローラ連携と描画統合
- `src/zstarview/ui/sky_worker.py`
  - 星空更新計算（天体計算 + 空色ディスク）をバックグラウンドで実行
  - `data_ready` シグナルで結果返却
- `src/zstarview/ui/cloud_controller.py`
  - 雲更新（ソース取得/レンダリング分離・例外正規化・running/pending制御・清掃トリガ）をバックグラウンドで実行
  - `cloud_started` / `cloud_ready` / `cloud_failed` シグナルでUIへ通知
- `src/zstarview/ui/composite.py`
  - 空色/雲合成、ハッチ処理、合成キャッシュ
- `src/zstarview/ui/cloud_state.py`
  - 雲画像・メタ・バナーの状態保持

### 2.5 クラウドデータ

- `src/zstarview/clouddisc/core.py`
  - 衛星選択、投影、サンプリング、画像化のオーケストレーション
- `src/zstarview/clouddisc/providers/goes.py`
- `src/zstarview/clouddisc/providers/hima.py`
- `src/zstarview/clouddisc/providers/_s3_io.py`
  - S3一覧/ダウンロードの共通I/Oヘルパー
  - タイムアウト/一般失敗を `TimeoutError` / `DownloadError` へ正規化

## 3. 主要データ構造

- `types.ViewerData`
  - 観測地点、タイムゾーン、都市名、視点中心
- `types.CelestialData`
  - 描画対象の天体データ一式
- `clouddisc.types.CloudMeta`
  - クラウド画像のデータソース情報

## 4. 処理フロー

### 4.1 起動フロー

1. `main()` が CLI とロガーを初期化。
2. `_startup_resolve_city()` で地点を決定し、次回用設定を保存。
3. 時刻オフセットと星カタログを読み込み。
4. `SkyWindow` を生成し、初回データ読み込み後に表示。

### 4.2 星空更新フロー

1. `SkyWindow.start_background_sky_data_update()` が `SkyDataWorker.update()` を呼ぶ。
2. `SkyDataWorker` がバックグラウンドスレッドで天体計算と空色ディスク生成。
3. `data_ready` シグナルで `SkyWindow._on_sky_data_calculated()` に反映。

### 4.3 雲更新フロー

1. `SkyWindow.start_background_cloud_update()` が `CloudController.update()` を呼ぶ。
2. `CloudController` は `SourceKey`（衛星・時刻スロット）でソース取得を管理し、`RenderKey`（`SourceKey` + 視点情報）で描画要求を管理する。
3. ソース取得は `fetch_source()`、描画は `render_from_source_with_coverage()` で分離実行する。
4. 描画結果は `request_id` で世代管理し、`latest-request-wins`（古い結果は破棄）を適用する。
5. `cloud_ready` で `image` + `missing_mask` + `coverage_ratio` を `SkyWindow` へ渡し、合成時に欠損領域へ薄い黄色ティントを重畳する。
6. `cloud_failed` 受信時はバナー表示に加え、雲画像と欠損マスクをクリアする。
7. 終了開始時（`aboutToQuit` / `closeEvent`）はタイマー停止に加えて `CloudController.shutdown()` を呼び、以後の更新受付を抑止。

## 5. スレッドモデル

- UIスレッド:
  - Qtイベントループ、描画、ウィジェット状態更新
- バックグラウンド:
  - 星空計算: `SkyDataWorker`
  - 雲更新: `SkyWindow` 管理スレッド
  - キャッシュ清掃: 雲更新側の補助スレッド

注意:
- UI操作（再描画要求を含む）はシグナル経由でUIスレッドに戻す。
- シャットダウン中はワーカー側 `emit` を抑止し、Qtオブジェクト破棄後のシグナル送出例外を回避する。

## 6. エラーモデル

- 起動系エラー:
  - `StartupAbortError` で起動シーケンスを中断
- クラウド系エラー:
  - `CloudDiscError` 派生（`TimeoutError`, `DownloadError`, `DataNotFoundError`, `RenderError`）
  - UIはバナー表示し、当該時点の雲レイヤーはクリアして非表示
  - 取得成功だが投影可能領域に欠損がある場合はエラー扱いにせず、部分描画 + 欠損ティントで継続表示
  - 更新ループは継続

## 7. 設定/永続化

- 設定ファイル:
  - `~/.config/zstarview/config.json`
- 保存情報:
  - 最終地点（都市形式または座標形式）
- 起動時表示オーバーライド（CLI）:
  - `--show-dso-initial true|false`
  - `--show-asterisms-initial true|false`
  - これらは起動時のみ有効で、設定ファイルへの永続化はしない。
- キャッシュ:
  - 衛星データ・エフェメリスはキャッシュディレクトリに配置

## 8. テスト戦略（現状と推奨）

- 現状:
  - 軽量な回帰テストを追加済み
    - 起動時地点保存経路のガード
    - 共通S3 I/Oの例外変換
- 推奨:
  - `astro.py` の純計算関数テストを拡張
  - `ui` はシグナル発火と状態遷移を中心にモックテスト
  - ネットワーク実アクセスを避け、providerはモック前提で検証

## 9. 機能設計（実装済み）

- 有名恒星ジャンプUI（仕様確定）
  - 目的:
    - 「代表的な有名恒星」をメニューから選択し、対象恒星が視野中心に入る向きへ即時移動できるようにする。
  - 初期スコープ:
    - 対象は `Name` が付与された恒星のうち `Vmag <= 2.0` を候補集合とする（現行データで約48件）。
    - 候補は赤緯で3分割して表示する:
      - 北天: `Dec >= +20°`
      - 赤道付近: `-20° < Dec < +20°`
      - 南天: `Dec <= -20°`
  - UI案:
    - ハンバーガーメニューに「Jump to Famous Star...」を追加。
    - 検索専用導線として、別メニュー「Search Famous Star...」を追加する（検索対象は固有名付き恒星の全件）。
    - ダイアログで領域タブ（北天/赤道付近/南天）+ 恒星リストを表示し、選択確定で視点を更新する。
    - ジャンプ確定後3秒間は、マウスホバー時と同じ見た目（円マーカー + 星名ラベル）で対象恒星を強調表示する。
  - 視点更新方針:
    - 恒星の `RA/Dec` を現在時刻・観測地点で `Alt/Az` へ変換し、`view_center = (alt, az)` へ設定する。
    - `alt` は既存仕様に合わせ `[0, 90]` にクランプする。
    - 更新後は既存の天球再計算経路（`request_sky_data_update`）を再利用する。
  - 受け入れ観点:
    - ダイアログで選んだ恒星が、確定後に画面中心近傍へ移動する。
    - 3カテゴリ分割で一覧操作が過密にならない（初期データで北12/赤道15/南21）。
    - 既存ショートカット（矢印/F11/Q/M）と競合しない。

- DSO（銀河・星団）カタログ生成基盤
  - OpenNGC を `pyongc`（開発依存）経由で取り込み、`src/zstarview/data/dso.csv` を生成するスクリプトを追加した。
  - 生成物は `Id, Name, Type, RAh, Dec, Vmag, MajorArcmin, MinorArcmin, PAdeg, SourceCatalog` を基本カラムとする。
  - 初期対象は `Type in {galaxy, open_cluster, globular_cluster}` かつ `Vmag <= 12.0`。
  - ランタイム依存には `pyongc` を追加せず、「生成して同梱」方式を維持する。

- DSO（銀河・星団）表示統合
  - `dso.csv` を起動時に読み込み、天球計算で `Alt/Az` へ変換して描画へ渡す。
  - 通常表示は「ネームドかつ形状情報あり」のみを薄い青系塗りで背景レイヤーに描画する。
  - ホバー判定は恒星/惑星と独立に実行し、DSOがヒットした場合は3x形状で強調表示する。
  - レイヤー順は「DSO本体 < 恒星/惑星本体」。ラベルは前景レイヤーで描き、ラベル同士は「DSOラベル < 恒星/惑星ラベル」。

- アステリウム表示（ホバー連動）
  - 目的:
    - 画面の恒星ホバーを起点に、その恒星を含む代表的アステリウムを文脈表示し、常時表示による過密を避ける。
  - データモデル:
    - `asterisms.py` に `Asterism(key, name, edges)` を定義し、`edges` は恒星カタログの `SourceId`（`HIP...`）同士の線分集合で記述する。
    - `ASTERISM_KEYS_BY_SOURCE_ID` で「恒星 -> 所属アステリウム一覧」を逆引きする。
  - 表示ルール:
    - 起動時の有効/無効はメニュー状態に従う（デフォルト有効、CLIで `--show-asterisms-initial` により上書き可能）。
    - 有効時は、マウスホバーで選択された恒星に `source_id` があり、かつアステリウム所属がある場合のみ表示する。
    - 所属数を `n` とし、表示対象はシステムクロック基準の 3 秒スロット `floor(unix_time/3) % n` で1件選ぶ。
    - 同時表示は常に1件のみ。選択されたアステリウムの線分と名称ラベルを描画する。
  - 描画スタイル:
    - 線分は既存オーバーレイと同様にアウトライン付きで描画し、テーマに応じた青系カラーを使用する。
    - 名称ラベルは DSO 強調表示と同系統の青色を使用する（light/day と night/black で色を調整）。
  - UI/タイマー:
    - ハンバーガーメニューの `Asterisms` トグルで有効/無効を切り替える。
    - DSO と同様に、起動時状態は CLI（`--show-*-initial`）で上書きできる。
    - `SkyWindow` に 1 秒周期のチェックタイマーを置き、tick ごとにホバー選択を再評価して必要時に再描画する。

## 10. 今後の設計課題

- `CloudController` 抽出後の責務整理は、必要時に再評価する（現時点では `SkyWindow` 側との分担は実用上バランスが取れている）。
- 描画パイプラインの入力データ契約（型/単位）は規約化済み。将来の変更は `docs/developer/input-data-contract.md` に従う。
- 仕様変更時に追従しやすいよう、CLI仕様と内部データ仕様の対応表を追加する。
- マウス操作のインタラクション拡張（検討中）:
  - 候補A: ドラッグで視点回転（地図アプリ風パン相当）
  - 候補B: クリック点を画面中心へ移動
  - UX比較観点: 誤操作率、既存ホバー/クリックとの競合、学習コスト、キー操作との一貫性
  - 方針: 当面は現行操作（キー中心）を維持し、必要性が高まった時点で段階導入する


## 11. 作業記録の扱い

- 実装履歴は [implementation-history.md](./implementation-history.md) に分離する。
- セッション単位の作業ログは `dev-notes/session-YYYY-MM-DD.md` に記録する。
