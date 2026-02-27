# zstarview 設計書（初版）

最終更新: 2026-02-27

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

### 2.3 描画

- `src/zstarview/render/draw.py`
  - 背景、星、惑星、ラベル、オーバーレイ描画
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
  - 雲更新（取得・例外正規化・running/pending制御・清掃トリガ）をバックグラウンドで実行
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
- `src/zstarview/clouddisc/providers/meteosat.py`
- `src/zstarview/clouddisc/providers/_s3_io.py`
  - S3一覧/ダウンロードの共通I/Oヘルパー
  - タイムアウト/一般失敗を `TimeoutError` / `DownloadError` へ正規化

補足（段階導入）:
- `providers/select.py` には将来プロバイダ向けの実験衛星定義を先行追加する場合がある。
- 実験衛星は `core.py` 側の provider 分岐へ接続されるまで `AUTO` 選択対象には含めない。
- `providers/meteosat.py` は PR3 で `core.py` に接続済み。現時点では `AUTO` には含めず、手動優先順位での選択対象とする。

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
2. `CloudController` が `CloudDisc.render_now()` を実行し、例外をUI向けステータスへ正規化。
3. `cloud_ready` / `cloud_failed` シグナルで `SkyWindow` に結果を返す。
4. `SkyWindow` が `CloudImageState` と合成キャッシュを更新し、再描画要求をUIスレッドへ配送。
5. 終了開始時（`aboutToQuit` / `closeEvent`）はタイマー停止に加えて `CloudController.shutdown()` を呼び、以後の更新受付を抑止。

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
  - UIはバナー表示し、更新ループは継続

## 7. 設定/永続化

- 設定ファイル:
  - `~/.config/zstarview/config.json`
- 保存情報:
  - 最終地点（都市形式または座標形式）
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

## 9. 今後の設計課題

- `CloudController` 抽出後の責務整理は、必要時に再評価する（現時点では `SkyWindow` 側との分担は実用上バランスが取れている）。
- 描画パイプラインの入力データ契約（型/単位）を明文化する。
- 仕様変更時に追従しやすいよう、CLI仕様と内部データ仕様の対応表を追加する。
- ダイナミックなタイムシフト操作（`-H` のGUI版）を追加する。
  - 目的:
    - 現在は `-H` が起動時固定値のため、実行中に時刻を前後できない。
    - GUI上で前進/巻き戻し/停止を行えるようにし、時刻探索を容易にする。
  - 仕様案:
    - キー操作中心で実装（メニュー操作には依存しない）。
    - 想定キーは `<` / `>`（実体は `Shift+,` / `Shift+.`）。キーボード配列差を考慮して予備キー併設も検討。
    - 押しっぱなし操作は `keyPressEvent`/`keyReleaseEvent` で押下状態を保持し、`QTimer` で連続更新する。
    - OSのキーリピート速度に依存しないため、操作感を一定化できる。
  - 内部設計案:
    - `SkyWindow` に時間制御状態を追加:
      - `time_shift_base: timedelta`
      - `play_rate: float`（正: 未来方向、負: 過去方向）
      - `last_tick_monotonic: float`
    - 実効オフセット:
      - `effective_delta_t = time_shift_base + elapsed * play_rate`
    - 星空更新は既存 `SkyDataWorker` 経路を再利用し、`delta_t` のみ動的更新する。
  - UI表示案:
    - オーバーレイへ `Shift` と `Rate` を追加表示（例: `Shift +02:30:00  Rate +10x`）。
  - 雲表示との整合:
    - 現在仕様どおり、実時間以外は雲を無効化する。
    - `effective_delta_t` が 0 近傍（しきい値例: ±30秒）のときのみ雲有効とし、境界チラつきを抑制する。
  - 段階導入:
    - PR1: 時間制御状態 + キー入力 + オーバーレイ表示
    - PR2: 速度プリセット/リセット導入と境界条件調整
    - PR3: 必要に応じてキー割当の設定化

## 10. 入力データ契約の明文化プラン（DONE / INPROGRESS）

状態:
- `DONE`: 1, 2, 3, 4, 5
- `INPROGRESS`: なし

目的:
- 角度・座標系・単位の取り違えを防止し、関数間契約を読み取れる状態にする。
- 既存挙動を壊さず、段階的に契約を強化する。

DONE:

1. 高: `alt/az` 引数順の統一（誤用バグ予防）
   - 実施:
     - `clouddisc/projectors/az.py` に `altaz_to_dir_ecef(alt_deg, az_deg, ...)` を導入。
     - 旧 `azalt_to_dir_ecef(...)` は互換ラッパーを経て削除済み。

2. 高: 型定義に単位・順序を反映
   - 実施:
     - `types.py` に `AltDeg/AzDeg/LatDeg/LonDeg/ViewCenterAltAz` を導入。
     - `ViewerData`, `PlanetBody` に `*_deg` 明示アクセサを追加（既存フィールド互換を維持）。

3. 中: `stars` 構造の契約を型化
   - 実施:
     - `types.py` に `StarsTable(TypedDict)` を追加。
     - `CelestialData.stars` と `calculate_visible_stars()` の戻り型を `StarsTable` へ変更。

4. 中: 正規化座標の仕様ドキュメントを実装へ一致
   - 実施:
     - `render/draw.py` の `altaz_to_normalized_xy_vectorized(...)` docstring を
       実装仕様（`90deg -> 1.0`、`1.0`超過あり得る）へ更新。

5. 低: パラメータ命名の単位接尾辞化
   - 実施:
     - `render/draw.py` で `alt/width/height/radius` 主要引数を
       `view_alt_deg/width_px/height_px/radius_px` へリネーム。

実施順:
- `1 -> 2 -> 4 -> 3 -> 5`
