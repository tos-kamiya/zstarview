# リファクタリング候補 引き継ぎ資料

作成日: 2026-07-19

## 目的

現行実装で単一ファイルに集中している責務を整理し、今後の機能追加や不具合調査で変更範囲を小さくするための候補を記録する。
これは実施済みの変更ではなく、次回以降の検討資料である。

## 優先度

### 優先度A: GUIウィンドウ本体の分割

対象: `src/zstarview/gui/window.py`（約3,094行）

現在は、ウィンドウの初期化、Qtアクション、イベント処理、設定反映、タイマー、ライフサイクル、表示状態の管理が一つのクラスに集中している。

候補:

- Qtアクションとメニュー構築を `window_actions.py` へ分離する。
- キーボード・マウス・視点変更などの入力処理を `window_input.py` へ分離する。
- 初期化・終了・バックグラウンドworker接続を `window_lifecycle.py` へ分離する。
- `SkyWindow` 本体には、依存オブジェクトの組み立てとQtイベントへの委譲だけを残す。

注意点:

- 大きなMixin再編を一度に行わない。
- 公開されている起動経路と `paintEvent()` の挙動を変えない。
- 分割後も `SkyWindow` の属性を直接参照する既存コードとの互換性を確認する。

### 優先度A: 描画ウィンドウMixinの責務分割

対象: `src/zstarview/gui/window_render.py`（約1,245行）

今回の恒星補間で、ベースフレーム、表示フレーム、恒星サーフェイス、HUD、volatile overlay、動的惑星が同じMixinに集まった。補間の不具合調査では、キャッシュ責務と描画責務を分けた方が追跡しやすい。

候補:

- `window_render_cache.py`: フレームキャッシュキー、キャッシュstamp、恒星サーフェイスキャッシュ。
- `window_render_frame.py`: `FrameContext` とベース／表示フレームの生成。
- `window_render_overlays.py`: 惑星、人工衛星、航空機、アステリズム、恒星ラベルなどの現在フレーム描画。
- `window_render_hud.py`: HUD、簡易表示、マウス位置に依存する表示。

第一段階では、現在の `SkyWindowRenderMixin` の公開メソッドを維持し、内部処理を補助Mixinまたは専用関数へ移すのが安全である。

### 優先度B: 更新処理Mixinの分割

対象: `src/zstarview/gui/window_updates.py`（約1,775行）

空、雲、衛星、航空機、台風、地形、水面、都市アウトライン、動的惑星、ステータス表示が同じファイルにある。スケジューラの変更が各レイヤーへ波及しやすい。

候補:

- `window_sky_updates.py`: sky worker、60秒スナップショット、動的惑星の2秒更新。
- `window_overlay_updates.py`: 雲、衛星、航空機、台風の取得・再投影。
- `window_terrain_updates.py`: 地形、水面、都市アウトライン。
- `window_status.py`: ステータス行と各レイヤーの状態表示。
- スケジューラ本体は `window_updates.py` に残し、各更新処理へ委譲する。

動的惑星更新は恒星補間と時間周期が異なるため、恒星補間ロジックと同じモジュールへ統合しない。

### 優先度B: 恒星描画モジュールの分割

対象: `src/zstarview/render/stars.py`（約1,309行）

恒星のラスター生成、描画キャッシュ、明るい恒星の下地、ダイヤモンド、ヒット判定、名前ラベルが混在している。

候補:

- `star_raster.py`: numpy canvas、矩形・ダイヤモンド・小星のラスター処理。
- `star_cache.py`: 恒星サーフェイスのキャッシュキーとQImageキャッシュ。
- `star_interaction.py`: 恒星・惑星のヒット判定と名前解決。
- `star_labels.py`: 可視ラベルの収集と描画。
- `stars.py`: 外部APIと描画経路の調整だけを残す。

恒星補間の3x3行列そのものは既に `render/star_interpolation.py` に分離されているため、そこへ恒星ラスター処理を戻さない。

### 優先度C: 描画パイプラインのレイヤー分割

対象: `src/zstarview/render/pipeline.py`（約1,016行）

描画順序の制御と、各レイヤーの具体的な描画処理が同じファイルにある。

候補:

- `pipeline_base.py`: ベースシーンと共通の描画順序。
- `pipeline_stars.py`: 恒星サーフェイス、明るい恒星、ラベル。
- `pipeline_overlays.py`: 惑星、衛星、航空機、台風など。
- `pipeline_hud.py`: 観測情報、簡易表示、ホバー表示。

ただし、まず `window_render.py` のキャッシュ責務を整理してから着手する。描画順序とキャッシュ境界を同時に変更すると、回帰原因の特定が難しくなる。

## テスト側の候補

- `tests/test_window_render_sync.py`（約6,575行）を、キャッシュ同期、表示フレーム、オーバーレイ、HUD、worker payload の領域別に分割する。
- `tests/test_planet_marker_style.py`（約1,776行）を、恒星・惑星・太陽／月・ラベルのテスト群に分ける。
- 分割時はテスト名と検証内容を変えず、失敗時の責務が明確になることを優先する。

## 推奨実施順

1. まず `window_render.py` のキャッシュ処理だけを抽出する。
2. 恒星描画のキャッシュ、ヒット判定、ラベルを分離する。
3. `window_updates.py` のスケジューラと各worker更新を分離する。
4. `window.py` の入力・アクション・ライフサイクルを分離する。
5. 最後に `pipeline.py` のレイヤー実装を分離する。

各段階で、公開メソッドの引数と `SkyWindow` の状態属性を維持し、対象テストを実行してから次へ進む。

## 分割しない方がよいもの

- `src/zstarview/render/star_interpolation.py`（約202行）は現在の責務に対して十分小さく、分割不要。
- `src/zstarview/gui/sky_worker.py`（約458行）はsnapshot計算とworker実行の境界が比較的明確で、現時点では分割優先度が低い。
- `CelestialData.star_time` は恒星スナップショットの基準時刻を示すデータ項目であり、独立モジュール化しない。

## 検証方針

- 各分割前後で対象テストを実行する。
- `QT_QPA_PLATFORM=offscreen` のGUIテストを含める。
- 恒星補間では、60秒境界、`-30`〜`+30`秒、75%補正、cache key、アステリズム・ラベル追随を確認する。
- Atlasとexport-imageが通常ビューアの分割変更から影響を受けないことを確認する。
- 分割だけの変更では、描画仕様や更新周期を同時に変更しない。

## 関連資料

- [`star-interpolation-handoff-2026-07-19.md`](star-interpolation-handoff-2026-07-19.md)
- [`docs/specification.md`](../docs/specification.md)
- [`docs/design/rendering-pipeline.md`](../docs/design/rendering-pipeline.md)
- [`docs/design/gui-screen-update-and-cache.md`](../docs/design/gui-screen-update-and-cache.md)
