# Handoff: インスタンス変数・メソッド引数のリファクタリング調査

作成日: 2026-09-04

## 概要

インスタンス変数、データクラス化できる属性、繰り返し渡されるメソッド引数を調査した。

## 最重要候補: SkyWindow

対象: src/zstarview/gui/window.py

SkyWindowCoreMixin.__init__ では、SkyWindowUserOptions と SkyWindowRuntimeOptions を受け取った後、多数の値を個別のインスタンス変数へ展開している。

混在している要素:

- 描画・表示設定
- オーバーレイの opacity と表示フラグ
- GUI の許可フラグと導出された capability flag
- 起動処理・更新スケジュール・検索状態
- Qt ウィジェット、コントローラー、描画キャッシュ

推奨分類:

1. UserOptions / RuntimeOptions: 初期化後に変化しない設定
2. 描画設定用の可変データクラス: opacity、表示切り替え、moon style、star visibility など
3. capability / derived state: toggle supported、requested enabled、opacity when enabled
4. 可変状態: 操作、検索、更新スケジュール、起動処理、描画結果・一時データ
5. SkyWindow 本体: Qt ウィジェット、コントローラー、サービス所有権、描画キャッシュ

既存 mixin は self.some_attribute を直接参照しているため、一度に全属性を移動せず、段階的に変更する。

SkyWindowState は既に操作、検索、更新予定、描画結果を広く保持している。起動フラグや refresh-due flag をすべてそこへ追加すると別の巨大な state container になるため、属性数や更新主体を見て StartupState、RefreshScheduleState などの小さな state object に分ける。SkyWindowState への移動そのものを目的にしない。

移行中は旧属性と新しい data class に同じ値を保持しない。新しい object を唯一の正とし、必要な場合のみ旧属性名を property で一時的に公開する。

## 高優先度: SkyDataWorker の計算要求

対象: src/zstarview/gui/sky_worker.py

update() から _run_update()、さらに compute_sky_snapshot() へ、多数の天体計算用引数が渡されている。

インスタンス変数へ移すのではなく、非同期処理用の frozen request dataclass を導入する。

候補:

- SkyComputationRequest
- SkyDiscRequest
- PlanetUpdateRequest

リクエスト値を self に保持すると、処理中に次の更新で値が上書きされる可能性がある。処理単位ごとに不変の request object を渡す。

### request の不変性に関する注意

`@dataclass(frozen=True)` が保証するのはフィールドの再代入防止までであり、list、dict、NumPy array、その他の可変 object の内部まで不変にはしない。非同期 request を構築するときは、各フィールドを次のように扱う。

- sequence は可能なら tuple に変換する
- mapping は読み取り専用 view または immutable value object に変換する
- NumPy array は所有権が不明な場合は copy を取り、必要に応じて `setflags(write=False)` を適用する
- catalog など意図的に共有する大きな object は、生存中に更新しない契約を型・ドキュメント・テストで明確にする

全データの無条件コピーは、大きな catalog や地形 array でメモリと遅延を増やす。「処理中に呼び出し側が更新し得るか」をフィールドごとに判定し、必要な値だけ snapshot 化する。

## 高優先度: CloudController の pending request

対象: src/zstarview/gui/cloud_controller.py

source/render の保留リクエストを dict[str, object] として保持している。

候補:

- CloudSourceRequest
- CloudRenderRequest

コントローラー自身が保持するのは、実行中・保留中 request、request key、latest request ID、worker 集合、shutdown 状態などの実行管理状態に限定する。

pending request は「次に実行したい入力」として保持し、request ID、source ID、source key などの実行時 metadata は active request へ昇格させる時点で付与する。ID の発行、latest ID の更新、active key の設定は同じ lock 内で原子的に行う。source/render のそれぞれに、queued request から active request へ変換する共通 helper を用意し、初回実行と pending 再実行で同じルールを使う。

## 中優先度: オーバーレイコントローラーの引数

対象: aircraft_controller.py、satellite_controller.py、terrain_controller.py、water_overlay_controller.py、urban_outline_controller.py、road_night_lights_controller.py、precipitation_controller.py

observer location、時刻、reason、描画条件などが update() から worker メソッドへ繰り返し渡されている。これらは毎回変わり得る要求データであり、安易にインスタンス変数へ移さない。長いシグネチャだけ controller ごとの request dataclass にまとめる。

## 低優先度: Cloud / Geo-satellite state の共通化

対象: src/zstarview/gui/cloud_state.py、src/zstarview/gui/geosatellite_state.py

画像、missing mask、grid、source key、render key、時刻、coverage、banner と set_result() 処理が重複している。共通の OverlayImageState または value object を検討するが、設定・request 整理の後に行う。

## 既存の有効な境界

- PreparedWindowCatalogs
- SkyWindowUserOptions
- SkyWindowRuntimeOptions
- SkyWindowState
- RenderInputs
- RenderSceneData
- RenderStyle

既存の境界を活用し、SkyWindow の flat な属性参照を減らす。

## 推奨実装順

1. SkyWindow の属性を設定・導出値・可変状態・外部オブジェクトに分類する
2. 表示設定データクラスを導入し、_render_style() の入力をそこから構築する
3. 起動・更新スケジュール状態を責務別の小さな state object へ段階的に移す
4. SkyDataWorker の巨大な計算引数を request dataclass 化する
5. CloudController の pending request 辞書を typed dataclass 化する
6. 他の overlay controller へ同じ整理を適用する
7. 最後に state の共通化を検討する

## 検証方針

各段階で関連テスト、GUI の offscreen テスト、ruff check、mypy、compileall、git diff --check を実行する。

非同期処理では、新しい request が実行中 request を上書きしないこと、stale result の破棄条件、shutdown 中の emit 抑制、pending request の再実行を確認する。

追加で、次をテストする。

- request 構築後に呼び出し元の list / array を変更しても、snapshot 化対象の request が変化しないこと
- pending request の最新の1件だけが実行され、active への昇格時に新しい request ID が付与されること
- source が差し替わった場合に、古い source ID / source key を持つ render result が破棄されること
- 移行中の互換 property と新しい state object が同じ値を参照し、二重管理になっていないこと

## 作業範囲

- 調査と候補整理のみ
- 実装変更なし
- 詳細な調査メモは dev-notes/session-2026-09-04.md に追記済み
