# GUIクラッシュ対策・並行処理再設計ハンドオフ

作成日: 2026-09-06

状態: 第1移行サイクル完了。診断基盤、Qt非依存IPC、cyclone取得のプロセス分離まで
実装済み。sky/render以降と長時間運用検証は未完了。
本資料の未完了項目は現在の仕様・設計を置き換えない。

## 1. 目的と診断上の前提

科学計算・データ変換中に子プロセスが異常終了しても、GUIが最後の正常な
表示を保ち、失敗した機能だけを停止・再起動できる構成にする。
同時に、GUI内の任意のワーカースレッドからQt関連オブジェクトを操作・解放する
経路を減らし、処理とデータの所有者を明確にする。

根本原因は未確定である。過去の説明・コメント・リリースノートにある
「NumPy等の競合を防いでクラッシュを解消」という断定は、検証済みの事実として
引き継がない。

- 1回目: GMNのキャッシュからdataclassを作るワーカーでGC中にSIGSEGV。
  別ワーカーはSkyfieldのtimescale読み込み中だった。
- 2回目: cycloneの座標列をPythonリストへ変換するワーカーでGC中にSIGSEGV。
  別ワーカーはexecutorの待機箇所にいた。科学計算同士の同時実行は確認できない。
- 両方ともGUIスレッドのPythonスタックは `render/text.py` の
  `_rect_overlap_count()` にあり、Qtの矩形操作を呼び出していた。
- GCを実行したのはGUIメインスレッドではなく、クラッシュしたワーカー。
  GC中という記載は、既存ヒープ破壊を検出したことの証明ではない。
  GC自身、拡張の走査・解放処理、オブジェクト寿命、以前の破壊などを区別できない。
- 04:09〜07:53の約3時間44分稼働はPython 3.14環境での記録。
  Python 3.13への変更後の安定性を示す記録ではない。
- 実行環境の番号がvenv 33から31に変わっている。各実行のアプリ版、依存版、
  GIL状態は別途確認する。修正版が実際にロードされたかも確認する。
- GIL無効版の既報と同一原因である根拠はない。ロック追加テストの成功も
  クラッシュ原因の特定や解消を証明しない。

プロセス分離後もGUI内のQt/Python由来のクラッシュは残り得る。
本提案は原因調査と並行して進める耐障害性の改善である。

## 2. 現行構成と移行対象

| 現行箇所 | 現行の責務・境界 | 移行先 |
| --- | --- | --- |
| `gui/application_services.py` | 2スレッドの共有executor、nativeロック、ephemeris | プロセス監督と要求管理。ephemerisは計算側が所有 |
| `gui/sky_worker.py` | 天体計算・空の画像生成、Qt Signalで結果通知 | 常駐sky/renderプロセス |
| `gui/meteor_controller.py` | 取得・キャッシュ・投影を共有ロック内で実行 | 取得プロセスと計算プロセス |
| `gui/tropical_cyclone_controller.py` | 取得・キャッシュ・大量の辞書化 | 取得プロセスで変換まで完結 |
| `gui/aircraft_controller.py`, `gui/satellite_controller.py` | 主にデータ取得 | 取得プロセス |
| `render/aircraft.py`, `render/satellites.py` | 描画に伴う座標計算 | 座標計算を計算側へ、GUIは完成済み描画データを使用 |
| `gui/cloud_controller.py`, `gui/geosatellite_controller.py` | 雲データの投影・画像化等 | rasterプロセス |
| `clouddisc/workers/cloud_source_worker.py` | 一回限りの子プロセス、JSON管理情報＋pickle結果 | 既存方式を段階的に新プロトコルへ接続 |
| terrain、urban、water各controller | 地理データ取得・変換 | 地理処理プロセス |
| `gui/precipitation_controller.py`, `gui/road_night_lights_controller.py` | 独自thread | 取得／地理処理プロセスへ統合 |
| 検索・起動時解決・月／太陽hover取得 | 共有executor利用 | 取得／計算へ分類。画像デコードの場所も確認 |
| GUIのpaintEvent、HUD、対話入力 | Qt描画・状態変更 | GUIメインスレッドに保持 |

一覧は入口の棚卸しである。実装前に全 `Thread`, `submit`, `QThread`,
プロセス起動箇所と、その呼び先を再検索する。export-image等のCLIは別入口として
互換性を確認し、GUI専用サービスへの依存を純粋計算関数に持ち込まない。

## 3. 提案するプロセス・スレッド構成

GUIが単一の監督役を所有し、要求の世代管理・受付上限・失敗処理を集約する。
各controllerはUI状態と要求内容を管理し、任意のPython関数をexecutorに渡す
方式から、型の決まったジョブを送る方式へ移行する。

- GUI: Qtイベントループ、入力、状態管理、準備済み画像・primitiveの描画。
  移行完了後はGUIプロセス内のPythonバックグラウンド計算スレッドをなくす。
  Qt内部スレッドまで無効にするという意味ではない。
- sky/render: 常駐、同時に1要求。カタログとephemerisを一度ロードする。
  天体・衛星・航空機・meteorの座標と、必要な空／星の画像を計算する。
- raster/地理処理: 重い雲・地形等を隔離する。初期移行では一回限りのプロセスを
  優先し、起動時間・キャッシュ効果を測ってから常駐化を判断する。
- 取得: HTTP・キャッシュ・JSON変換をGUIから隔離する。長いダウンロードが
  sky更新を妨げない独立枠とし、実行数・待機数に上限を設ける。

全controllerごとの常駐プロセスは作らない。開始時の候補はsky 1枠と
バックグラウンド重処理／取得合計2枠。既存Cloud子プロセスも上限に含め、
CPU・メモリ測定で調整する。多段の子プロセス生成は原則避ける。

起動は `sys.executable -m ...` を使い、GUIを起動した後のforkによる状態継承は
行わない。Qtオブジェクトはプロセスをまたいで渡さない。workerが画像生成に
Qtを使う場合はworker内で初期化・終了を完結し、offscreen環境を検証する。

GUI側の監督にはQProcessとイベント通知を第一候補とする。stdoutには小さな
通知のみを流し、巨大な結果や長時間待機、重い読み込みをGUIコールバックへ
持ち込まない。具体的APIは実装開始時に対応PySide6版で確認する。

## 4. IPCとデータ所有権

制御はバージョン付きJSON、配列・画像は専用作業ディレクトリ内の成果物とする。
初期段階でMessagePackや共有メモリへの依存は追加しない。

要求に含める項目:

- protocol_version、session_id、worker_epoch、request_id、job_kind
- layer_generation、view_generation、入力データのrevision
- UTC表示時刻、観測地・標高、視線、画面サイズ・device pixel ratio、必要な設定
- 親が管理する期限とworker内の処理時間予算
- 入力成果物の相対パス、schema、サイズ、配列shape/dtype等

結果は同じ識別情報に加えてstatus、成果物一覧、データ時刻、coverage、
処理時間、分類済みエラーを持つ。異常終了時は結果JSONがなくても親が
終了状態から失敗を確定できるようにする。

| データ | 初期形式 | 所有・転送方針 |
| --- | --- | --- |
| 少数の天体・ステータス | JSON | 度・距離の単位、UTC、nullを明文化。NaN/Infinityを送らない |
| 大きな座標列・polygon | 数値型のNPY＋JSON管理情報 | object dtype禁止、pickle禁止、offset配列でring等の境界を表す |
| カタログ・ephemeris | バージョン付き入力パス | workerがロードし保持。毎フレーム転送しない |
| 完成画像 | 固定形式RGBAバイト列＋幅・高さ・stride | 圧縮コストを避ける。GUIはサイズ検証後に所有コピーを作る |
| 障害ログ | stderr／worker.log | 起動情報、終了コード、要求情報と紐付け、保存量を制限 |

高密度配列は可能ならworker内で最終描画データまで縮約する。GUIにNPYを
読み込ませる必要がある領域は残存するNumPy利用として明示し、障害境界の
達成度に含める。PNGはデバッグ出力や小さな低頻度画像の候補に留める。

成果物は一時名へ書き、close後にatomic renameし、最後にmanifestを確定する。
親はschema、ID、世代、許可サイズ、寸法、パスの作業領域内包含を検証する。
古い結果はUIへ反映せず廃棄する。新しい画像はロード完了後に一括交換する。

GUIが外部バッファからQImageを構築する場合、コピー完了までバッファを保持し、
独立した所有コピーになってから成果物を削除する。将来ゼロコピーにする場合は
参照寿命・release通知・クラッシュ時回収を別途設計する。

既存Cloud方式はpickleを使用しているため、そのまま汎用IPC標準にはしない。
移行中は内部の信頼済み成果物だけを扱う既存adapterとして保持し、後段で置換する。

## 5. 要求管理、復旧、終了

- 同種の視線更新は「実行中1＋最新待機1」。古い待機要求を置換する。
  取得結果は視線と独立したrevisionで保持し、視線変更のたびに再取得しない。
- 対話中は直近フレームや既存の簡易表示を使い、操作停止後に最新要求を優先する。
  子プロセス待ちでpaintEventを止めない。
- sky全体の結果は同一要求の時刻・視線として反映し、異なる世代の星と空を
  混在させない。独立レイヤーには個別世代とデータ時刻を持たせる。
- deadline超過、異常終了、壊れた結果、通常の取得失敗を区別する。
  最終正常結果を保持し、stale／unavailableをUIに表示する。
- 再起動はバックオフ付き・回数上限付き。例として5分以内に3回異常終了した
  レイヤーは自動再試行を止め、手動再試行可能にする。数値は運用で調整する。
- 通常のキャンセルは古い結果を無視する。長時間処理は協調キャンセルし、
  応答しないworkerは期限後にterminate、猶予後にkillする。
- 終了時は新規受付停止、待機破棄、子プロセス終了、参照解放後の成果物回収。
  GUIスレッドで無期限waitしない。親終了時の孤児化も対応OSで検証する。
- ログと失敗成果物の保存には容量・件数制限を設ける。PIDだけで所有を判定せず、
  session IDで他の起動中アプリの作業領域を削除しない。

## 6. 移行手順と実装担当者の成果物

1. 診断情報を整備する。実行中アプリのPython・依存版・GIL状態・app revisionを
   採取し、可能ならcore dumpのCスタックをシンボル化する。Qtラッパーの寿命、
   GUI外のGC・解放との関係も調査する。GC無効化を恒久対策にしない。
2. Qtに依存しないIPC schemaとプロセス監督を独立モジュールとして実装する。
   テスト用workerの異常終了・timeout・壊れた結果で親の復旧を検証する。
3. cyclone取得・辞書化を最初の小さな移行単位とする。今回の原因箇所と断定せず、
   小さいインターフェースで復旧・世代管理を実証する目的で選ぶ。
4. 常駐sky workerへ天体計算を移し、カタログ・ephemerisの所有を移す。
   衛星・航空機の描画内計算も移し、meteorを接続する。
   時刻・視線の契約を保ち、起動時の `ApplicationServices` のastro import等も整理する。
5. 空・星・雲画像生成、地理処理を移行する。既存Cloud subprocessを監督へ統合し、
   pickle成果物を置換する。画像コピー・遅延・メモリ量を測定する。
6. 降水・道路夜間表示の独自thread、検索・hover・起動処理を移行する。
   全てのスレッド入口を再監査してから旧executorを廃止する。
7. `native_work_lock` は旧経路が残る間は保持する。GUI側でそのロックを待つ
   設計を追加しない。移行完了後、残存利用を検査し、AGENTS.mdの指針を
   「プロセス境界・所有権・GUI外計算」に更新する。

各段階は独立した変更・比較可能な実行経路とする。移行中の切替は内部設定で
よいが、worker失敗時に旧GUI内計算へ自動フォールバックしない。
恒久的な二重実装を残さず、各段階の合格後に旧経路を削除する。

## 7. 検証・受け入れ条件

- 異常終了を注入してもGUIが操作でき、直近正常フレームが維持される。
  制限付き再起動で最新要求を回復できる。親への影響が残ればその境界を明記する。
- timeout、途中ファイル、未知schema、過大shape、古い世代、逆順完了、
  子の再起動前の結果をテストする。結果通知と終了通知の順序差も扱う。
- 起動・リサイズ・視線変更・時刻変更・機能無効化・終了中の競合を検証する。
  誤った時刻／視線の画像を一瞬でも反映しない。
- 既存の座標・時刻・投影の純粋ロジックテストを再利用し、固定入力の数値差と
  画像差を比較する。Qtを必要とする統合試験は通常の単体試験と分離する。
- 処理時間だけでなく、入力から表示までのp50/p95、GUI停止時間、親と子の合計RSS、
  転送バイト数、ディスク使用量、worker起動時間を現行版と比較する。
  性能予算は移行前に実測して決める。
- 報告時のHapeville起動条件、周期画像保存有効でPython 3.13/3.14を比較する。
  再インストール時の依存版の差を記録し、まず各8〜12時間を複数回実施する。
  無クラッシュを完全解消の証明とは扱わない。
- 正常終了と強制終了の後に、孤児プロセス・参照中成果物の削除・無制限の
  作業ディレクトリ増加がないことを確認する。

実装が確定した段階で `docs/design.md` に実構成を、`docs/specification.md` に
失敗・stale表示・復旧のユーザー動作を、release-notesに出荷済み変更を反映する。
本資料の未検証案を現行仕様へ転記しない。

## 8. 2026-09-06時点の実施状況

### 8.1 目的と診断上の前提

- 完了: 根本原因を断定せず、プロセス分離を耐障害性向上と原因調査の補助として扱う方針を維持した。
- 完了: 起動時および `zstarview-diagnose-runtime` で、app version/revision、Python、GIL状態、OS、PID、session ID、主要依存版をASCII安全に採取できるようにした。
- 未完了: 実運用中の全依存版・GIL状態・app revisionを長時間実行ごとに収集する運用、core dumpのCスタック解析、Python 3.13/3.14の比較は未実施。

### 8.2 現行構成と移行対象

- 完了: Qt非依存の `zstarview.processes` protocol/supervisorを追加した。
- 完了: cycloneの取得、キャッシュ処理、polygonの辞書化を専用subprocessへ移し、旧GUI内 `_run_update()` 経路を削除した。
- 未完了: sky、雲、地理処理、衛星、航空機、meteor、terrain、urban、water、precipitation、road night lights、検索、hoverの棚卸しと移行は残っている。

### 8.3 プロセス・スレッド構成

- 完了: cycloneはGUIの共有executorではなく、`sys.executable -m zstarview.tropical_cyclones.worker` で起動する。
- 完了: GUIはQTimerからsupervisorをpollし、timeout、terminate、kill、限定再起動を管理する。
- 未完了: 常駐sky/render worker、取得／raster／地理処理の実プロセス枠、QProcessへの移行は未実施。現状はPopen＋QTimer方式である。
- 未完了: GUI内の他のPythonバックグラウンドスレッドと共有executorは残っている。

### 8.4 IPCとデータ所有権

- 完了: versioned JSON request/result、session ID、worker epoch、request ID、layer/view generation、input revisionを実装した。
- 完了: manifestのschema、ID、世代、作業領域内パス、成果物サイズを検証し、JSON payloadをGUIが読み込み後にartifactを解放する。
- 完了: 成果物のatomic rename、件数・容量上限、session単位の回収を実装した。
- 未完了: NPY＋offset配列、完成RGBA画像、カタログ／ephemerisのworker常駐所有はsky等の移行時に実装する。
- 未完了: 既存cloud subprocessのpickle成果物は未置換であり、汎用IPC標準にはしていない。

### 8.5 要求管理、復旧、終了

- 完了: cycloneは実行中1件＋最新待機1件とし、古い結果を要求IDで破棄する。
- 完了: worker exit／timeoutには最大2回のbackoff付き再起動を行い、通常の取得失敗・manifest不正は自動再試行しない。
- 完了: 終了時の新規受付停止、worker停止、artifact解放、最後の正常結果を維持する基本経路を実装した。
- 未完了: プロセスグループ単位の終了、GUI強制終了時の孤児プロセス回収、実OSでの孤児成果物検証は未実施。
- 未完了: 失敗ログの保存量制限、再起動回数の運用調整、全レイヤー共通の手動再試行UIは未実施。

### 8.6 移行手順

| 手順 | 状態 | 備考 |
| --- | --- | --- |
| 1. 診断情報整備 | 部分完了 | runtime diagnosticsは実装済み。core解析・長時間比較は未完了。 |
| 2. IPC schema／process supervisor | 完了 | 異常終了、timeout、壊れた結果、再起動、成果物保持をテスト済み。 |
| 3. cyclone取得／辞書化 | 完了 | subprocess経路へ移行し、旧GUI worker経路を削除済み。 |
| 4. 常駐sky worker | 未着手 | 天体計算、ephemeris、衛星、航空機、meteorを含む。 |
| 5. 空・星・雲画像／地理処理 | 未着手 | 既存cloud pickle adapterも未置換。 |
| 6. 独自thread／検索／hover移行 | 未着手 | 旧executor廃止前の全入口再監査も未実施。 |
| 7. `native_work_lock` 整理 | 未着手 | 旧経路が残るためlockは現状維持。 |

### 8.7 検証・受入条件

- 完了: 全テスト、IPC異常終了／timeout／壊れた結果／再起動／artifact保持、cycloneのQCoreApplication統合経路を検証した。
- 完了: 全体pytestは2026-09-06時点で `1667 passed`、Ruffと差分チェックも通過した。
- 未完了: GUIの実クラッシュ注入、直近正常フレームを保持したままの再起動復旧、起動・リサイズ・時刻変更・終了競合の網羅検証は未完了。
- 未完了: p50/p95、GUI停止時間、RSS、転送量、ディスク使用量、起動時間の現行版比較は未実施。
- 未完了: Hapeville条件でのPython 3.13/3.14各8〜12時間の複数回比較、孤児プロセス・成果物の実機検証は未実施。

### 8.8 文書・出荷扱い

- 完了: `docs/design.md` に現時点のプロセス構成と診断方針を反映した。
- 未完了: `docs/specification.md` と `release-notes.md` への反映は、ユーザー向け出荷範囲と運用検証が確定していないため保留している。
