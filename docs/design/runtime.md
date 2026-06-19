# ランタイム、フロー、キャッシュ

この文書は、3 つのアプリケーション入口で共通する実行時の振る舞いをまとめる。

`zstarview-debug` は `zstarview` と同じ GUI 起動フローを使うが、console script 版として配布されるため、Windows では terminal から起動して startup log を追いやすい。

## スレッド

- GUI 更新は Qt の UI スレッドで行う。
- 長時間の I/O と重い計算はバックグラウンド worker に逃がす。
- 視点変更や再取得は latest-request-wins で扱う。
- worker が busy の間に古い周期 tick を溜め込まない。

## GUI 状態更新

この節では、GUI のイベント、タスク、キュー、worker スレッドが、どの順で状態更新に結び付くかを定義する。

### 入力イベント

- マウス、キーボード、メニュー、ダイアログ確定などの入力は UI スレッドで受ける。
- 入力イベントは、その場で重い計算をせず、更新要求または操作コマンドへ正規化する。
- 視点変更、トグル切替、検索確定、再取得要求は、まず UI state を更新するための意図として扱う。

### task と queue

- UI スレッドは、更新要求を controller 単位の task に変換して worker pool に送る。
- 共有 worker pool は、同時実行を避けるため基本的に 1 件ずつ task を処理する。
- 各 task は request id か同等の識別子を持ち、古い task の結果で UI state を巻き戻さない。
- 定期更新の tick は、worker busy 中は溜め込まず、次回 idle で 1 件だけ進める。

### worker スレッド

- worker は、天体計算、雲投影、地形地平線、都市アウトライン、航空機、人工衛星、検索候補の解決などの重い処理を担う。
- worker は UI widget を直接触らず、結果データだけを返す。
- worker の結果は、完了時に UI スレッドへ signal または同等の queued handoff で戻す。

### UI state の反映

- UI スレッドは、worker の結果を受けて `SkyWindowState` や各 overlay state を更新する。
- 更新後は必要な frame cache を無効化し、再描画要求を出す。
- `latest-request-wins` の条件により、古い結果は state へ反映せず破棄してよい。

### shutdown と cancellation

- window teardown 開始後は、新規 task の投入を止める。
- 進行中 task は可能ならキャンセルし、破棄済み UI への更新通知を避ける。
- 終了処理は、worker の完了通知よりも UI の安全性を優先する。

## 処理フロー

- 起動フロー: 引数解析、地点解決、時刻正規化、保存状態復元を行い、その後に GUI かレンダリング系へ進む。
- fast-mode フロー: 視点変更のたびに軽い再投影を行い、応答性を優先する。
- normal-mode フロー: 時刻、キャッシュ状態、オーバーレイ準備状況が変わったら全体を再計算する。
- export フロー: 必要なレイヤーを待ち、1 回だけ描画して最終画像を書き出す。

### 起動時のスプラッシュ warm-up

- 観測地点と観測高度が確定した後、通常 GUI を表示する前に、静的レイヤーの warm-up をスプラッシュ表示中へ寄せてよい。
- warm-up の対象は、DEM、terrain horizon、night light の alpha grid と ridge glow、水面、都市アウトラインなど、起動後に変化しない補助レイヤーである。
- この段階では、既存 cache の読み込み、欠損時のダウンロード、派生 cache の生成、初回 frame のための描画準備までを行ってよい。
- 雲、Geo-satellite 雲、航空機、人工衛星、台風・サイクロンのような時間依存レイヤーは、通常 GUI の表示後に別 task として開始してよい。
- startup log は、この warm-up の進行と失敗を利用者へ見せるための可視化手段として扱ってよい。

## エラーと復帰

- 失敗はできるだけ利用者に見える状態で返す。
- 補助オーバーレイの失敗で、本体の星空まで消さない。
- 起動失敗と、後段の補助レイヤー失敗を分けて扱う。
- download / projection / partial data / unavailable の状態文言を表示に出す。

## ログと診断出力

- ログは診断用であり、UI 表示そのものとは分けて扱う。
- terminal、console、log、CLI help、exception message、subprocess の stdout/stderr に出る文字列は ASCII-only を原則とする。
- UI 専用文字列は、画面表示のために必要な場合だけ非 ASCII を許容する。
- 外部入力や外部データに含まれる非 ASCII 文字は、必要なら escape するか要約してログに残す。
- ユーザー向けの日本語説明は、ログではなく UI や仕様文書へ出す。

## 設定の保存と復元

- `zstarview` と `zstarview-export-image` は共有設定を `~/.config/zstarview/config.json` から読む。
- `zstarview-gui` は起動前ダイアログ用の launch profile を別ファイルとして読み書きする。
- `zstarview` と `zstarview-export-image` は、GUI 専用の launch profile を既定値として上書きしない。
- 起動時は、CLI 引数、保存済み設定、既定値の順で初期状態を決める。
- 保存は、UI の確定タイミングか終了時に行い、途中の transient state は残さない。
- 保存対象は、前回地点、window geometry、dialog 初期値、必要なら structured place payload などの永続化に適した値に限る。

## キャッシュ方針

- キャッシュは寿命で分けて扱う。
- 短命キャッシュは、航空機、人工衛星、短時間の検索結果、直近の再描画結果などに使う。
- 準静的キャッシュは、地形、建物、DEM、夜間光などのように、データ元で更新されるが数日から数週間は再利用してよいものに使う。
- 長寿命キャッシュは、更新頻度が低く、再取得のタイミングを明示的に管理したい外部データに使う。
- 期限判定の一次情報は mtime ではなく、メタデータに寄せる。
- 失敗しても読める stale cache は、可能なら使い続ける。
- 壊れた cache は stale ではなく invalid として破棄してよい。
- 準静的キャッシュは TTL 超過でも即削除せず、再取得成功までは古い結果をフォールバック利用してよい。
- 短命キャッシュは、鮮度が落ちたら次回更新時に再取得候補へ回す。

## 状態との関係

状態の定義と構造は [data-model.md](data-model.md) に分離している。
runtime 側では、各状態がどのタイミングで更新されるか、どのように復帰するかだけを扱う。
