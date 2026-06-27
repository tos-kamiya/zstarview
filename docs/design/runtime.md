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
- ridge glow は、まず粗い全域モデルを作り、その後に主稜線近傍だけ高解像度の局所モデルを重ねる二段構成としてよい。局所モデルの外側は欠損点を外挿で埋めず、粗い全域モデルを背景として使ってよい。粗密の境界は滑らかにつないで、主稜線付近の段差を抑えてよい。描画段階では terrain mask を適用して見える部分だけを残してよい。
- この段階では、既存 cache の読み込み、欠損時のダウンロード、派生 cache の生成、初回 frame のための描画準備までを行ってよい。
- 雲、Geo-satellite 雲、航空機、人工衛星、台風・サイクロンのような時間依存レイヤーは、通常 GUI の表示後に別 task として開始してよい。
- 雲と Geo-satellite については、既存の raw 画像 cache があれば再ダウンロードを省略してよいが、最終的な repaint は投影完了側の結果で発生させてよい。
- startup log は、この warm-up の進行と失敗を利用者へ見せるための可視化手段として扱ってよい。

### `zstarview-export-image` の取得順

- `zstarview-export-image` は単発処理なので、GUI の splash warm-up と同じ順序に固定しなくてよい。
- 1 枚の画像の鮮度を上げたい場合は、寿命の短いレイヤーを先に開始してよい。特に aircraft と satellites は、cloud より先に取得を開始してよい。
- cloud は GUI と同様に独立した取得経路として扱ってよい。ダウンロードだけを先行させ、取得後に (altitude, azimuth) グリッドへの取り込みを行い、最終的なカメラ依存の描画は後段で行ってよい。
- その後に、準静的な terrain、urban outline、water、night light / ridge glow をまとめて取得してよい。
- 重要なのは固定順序よりも、単一の基準時刻を使って最終描画の整合性を保つことと、依存関係のある結果だけを最後の合成時点で揃えることである。

### PNG メタデータの埋め込み

- `zstarview-export-image` が PNG に保存する場合は、書き出し直前に画像メタデータをまとめて付与してよい。
- メタデータは、描画済み画像から逆算するのではなく、同じ `viewer_data` / `celestial_data` / `style` / `scene` から生成する。
- 少なくともアプリのバージョン番号は別キーで保持し、HUD 相当の表示情報は機械可読な構造体として 1 つにまとめてよい。
- HUD 相当情報の基底は、`render.background.format_overlay_info_lines()` が返す固定表示行とする。ここには地点、時刻、視線方向、`vmag` 上限を含めてよい。
- export-image は headless なので、マウスホバーに依存する一時的な HUD 情報は省いてよい。
- `--place` は、ユーザー入力の生クエリと正規化後の解決結果を別フィールドで保持してよい。解決結果には、表示名、座標、タイムゾーン、ソース名などを含めてよい。
- `--search` の選択結果や cloud coverage ratio のような export-image 固有の補助情報は、HUD 相当の構造体に追加フィールドとして載せてよい。
- PNG への実際の書き込みは Qt の text chunk 機能で行ってよい。キー名は将来の互換性を考えて `zstarview.*` の名前空間に寄せる。
- stdout 出力や SIXEL 出力は同じ内部メタデータ生成経路を共有してよいが、外部に現れるのは画像バイト列のみである。

#### 正規形式

PNG text chunk には、正規フォーマット `zstarview.export-image-metadata.v1` を 1 つの JSON payload として入れる。

```json
{
  "schema": "zstarview.export-image-metadata.v1",
  "version": "1.32.6",
  "hud": {
    "lines": [
      "Matsue, Shimane, Japan",
      "2026-06-25 02:50:00 JST",
      "Alt 35 deg  Az 120 deg (ESE)",
      "Vmag limit 6.0"
    ],
    "view": {
      "city_name": "Matsue",
      "lat_deg": 35.47,
      "lon_deg": 133.05,
      "view_center_alt_deg": 35.0,
      "view_center_az_deg": 120.0,
      "vmag_limit": 6.0
    }
  },
  "place": {
    "query": "Matsue",
    "resolved": {
      "display_name": "Matsue, Shimane, Japan",
      "lat_deg": 35.47,
      "lon_deg": 133.05,
      "timezone_name": "Asia/Tokyo",
      "source": "nominatim"
    }
  },
  "extra": {
    "search_target": null,
    "cloud_coverage_ratio": 0.77
  }
}
```

- `schema` は固定値 `zstarview.export-image-metadata.v1` とする。
- `version` はアプリ本体の `__version__` を入れる。
- `hud.lines` は画面上に出る固定テキスト行そのものを保持する。
- `hud.view` は再利用しやすい正規化済み数値を保持する。
- `place.query` は `--place` の生入力であり、`place.resolved` は解決後の正規化結果である。
- `extra` は export-image 固有の任意拡張領域として使う。
- 実装は `iTXt` に JSON を 1 本入れることを正規の書式とし、複数 chunk に分割する実装は許容しない。

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
