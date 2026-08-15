# GUI 画面更新とキャッシュ

この文書は、`zstarview` GUI における「いつ画面を更新するか」と「何をキャッシュして再利用するか」をまとめる。

対象は主に `src/zstarview/gui/window.py`、`src/zstarview/gui/window_updates.py`、
`src/zstarview/gui/window_update_*.py`、`src/zstarview/gui/window_render.py`、
`src/zstarview/gui/composite.py` である。
更新ハンドラの公開入口は `window_updates.py` のままにし、HUD 状態、天空データ、
雲、その他オーバーレイの mixin は `window_update_status.py`、
`window_update_sky.py`、`window_update_cloud.py`、`window_update_overlays.py`
に分ける。

## 1. 基本方針

- UI スレッドは入力受付、状態反映、再描画要求だけを担う。
- 重い計算や外部 I/O は worker へ逃がす。
- 画面更新は「状態が変わったので再描画を要求する」と「背景処理の結果が届いたので再描画する」の 2 系統で進む。
- キャッシュは 1 枚の画像だけでなく、レイヤー単位の中間結果にも置く。
- 古い結果は `render_generation` や現在の視点条件で弾き、最新状態だけを残す。
- 更新経路はできるだけ少ない種類にまとめる。
- 可能であれば、1 回の描画に使う sky / cloud / terrain / stars / overlays の入力を同じ世代としてまとめて publish する。
- ただし atomic 化が単純化と衝突する場合は、単純化を優先し、atomic 化は optional の将来改善として扱う。

## 2. 画面更新の流れ

GUI の更新はおおむね次の順で進む。

1. ユーザー操作や定期 tick が入る。
2. `window_updates.py` 経由の各ハンドラが UI state を更新する。
3. 必要なら `_compositor.invalidate()` で合成キャッシュを無効化する。
4. `request_client_update()` で Qt の再描画を予約する。
5. `paintEvent()` が走り、`window_render.py` のキャッシュ付き描画に入る。
6. `QImage` の再利用が可能なら既存キャッシュを返し、条件が変わっていれば再生成する。

このため、画面のちらつきを避けたい場合でも、状態更新後に必ず即時 repaint するのではなく、無効化と再描画要求を分けて扱う。

## 3. どこで更新が起こるか

### 入力イベント

- マウス、キーボード、メニュー、ダイアログ確定は UI スレッドで受ける。
- それらはその場で重い処理をせず、`SkyWindowState` の更新と background task の起動に変換する。

### 背景処理の完了

- sky data、cloud、terrain horizon、water overlay、urban outline、satellite、aircraft、tropical cyclone などの worker 結果は、完了時に UI state へ反映される。
- 反映後は必要なキャッシュを無効化し、`request_client_update()` を呼ぶ。
- 可能なら、これらの結果は個別に即時反映せず、同じ描画世代に属する入力を揃えてからまとめて publish してよい。

### サイズ変更

- client widget の resize では `self._disc_generation` を進め、古い描画結果を捨てる。
- sky disc、cloud image、missing mask などの古い描画資産を明示的に消し、`_compositor.invalidate()` を呼ぶ。
- resize 中は view 更新を先に進めず、新しいサイズで再計算し直す。
- `Square Window` はワンショットのリサイズ処理として扱う。利用者がメニュー項目を選んだ時点で、クライアント領域の短辺に合わせて 1 回だけ正方形へ補正してよい。
- 補正処理の成立条件は、その呼び出しで `square_client_area()` を実行したかどうかだけにする。継続監視や pending 状態は持たない。
- 補正後は通常の resize と同じ更新・再配置経路に戻してよい。

### 視点変更

- viewport interaction mode 中は重い再計算を抑え、fast-mode の軽量描画を使う。
- 交互に発生する視点変更では、古い cloud / satellite / aircraft の投影結果を残し続けない。
- ただし表示上の buffer を一時保持する場合は、`preserve_cloud_buffers` のように明示的に分ける。

### 単純化した更新ポリシー

以下のように更新イベントを絞ると、描画入力の食い違いを減らしやすい。

- 常時周期処理は、700msのポーリングではなく、カレンダー秒境界に同期する1秒tickへまとめる。繰り返しタイマーの遅延を累積させず、各tick後に次の秒境界までのdelayを再計算するsingle-shot方式を基本とする。
- tickの識別には現在のカレンダー秒bucketを使う。イベントループの遅延やサスペンドで秒が飛んだ場合は最新bucketへ直接進み、欠落したtickを再生しない。
- フレームキャッシュにかかる描画リクエストは、ユーザー操作、視線変更、ウィンドウサイズ変更、偶数秒の表示更新、1 分ごとの再描画に限定する。
- 航空機と人工衛星の2秒周期更新は、フレームキャッシュとは別の短命レイヤーとして扱い、偶数秒tickで同時にpublishする。
- 恒星補間では、60秒ごとのsky snapshot更新と、偶数秒tickでの恒星行列更新を分離する。行列更新は短命レイヤーの更新機会を利用するが、sky workerの重いスナップショット計算を再実行してはならない。
- 惑星、太陽、月は次の偶数秒tickを対象時刻として専用workerで先行計算する。結果は対象時刻とgenerationを持つ準備済みフレームとして保持し、対象tickでのみpublishする。空色ディスク、雲の色、夜間光のキャッシュはこの更新で無効化しない。
- 惑星workerの完了signalは結果を保存するだけとし、それ単独ではrepaintを要求しない。計算が対象tickに間に合わない場合は前回の惑星位置を維持し、遅れて届いた期限切れ結果は破棄する。
- 偶数秒tickでは、準備済み惑星位置、航空機、人工衛星、瞬き、恒星補間行列を同じ対象時刻で更新した後、`request_client_update()`を1回だけ呼ぶ。
- 恒星補間の行列が変わった場合は、重い base scene cache 全体を破棄せず、恒星サーフェイスと明るい恒星の表示状態だけを現在の行列で再合成する。
- 検索ターゲット、マウスホバー、ラベル、HUD は 2 秒周期レイヤーとは別の、即応 UI オーバーレイとして扱う。
- これらの UI オーバーレイは、マウス位置、検索ターゲット、選択状態、視線変更に対して直ちに再描画されてよい。
- 視線変更と resize は fast-mode への遷移を伴う。
- fast-modeの復帰はカレンダー秒tickへ統合せず、最後のviewport操作から1000msの独立したsingle-shot settleタイマーで扱う。再びfast-mode遷移が起きたら既存タイマーを破棄して新しいタイマーに差し替える。
- cloud は 10 分ごとに source を再取得し、alt/az の 2D サーフェイスへの変換が終わったら静かに差し替える。
- cloud の差し替えは単独で repaint を起こさず、次の 1 分再描画で自然に使われるようにしてよい。

このポリシーは、画面を更新するタイミングを少数の基準に寄せることで、視線変更と背景更新が競合する余地を減らすためのものである。

## 4. 主要なキャッシュ層

### 4.1 ウィンドウ内のフレームキャッシュ

`src/zstarview/gui/window_render.py` では、次の 4 層を持つ。

- `_frame_cache_key` / `_frame_cache_image`
  - base scene の描画結果を保持する。
  - sky disc、terrain、water、urban outline などの条件が変わると無効になる。
- `_present_frame_cache_key` / `_present_frame_cache_image`
  - base scene に、航空機、人工衛星、台風・サイクロンなどの通常オーバーレイを重ねた表示を保持する。
  - この層は、volatile な HUD / status line / mouse hover / search marker を含めない。
  - ラベル候補はこの段階で集約・保持するが、hover や persistent search などの一時表示は paint 時に重ねる。
- `_display_frame_cache_key` / `_display_frame_cache_image`
  - display tone curveが有効な場合だけ、補正済みのpresent frameを保持する。
  - キーには元のpresent frameの画像識別子、`BLACK,WHITE`、LUT生成規則のversionを
    含める。
  - tone curveが`off`へ変わったら両方を`None`にし、present frameを直接表示する。
  - HUD、status line、mouse hover、search markerはこの画像へ焼き込まない。
- `_fast_frame_base_cache_key` / `_fast_frame_base_cache_image`
  - viewport interaction 中に使う縮小版 base scene を保持する。
- `_fast_frame_cache_key` / `_fast_frame_cache_image`
  - fast-mode の base scene に、軽量ガイドや fast overlay を重ねた表示を保持する。
  - 通常表示と同じく、HUD / status line / mouse hover / search marker は含めない。

これらは `frame_key` が一致する限り再利用される。`frame_key` には次のような条件が入る。

- ウィンドウサイズ
- 視点中心
- 時刻
- テーマ
- 補助レイヤーの ON/OFF
- `QImage.cacheKey()` 相当の画像識別子
- cloud / terrain / water / urban / night light の中間状態

通常のフレームキャッシュキーは、HUD の表示位置、status line の文言、マウス位置、
hover 対象、search highlight のような volatile UI 状態を含めない。これらは
`paintEvent()` の最後に `_draw_volatile_overlay_layers()` で毎回描画する。つまり、
重い base/present frame は安定した描画入力だけで再利用し、即応 UI オーバーレイは
キャッシュを汚さずに現在状態を反映する。

display tone curveが有効な場合は、present frameの生成後に表示用キャッシュを
解決してからウィンドウへ描画し、その後でvolatile UI overlayを直接重ねる。
present frameが変わらず、tone curve設定も同じなら、hoverやstatus lineの更新だけで
全画面LUT変換を再実行しない。tone curveが無効な通常環境では、表示用キャッシュを
割り当てない。

### 4.2 合成キャッシュ

`src/zstarview/gui/composite.py` の `SkyCompositorCache` は、sky image と cloud / terrain / glow の合成をまとめて再利用する。

- `invalidate()` で合成結果と glow mask の両方を捨てる。
- `_glow_mask_cache` は夜間光の通常成分を保持する。
- `_edge_glow_mask_cache` は ridge glow 用の別マスクを保持する。
- `draw()` は合成キー `comp_key` を作り、同じ条件なら前回の合成結果を使う。

ここで重要なのは、合成キャッシュは「1 枚の完成画像」だけでなく、「夜間光のマスク」や「境界付近の補助マスク」も再利用する点である。

### 4.2.1 恒星補間用のレイヤーキャッシュ

恒星補間を導入する場合は、恒星を完成済みの base frame に埋め込んだまま変換しない。
少なくとも次の状態を別々に管理する。

- `dark_star_surface`
  - `vmag > 4.0` の暗い恒星を含む透明サーフェイス。
  - 60 秒スナップショットの生成時に更新し、補間中は画像を再生成しない。
- `bright_star_state`
  - `vmag <= 4.0` の恒星位置、サイズ、色、下敷き、ひし形描画に必要な入力。
  - 現在の補間行列で画面位置へ変換し、合成時に直接描画する。
- `star_interpolation_transform`
  - 60秒区間の中間時刻から現在時刻まで（`-30`〜`+30`秒）の 3x3 行列と、その算出に使った表示範囲・対応点条件。
  - 行列の変更だけでは sky disc、terrain、cloud、urban などの重いキャッシュを無効化しない。

合成順は、空色ディスクの後に、恒星下敷きを `SourceOver` で置き、その後に暗い恒星サーフェイスと
明るい恒星本体・ひし形を適切な合成モードで重ねる。暗色下敷きを加算用画像へ混ぜてから
`Plus` で合成すると黒が暗くならないため、下敷きと発光成分は少なくとも合成モード上は分離する。

恒星補間中の再描画要求では、`request_client_update()` を使って Qt の repaint を予約する。
これは worker の sky snapshot 完了通知とは別の経路であり、Qt のイベントループによる要求の
集約を許容する。現在の `present-frame` cache を再利用する場合は、補間行列またはその世代を
cache key に含め、古い行列で作った表示画像を返さないようにする。

### 4.3 レイヤー状態キャッシュ

各 overlay state は、描画に必要な最新データを保持する。

- `CloudImageState`
  - cloud image、missing mask、alt/az grid、coverage ratio、banner を保持する。
- `TerrainHorizonState`
  - 地形地平線の profile、二次稜線、サンプル列、ground elevation を保持する。
- `WaterOverlayState`
  - sea-level / DEM の点列、現在の active dots、banner を保持する。
- `SatelliteState`
  - 軌道要素由来の records と更新時刻を保持する。
- `AircraftState`
  - 短命な航空機 snapshots と時刻を保持する。
- `TropicalCycloneState`
  - ストーム snapshot 群、更新期限、ソース情報を保持する。

これらの state は、次の描画時に再利用される「入力キャッシュ」であり、毎回ゼロから描画するための中間結果ではない。

### 4.4 原子的 publish を狙う入力束ね

将来的に atomic 化を採る場合は、次の考え方が自然である。

- 1 回の描画で使う入力を `render_generation` あるいは同等の単一世代 ID に束ねる。
- sky disc、cloud projection、terrain horizon、water overlay、urban outline、star table のような base frame 入力は、その世代で揃った時点でまとめて publish する。
- 検索ターゲット、ホバー、ラベルのような即応 UI オーバーレイは、base frame と別の publish 単位としてよい。
- publish 後の repaint は、その世代の入力が揃ったことを前提に 1 回だけ起こす。

ただし、この束ね方が更新ポリシーの単純化を難しくするなら、無理に導入しない。  
その場合は、単純化した更新トリガーを優先し、atomic 化は保留でよい。

## 5. cloud の更新と部分再描画

cloud はこのアプリで最も「2 段階」に分かれた更新対象である。

1. source 取得
2. 現在の視点への projection

`window_updates.py` では、source が届いた時点と projection が届いた時点を別イベントとして扱う。

- source が届いたら `cloud_state.set_source_ready(...)` を呼ぶ。
- projection が届いたら `cloud_state.set_result(...)` を呼ぶ。
- source だけある状態では、再投影のために `_request_cloud_projection_update()` を呼んでよい。
- source ready 自体は repaint の直接トリガーにしない。GUI は source 更新を state と投影要求に反映し、実際の再描画は projection 完了時または別の更新契機に任せてよい。
- すでに別の視点条件に変わっている結果は `render_generation` を見て破棄する。

Geo-satellite も同じ 2 段階を使い、source ready は raw 画像の復元や投影準備に留め、最終的な再描画は ready 側で起こしてよい。

この分離により、通信待ちと投影計算を別のキャッシュ単位で扱える。

## 6. 起動時 warm-up とキャッシュ

起動時は splash 表示中に静的レイヤーを warm-up してよい。

- DEM
- terrain horizon
- night light / ridge glow
- water overlay
- urban outline

ここでの狙いは、最初の本描画時に必要な重いキャッシュを先に埋めておくことにある。

一方で、雲、Geo-satellite 雲、航空機、人工衛星、台風・サイクロンのような動的レイヤーは、GUI 表示後に遅延開始してよい。

## 7. 失敗時の扱い

- 補助レイヤーの失敗で本体の空表示まで消さない。
- 失敗しても stale な cache が読めるなら、そのまま表示に使ってよい。
- 壊れた cache は invalid とみなし、再生成へ回す。
- 画面上の status line は、取得中、投影中、失敗、部分欠損を区別してよい。

## 8. 読み方の目安

この文書だけで全体像を追うより、次の順で読むと把握しやすい。

1. `docs/design/runtime.md`
2. `docs/design/gui-screen-update-and-cache.md`
3. `src/zstarview/gui/window_updates.py`（公開入口とスケジューラ）
4. `src/zstarview/gui/window_update_sky.py` / `window_update_cloud.py` / `window_update_overlays.py` / `window_update_status.py`
5. `src/zstarview/gui/window_render.py`
6. `src/zstarview/gui/composite.py`

## 9. 関連ファイル

- [`src/zstarview/gui/window.py`](../../src/zstarview/gui/window.py)
- [`src/zstarview/gui/window_updates.py`](../../src/zstarview/gui/window_updates.py)
- [`src/zstarview/gui/window_update_status.py`](../../src/zstarview/gui/window_update_status.py)
- [`src/zstarview/gui/window_update_sky.py`](../../src/zstarview/gui/window_update_sky.py)
- [`src/zstarview/gui/window_update_cloud.py`](../../src/zstarview/gui/window_update_cloud.py)
- [`src/zstarview/gui/window_update_overlays.py`](../../src/zstarview/gui/window_update_overlays.py)
- [`src/zstarview/gui/window_render.py`](../../src/zstarview/gui/window_render.py)
- [`src/zstarview/gui/composite.py`](../../src/zstarview/gui/composite.py)
- [`src/zstarview/gui/window_state.py`](../../src/zstarview/gui/window_state.py)
- [`src/zstarview/gui/cloud_state.py`](../../src/zstarview/gui/cloud_state.py)
- [`src/zstarview/gui/terrain_state.py`](../../src/zstarview/gui/terrain_state.py)
- [`src/zstarview/gui/water_overlay_state.py`](../../src/zstarview/gui/water_overlay_state.py)
