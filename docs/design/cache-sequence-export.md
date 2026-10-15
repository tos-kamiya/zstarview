# キャッシュからの毎分画像出力（設計案）

作成: 2026-10-15。`--list-cache-ranges` は実装済み。画像出力・HUD・再開は未実装。
通常GUIの仕様変更は含まない。

## 目的と再現の意味

蓄積済みキャッシュから、指定期間の各分に対応する空を再描画し、PNG列として保存する。
地点、観測高度、画角、方位、画像サイズ、明るさなどは出力時に指定して全期間に適用する。
星・惑星・太陽・月は各フレーム時刻で再計算し、雲・航空機・人工衛星はその時刻までに
取得済みだった履歴を選ぶ。当時のウィンドウ、操作、ロード途中の表示を復元するものではない。

通常GUIの過去時刻表示は最新のデータレイヤーを取得し続ける。その入口や
`zstarview-export-image` を1分ごとに起動する方式では、この目的を満たさない。
出力専用の履歴選択経路を設け、レイヤーのネットワーク取得とライブcontrollerを使わない。
地点解決は通常CLIと同じ入力を受け付け、必要な場合のみNominatim検索を行う。

## 過去の決定との関係

- `dev-notes/plan-24h-replay-cache.md`: 初期案。毎分scene manifest、表示設定の記録、
  `zstarview-export-replay` によるPNG列出力を提案したが、画面記録部分は廃案。
- `dev-notes/plan-historical-layer-replay.md` および `session-2026-10-11.md`:
  表示設定を再生時に指定する方針へ変更し、履歴selectorを実装。その後、通常GUIは
  過去の天体表示中も現在のデータ取得を続ける方針へ変更した。
- `session-2026-10-13.md`: 保持期間は固定24時間から設定可能な1〜720時間へ変更。

今回は専用の連続出力を追加する。毎分の画面記録、session ID、新しい取得用アーカイブは
追加しない。以下のmanifestは出力処理を再開・検証するためのもので、GUI記録ではない。
現行動作の参照元は `cache-history-replay.md` と実装とし、古い計画の完了記載だけで
現在のGUI動作を判断しない。

## CLI案

入口は `zstarview-export-replay`、実体は `cli/export_replay.py` とする。
初期段階から時間間隔は60秒固定とし、動画エンコードは対象外とする。

```sh
zstarview-export-replay \
  -p "Matsue Station" -A 30 -Z S -a 0.6 \
  --height-add-m 1.7 --cloud-opacity 0.5 \
  --image-size 1920,1080 \
  --start 2026-10-14T00:00:00+09:00 \
  --end 2026-10-15T00:00:00+09:00 \
  --output-dir frames
```

- 地点は `zstarview` と同じ `-p / --place QUERY` または位置引数で指定する。
  位置引数は都市、緯度経度、対応する地図URL、タワー・山の名前など、既存resolverの
  入力形式と優先順位に従う。`--place-countrycode` と `--place-lang` も共通化する。
  地点省略時は通常CLIと同じ前回地点を読み、実効地点を表示する。保存設定は変更しない。
- 表示設定はCLIで直接指定する。JSON設定ファイルを必須にせず、初版では
  `--view-config` を設けない。共通parserから対応引数を組み立て、値域・別名・既定値を共有する。
  解決済み地点・高度・設定はmanifestへ固定し、再開時には再検索・再解決しない。

| 用途 | 引数 | 意味 |
| --- | --- | --- |
| 中心高度 | `-A / --view-center-alt` | 通常CLIと同じ高度角 |
| 中心方位 | `-Z / --view-center-az` | 通常CLIと同じ角度・方位名 |
| 航空機 | `-a / --aircraft-opacity` | 航空機の不透明度。高度指定には使わない |
| 建物上から観測 | `-B / --use-building-top` | 建物上端を観測基準とする |
| 追加高さ | `--height-add-m` | 基準面からの追加高さ。`--observer-height-m` も同じ互換別名 |
| 画角 | `--edge-fov-deg`, `--content-fov-deg` | 通常CLIと同じ投影設定 |
| レイヤー | `--cloud-opacity`, `--satellite-opacity` 等 | 対応する各レイヤーの既存opacity引数と短縮名 |
| 画像サイズ | `--image-size WIDTH,HEIGHT` | 単発exportと同じ形式・既定値でピクセル寸法を指定 |
| 保存先 | `-o / --output-dir DIRECTORY` | PNG列とジョブmetadataの出力directory |

`-B` と追加高さの併用、タワー名指定の高さ計算も通常CLIと同じ意味にする。
建物上端や地面標高はローカルデータから解決し、不足時にネットワーク取得せずエラーとする。
GUIの位置座標や `restore` を含む `--window-geometry` は初版では受け付けず、
サイズは `--image-size` に統一する。背景・雲・地形・水面・都市・夜間光などの
対応するopacityと、表示モード・描画設定も共通の引数定義を使う。
未対応の動的レイヤーは明示的な正のopacity指定をエラーとし、未指定時は0へ解決して
helpに示す。それ以外の対応設定の既定値は通常CLI/単発exportと揃える。

- `--start` / `--end`: タイムゾーン付きISO 8601。秒・小数秒は0を要求し、暗黙に丸めない。
  内部はUTC、区間は開始を含み終了を含まない。上記は1440枚。終了が開始以下ならエラー。
- `--cache-dir`: 履歴と静的データの入力root。省略時は通常のcache root。
  入力rootと出力directoryは分離し、出力を入力root配下に置くことを拒否する。
- `--plan-only`: データ選択と不足一覧の出力まで行い、描画しない。
- `--allow-partial-data`: レイヤー不足を許容して全時刻を出力する。既定は不足時に事前検査で停止。
- `--resume`: 確定済みplanから未完了画像だけを生成する。既存画像の無条件上書きは行わない。

端末出力とエラー文はASCIIとする。既存の単発exportの全オプションをそのまま
転送する方式は採らず、オフラインで成立する表示設定だけを受け付ける。

## 地点ごとのキャッシュ時刻範囲

この節の範囲照会は実装済み。ユーザーの指定によりメタデータのみを確認する。
雲・画像・静的source・暦の本体は読み込まず、内容検証も行わない。
`common_ranges` は動的レイヤーのメタデータ上の共通候補区間とし、暦の対応期間を
検証しないため `ephemeris_ranges` / `celestial_common_ranges` は `null` とする。
`verification_level=metadata-only`、`payloads_verified=false` を返し、描画成功は保証しない。

同じコマンドの `--list-cache-ranges` で確認する。描画、Qt初期化、出力directory作成は
行わず、`--start` / `--end` / `--output-dir` は不要とする。開始・終了指定があれば
その区間へ一覧を絞り込む。`--json` は機械処理用出力を選ぶ。

```sh
zstarview-export-replay -p "Matsue Station" --list-cache-ranges
zstarview-export-replay matsue --list-cache-ranges --json
zstarview-export-replay '35.46,133.06' -a 0 --list-cache-ranges
```

一覧は解決済み座標と地点名を表示し、地点に合う履歴が残っている各動的レイヤーについて、
件数、source時刻の最小・最大、取得時刻の最小・最大、毎分の再生可能区間を出す。
opacityが0のレイヤーも保存状況は表示するが、再生可能な共通区間の計算からは除外する。
該当履歴がなければ明示し、単なる最古・最新の間を連続再生可能とは表示しない。

再生可能区間は画像出力と同じselectorで判定する。UTCの分境界において、有効期限、
取得時刻、source時刻、scopeを満たす候補がある分を連結し、半開区間 `[start, end)` と
欠損区間を表示する。たとえば10:00、10:01だけ再生可能なら `[10:00, 10:02)` とする。
航空機ゼロ機の有効snapshotも可能な分に数える。レコードの有効期限で最終source時刻より
後まで再生できる場合があるため、source時刻範囲と再生可能区間は別々に示す。

全有効レイヤーの区間の積集合を「共通の再生可能区間」とし、必要な衛星groupもすべて
満たすことを条件にする。ローカル暦・静的データの内容は検証しない。
`-B` の建物上端をメタデータから確定できない場合は高さを未確認とし、scopeを推測しない。
完全出力可能とは表示しない。
動的レイヤーをすべて無効にした場合は、共通のキャッシュ区間は適用外と表示する。
空の積集合は、異なるレイヤーの取得期間が重ならない場合と履歴自体がない場合を区別する。

native雲は既存のobserver座標許容差、航空機はbbox、Horizonsはobserver scopeで
一致判定する。Geo-satelliteなど地点共通の入力も、その地点がsourceの対応範囲か検証する。
地点文字列の一致や最寄り地点の起動metadataだけで履歴を選ばない。
`-B` や高さによってsource scopeが変わる場合は、その解決済みscopeも適用する。

一覧は履歴schema、時刻、scope、参照pathとファイルの存在に基づく。digestは記録された
値の形式だけを確認し、本体のハッシュは計算しない。航空機・人工衛星では既存JSONに
含まれる時刻・bbox・groupなどを参照し、描画用の要素オブジェクトは生成しない。
表示はUTCのoffset付き時刻とし、JSONはレイヤー別の状態、区間配列、欠損理由、確認時刻を持つ。
本体が破損していても候補に含まれ得るため、実際の出力開始時に完全検証する。
該当なしの照会自体は正常終了0とし、不正入力・地点解決失敗・走査失敗は1とする。

`-p` の地点解決は既存place cache、同じquery条件の起動metadata、通常の検索の順で
再利用できるようにする。起動metadataが複数の異なる座標を指す場合は自動選択せず、
通常resolverへ戻す。外部検索を行うのは地点解決だけで、範囲照会や画像出力の
気象・航空機・人工衛星データを新たに取得しない。解決済み座標を指定すれば検索は不要。

## レイヤー選択

フレーム時刻を `t` とし、動的データは原則として取得時刻とsource時刻の両方が
`t` 以下のものに限定する。未来の取得結果を遡って使わない。
有効期限は `t` から測り、実行日の現在時刻や保持期間によって候補を除外しない。
保持期間は削除方針であり、残存データの再生可能期間とは分離する。

| レイヤー | 入力と選択 | 初版での扱い |
| --- | --- | --- |
| 天体 | ローカルcatalogと暦、フレーム時刻 | 毎分再計算。暦不足・範囲外は事前エラー |
| GOES/Himawari雲 | `source_history` の完成済みCloudSourceDataとsidecar | 観測地点一致。取得時刻順に切替え、次の取得成功まで保持。末尾は最後の取得から30分 |
| Geo-satellite雲 | `raw/history` のPNGとsidecar | 種別・area・size一致。取得時刻順に切替え、次の取得成功まで保持。末尾は最後の取得から60分 |
| 航空機 | bbox一致の履歴snapshot | 最新取得snapshot。各機のlast_contactが未来でなく10分以内 |
| 人工衛星 | group別要素、必要な場合observer scope | 取得・epochが未来でなくgroup有効期間内。epoch、取得時刻の順で最新を選び、tへ伝播 |
| 地形・都市・水面・背景画像 | ローカルの対応するsource/派生cache | 地点・範囲・schemaを検証して使う。必要な描画変換は出力側で実行 |
| 降水・台風等 | 毎分再現に使える履歴契約が未整備 | 初版は非対応。有効指定はエラーとし、現在値を代用しない |

雲は設定でnativeまたはGeo-satelliteを明示し、欠損時に別providerへ自動切替しない。
同時刻の候補は取得時刻と安定したrecord ID/pathで順序を確定する。
雲では最新の取得結果を優先し、後から取得した観測時刻の古いsourceもGUIと同じく
その取得時点で切り替える。観測時刻が取得時刻より未来の不整合recordは除外する。
後続の取得までの保持は観測鮮度の保証ではなく、既存画面の継続表示を近似する規則である。
起動・終了を記録していない履歴からGUIの停止期間は識別しない。
通常GUIや既存の履歴readerの鮮度条件は変更せず、replayのselectorに限定して適用する。
sourceが同じ間は雲分布を保持するが、太陽位置や照明はフレーム時刻で計算する。
航空機の位置更新・外挿は既存描画のモデルに従い、実際の毎分観測があるとは表現しない。
Horizonsなど時刻固定の記録は、要素として任意時刻へ伝播できると仮定しない。
対応時刻・observer scopeを満たさない場合は不足にする。

航空機ゼロ機の成功snapshotは `empty` として有効に扱う。現行の
`load_aircraft_cache_for_time` は空のsnapshotを候補にしないため、そのまま再利用すると
「飛行機なし」と「データなし」を混同する。専用adapterで区別する。

## 構成と処理順

`replay/` にQt非依存の `models.py`、`inventory.py`、`selection.py`、
`plan.py` を追加し、レイヤー固有の読込adapterを分ける。既存の起動metadata処理と共存する。
`cli/export_replay.py` はCLI、Qt初期化、連続描画と進捗表示を担当する。

1. 共通CLI引数を検証して地点を解決し、表示設定と対象時刻列を確定する。
2. 履歴を一度走査して、時刻・scope・provider別に索引化する。毎分の再帰走査はしない。
3. 各時刻のレコード参照を選び、静的依存を含めたplanと不足一覧を作る。
4. 選択ファイルを読込・digest検証し、出力directory内の `inputs/` へ一時保存して
   atomic renameする。内容digestで重複排除する。以後はこの固定入力だけを使う。
5. planを確定し、単一のQt application、catalog、静的sceneを初期化する。
6. 時刻順に天体計算、履歴読込、時刻依存の照明・投影、scene構築、PNG保存を行う。
7. 各PNGと処理結果を確定し、完了・不足・失敗件数と終了コードを出す。

固定入力は新たな収集用履歴ではなく、出力ジョブの依存物である。長時間出力中に
通常GUIのcleanupが元cacheを消しても継続できるようにする。hard linkに依存せず、
読みながらdigestを計算してコピーする。収集中に元ファイルが消えた場合は、その参照を
不足として扱い、別の履歴を勝手に選び直さない。静的依存も固定対象に含める。
必要容量を事前に見積もり、保存中の容量不足はジョブを停止する。

## HUDと操作用UI

画像出力は `zstarview-export-image` に揃え、既定ではHUDとステータスを表示しない。
ハンバーガーメニュー、サイズ変更グリップ、ウィンドウ枠、マウスやhover表示も出力しない。
これらの操作用UIは `--hud all` にも含めない。

`--hud ITEMS` はカンマ区切りで表示項目を指定する。既定は `none`。

| 項目 | 表示内容 |
| --- | --- |
| `datetime` | フレームの日時と表示タイムゾーン |
| `location` | 観測地点名 |
| `coordinates` | 観測地点の緯度・経度 |
| `height` | 地面標高、観測基準の高さ、追加高さ（既存HUDの高さ表示に準拠） |
| `view-direction` | 視線の高度角・方位角 |
| `time-of-day` | 昼夜を示す既存マーカー |
| `none` | 全項目を非表示 |
| `all` | 上記の表示項目すべてを表示 |

限界等級は現行HUDで表示されていないため、`magnitude-limit` は設けない。
日時表示は `--hud datetime` に統一し、`--include-timestamp` は設けない。
未知の項目、空の項目、`none` / `all` と他項目の混在は入力エラーとする。
項目の列挙順は配置順に影響させず、既存HUDの情報順を基本とした固定順で表示し、
非表示項目のスペースは詰める。重複項目は一つにまとめる。

配置は既存CLIの `--observation-info` に揃え、`top`（左上）と `bottom`（左下）のみ
受け付ける。既定は通常CLIと同じ `bottom`。replayでは `auto` と `off` は入力エラーとし、
helpにも許容値を `top|bottom` と表示する。表示・非表示は `--hud` に一本化する。
`--hud-position` は設けない。`time-of-day` も選択した位置のHUD領域内で文字と
重ならない位置に置く。マウス位置による自動移動やフェードは行わない。
`datetime` と `time-of-day` は実行時刻ではなく各フレームの対象時刻から生成する。

```sh
zstarview-export-replay matsue \
  --start 2026-10-14T00:00:00+09:00 \
  --end 2026-10-15T00:00:00+09:00 \
  --hud datetime,location --observation-info bottom \
  --image-size 1920,1080 -o frames
```

方位グリッドは既存と同じ `--include-direction-grid` で独立に指定する。
天体名ラベルも既存の描画設定に従い、`--hud none` で消さない。
解決済みのHUD項目と位置もmanifestの表示設定に保存し、再開時に維持する。

## 描画の再利用と実行モデル

描画は既存の `FrameContext` / `RenderSceneData` / `RenderStyle` と
`export_image_render.py` の共通パイプラインを使う。画面枠、マウス、メニュー、
一時的なstatusやUIフェードは対象外。HUDは上記の明示指定に従って合成する。
毎フレームsceneを明示的に組み立て、前フレームのレイヤーが欠損時に残らないようにする。

Qt描画は主スレッドで行い、初版はフレームを直列処理する。Astropy/Skyfield/NumPyの
native workは `ApplicationServices.native_work_lock` で直列化する。
雲workerを使う場合も前の処理の完了を待つ。照明画像はsource keyだけでは再利用せず、
時刻・太陽方向・地点・画角・解像度・表示設定を含めて管理する。

## 既存コードの再利用時に必要な分離

- `clouddisc/history_cache.py`: metadataとartifactの検証は再利用可能。
  現行loaderは各呼出で全件を走査し、参照pathや選択理由を返さないため、索引と
  provenanceを返す純粋な選択処理を抽出する。pickleはアプリ自身の信頼済みcacheに限定し、
  digest検証を任意pickleの安全性保証とは扱わない。型/schema不互換は明示する。
- `geosatellite/cache.py`: 現行の履歴readerには現在時刻による保持期限判定と
  legacy履歴の書込がある。専用の読取専用adapterで取得時刻も検査する。
  同じ撮影slotへの再取得でsidecarが更新される可能性があり、最初の取得時刻を
  復元できない場合は現存sidecarの取得時刻を採用する。
- `geosatellite/pipeline.py`: `historical=True` や `raw_png` 指定でも入口で
  purge/cleanupが走る。投影・雲source生成を副作用のない関数へ分離して使う。
- `aircraft/cache.py` / `satellites/cache.py`: decoderとscope/有効期限規則を再利用し、
  fetchやライブの時刻判定を経由しないadapterを置く。
- `cli/export_image.py` / `export_image_inputs.py` / `export_image_layers.py`:
  mainやfetch付き初期化は呼ばない。catalog/設定の準備、暦検証、静的データ読込、
  scene構築から必要な部分だけを共通化し、変更範囲を小さく保つ。

## 出力、欠損、再開

```text
frames/
  manifest.json
  inputs/<sha256>...
  frames/20261013T150000Z.png
  frames/20261013T150000Z.json
  frames/20261013T150100Z.png
  frames/20261013T150100Z.json
```

manifestにはschema、アプリ/renderer版、正規化した設定とdigest、対象時刻列、
フレームごとの入力参照を保存する。入力参照は相対path、digest、source時刻、取得時刻、
scope、選択規則を含む。各フレームJSONは画像digest、plan digest、各レイヤーの
`ok / empty / disabled / missing / corrupt / incompatible / unsupported` と理由を持つ。
元cacheの絶対pathを再開の必須条件にしない。

既定では有効レイヤーの不足を事前報告してPNG生成前に停止する。
`--allow-partial-data` の場合は不足レイヤーを空にして毎分のPNGを作り、JSONで明示する。
描画自体の例外は成功扱いにせず、そのフレームを失敗として残す。
終了コードは0=全対象成功、2=不足を含む完了、1=入力エラー/事前検査不合格/処理失敗。

PNGは一時ファイルからrenameし、フレームJSONを最後に確定する。
再開時はmanifestと固定入力、PNGのdigestが一致する完成フレームだけをスキップする。
PNGだけ残った場合やdigest不一致は黙って上書きせず、不整合を報告する。
初版は同一アプリ/renderer版での再開を要求し、異なる環境でのpixel完全一致は保証しない。
既存出力directoryは `--resume` なしでは拒否する。

## 実装順序と検証

1. 共通CLI引数、時刻列、解決済み表示設定の型、索引と `--list-cache-ranges` を実装する。
2. 雲・航空機・人工衛星の読取専用adapterと静的依存の準備を実装する。
3. 固定入力、単一フレームの共通描画、毎分ループとmanifestを実装する。
4. 部分出力と再開を実装し、小さい実cache期間で実地確認する。
5. 実装後に利用者向け仕様、コマンド例、release notesを更新する。

テストはUTC境界・1440枚の時間列、未来の取得排除、scope不一致、鮮度境界、
空snapshot、壊れたdigest、非対応schema、欠損時の残像防止、再開時の整合性を対象とする。
CLI互換性は地点入力各種、`-A` / `-Z` / `-a`、`-B` と追加高さ、opacity、画像寸法を確認する。
HUDは既定で非表示、各項目の独立表示、all/noneと不正指定、上下の配置、auto/offの拒否、対象時刻の使用、
非表示項目の余白除去、操作用UIの非出力を確認する。
範囲照会は途中欠損、積集合が空、scope不一致、全無効、時刻境界と実際のplanの一致を検証する。
描画テストはQt offscreenで少数フレームに限定する。地点解決後のネットワーク呼出と元cacheへの
書込/削除を失敗させた状態で全経路が完走すること、固定後に元cacheがなくても
継続できることを確認する。雲のsourceが変わらない隣接分でも天体と照明が更新されることを確認する。
24時間分の性能・出力容量は実測して記録し、設計段階では所要時間を保証しない。
