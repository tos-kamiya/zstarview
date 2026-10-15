# 実装プラン: replayのキャッシュ範囲照会

作成: 2026-10-15
状態: 範囲照会を実装済み。性能問題の修正としてメタデータのみを確認する方式へ変更。
以下の初期プランの本体digest・型・軌道要素・静的source・暦の内容検証は照会では実行しない。
これらは画像生成前の完全検証へ延期し、一覧の候補区間と本体検証済みの区間を区別する。
雲の選択規則も更新済み: 取得時刻順で次の取得成功まで保持し、末尾のみ最後の取得から
native 30分/Geo-satellite 60分で区切る。下記初期プランのsource基準の雲有効期限は廃止。
検証: 全体1818件、最終変更後の関連43件が成功。新規コードのmypy、Ruff、diff checkが成功。
設計: `docs/design/cache-sequence-export.md`

## 今回の到達点

`zstarview-export-replay LOCATION --list-cache-ranges` を実用可能にする。
`-p`、位置引数、必要な観測設定、レイヤーopacity、任意の期間絞り込み、JSON出力を
受け付け、その地点の履歴と連続再生可能区間を示す。
正式オプション名は複数形 `--list-cache-ranges` とする。

画像生成、HUD描画、固定入力コピー、manifest、再開、動画出力は次段階。
大規模な事前リファクタは行わず、照会に必要な処理を追加し、副作用を避けるために
必要な箇所だけ既存コードから分離する。通常GUIと単発exportの動作を維持する。

## 1. 照会の契約とデータ型

新規 `replay/models.py` に以下のQt非依存型を置く。

- `ResolvedReplayLocation`: 座標、timezone、地面標高、観測基準高、追加高さ、解決元。
- `HistoryRecord`: レイヤー、provider/group、scope、source/取得時刻、参照path、
  schemaとdigest、読込検証状態。レイヤー固有payloadはadapterが扱う。
- `LayerSelection`: 選択record、ok/empty/missing/corrupt/incompatible/unsupported、理由。
- `MinuteRange`: UTC分境界の半開区間。
- `CacheRangeReport`: schema version、検証時刻、地点、実効設定、レイヤー別件数・
  source/取得時刻範囲・再生可能区間、共通区間、静的依存の状態、診断件数。

一覧の「件数」は検証済みかつscope一致のレコード数とする。破損・互換性不一致・
scope不一致の件数は別に出し、検証済み件数へ混ぜない。
時刻はtimezone付き入力からUTCへ正規化する。分境界の時刻列と区間連結を純粋関数にする。

完了条件: JSONの構造と状態の意味が確定し、後の画像生成も同じ選択結果を利用できる。

## 2. 地点と観測scopeの解決

`replay/location.py` を追加し、通常のresolverを利用する薄い入口を設ける。
位置引数の都市/座標/URL/山/タワー、`-p` と検索条件、前回地点を既存規則で解釈する。

- `resolve_launch_location` の保存callbackを無操作にし、照会による前回地点更新を防ぐ。
- place cacheと同じ検索条件のlaunch metadataを検査する。複数座標の候補は
  勝手に選ばず通常resolverへ戻す。検索は地点解決に限り許容する。
- 標高・`-B` は取得を伴う経路を使わず、ローカル読込で解決する。
  `location_resolver/resolve.py` の標高・建物上端処理を調べ、必要なら
  読取専用依存を注入できる小さい境界を追加する。GUI側の既定は変更しない。
- `--cache-dir` を指定したら、scope解決に使うDEM/建物/launch履歴もそのrootから読む。
  内部のglobal pathを一時書換えする方式は使わない。
- 高さ不足は地点の平面座標まで失敗させず、観測条件の不足として報告する。
  高さが必要なscopeだけを未確定にし、完全再生可能とは表示しない。

完了条件: 地点名と座標を通常CLIと同じ意味で解決でき、レイヤー取得と設定保存が走らない。

## 3. 履歴inventoryと読取専用adapter

`replay/inventory.py` と `replay/adapters/` を追加する。各rootを一度走査し、
scope/groupと時刻で索引化する。既存decoder・schema・定数を可能な範囲で再利用する。

| adapter | 元の実装 | 照会に必要な処理 |
| --- | --- | --- |
| native雲 | `clouddisc/history_cache.py` | sidecar、path包含、digest、型とobserver検証。source_historyを対象とする |
| Geo-satellite | `geosatellite/cache.py` | PNG/sidecarの検証とarea/size/kind一致。取得時刻も確認 |
| 航空機 | `aircraft/cache.py` | bbox一致、空の成功snapshotとmissの区別、各機の時刻検査 |
| 人工衛星 | `satellites/cache.py` | group/observer scope、epoch、取得時刻、group別有効期限の検査 |

履歴を優先し、互換latestを使う場合も同じ時刻・scope検証と重複排除を行う。
取得失敗だけのmetadataや部分ファイルは成功レコードに数えない。
現在時刻によるretention判定、cleanup、legacy archive書込、fetchを呼ばない。
Geo-satelliteの既存reader/pipelineをそのまま呼ばない。
artifact検証は一つにつき一度にまとめ、巨大payloadを全件メモリに保持しない。
native雲pickleは信頼済みアプリcacheのみを対象とし、digestと型・schemaを検査する。

完了条件: cacheを読取専用にしても照会でき、削除・移行・ネットワーク取得が一切走らない。

## 4. selectorと区間計算

`replay/selection.py` にフレーム時刻を引数とするselectorを追加し、
`replay/coverage.py` で毎分の結果を区間へまとめる。画像生成も同じselectorを使う。

- source時刻と取得時刻が対象時刻以下であることを要求する。
- native雲30分、Geo-satellite60分、航空機last_contact10分、衛星group別有効期限を
  既存規則に沿って適用する。境界は各定数と現行判定を確認して明示テストする。
- 航空機は最新取得snapshotを選び、全機が期限外ならmissing、元から空の成功ならempty。
  古いsnapshotの機体を継ぎ足さない。
- Horizons等は任意時刻へ伝播可能な要素と区別し、実際の記録時刻の適用条件を調べて
  明示する。条件が満たせない分はmissingとし、時刻固定記録から連続区間を捏造しない。
- 同順位候補は取得時刻と安定したID/pathで順序を確定する。
- opacity 0のレイヤーも保存状況を列挙し、共通区間からのみ除外する。
- 共通区間は全有効レイヤーと必要な全衛星groupの積集合。全無効は適用外とする。

期間未指定時はscope一致レコードの最も早い適用開始から最も遅い有限の有効期限までを
候補範囲にする。保持期間や実行時刻を上限にしない。--start/--endがあれば交差させる。
source時刻が未来でなくても取得時刻以前には適用しないため、開始は両者を考慮する。
有効分10:00と10:01は `[10:00,10:02)` にまとめる。欠損をまたいで連結しない。
毎分全履歴を走査せず、索引と時刻cursor/二分探索で候補を絞る。

完了条件: 最小/最大source時刻と再生可能区間が区別され、欠損を含む履歴でも正しく列挙できる。

## 5. 静的依存と暦の事前検査

`replay/dependencies.py` で有効な地形・都市・水面・背景などに必要なローカル入力を調べる。
取得付き `export_image_layers.py` 関数を呼ばず、既存のcache検証関数を再利用する。
sourceがあり派生生成で対応可能な場合と、source自体の欠損を区別する。
照会では派生生成を行わない。必要な派生生成は将来の出力側の責務とする。

暦はローカルファイルから読み、対象天体に必要なcoverageを検査する。
`astro.load_ephemeris()` はSkyfield loader経由のため、単に呼ぶだけでオフラインと
仮定しない。ファイル指定の読込を使い、fallback downloadを禁止する。
native処理は `gui/application_services.py` の共有lockに従って直列化する。
Qt applicationやGUI windowは生成しない。

動的レイヤーの共通区間に暦coverageを交差させる。静的不足は別項目に表示し、
「動的データは揃う」と「全依存が揃い出力可能」を区別する。

完了条件: 暦・静的データ不足を隠さず、追加ダウンロードなしで説明できる。

## 6. CLIと出力の接続

`cli/export_replay.py`、専用引数定義、`pyproject.toml` のscript登録を追加する。
共通のlocation/observing/opacity引数定義と検証を再利用し、全export parserを流用しない。
共通関数の変更が必要なら任意引数追加などで既存呼出の意味を維持する。

- 今回の必須モードは `--list-cache-ranges`。指定なしの画像生成は未対応と明示する。
- `--json`、`--cache-dir`、任意の `--start` / `--end` を受け付ける。
  開始・終了は片方だけでも範囲絞り込みに使える。両方ある場合はstart < endを要求する。
- `-A` / `-Z` / `-a`、`-B`、追加高さ、対応レイヤーのopacityを共有する。
  画角や画像寸法が入力範囲に影響する場合も同じ設定を保持する。
- HUDや画像保存・再開のオプションは今回先行実装しない。
- 未対応動的レイヤーは未指定なら無効、明示的な正のopacityはエラー。
- ASCIIの表形式と、schema version付きJSONを作る。Unicode地点名は端末ではescapeし、
  JSONはensure_ascii=Trueとする。JSON実行時のログはstderrへ出し、stdoutを汚さない。
- 該当履歴なし・破損recordの除外・共通区間なしは結果の状態として出し終了0。
  引数エラー、地点解決不能、権限不足などの走査不能は終了1とする。
  argparse標準の終了2との違いもCLI側で統一する。
- 照会中に消えたファイルは診断へ記録する。directory全体を読めない場合まで
  「履歴なし」に置き換えない。表示は照会時点の結果とする。

完了条件: 通常テキスト/JSONで同じ情報を取得でき、出力directoryや設定を作成しない。

## 7. 検証と仕上げ

意味のあるfixtureを使い、以下を段階ごとに検証する。

- CLI: 場所指定、既存別名・値域、-aの意味、期間入力、モード制約と終了コード。
- selector: 未来取得、scope不一致、鮮度境界、空snapshot、破損、schema不一致、
  同時刻候補、Horizons、期限より古いが残存するcache。
- coverage: 孤立した1分、途中欠損、隣接区間連結、空の積集合、全無効、片側期間指定。
- 実行全体: Qt未初期化、ネットワークfetch/cleanup/設定保存を禁止したfixtureで照会が完走。
  地点検索は別テストでmockし、検索だけが許容されることを確認する。
- JSON: stdoutを直接decodeし、UTC半開区間、件数、diagnosticsが表形式と一致する。
- 実cache: 利用可能な地点を一つ選び、読取専用で範囲と不足理由を確認する。
  データがない場合はfixtureで代替し、実cache検証未実施を明記する。

新規の照会テストと関連するCLI/location/cache regressionを先に実行する。
共有parser/resolverを変更した場合は該当既存テストも実行し、Ruff・mypy・diff checkを行う。
変更範囲に応じて既存の全体テストも実行し、結果をセッションログに記録する。
仕様書に実装済みの範囲照会だけを追加し、画像生成設計は未実装のまま維持する。
バージョン変更と公開・pushは今回の範囲に含めない。

## 完成の判定

地点指定で履歴保存範囲と欠損を確認でき、同じselectorで後続の毎分出力対象を決められる。
取得済み履歴があっても入力不整合や依存不足なら、その理由を確認できる。
通常GUIと単発exportの既存動作を維持し、照会がレイヤー取得やcache書換えを起こさない。

実装コミットは、(1) inventory/selector、(2) 地点/coverage/依存検査、
(3) CLI/出力/検証・仕様更新を目安にする。機械的に分けず、各コミットで整合性を保つ。
