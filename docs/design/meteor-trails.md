# GMNメテオ軌跡

この文書は、Global Meteor Network（GMN）の高レベルtrajectory summaryを取得し、
観測時刻に観測地点から見えた発光区間をAlt/Azへ固定し、GUIへ描画するまでの内部設計をまとめる。

## モジュール責務

- `meteors/client.py`
  - GMN公式配布元へのHTTP GET
  - `zstarview/<version> (+gmn-meteors)` User-Agentの付与
- `meteors/parser.py`
  - 公式daily directory indexから通常の日次ファイル名を抽出
  - semicolon区切りtrajectory summaryの必須列を正規化
  - 不正な個別レコードをスキップ
- `meteors/repository.py`
  - 表示窓を覆う配布ファイルの発見、取得、キャッシュ
  - exact UTC window filter、軌跡IDによる重複排除、部分結果の返却
- `meteors/projection.py`
  - 500km候補抽出
  - WGS84地理座標からECEFへの変換
  - 500km候補の各端点を発生時刻の観測地点基準Alt/Azへ変換
- `meteors/service.py`
  - 表示時刻以前の最新観測を終端とする24時間窓を決め、最新100軌跡へ制限してprojectionを接続
- `gui/meteor_controller.py`, `gui/meteor_state.py`
  - `Meteor trails`を有効にした後だけworker poolで取得・変換を実行
  - 成功後は1時間ごと、失敗後は10分後に再試行
- `render/meteors.py`
  - 発生時刻に保存したAlt/Az端点を画面座標へ投影して細線として描画
  - 線の開始端に小型の観測時刻差ラベルを直接描画する（重なり回避なし）
  - 表示時刻から24時間は100%、そこから72時間で30%まで線形減衰、以降は30%

## 表示窓の決定

GMNの高レベルデータは定期更新されるが、個々の軌跡は観測直後に登録されるとは限らない。
そのため、`display_time - 24h`から`display_time`までを固定的に要求窓にはしない。

repositoryで表示時刻以前のレコードを探索し、その中の最新`Beginning UTC Time`を
`window_end`とする。`window_start = window_end - 24h`とし、両端を含むレコードだけを
最終結果へ残す。過去表示でも`display_time`より後の観測は候補に含めない。利用可能な
最新時刻を見つけるための探索範囲は、GMNの通常の公開遅延を覆えるよう複数のdaily fileへ
段階的に広げてよいが、無制限には遡らない。

フェードの年齢基準には`window_end`ではなく`display_time`を使う。不透明度は経過24時間
まで1.0、そこから72時間かけて0.3まで線形減衰し、それ以降は`0.3`のままとする。
GMNの登録は最短でもおおよそ1日遅れるため、24時間の据え置きがないと基準値に届かない。
公開遅延のある古い窓ほど薄くなるが、表示対象の軌跡は消さない。
各端点の位置は各観測の実際の開始UTC時刻に観測地点から見えたAlt/Azとして保持し、
表示時刻が変わっても再投影しない。

ステータス用には`display_time - window_start`と`display_time - window_end`を時間へ変換する。
古い側は切り上げ、新しい側は切り捨て、`M {count}, {oldest}-{newest}h ago`としてASCIIだけで
組み立てる。例えば54時間前から30時間前の0件窓は`M 0, 54-30h ago`となる。

表示する軌跡は、選択窓内の`Beginning UTC Time`が新しい順に最大100件とする。100件を超える
場合は古い軌跡から表示対象外とする。窓の時間範囲は表示対象を制限する前の選択窓に基づき、
フェードは描画時の表示時刻からの経過時間だけを使う。

各線の開始端には表示時刻との差を`-32h`のような小型ラベルで描画する。これは通常のラベル
候補・重なり回避経路を使わず、軌跡ごとに直接配置する補助情報とする。

## GMN配布単位

GMNの`daily/`配下にある通常ファイルはUTC 00:00区切りではない。ファイル名は
`traj_summary_YYYYMMDD_solrange_A-B.txt`で、太陽黄経1度分のおおむね24時間を収録する。
同じ暦日に複数ファイルが始まる場合や、ファイル名の日付が1日飛ぶ場合がある。

そのため日付から`solrange`を推定してURLを組み立てない。daily directory indexを解析し、
要求窓の開始日・終了日に前後1日を加えた範囲の通常ファイルを候補とする。候補を読み込んだ
後、各レコードの`Beginning UTC Time`で要求窓へ厳密に絞り込む。`latest_daily`と
`yesterday`は可変内容の別名なので、永続キャッシュの識別子には使用しない。

## パースと入力検証

trajectory summaryはsemicolon区切りで、コメント行は`#`から始まる。表示用ヘッダーは
名称行と単位行に分かれ、同名のsigma列を多数含むため、GMNの公開カラム仕様に対応する
固定位置で必須値を読む。

必須値は軌跡ID、開始UTC時刻、始点・終点の緯度、経度、WGS84楕円体高である。緯度、
経度、高度の値域、数値のfinite性を検証し、必須値が不正なレコードだけをスキップする。
Duration、Peak AbsMag、Vinit、IAU shower codeは任意値として保持する。`...` shower codeは
sporadicを表すため内部では`None`へ正規化する。

## キャッシュ

キャッシュルートは`CACHE_PATH/meteors/gmn/schema-1/`とする。

```text
daily-index.json
daily/
  traj_summary_YYYYMMDD_solrange_A-B.txt
  traj_summary_YYYYMMDD_solrange_A-B.txt.json
```

directory indexと各配布ファイルは、対応するmetadataの`fetched_at_utc`から6時間をfreshとする。
ファイル名の代表日や観測レコードの時刻による鮮度判定の例外は設けない。通信失敗時は、
SHA-256、schema、ファイル名を検証できる既存キャッシュをstaleとして使用する。配布ファイルの書き込みは
同じディレクトリの一時ファイルからatomic renameし、メタデータには取得時刻とSHA-256を
保存する。ハッシュ検証時は配布ファイルの改行コードを変換せず、取得時と同じテキストを
使用する。

`fetched_at_utc`は、そのファイルをGMNから正常に取得できた最後の時刻である。取得内容が
前回と同じ場合も正常取得なら更新し、通信失敗時は更新しない。GUIは有効化中におおむね
1時間ごとに更新確認を行うが、実際の再取得は各ファイルの6時間TTLが切れた場合だけ行う。

1ファイルの取得失敗は窓全体の失敗へ昇格させない。取得できたファイルだけを解析し、
利用不能なファイル名を結果へ保持する。directory index自体をオンライン・キャッシュの
どちらからも読めない場合だけ、ファイル発見不能として呼び出し元へ例外を返す。

## 地理的選別と観測時Alt/Az変換

候補抽出は観測地点と始点、終点、球面上の中点の大円距離を計算し、いずれかが500km
以内なら後段へ渡す。これは対象地域を限定するための距離選別であり、地平線による可視性
判定や線分クリップは行わない。

始点、終点、観測地点はAstropyの`EarthLocation`でWGS84 ECEFへ変換する。観測地点の
測地緯度・経度からENU基底を作り、各端点の観測地点相対ベクトルをAlt/Azへ変換する。
端点が観測時刻の水平線下にある場合も、500km以内の候補であれば除外しない。

## 観測時Alt/Az固定

クリップ後の各端点を観測地点基準のENUへ変換して、流星の開始UTC時刻におけるAlt/Azを
求める。`CelestialMeteorTrail`には始点・終点のAlt/Azを保存し、ICRS RA/Decへ固定しない。

後続の描画段階では、保存済みの観測時Alt/Azをそのまま画面座標へ投影する。これにより、
表示時刻が観測時刻からずれても、軌跡は「その観測地点から発生当時に見えた方向」に残る。
地球の自転によって歴史的な軌跡を現在の恒星背景へ追従させることは、このレイヤーの目的と
しない。画面座標、色、不透明度、表示時刻基準の線形フェードはこの段階のモデルへ含めない。
