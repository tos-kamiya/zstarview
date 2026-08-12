# GMNメテオ軌跡

この文書は、Global Meteor Network（GMN）の高レベルtrajectory summaryを取得し、
観測地点から見えた発光区間を天球座標へ固定するまでの内部設計をまとめる。GUI状態、
再投影、描画は後続段階で扱う。

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
  - 1000km候補抽出
  - WGS84地理座標からECEFへの変換
  - 発生時刻の幾何学的地平線での線分クリップ
  - 可視区間のAlt/AzからICRS RA/Decへの固定
- `meteors/service.py`
  - 表示時刻から24時間窓を作り、repositoryとprojectionを接続

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

directory indexと直近2日以内の配布ファイルは取得から6時間をfreshとする。それより前の
取得済み配布ファイルは完了済み履歴として長期再利用する。通信失敗時は、SHA-256、schema、
ファイル名を検証できる既存キャッシュをstaleとして使用する。配布ファイルの書き込みは
同じディレクトリの一時ファイルからatomic renameし、メタデータには取得時刻とSHA-256を
保存する。ハッシュ検証時は配布ファイルの改行コードを変換せず、取得時と同じテキストを
使用する。

1ファイルの取得失敗は窓全体の失敗へ昇格させない。取得できたファイルだけを解析し、
利用不能なファイル名を結果へ保持する。directory index自体をオンライン・キャッシュの
どちらからも読めない場合だけ、ファイル発見不能として呼び出し元へ例外を返す。

## 地理的選別と地平線クリップ

候補抽出は観測地点と始点、終点、球面上の中点の大円距離を計算し、いずれかが1000km
以内なら後段へ渡す。これは処理量削減であり可視性判定ではない。

始点、終点、観測地点はAstropyの`EarthLocation`でWGS84 ECEFへ変換する。観測地点の
測地緯度・経度からENU基底を作り、各端点の観測地点相対ベクトルをUp軸へ射影する。
両端が負なら発光区間全体を地平線下として除外する。一方だけが負なら、ECEF線分上で
Up成分が0になる点を線形補間し、可視側だけを残す。

## 天球座標固定

クリップ後の各端点を観測地点基準のENUへ変換してAlt/Azを求める。そのAlt/Azを流星の
開始UTC時刻と観測地点を持つAstropy `AltAz` frameとして構築し、ICRSへ変換する。
`CelestialMeteorTrail`には始点・終点のRA/Decを保存する。

このモデルは発生時刻の地理的Alt/Azを保存しない。後続の描画段階では表示時刻のAlt/Azへ
再投影することで、恒星背景と同じ天球固定の動きを得る。画面座標、色、不透明度、
18時間固定＋6時間線形フェードはこの段階のモデルへ含めない。
