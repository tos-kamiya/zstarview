# 24時間キャッシュと過去時刻表示

最終更新: 2026-10-11

この文書は、雲・航空機・人工天体の24時間キャッシュ保持と、過去の天体表示中に各データレイヤーを現在値で更新し続ける動作を定める。完全な画面録画や起動状態の復元は対象にしない。

## 目的と範囲

- 通常起動中に取得したデータを24時間キャッシュへ蓄積する。
- 過去24時間以内の天体表示中も、雲、航空機、人工天体は通常経路で現在時刻のデータを取得・更新して表示する。
- データの取得時刻、観測時刻、対象範囲を記録し、通常の鮮度・scope判定に使う。
- 画面PNG、毎分のscene manifest、ウィンドウ状態や表示スタイルの復元は要件に含めない。
- プロセスsession IDはキャッシュ選択に使わない。

GUIでは既存の `--datetime`、`--hours`、`--days` で時刻を選択する。GUI以外の呼び出し元とheadless画像出力は、従来の時刻シフト時のレイヤー可否を維持する。

履歴が存在するのはこの機能を含むバージョンで取得が成功した後に限る。導入前のデータは復元できない。画面上の星・惑星など既存の天体表示は、選択時刻に対する通常の計算を行う。

## 時刻と表示モデル

- 要求時刻はUTCで扱う。UIが現地時刻を入力・表示する場合はタイムゾーンでUTCへ変換する。
- 過去表示中も指定時刻は天体計算に使う。`--hours` の相対指定は現在時刻に対するオフセットなので、表示対象時刻も実時間とともに進む。
- 雲、航空機、人工天体のデータ要求には現在UTC時刻を使い、通常のネットワーク取得、cache hit、rate limit、retry、stale fallbackを継続する。指定時刻より新しいsource dataも利用する。
- 各sourceの既存の鮮度条件とscope判定を維持する。取得・cache利用に失敗した場合は通常時と同じく当該レイヤーを `unavailable` または `cache` とし、他レイヤーは継続する。
- レイヤーのsource時刻・取得時刻が天体表示の指定時刻と異なることは許容する。これは過去の気象・航空機・人工衛星状態を再現するものではない。
- データ取得・雲画像生成途中の部分成果物を履歴として公開しない。雲は完成済みsource/gridだけを選ぶ。

## 地点解決の起動メタデータ

`--place` を指定した起動では、入力キーワードと解決結果を一つの起動メタデータJSONへ記録する。これは画面状態のmanifestではなく、地点の来歴とキャッシュscopeを確認するための記録である。`--place` 以外の起動方法に同じ記録を必須とはしない。

保存先は `CACHE_PATH/replay/v1/launches/YYYY/MM/DD/HH/mm/<launch-id>.json` とする。`launch-id` は同じ分内の複数起動の衝突を避けるUUIDとし、プロセスsession IDには依存しない。時刻階層とメタデータの起動時刻はUTCを使う。

```json
{
  "schema_version": 1,
  "launch_id": "uuid",
  "started_at_utc": "2026-10-11T03:04:05Z",
  "location": {
    "input": {
      "option": "--place",
      "query": "松江駅",
      "countrycode": null,
      "language": "en"
    },
    "resolved": {
      "display_name": "Matsue Station, Matsue, Shimane, Japan",
      "latitude_deg": 35.4641778,
      "longitude_deg": 133.0628539,
      "observer_height_m": 1.7,
      "timezone": "Asia/Tokyo",
      "source": "nominatim"
    }
  }
}
```

要件:

- `query` はCLIに渡された `--place` の文字列を保持する。検索用に正規化した値で置き換えない。
- `latitude_deg` と `longitude_deg` は実際の地点解決結果を保持し、有限値と緯度・経度範囲を検証する。
- `display_name`、`observer_height_m`、`timezone`、`source` は解決結果で取得できる場合に保存する。高度は場所の標高と観測者の追加高さを混同しないよう、実効観測高度として定義する。
- `countrycode` と `language` は指定値または検索時に適用した値を記録する。未指定値は `null` または仕様上の既定値を一貫して保存する。
- `source` は地点解決種別（例: `nominatim`）を記録する。
- メタデータは地点解決成功後に書き、JSONを一時ファイルへ完全に書き込んでからatomic renameする。書き込み失敗は通常起動を妨げず、ログへ記録する。
- メタデータは起動時刻から24時間で削除する。個人の入力文字列を含むため、通常キャッシュと同じローカルキャッシュ領域に置き、外部送信しない。

## レイヤーごとの履歴データと選択

全履歴レコードはUTCの `fetched_at_utc` を持つ。元データが表す時刻は別フィールドとして保持し、取得時刻から推定しない。履歴ファイル名またはディレクトリの取得時刻とJSON内部時刻が矛盾するレコードは利用しない。

### 雲

- 既存のGOES/Himawari配信元時刻入りsource cacheを履歴の基本とし、既存の24時間cleanupを継続する。
- レコードには少なくともsource、source/product時刻、observer緯度経度、source key、取得完了状態と参照assetを結び付ける。
- 要求時刻に対応する完成済みsource一式を選び、通常の投影・表示設定で描画する。異なるsource時刻やobserver scopeのタイルを一つの履歴画像として混成しない。
- 雲sourceの取得途中・projection途中の成果物は選択対象にしない。

#### Geo-satellite

- MET Norway Geo-satelliteの元画像PNGを、UTCの `captured_at_utc` と画像種別ごとの履歴として保存する。取得時刻 `fetched_at_utc` は撮影時刻と区別してsidecarに記録する。
- 保存先は `GEOSATELLITE_CACHE_ROOT/raw/history/YYYY/MM/DD/HH/mm/<kind>.png` とし、同じディレクトリに `<kind>.json` sidecarを置く。時刻階層は撮影時刻を使うため、同じ画像スロット・種別の重複取得は同一レコードへ集約する。
- sidecarにはschema version、画像種別、撮影時刻、取得時刻、area/size、source URL、画像の相対pathとSHA-256を記録する。読込時はcache root内のpath、digest、画像種別、時刻を検証する。
- 過去の天体表示中も通常のlatest-available解決、source cache、ネットワーク取得を使い、現在時刻の最新画像を選ぶ。24時間のraw-image履歴は引き続き保持するが、過去再生用の画像選択には使わない。
- 元画像は観測地点や画角に依存しないため、observer scopeやカメラ設定をcache keyに含めない。過去のraw PNGから現在の表示地点・画角に合わせて派生データを再構成する。
- 履歴の期限は保存した取得時刻から24時間とする。ライブ経路は最新の利用可能時刻を引き続き選ぶ。既存利用者との互換性のため `raw/latest.png` と `raw/latest.json` も更新するが、過去時刻の選択には時刻別履歴を使う。

### 航空機

- 5分程度の通常更新で取得成功した各snapshotを不変履歴として保存する。bboxは observer latitude/longitude から得た範囲であり、選択キーに含める。
- snapshotは `fetched_at_utc` とOpenSkyの状態時刻（例: `last_contact`）を別に持つ。
- 過去の天体表示中も通常のfetch/cache経路で現在のsnapshotを更新し、現在のobserver bboxに一致する最新のsnapshotを使う。
- snapshotは従来どおり現在時刻に対する鮮度・rate-limit・stale fallbackを適用する。航空機の状態は選択中の天体表示時刻の過去状態を意味しない。

### 人工天体

- 各成功取得の軌道要素または観測記録を取得時刻階層へ保存する。既存の最新値ファイルはライブ表示との互換用に維持してよい。
- TLE/OMM等は通常の現在時刻に基づく取得・cache freshness規則を使い、取得済みの最新要素から選択中の天体表示時刻へ伝播して描画する。
- Horizons等のsourceも通常経路で現在時刻を渡して更新する。過去表示中も既存のobserver scopeを維持する。
- 要素取得に失敗した場合は人工天体レイヤーへ通常のstale fallbackまたは `unavailable` を適用する。

## データ構造と保存

現在の航空機・人工天体履歴配置を継続し、レコードに再生選択に必要なscopeと時刻情報を保持する。例:

```text
AIRCRAFT_CACHE_ROOT/history/YYYY/MM/DD/HH/mm/<request-id>-<bbox>.json
SATELLITE_CACHE_ROOT/history/YYYY/MM/DD/HH/mm/<request-id>-<group>-<scope>.json
GEOSATELLITE_CACHE_ROOT/raw/history/YYYY/MM/DD/HH/mm/<kind>.png
GEOSATELLITE_CACHE_ROOT/raw/history/YYYY/MM/DD/HH/mm/<kind>.json
```

概念上の共通envelope:

```json
{
  "schema_version": 1,
  "kind": "aircraft_snapshot",
  "fetched_at_utc": "2026-10-11T03:05:00Z",
  "valid_time_utc": "2026-10-11T03:04:42Z",
  "scope": {
    "observer_latitude_deg": 35.4641778,
    "observer_longitude_deg": 133.0628539,
    "bbox": [132.8, 35.2, 133.3, 35.7]
  },
  "payload": {}
}
```

`payload` はsource固有schemaを維持する。雲は各既存source assetとgrid manifest、航空機は状態snapshot、人工天体はelement setまたは時刻固定の計算結果を参照できる形にする。大きな雲assetをJSONへ埋め込まず、相対pathと必要に応じてSHA-256を記録する。キャッシュ参照はcache root内に限定する。

履歴書き込みは一時ファイル、schema検証、atomic renameの順で確定する。最新値ファイル、航空機の10分stale fallback、人工天体の現在表示でのfresh/stale条件は履歴保持仕様によって変更しない。

## 期限切れデータの削除

- レイヤー履歴と起動メタデータはUTCの取得・起動時刻から24時間経過後に期限切れとする。
- cleanupは24時間を超えたレコードと空の履歴ディレクトリを削除する。現在の通常表示が読み取り中のファイルはその読み取りが終わるまで削除しない。
- 履歴の保存期間は、航空機のライブstale fallback期間、人工天体要素の表示有効期間とは独立させる。
- 削除・書き込みに失敗してもGUI表示を停止させない。失敗はログへ記録する。

## 実装上の注記

- 雲source履歴は完成済みの正規化 `CloudSourceData` artifactとJSON sidecarを、CloudDisc cache root内の取得時刻階層へ保存する。sidecarにはschema、取得時刻、source時刻/key、observer座標、artifact相対path、SHA-256を記録する。artifactはsidecar公開後にのみ再生候補となり、読込時にcache root包含、digest、型、observer scope、時刻を検証する。
- 航空機・人工天体の履歴は既存の `history/YYYY/MM/DD/HH/mm/` 配置を読む。航空機はmatching bbox、人工天体はgroupと必要なHorizons observer scopeを選択キーにする。
- 過去レイヤーの読み込み開始時に前の時刻のレイヤー状態をクリアする。cache missや検証失敗で、前の画面または現在値を残さない。
- 実装計画と検証範囲は `dev-notes/plan-historical-layer-replay.md` を参照する。
