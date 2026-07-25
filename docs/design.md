# zstarview 設計書

最終更新: 2026-07-15

この文書は、`zstarview` の内部設計の入口である。
`docs/design/` 以下に、責務ごとに分割した詳細文書を置く。

`zstarview` には次の 3 つのアプリケーション入口がある。

- `zstarview`
  - CLI 引数で初期状態を与えて GUI を起動する
- `zstarview-gui`
  - 起動前ダイアログを先に開いてから GUI を起動する
- `zstarview-export-image`
  - GUI を起動せず、1 枚の画像を headless で書き出す

補助的に `zstarview-debug` という console script 版の GUI ランチャーもある。
これは `zstarview` と同じ `main()` を呼ぶが、Windows では terminal を伴って起動しやすく、起動時ログを見たい診断用途に向く。

補助的に `zstarview-install-overturemaps-exe-cli` という staging 用の console script も置いてよい。
これはダウンロード済みの `overturemaps` Windows 実行ファイルを `CACHE_PATH/overturemaps.exe` にコピーするだけの単純な入口として扱う。
`overturemaps` 本体の探索やダウンロードは行わず、ファイル配置だけを責務にする。

これら 3 つは、地点解決、時刻解釈、描画、キャッシュ、外部データ取得の核心を共有する。
差分は「どの入口から始まるか」「対話 UI を持つか」「1 枚の画像で終わるか」にある。

### 夜間光データの配布

1.35.0 では、夜間光の元データを EOG の 2025 年次 VIIRS Nighttime Lights
VNL v2.2 `average_masked` とする。元配布物をリリース時の前処理でGeoTIFF
に変換し、PyPIパッケージには含めず、GitHub Release `night-lights-2025`の
リリースアセットとして配布する。実行時の固定入口は次のmanifest URLとする。

`https://github.com/tos-kamiya/zstarview/releases/download/night-lights-2025/manifest.json`

実行時はmanifestを取得してタイル一覧とSHA-256を確定し、観測地点から最大
夜間光距離までのレイが通過しうるタイルだけを取得する。ダウンロードした
ファイルは一時ファイルへ書き込み、SHA-256、GeoTIFFのCRS・解像度・寸法・
範囲を検証した後にキャッシュへatomic renameする。GitHub Pagesやアプリ作者の
Webサイトを実行時の必須依存にはしない。

manifest・各タイル・NOTICEは同じReleaseのAssetsとして管理する。manifestには
データセット識別子、データ年、各Assetの相対ファイル名・URL・SHA-256、元データ
URL・元ファイルのSHA-256、変換スクリプトのGitコミットを記録する。

配布アセットには、データ年・製品名・変換方法・SHA-256 checksum を記録した
manifest または README を添付する。EOG VNL は CC BY 4.0 の条件に従い、
EOG の帰属、GeoTIFF への変換を行った旨、関連する出典論文を配布物と
アプリのクレジットに記載する。

白背景で雲、明るい恒星、惑星、航空機、人工衛星、都市アウトライン、地形、水面、Earth guide を読む別入口アプリは、`zstarview` 本体の表示モードではなく別アプリ入口として設計する。
詳細は [atlas.md](design/atlas.md) に分離する。

### OvertureMaps exe staging

`import_overture_buildings.py` は `overturemaps` CLI を呼び出して Overture building data を取得する。現行の呼び出し経路は `--overturemaps-bin` を受け取り、`shutil.which(...)` で解決した実行ファイルか、明示されたパスをそのまま使う。

`zstarview-install-overturemaps-exe-cli` は、この経路に対して手元の Windows 版 exe を用意するための補助である。
目的は、GitHub Releases から取得した `v1.0.1` 以上の `*-windows-x86_64.exe` のような資産を、固定名 `overturemaps.exe` として `CACHE_PATH` 配下に staging することにある。
コピー後の実行パスを `--overturemaps-bin` に渡せば、既存の import パイプラインや GUI 側の呼び出しと接続できる。
この helper はリリース資産名の suffix を解釈しないため、バージョン番号やアーキテクチャ名は destination には残らない。
将来もし自動探索を足すなら、`CACHE_PATH/overturemaps.exe` を first-class に見る lookup を別途設計してもよいが、この helper 自体には含めない。

実行ファイルの解決順は、呼び出し側が `--overturemaps-bin` を明示した場合を最優先にし、次に `CACHE_PATH/overturemaps.exe` の staging 済みファイルを見て、最後に `PATH` 上の `overturemaps` を探索する形が自然である。
この順序なら、利用者が明示したパスを壊さず、Arm64 Windows の回避策として staging を使え、従来のインストール形態もそのまま残せる。

### PLATEAU building preparation

PLATEAUは、通常起動時の暗黙のダウンロードソースではなく、利用者が明示的に準備する任意の高密度ソースとして扱う。
準備入口は `zstarview-download-plateau-buildings` とし、APIからCityGMLを取得し、既存のデータソース非依存なderived building tile形式へ変換して `CACHE_PATH/plateau_buildings` 配下へ保存する。

CityGMLの取得・展開・変換はこの準備CLIプロセスの責務であり、通常のGUI workerや描画スレッドに持ち込まない。
CLIは大容量データを扱うため、取得前に対象データセット、年度、推定サイズ、必要空き容量を確認し、利用者の明示確認後に処理を開始する。
ZIPと展開ディレクトリは作業ディレクトリに置き、変換が正常完了した場合だけ完成済みキャッシュへ反映する。失敗時は部分キャッシュを有効化せず、可能なら一時データを削除する。

PLATEAU派生キャッシュの完了判定はディレクトリの存在だけに依存しない。`cache_meta.json` のメタデータスキーマ、派生タイルスキーマ、データソース、変換完了状態と `bldg/tile_index.json`、参照先タイルを検証し、観測地点と検索半径を含むキャッシュだけを選択する。

準備CLIの`--list`は、同じキャッシュ検証を使う読み取り専用の一覧経路である。通常出力は自治体コード・データ年度・保存パスの1行形式とし、`--city-code`で絞り込める。`--jsonl`では`cache_meta.json`の詳細に保存パスを加え、1キャッシュ1行で出力する。カタログAPIには接続せず、無効なキャッシュは一覧から除外する。
取得済みの元CityGMLを残すことは既定にせず、検査や再変換が必要な場合だけ明示オプションで保持する。
準備CLIを再実行した場合は、PLATEAUカタログの整備年度、登録年度、仕様、建物ファイル数、建物ファイルサイズ合計を既存の`cache_meta.json`と照合する。全項目が一致するキャッシュは再利用し、いずれかが異なる場合は利用者の確認後に再ダウンロードして置き換える。通常の`zstarview`起動時にはこのカタログ照合を行わない。

PLATEAU建物形状はLOD0/LOD1を基本のフォールバックとし、利用可能なLOD2については`RoofSurface`だけを追加利用する段階方式とする。LOD2の壁面全体やLOD3/LOD4は、空を背景とした建物シルエット表示に対してデータ量と処理量が過大になるため対象外とする。PLATEAUのCityGML APIが提供するLOD別地物数を、対象データの詳細度確認に利用してよい（[CityGML API](https://docs.plateauview.mlit.go.jp/api/rest/operations/datacatalogcitygmlconditions/)）。

派生タイルの建物レコードは、従来の外周リング・高さ・範囲を後方互換の基本情報として保持する。LOD1では、壁面を多数の線として直接描画せず、単純立体の屋根面を観測地点から見たシルエットへ変換し、既存の屋根ライン表示と同じ描画経路へ渡す。LOD2では、CityGMLの`RoofSurface`ポリゴンと各頂点の標高を保存し、各屋根面の外周を観測地点基準のAlt/Azへ投影する。同一建物内で標高差が小さい面群は、投影後に共有辺を相殺して外周リングへ連結し、完全包含される内側の面や穴の境界は描画対象から除外する。標高差の許容値は観測地点を0mとし、距離に比例して増加させ、0kmで-1m、1kmで0m、4kmで3m、4km以上では3mに固定する。許容値が負の場合は面の結合を行わず、個別のRoofSurfaceをそのまま描画する。壁面内部や底面は描画しない。表示密度を抑えるため、LOD2屋根面の描画は当面40m以上の建物に限定し、それ以外はLOD1/LOD0へフォールバックする。

LOD1の初期描画は、箱の上面を屋根ポリゴンとして扱う。上面の各頂点を観測地点基準の東・北・上座標へ変換し、Alt/Az座標上へ投影したポリゴン境界を、既存の建物外周ラインと同じ`UrbanOutlinePolyline`経路へ渡す。壁面内部や底面は描画対象にしない。方位角の0度/360度境界、視野境界、検索半径境界の分割、エッジのサンプリング、細い輪郭の簡略化は既存のOverture Maps用処理を再利用する。これにより、LOD1の立体情報を使いながら画面上の線の密度を増やさない。

形状情報を追加する派生タイルでは、`DERIVED_TILE_SCHEMA_VERSION`を3へ更新し、`geometry_lod`、`geometry_mode`、LOD別建物数、`roof_surfaces`、詳細サーフェス数をメタデータまたはタイルレコードへ保存する。スキーマ1/2のキャッシュは新しい形状方式に対応済みとはみなさず、準備CLIの再実行時に再生成する。キャッシュ更新判定では、既存の整備年度・登録年度・仕様・ファイル数・ファイルサイズの照合を引き続き適用する。

形状方式の実験対象は、東京タワーを含む港区の自治体コード`13103`とする。`13103`のCityGMLカタログでLOD2の提供状況を確認し、`RoofSurface`が存在する建物について派生タイルを生成する。東京タワー周辺では、Overture Mapsの従来表示、PLATEAUのLOD0、LOD1、LOD2屋根面シルエットを比較し、輪郭の再現性、画面の混雑、変換時間、キャッシュ容量、起動後の描画時間を評価する。

### Building source selection

建物ソースの選択はGUI側のcontrollerが行い、各ソースの取得・変換処理とは分離する。
選択順は、観測地点に対応する有効なPLATEAU派生キャッシュ、通常のOverture派生キャッシュ、最後にOvertureのオンデマンド取得である。PLATEAUが選択された場合は、近距離建物と遠距離skyscraper補助のOverture取得を行わず、PLATEAU派生タイルだけを共通の建物読込・輪郭計算経路へ渡す。
PLATEAU対象地域かどうかだけではPLATEAU取得を開始せず、PLATEAU準備済みキャッシュの有無を優先条件とする。

PLATEAUキャッシュが有効な地点では、Overtureデータを無条件に重ねない。これにより同一建物の二重表示を避け、PLATEAUを選択した利用者のデータソース意図を保持する。
PLATEAUキャッシュがない場合は、対象地域内であっても既存のOverture経路を使い、起動時の大容量ダウンロードを避ける。

### 外部 API の識別

外部 HTTP API へのリクエストは、`zstarview/<app-version> (+service)` 形式の識別可能な `User-Agent` を持つ。  
共通基底は現行版 `zstarview/1.32.11` で、各クライアントは短い suffix を足して区別する。
将来の版更新時は、`zstarview/<current-version>` の基底部分だけを差し替えればよい。

- `build_user_agent("water-overlay")` -> `zstarview/1.32.11 (+water-overlay)`
- `build_user_agent("nominatim")` -> `zstarview/1.32.11 (+nominatim)`
- `build_user_agent("night-lights")` -> `zstarview/1.32.11 (+night-lights)`
- `build_user_agent("overture-release")` -> `zstarview/1.32.11 (+overture-release)`
- `build_user_agent("geosatellite")` -> `zstarview/1.32.11 (+geosatellite)`
- `build_user_agent("tropical-cyclone")` -> `zstarview/1.32.11 (+tropical-cyclone)`
- `build_user_agent("ip-api")` -> `zstarview/1.32.11 (+ip-api)`
- `build_user_agent("opensky")` -> `zstarview/1.32.11 (+opensky)`
- `build_user_agent("satellites-celestrak")` -> `zstarview/1.32.11 (+satellites-celestrak)`
- `build_user_agent("satellites-horizons")` -> `zstarview/1.32.11 (+satellites-horizons)`
- `build_user_agent("satellites-wheretheiss")` -> `zstarview/1.32.11 (+satellites-wheretheiss)`
- `build_user_agent("copernicus-dem")` -> `zstarview/1.32.11 (+copernicus-dem)`
- `build_user_agent("s3")` -> `zstarview/1.32.11 (+s3)`
- `build_user_agent("skyfield-loader")` -> `zstarview/1.32.11 (+skyfield-loader)`

この方針は、サービス運営側のトラフィック識別を助けつつ、障害調査でどの経路が使われたかを追いやすくする。  
仕様上の公開一覧は `docs/specification.md` に置き、実装の増減があってもそこへ反映する。

- 夜間光は GeoTIFF 由来の base layer と、DEM と距離だけで作る edge glow layer に分ける。base layer は `night_light_opacity`、edge glow layer は `ridge_glow_opacity` で別々に調整する。
- base layer は副稜線レイヤー配列の順序をそのまま使う。`0` 番は最初の副稜線であり、主稜線は入力しない。
- edge glow は夜間光の色を借りず、`GlowMask` の固定 tint で描画する。ridge glow はまず粗い全域モデルを作り、その上に主稜線近傍だけ高解像度の補助モデルを重ねる二段構成としてよい。細密モデルの外側は外挿で埋めず、粗いモデルを背景値として使ってよい。
- `night_light` と `ridge glow` は別の強度成分として扱うが、どちらも最終的には同じ `GlowMask` 系へ折りたたんで描画する。night light は従来どおり直接描画し、ridge glow は粗い基盤モデルと局所的な補助モデルを合成してから描画時にマスクする。粗密の切り替えは硬く切らず、境界では滑らかにつないでよい。
- ただし glow の可視化は cloud のダウンロード完了と結びつけず、sky snapshot 完了時点で base sky/glow を先に再描画してよい。cloud overlay は別トリガーで後から重ねてよい。
- cloud overlay の描画色は昼間は白を基準にし、夜間は太陽高度だけを使って青灰色へ補間する。局所的な sky-disc 色から雲色を再推定しないことで、夕焼け域の過剰な色付きとキャッシュ複雑化を避ける。
- `zstarview-export-image` は単発処理なので、GUI の splash warm-up と同じ順序に固定しなくてよい。短命な aircraft / satellites を先に開始し、cloud は独立した取得経路として早めに走らせ、取得後に (alt, az) グリッドへの取り込みを行い、その後に terrain / urban / water / night light をまとめて揃える設計を許容してよい。

night light の有効条件は terrain horizon の生成結果の有無に合わせる。terrain horizon がまだない間は夜間光の alpha grid を作らず、terrain horizon が用意できた時点で 1 回だけ alpha grid を生成して保持する。以後は同じ terrain 条件ではその grid を使い回し、terrain horizon が再計算されたときだけ night light 側も再生成する。

### Water polygon and boundary overlay

水面の取得経路は、Overpass API の `natural=water`、`waterway=riverbank`、および対応する multipolygon を `WaterPolygonFootprint` として保持する。各フットプリントは外周リング・内周リング・水面種別・タグを持つ。`waterway=river` などの中心線だけの要素は、閉じた水域ポリゴンへ復元できない場合、水面フットプリントとして扱わない。

取得後は観測地点中心のローカル座標へ変換し、距離に応じたグリッドでリングを簡略化する。簡略化済みのフットプリントをキャッシュおよび後続処理の共通入力とし、取得元の高密度な頂点列を描画経路へ渡さない。これにより、水面ドットと境界線を追加しても、元データの頂点数に比例した描画負荷を避ける。

水面表示は同じ簡略化済みフットプリントから二つの出力を作る。

- 水面ドット: 観測地点からのレイを走査し、ポリゴン内部に入ったサンプルを `WaterOverlayPoint` に変換する。既存の距離ベースの間引き、透明度、カテゴリ色を適用する。
- 境界折れ線: 外周リング（必要に応じて内周リング）を観測地点基準の Alt/Az へ投影し、表示可能な隣接頂点を `WaterOverlayPolyline` 相当のフラグメントへまとめる。ポリゴンを面として塗りつぶすのではなく、低密度の輪郭線として描画する。

境界折れ線の生成では、まずローカル XY 上で水面走査と同じ最大距離の円にリングをクリップする。次に、距離円との交点を補間したリングを Alt/Az へ変換し、前方半球と表示 FOV の判定を行う。方位角の 0 度/360 度境界、視野境界、地平線下の表示限界をまたぐ部分は分割し、画面外の線分を描画しない。これは水面ドットと境界線で取得範囲・表示範囲の差が生じることを防ぐためである。

境界線は水面ドットの上に描画し、ポリゴン内部のドットは従来どおり残す。境界線の太さと透明度は距離に応じてドットと同等またはそれより弱くし、遠距離の川・湖が輪郭として読めることと、都市アウトラインや地形線が過密にならないことを両立する。この境界折れ線経路は本設計で定義する後続実装であり、既存の水面ドット経路との互換性を維持する。

## 文書構成

- [overview.md](design/overview.md)
  - 設計方針、全体アーキテクチャ、3 つのアプリの共通前提
- [runtime.md](design/runtime.md)
  - スレッドモデル、GUI 状態更新、処理フロー、エラー処理、キャッシュ方針
- [gui-screen-update-and-cache.md](design/gui-screen-update-and-cache.md)
  - GUI の更新トリガー、再描画の流れ、フレーム/合成/状態キャッシュの整理
- [data-model.md](design/data-model.md)
  - 主要データ構造、状態オブジェクト、アプリ間で共有する scene/state の境界
- [rendering-pipeline.md](design/rendering-pipeline.md)
  - 描画パイプライン、オーバーレイ合成、ラベル、外部依存
- [atlas.md](design/atlas.md)
  - Atlas の入口、表示プロファイル、白背景向け描画設計
- [legacy-archive.md](design/legacy-archive.md)
  - 分割前の単一ファイル版の記録。参照用の履歴として残す

利用者向けの機能説明は [specification.md](specification.md) を参照する。

## 起動時のレイヤー分担

- スプラッシュ表示中は、観測地点と観測高度が確定したあとに変化しない静的レイヤーを先行して warm-up してよい。
- この静的レイヤーには、DEM、地形地平線、夜間光と ridge glow、水面、都市アウトラインを含めてよい。
- これらの準静的レイヤーは、既存キャッシュの再利用、欠損時のダウンロード、初回描画用キャッシュの作成までを含めてよい。
- 雲、Geo-satellite 雲、航空機、人工衛星、台風・サイクロンのような動的レイヤーは、通常 GUI が見えてから遅延させてよい。
