# zstarview 設計書

最終更新: 2026-07-08

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
