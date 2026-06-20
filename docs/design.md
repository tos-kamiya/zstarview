# zstarview 設計書

最終更新: 2026-06-20

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

これら 3 つは、地点解決、時刻解釈、描画、キャッシュ、外部データ取得の核心を共有する。
差分は「どの入口から始まるか」「対話 UI を持つか」「1 枚の画像で終わるか」にある。

### 外部 API の識別

外部 HTTP API へのリクエストは、`zstarview/<app-version> (+service)` 形式の識別可能な `User-Agent` を持つ。  
共通基底は現行版 `zstarview/1.31.2` で、各クライアントは短い suffix を足して区別する。
将来の版更新時は、`zstarview/<current-version>` の基底部分だけを差し替えればよい。

- `build_user_agent("water-overlay")` -> `zstarview/1.31.2 (+water-overlay)`
- `build_user_agent("nominatim")` -> `zstarview/1.31.2 (+nominatim)`
- `build_user_agent("night-lights")` -> `zstarview/1.31.2 (+night-lights)`
- `build_user_agent("overture-release")` -> `zstarview/1.31.2 (+overture-release)`
- `build_user_agent("geosatellite")` -> `zstarview/1.31.2 (+geosatellite)`
- `build_user_agent("tropical-cyclone")` -> `zstarview/1.31.2 (+tropical-cyclone)`
- `build_user_agent("ip-api")` -> `zstarview/1.31.2 (+ip-api)`
- `build_user_agent("opensky")` -> `zstarview/1.31.2 (+opensky)`
- `build_user_agent("satellites-celestrak")` -> `zstarview/1.31.2 (+satellites-celestrak)`
- `build_user_agent("satellites-horizons")` -> `zstarview/1.31.2 (+satellites-horizons)`
- `build_user_agent("satellites-wheretheiss")` -> `zstarview/1.31.2 (+satellites-wheretheiss)`
- `build_user_agent("copernicus-dem")` -> `zstarview/1.31.2 (+copernicus-dem)`
- `build_user_agent("s3")` -> `zstarview/1.31.2 (+s3)`
- `build_user_agent("skyfield-loader")` -> `zstarview/1.31.2 (+skyfield-loader)`

この方針は、サービス運営側のトラフィック識別を助けつつ、障害調査でどの経路が使われたかを追いやすくする。  
仕様上の公開一覧は `docs/specification.md` に置き、実装の増減があってもそこへ反映する。

- 夜間光は GeoTIFF 由来の base layer と、DEM と距離だけで作る edge glow layer に分ける。base layer は `night_light_opacity`、edge glow layer は `ridge_glow_opacity` で別々に調整する。
- base layer は副稜線レイヤー配列の順序をそのまま使う。`0` 番は最初の副稜線であり、主稜線は入力しない。
- edge glow は夜間光の色を借りず、`GlowMask` の固定 tint で描画する。ridge glow は alpha grid の段階では主稜線で強く切らず、主稜線の下側も含めた中間強度を保持しておき、描画時に terrain mask で最終クリップしてよい。
- `night_light` と `ridge glow` は別の強度成分として扱うが、どちらも最終的には同じ `GlowMask` 系へ折りたたんで描画する。night light は従来どおり直接描画し、ridge glow は未マスクの中間 grid を経由してから描画時にマスクする。
- ただし glow の可視化は cloud のダウンロード完了と結びつけず、sky snapshot 完了時点で base sky/glow を先に再描画してよい。cloud overlay は別トリガーで後から重ねてよい。
- `zstarview-export-image` は単発処理なので、GUI の splash warm-up と同じ順序に固定しなくてよい。短命な aircraft / satellites を先に開始し、cloud は独立した取得経路として早めに走らせ、その後に terrain / urban / water / night light をまとめて揃える設計を許容してよい。

night light の有効条件は terrain horizon の生成結果の有無に合わせる。terrain horizon がまだない間は夜間光の alpha grid を作らず、terrain horizon が用意できた時点で 1 回だけ alpha grid を生成して保持する。以後は同じ terrain 条件ではその grid を使い回し、terrain horizon が再計算されたときだけ night light 側も再生成する。

## 文書構成

- [overview.md](design/overview.md)
  - 設計方針、全体アーキテクチャ、3 つのアプリの共通前提
- [runtime.md](design/runtime.md)
  - スレッドモデル、GUI 状態更新、処理フロー、エラー処理、キャッシュ方針
- [data-model.md](design/data-model.md)
  - 主要データ構造、状態オブジェクト、アプリ間で共有する scene/state の境界
- [rendering-pipeline.md](design/rendering-pipeline.md)
  - 描画パイプライン、オーバーレイ合成、ラベル、外部依存
- [legacy-archive.md](design/legacy-archive.md)
  - 分割前の単一ファイル版の記録。参照用の履歴として残す

利用者向けの機能説明は [specification.md](specification.md) を参照する。

## 起動時のレイヤー分担

- スプラッシュ表示中は、観測地点と観測高度が確定したあとに変化しない静的レイヤーを先行して warm-up してよい。
- この静的レイヤーには、DEM、地形地平線、夜間光と ridge glow、水面、都市アウトラインを含めてよい。
- これらの準静的レイヤーは、既存キャッシュの再利用、欠損時のダウンロード、初回描画用キャッシュの作成までを含めてよい。
- 雲、Geo-satellite 雲、航空機、人工衛星、台風・サイクロンのような動的レイヤーは、通常 GUI が見えてから遅延させてよい。
