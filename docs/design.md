# zstarview 設計書

最終更新: 2026-06-16

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

夜間光レイヤーは副稜線レイヤー配列の順序をそのまま使う。`0` 番は最初の副稜線であり、主稜線は入力しない。主稜線側の sky glow は別のパスで描画し、街灯側の opacity とは独立に制御する。
街灯の帯は、`az` ごとに近い距離からレイを進めて光を積算する。`az` は 2° 刻みの試行でよく、距離方向は 1 km 刻みで十分である。各レイでは、手前の稜線や地平線により柱の下端が持ち上がることがあるため、可視部分だけを積算して奥へ進む。
sky glow の色は、raw な sky disc の地平線付近サンプルを基準にしつつ、night lights の基準色を混ぜた補助色として決める。夜間で sky disc がほぼ黒に落ちる場合でも、sky glow が消えすぎないようにするための扱いである。
sky glow の描画処理は `src/zstarview/render/ridge_glow.py` 側に分離し、夜間光本体とは別の責務として扱う。
将来的には、街灯の帯表現をフル解像度のポリゴンから切り離し、`GlowMask` へ rasterize してから合成する。`GlowMask` は `float32` の alpha だけを持つ。街灯側は、低解像度の画面グリッドを逆投影して `alt` / `az` を求め、night-light profile の azimuth サンプルを補間した連続 alpha field として描く。地平線より上だけを固定高さの減衰で残し、最終的な可視化は密度場として扱う。
この方式では、幾何の輪郭を作るシャープ層と、にじみを作るグロー層を分ける。前者は現在の稜線や帯ポリゴンを担い、後者は夜間光の柱状寄与をまとめて保持する。最終的な見た目調整は、低解像度の alpha マスクに対する平滑化と、固定色での tint で行う。街灯の描画は ridge glow と分離して扱い、まずは街灯のみを置き換える。

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
