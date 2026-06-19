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

夜間光レイヤーは副稜線レイヤー配列の順序をそのまま使う。`0` 番は最初の副稜線であり、主稜線は入力しない。ridge glow は夜間光とは独立した補助レイヤーとして扱い、主稜線のうち遠方サンプルのみを薄く持ち上げる。
街灯の帯は、`az` ごとに近い距離からレイを進めて光を積算する。`az` は 2° 刻みの試行でよく、距離方向は 1 km 刻みで十分である。各レイでは、手前の稜線や地平線により柱の下端が持ち上がることがあるため、可視部分だけを積算して奥へ進む。
ridge glow の色は sky disc の色と夜間光の色を直接混ぜるのではなく、`GlowMask` の固定 tint で決める。実装上は alpha-only の低解像度マスクへ畳み込み、復元時に一度だけ tint して合成する。
現在の実装では、GlowMask を night light と ridge glow の 2 経路で共有する。night light 側は低解像度の画面グリッドを逆投影して `alt` / `az` を求め、night-light profile の azimuth サンプルを補間した連続 alpha field として描く。ridge glow 側は主稜線の遠方サンプルだけを使い、主稜線 altitude をマスクにした小さな追加 alpha を積算する。どちらも最終的な可視化は密度場として扱う。
この方式では、幾何の輪郭を作るシャープ層と、にじみを作るグロー層を分ける。前者は現在の稜線や帯ポリゴンを担い、後者は夜間光の柱状寄与と ridge glow の微小な補助分をまとめて保持する。見た目の色は固定 tint に任せ、sky disc との混色は行わない。

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
