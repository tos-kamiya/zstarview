# 設計概要

この文書は、次の 3 つのアプリケーション入口で共通する設計の土台をまとめる。

- `zstarview`
- `zstarview-gui`
- `zstarview-export-image`
- `zstarview-debug`
  - `zstarview` と同じ GUI 起動処理を console script として提供する診断用入口である

## 共通モデル

3 つの入口は、次の領域で同じドメインモデルを共有する。

- 地点解決
- 時刻解析とタイムゾーン処理
- 天球幾何と天体計算
- 星、惑星、ラベル、ガイド、補助レイヤーの描画規則
- キャッシュの扱いと長寿命外部データ

違いはオーケストレーションにある。

- `zstarview` は CLI で与えた状態から GUI を直接起動する。
- `zstarview-gui` は起動前ダイアログを先に開き、より豊かな launch profile を保存する。
- `zstarview-export-image` は headless で動作し、必要なデータを待って 1 枚だけ画像を書き出す。
- `zstarview-debug` は `zstarview` と同じ処理だが、Windows で terminal から起動しやすいよう console script として配布する。

## アプリケーション入口

### `zstarview`

- CLI 引数を受け取り、そのままメイン GUI を起動する。
- 地点、時刻、検索、描画の引数モデルは、可能な限り export 側と共有する。
- 共有設定は `~/.config/zstarview/config.json` に保存する。

### `zstarview-gui`

- 起動前ダイアログを先に開く。
- ダイアログでの選択結果を保持するため、より豊かな launch profile を使う。
- `zstarview` / `zstarview-export-image` の legacy 設定ファイルとは独立して保存する。

### `zstarview-export-image`

- headless で動作し、1 枚の画像を生成したら終了する。
- `zstarview` と共通の引数定義を使える範囲では共有する。
- 出力先、画像サイズ、タイムアウト、部分データ許可、sixel 出力など、画像出力専用オプションを追加する。

### 共通の CLI と解決ルール

- `--place`、`--place-countrycode`、`--place-lang` はオンラインの place 検索を使う。
- `--search` と `--search-keep-marker` は GUI の検索と export-image の検索で共有する。
- `--list` は export-image 専用とする。
- tower / mountain viewpoint のデータセット問い合わせは即時終了パスとして扱い、GUI 本体へは入らない。

### 利用者向けのアプリケーション像

文書上は、この製品を 1 つの単一ランチャーではなく、次の 3 系統として扱う。

- CLI 駆動の GUI ランチャー
- ダイアログ先行の GUI ランチャー
- headless の画像出力器

この切り分けは、インストール説明、利用例、エラー処理、設定永続化の説明に直接関係する。

## 設計方針

- UI の仕事と長時間処理・ネットワーク I/O を分離する。
- 入口ごとの分岐より、共有ドメインヘルパーを優先する。
- 失敗は可能な限り利用者に見える状態へ正規化する。
- 描画パイプラインは GUI と export で再利用できるようにする。
- 更新頻度の低い補助レイヤーは独立パイプラインとして扱う。
- ターミナル、ログ、サブプロセス向け文言は、UI 表示が必要な場合を除き ASCII を原則とする。

## 上位レイヤー

- 起動と設定
- ドメイン計算
- 描画
- UI オーケストレーション
- 外部データ取得とキャッシュ管理

以降の `docs/design/` は、この上位レイヤーを各責務ごとに掘り下げる。

## ソース構造の概観

`src/zstarview/` 配下は、責務ごとに大きく次のまとまりに分かれる。

- `cli/`
  - `zstarview` と `zstarview-export-image` の引数解析、即時終了コマンド、出力系 CLI
- `gui/`
  - ウィンドウ本体、イベント処理、worker、状態反映、起動前ダイアログ
- `render/`
  - 共通描画パイプライン、ラベル、ガイド、背景、オーバーレイ合成
- `search/`
  - 星、アステリズム、place、JPL 小天体、衛星検索の共通解決層
- `location_resolver/`
  - 都市名、座標、Nominatim、`auto`、タワー / 山の地点解決
- `clouddisc/`
  - 雲ソースの取得、射影、キャッシュ、サンプリング
- `satellites/`
  - 人工衛星の取得、正規化、追跡、キャッシュ
- `aircraft/`
  - OpenSky 航空機の取得、正規化、短時間予測、キャッシュ
- `terrain/`
  - DEM と地形地平線の計算
- `geosatellite/`
  - 実験的 Geo-satellite 雲経路の取得と投影
- `tropical_cyclones/`
  - 台風・サイクロン補助レイヤーの取得と正規化
- `data/`
  - 同梱資産、生成スクリプト、データ変換補助
- `utils/`
  - 共通ユーティリティ、表示整形、時刻・位置変換

この粒度は設計の入口としての概観であり、個別モジュールの責務は各詳細文書とソースコードで追う。
