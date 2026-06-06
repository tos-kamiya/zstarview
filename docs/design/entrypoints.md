# アプリケーション入口

この文書では、共有コアの上に載る 3 つのアプリケーション入口を整理する。

## `zstarview`

- CLI 引数を受け取り、そのままメイン GUI を起動する。
- 地点、時刻、検索、描画の引数モデルは、可能な限り export 側と共有する。
- 共有設定は `~/.config/zstarview/config.json` に保存する。

## `zstarview-gui`

- 起動前ダイアログを先に開く。
- ダイアログでの選択結果を保持するため、より豊かな launch profile を使う。
- `zstarview` / `zstarview-export-image` の legacy 設定ファイルとは独立して保存する。

## `zstarview-export-image`

- headless で動作し、1 枚の画像を生成したら終了する。
- `zstarview` と共通の引数定義を使える範囲では共有する。
- 出力先、画像サイズ、タイムアウト、部分データ許可、sixel 出力など、画像出力専用オプションを追加する。

## 共通の CLI と解決ルール

- `--place`、`--place-countrycode`、`--place-lang` はオンラインの place 検索を使う。
- `--search` と `--search-keep-marker` は GUI の検索と export-image の検索で共有する。
- `--list` は export-image 専用とする。
- tower / mountain viewpoint のデータセット問い合わせは即時終了パスとして扱い、GUI 本体へは入らない。

## 利用者向けのアプリケーション像

文書上は、この製品を 1 つの単一ランチャーではなく、次の 3 系統として扱う。

- CLI 駆動の GUI ランチャー
- ダイアログ先行の GUI ランチャー
- headless の画像出力器

この切り分けは、インストール説明、利用例、エラー処理、設定永続化の説明に直接関係する。
