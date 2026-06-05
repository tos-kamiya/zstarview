# 起動時の対象検索

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `--search QUERY` | 起動時に対象を解決します。GUI の `対象検索...` と同じ local-first 検索を使い、ISS はアプリ側の現在位置を優先し、ISS として認識されたのに現在位置を取得できない場合はエラーにして JPL へはフォールバックしません。JPL-backed spacecraft は JPL 経路を使い、CLI / export-image 側の JPL 解決は major body / small body を含めてもよいです。`=` がなければラベルまたは ID で検索します（例: `Ceres`, `2000001`）。 | |
| `--search label=QUERY` | ラベルだけを検索します（例: `label=Ceres`）。 | |
| `--search id=QUERY` | ID だけを検索します（例: `id=2000001`）。 | |
| `--search-keep-marker` | 選択対象をマーカーとラベル付きで継続表示します。 | |
| `--list` | `zstarview-export-image` 専用です。候補を表示して終了します。 | |

`-A` または `-Z` が指定されている場合は、その軸を固定し、未指定側だけを検索結果で補います。検索結果の高度角は `-5°` にクランプされます。
