# 地面

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `-d`, `--terrain-horizon-opacity OPACITY` | 地形地平線ポリラインの不透明度を指定します（0.0〜1.0）。0.0 で DEM ダウンロード・地形地平線計算・描画・地球ガイドを無効化します。※4 | `0.003` |
| `-e`, `--earth-guide-opacity OPACITY` | 地平線下の地球ガイド線レイヤーの不透明度を指定します（0.0〜1.0）。0.0 で地球ガイド描画を無効化します。※4 | `0.028` |
| `--water-surface-opacity OPACITY` | 水面ドットの不透明度を指定します（0.0〜1.0）。0.0 で水面データの取得と描画を無効化します。※5 | `0.4` |
| `--night-light-opacity OPACITY` | NASA 夜間光の街灯部分の不透明度を指定します（0.0〜1.0）。0.0 で Black Marble のダウンロードと街灯描画を無効化します。 | `0.22` |
| `--ridge-glow-opacity OPACITY` | 夜間光プロファイル由来の ridge glow レイヤーの不透明度を指定します（0.0〜1.0）。0.0 で無効化します。 | `0.08` |
| `--road-light-opacity OPACITY` | OSM道路光レイヤーの不透明度を指定します（0.0〜1.0）。0.0 で無効化します。 | `0.12` |
| `--road-light-max-candidates N` | 重い処理の前に道路 `way` 候補を最大 `N` 件まで残します。0 でレイヤーを無効化します。 | `5000` |
| `-u`, `--urban-outline-opacity OPACITY` | 都市アウトラインの不透明度を指定します（0.0〜1.0）。0.0 で無効化します。 | `0.2` |
| `-r`, `--urban-outline-radius-km RADIUS_KM` | 観測地点からこの半径内の建物を都市アウトラインとして取得・描画します。 | `2.5` |
| `--urban-outline-skyscraper-radius-km RADIUS_KM` | 遠距離の高層建築補助レイヤーの外側半径です。0 で探索を無効化し、それ以外は `--urban-outline-radius-km` 以上にします。 | `60.0` |
| `--urban-outline-max-candidates N` | 重いリングサンプリングの前に都市アウトライン候補を最大 `N` 件まで残します。0 でレイヤーを無効化します。 | `5000` |
| `--urban-outline-skyscraper-only` | 通常の近距離都市アウトラインを省略し、遠距離の高層建築補助レイヤーだけを描画します。 | 無効 |
| `--urban-outline-download-timeout-seconds SECONDS` | Overture都市アウトラインの各ダウンロードサブプロセスを待つ最大時間を指定します。 | `120.0` |

<a id="about-water-surface"></a>

#### 注釈

※4 地形地平線表示は初回利用時に Copernicus DEM タイルをダウンロードし、以後はローカルキャッシュを再利用します。地球ガイドは独立したレイヤーです。

※5 水面は海側のローカル sea-mask タイルと、Overpass API 経由の OpenStreetMap 内陸水域データを使います。取得データは観測地点ごとにキャッシュされます。
