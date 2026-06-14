# オーバーレイ

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `-c`, `--cloud-opacity CLOUD_OPACITY` | 雲の不透明度を指定します（0.0〜1.0）。0.0 で、そのセッション中の雲描画を無効化します。`--geo-satellite true` を有効にしていても同様です。※2 | `0.07` |
| `--cloud-stripe MODE[,COUNT[,WIDTH]]` | 雲ストライプの方式を指定します。`width` は中心対称のストライプを描き、雲量に応じて見かけの線幅を変えます。`alpha` は線幅を固定したまま線の alpha を変えます。`COUNT` は既定の 600x600 星レンダリング面でのストライプ密度として扱い、実際の描画時には星レイヤーの縮小レンダリング面サイズに合わせてスケールします。`width` は `width,50,0.85`、`alpha` は `alpha,50,0.25` に展開されます。count または width を `0` にすると雲描画を無効化します。 | `width,50,0.85` |
| `--cloud-missing-tint-opacity OPACITY` | 雲欠損領域を示す黄色の濃さを指定します（0.0〜1.0）。 | `0.176` |
| `--night-light-opacity OPACITY` | NASA 夜間光オーバーレイのうち、街灯部分の不透明度を指定します（0.0〜1.0）。0.0 で、起動中の Black Marble のダウンロードと街灯描画を無効化します。 | `0.12` |
| `--ridge-glow-opacity OPACITY` | 地形地平線の上に描く ridge glow の不透明度を指定します（0.0〜1.0）。0.0 で、街灯レイヤーとは独立して ridge glow を無効化します。 | `0.5` |
| `--water-surface-opacity OPACITY` | 水面ドットの不透明度を指定します（0.0〜1.0）。0.0 で、起動中の水面データの取得とドット描画を無効化します。※5 | `0.4` |
| `-a`, `--aircraft-opacity OPACITY` | 航空機オーバーレイの不透明度を指定します（0.0〜1.0）。0.0 で、起動中の航空機問い合わせと描画を無効化します。 | `0.5` |
| `--satellite-opacity OPACITY` | 人工衛星オーバーレイの不透明度を指定します（0.0〜1.0）。0.0 で、起動中の軌道要素取得と描画を無効化します。 | `0.5` |
| `--tropical-cyclone-opacity OPACITY` | 台風・サイクロンオーバーレイの不透明度を指定します（0.0〜1.0）。0.0 で、台風 API の取得と描画を無効化します。 | `0.7` |
| `--show-guidelines-initial true\|false` | 起動時にガイドライン表示を有効にするかを指定します。対象は幾何学的地平線、天の赤道、黄道、never-rises 円、方位ラベル、天頂マーカー、天の極マーカーです。 | `show` |
| `-d`, `--terrain-horizon-opacity OPACITY` | 地形地平線ポリラインの不透明度を指定します（0.0〜1.0）。0.0 で DEM ダウンロード・地形地平線計算・描画・地球ガイドを無効化します。※4 | `0.003` |
| `-e`, `--earth-guide-opacity OPACITY` | 地平線下の地球ガイド線レイヤーの不透明度を指定します（0.0〜1.0）。0.0 で、起動中の地球ガイド描画を無効化します。※4 | `0.028` |
| `--ground-tint-opacity OPACITY` | 幾何学的地平線または地形地平線より下の地面色塗りの強さを指定します（0.0〜1.0）。 | `0.1` |
| `-u`, `--urban-outline-opacity OPACITY` | 都市アウトライン重ね表示の不透明度を指定します（0.0〜1.0）。0.0 で起動中は表示を無効化します。 | `0.2` |
| `-r`, `--urban-outline-radius-km RADIUS_KM` | 観測地点からこの半径内の建物を都市アウトラインとして取得・描画します。 | `2.5` |
| `--urban-outline-skyscraper-radius-km RADIUS_KM` | 遠距離の高層建築補助レイヤーの外側半径です。`0` を指定すると起動中は skyscraper tile 探索を無効化します。それ以外の値は `--urban-outline-radius-km` 以上でなければなりません。 | `60.0` |
| `-b`, `--urban-outline-min-building-height-m METERS` | 非推奨。都市アウトラインの取得/キャッシュ時にこの高さ未満の建物を除外します。性能調整は `--urban-outline-max-candidates` を使ってください。 | `0.0` |
| `--urban-outline-max-candidates N` | 重いリングサンプリングの前に都市アウトライン候補を最大 `N` 件まで残します。都市アウトライン間引きの主な性能調整用オプションです。`0` でレイヤーを実質無効化できます。 | `5000` |
| `--urban-outline-feature-type {both,building}` | 都市アウトラインに含まれるデータのうち、表示に使うものを選びます。`both` は `building` と `building_part` を組み合わせ、part がある場合はそちらを優先します。 | `both` |

<a id="about-water-surface"></a>

#### 脚注

※2 雲の描画は気象衛星（**Himawari** / **NOAA GOES**）の赤外線データを公開 S3 バケットから取得して行います。ネットワーク関連の注意や回避策は「トラブルシューティング」を参照してください。Geo-satellite を有効にしていても、`-c 0` はユーザーが手動で再有効化するまで雲描画を無効のままにします。

※5 水面の描画は 2 つの経路を使います。海側のドットは [OSM Water Polygons](https://osmdata.openstreetmap.de/data/water-polygons.html) 由来のローカル sea-mask タイルから、川・湖・池などの内陸水域のドットは [Overpass API](https://overpass-api.de/) 経由で取得した OpenStreetMap データから生成します。海岸沿いの少し高い観測地点では、1 回の Overpass 問い合わせで 9MB 前後の内陸水域データを取得することがあります。取得した生データは観測地点ごとにキャッシュされますが、観測地点を何度も変えて実行すると、ダウンロード量が積み上がり、公開 Overpass インスタンスの目安である 1GB/日の安全ガイドラインに近づくことがあります。

※4 地形地平線表示は初回利用時に Copernicus DEM タイルをダウンロードし、以後はローカルキャッシュを再利用します。有効時はディスク内の地面/空の塗り分け境界にも地形プロファイルを使い、地球ガイドも同じ地面トーンの色で描画されます。
