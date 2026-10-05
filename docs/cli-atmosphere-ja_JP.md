# 大気

## 大気

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `-S`, `--sky-opacity SKY_OPACITY` | 空色ディスクの総合表示強度を指定します（0.0〜1.0）。値が大きいほど、彩度の比例増加を抑えて輝度を優先します。0.0では空色ディスクと明るい天体の暗色下敷きを無効化します。 | `0.3` |
| `--sky-disc-altaz-rings {off,dimalt,altaz}` | 常時表示の空ディスク方位/高度オーバーレイです。`dimalt` は控えめな高度リング、`altaz` はフルグリッドを表示します。 | `dimalt` |
| `--sky-disc-altaz-rings-hover {off,dimalt,altaz}` | ホバー時の空ディスク方位/高度オーバーレイです。意味は上記と同じです。 | `altaz` |
| `-c`, `--cloud-opacity CLOUD_OPACITY` | 雲の不透明度を指定します（0.0〜1.0）。0.0 で、そのセッション中の雲描画を無効化します。`--geo-satellite true` を有効にしていても同様です。夜間は雲の視認性を保つため、太陽高度に応じて実効値を最大30%まで滑らかに持ち上げます。※2 | `0.3` |
| `--cloud-mode {voxel,shell}` | 3次元ボクセル内で散乱・透過を計算する方式、または高度別の雲シェルを画面へ投影する方式を選びます。shell 方式では`--cloud-stripe`で見た目を選びます。 | `voxel` |
| `--geo-satellite true\|false` | 対応する Europe workflow band 内で、実験中の MET Norway 赤外画像を使います。voxel と shell の両モードに対応し、voxel では表示画像から1〜9 kmへの配分を推定します。 | `false` |
| `--cloud-stripe MODE[,COUNT[,WIDTH]]` | 既定値は`--cloud-mode`に応じて変わり、voxelでは`cutout,16`、shellでは`halftone2,30,1.7`です。shell 方式では雲の見た目を選びます。`halftone2` はシェル別ドット、`halftone` は以前の全シェル集約表現です。`width` は雲量に応じて線幅を連続的に変え、`width-quantized` は5段階で線幅を変え、`alpha` は線幅を固定して alpha を変えます。voxel 方式では`--cloud-stripe cutout[,COUNT]`で疎な透明線の本数を指定できます。線は右下がり45度で、線幅は513x513の雲画像上で1pxです。`COUNT=0`はcutoutだけを無効にします。cutoutはvoxel専用です。 | voxel: `cutout,16`; shell: `halftone2,30,1.7` |
| `--cloud-missing-tint-opacity OPACITY` | 雲欠損領域を示す黄色の濃さを指定します（0.0〜1.0）。 | `0.176` |
| `-P`, `--precipitation-opacity OPACITY` | 任意で有効化するOpen-Meteoモデル予報降水の雨線について、不透明度を指定します（0.0〜1.0）。現在時刻に最も近い15分予報区間の中心を選択します。ネイティブな15分モデルの対象外地域では、時間予報から補間される場合があります。正の値を指定する場合は、非商用Free API利用規約への初回同意が必要です。 | `0.0` |
| `--tropical-cyclone-opacity OPACITY` | 台風・サイクロンオーバーレイの不透明度を指定します（0.0〜1.0）。0.0 で、台風 API の取得と描画を無効化します。時刻をずらした表示では自動的に非表示になります。 | `0.7` |
| `-a`, `--aircraft-opacity OPACITY` | 航空機オーバーレイの不透明度を指定します（0.0〜1.0）。0.0 で、起動中の航空機問い合わせと描画を無効化します。 | `0.0` |
| `--satellite-opacity OPACITY` | 人工衛星オーバーレイの不透明度を指定します（0.0〜1.0）。0.0 で、起動中の軌道要素取得と描画を無効化します。 | `0.7` |
| `--meteor-trails-opacity OPACITY` | GMNメテオ軌跡の不透明度を指定します（0.0〜1.0）。0.0ではその起動中の取得・描画・メニューからの再有効化を無効にします。 | `0.5` |
| `--meteor-trails-max-candidates N` | 地理的選別後に表示するGMNメテオ軌跡を、新しいものから最大 `N` 件に制限します。`0` で無制限にします。 | `150` |
#### 脚注

※2 雲の描画は気象衛星（**Himawari** / **NOAA GOES**）の赤外線データを公開 S3 バケットから取得して行います。ネットワーク関連の注意や回避策は「トラブルシューティング」を参照してください。Geo-satellite を有効にしていても、`-c 0` はユーザーが手動で再有効化するまで雲描画を無効のままにします。
