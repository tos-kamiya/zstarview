# 衛星雲の夜間照明を画像で比較する

`scripts/prototype_cloud_cell_shading.py` は、実際の衛星雲データと
zstarview の voxel 描画処理を使う比較用スクリプトです。月の位置・月相は
指定地点・時刻の天体計算から、地表光は既存の夜間光タイルから取得します。

## 4通りをまとめて比較する

リポジトリのルートで、開発環境の Python を使って実行します。

```bash
.venv/bin/python scripts/prototype_cloud_cell_shading.py Tokyo \
  --cloud-only --compare-night-lighting \
  --view-alt 45 --view-az 180 --size 800x600 \
  --output /tmp/cloud-night-comparison.png
```

`--datetime` を省略すると実行時の時刻を使います。地点は都市名のほか、
`@35.68,139.69` のような緯度・経度も指定できます。

以下の画像と、それぞれに対応する `.json` を保存します。

| ファイル末尾 | 照明 |
| --- | --- |
| `-ambient.png` | 弱い環境光のみ |
| `-moon.png` | 環境光＋月光 |
| `-ground.png` | 環境光＋直下の地表光 |
| `-combined.png` | 環境光＋月光＋直下の地表光 |
| 指定した出力名 | ラベル付きの2×2比較画像 |

並び順は左上が環境光、右上が月光、左下が地表光、右下が両方です。
すべて同じ衛星データ、背景、地点、時刻、視点、雲量補正、不透明度を使います。
地表光のサンプリング結果も比較中は使い回します。比較画像に並べる際に
画素の明るさや画像サイズを変更しません。

## 1枚だけ生成する

```bash
.venv/bin/python scripts/prototype_cloud_cell_shading.py Tokyo \
  --cloud-only --night-lighting combined \
  --output /tmp/cloud-night.png
```

`--night-lighting` は `ambient`、`moon`、`ground`、`combined` から選びます。
省略時は `combined` です。比較モードとは併用しません。

`--cloud-only` は黒い背景に雲だけを描きます。代わりに `--sky-color-only`
を指定すると空色を背景にします。両方省略すると従来どおり空・星・天体の
出力画像を背景に使います。背景の月が画面外でも、照明には月を使います。

`--datetime '2026-10-10T12:00:00Z'` のように UTC の ISO 時刻を指定できます。
衛星雲が取得可能な時刻を選んでください。衛星の観測時刻は描画時刻と
異なることがあるため、両方を JSON に記録します。

## 月の照明条件だけを上書きする

新月付近などで月光が弱いときは、仮想的な月を照明に使えます。
以下は現在の衛星雲に、高度45度・南の方位・満月の照明を加える例です。

```bash
.venv/bin/python scripts/prototype_cloud_cell_shading.py Tokyo \
  --cloud-only --compare-night-lighting \
  --moon-alt 45 --moon-az 180 --moon-phase 180 \
  --size 800x600 --output /tmp/cloud-full-moon-comparison.png
```

| オプション | 意味・範囲 |
| --- | --- |
| `--moon-alt DEG` | 月の高度。−90〜90度。地平線以下では月光を使いません。 |
| `--moon-az DEG` | 月の方位。0〜360度。北0、東90、南180、西270。360は0として扱います。 |
| `--moon-phase DEG` | 月と太陽の離角。0は新月、90は半月、180は満月です。 |

各項目は個別にも指定でき、省略した項目には実際の月の値を使います。
単独出力の `--night-lighting moon` や `combined` でも利用できます。
比較では月光を含む2枚だけに上書きした月の条件が反映されます。
太陽、地表光、雲の取得時刻、背景に描かれる月の位置・見た目は変えません。
JSON の `moon` に照明に使う値、`moon.actual` に実際の月の値、
`moon.override_applied` に上書きの有無を記録します。
昼夜の切り替えは実際の太陽高度に従うため、昼間は仮想月を指定しても月光は出ません。
これらのオプションは voxel 方式の比較専用で、アプリ本体の設定には影響しません。

## 比較の読み方と制約

- 太陽高度が−4度以上では夜間照明を使わず、4枚が同じになることがあります。
  −4〜−9度では昼光と夜間照明を滑らかに混ぜます。
- 月が地平線以下、または新月に近い場合は月光の差が小さくなります。
  端末と JSON に月の高度・太陽との離角を記録します。
- 地表光のデータ取得に失敗した領域は寄与が0になります。警告を端末へ出力します。
  JSON の `grid.ground_light_max` が0でも、元の地表が暗い場合と欠損の場合の両方があります。
- `ambient` は現在の弱い環境光を使います。夜間照明を追加する前の
  雲量依存の青灰色を再現するモードではありません。
- 雲量から引く値の既定値は、現在のアプリと同じ `0.09` です。
  `--cloud-amount-subtract` で変更できますが、比較する画像間では同じ値を使います。
- 照明計算はアプリと共有します。試作スクリプトは指定サイズの画素を直接描画し、
  cutout線は加えません。不透明度は視線方向の減衰に適用する従来の試作仕様です。
  アプリの513×513雲画像と描画時の不透明度適用を画素単位で再現するものではありません。
- 対象は native CloudDisc（Himawari/GOES）経路です。Europe の Geo-satellite 入力は
  このスクリプトの対象外です。`--cloud-model altaz` の旧方式では夜間照明比較を使えません。
- 衛星データ・天体暦・夜間光タイルが未取得の場合はネットワークアクセスが必要です。
  `--timeout` は衛星データ取得workerの期限であり、スクリプト全体の期限ではありません。
