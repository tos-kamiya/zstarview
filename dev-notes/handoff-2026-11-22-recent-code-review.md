# 直近の描画と CLI のコードレビュー引き継ぎ

2026-11-22。対象は `88536ff6..b663875a` のソース差分（1.64.1〜1.64.3）。雲のマーク、地球ガイドの都市光、関連する CLI と描画経路を確認した。以下は実装候補で、製品コードへの適用は未着手。

## 1. 雲の斜線マスクを対角線ごとに計算する

`apply_cloud_voxel_cutout()` は、斜線ごとに画像全体で距離とマスクを計算する。判定は中心補正後の `x - y` のみに依存するため、`height + width - 1` 個の対角線について一次元で判定し、その結果を画素へ展開できる。変更先は [cloud_voxel_style.py](../src/zstarview/cloud_voxel_style.py#L125)。

使い捨ての試作では、513×513 の不透明な uint8 RGBA 画像に既定の 16 本を描く処理が **8.25 ms から 1.19 ms** になった。入力生成を計測区間から外し、各実装を 30 回ずつ 3 セット実行して、最短セットの 1 回あたりの時間を比較した。これは `cutout` 単体の値で、既定の `holes` やアプリ全体の描画時間への効果は計測していない。

試作は、座標差を `np.arange(1 - height, width, dtype=np.float32) - (width - height) / 2` で列挙し、既存と同じオフセットと線幅で判定する。画素 `(x, y)` の表の添字は `x - y + height - 1`。現行実装と画素単位で一致した条件は、高さ×幅が 12×9、9×12、513×513、1025×1025、513×512、3×1、本数がそれぞれ 1、16、100。実装時は本数 0 の早期リターン、入力の非破壊、RGB の維持も確認する。

## 2. 地球ガイドで太陽高度の変換結果を共有する

地表の明るさと都市光の重みは、どちらも同じ `sun_strength` 配列をクリップし、`arcsin` と度数変換を適用する。塗りつぶし経路で度単位の太陽高度を一度だけ求めて共有すれば、変換と中間配列を減らせる。既存関数の入力契約と、都市光が −8 度から −14 度にかけて現れる挙動を維持する。昼・薄明・夜の描画結果を比較し、描画時間への効果を測る。変更先は [earth_guide.py](../src/zstarview/render/earth_guide.py#L438)。

## 3. 雲のマーク設定を CLI パーサーで型にする

`_parse_cloud_stripe()` が返すタプルの第 3 要素は常に `1.0` で、`CloudVoxelMarkStyle.from_cli()` が破棄する。パーサーと既定値を `CloudVoxelMarkStyle` に揃え、GUI・画像出力・リプレイへの受け渡しを確認する。公開書式の `holes[,RADIUS]` と `cutout[,COUNT]`、値 0 による無効化、保存済みプロファイルの読み込みを維持する。関連箇所は [args.py](../src/zstarview/cli/args.py#L104) と [cloud_voxel_style.py](../src/zstarview/cloud_voxel_style.py#L32)。

## レビュー時の検証

リポジトリのルートで以下を実行し、135 件すべて成功した。

```bash
QT_QPA_PLATFORM=offscreen .venv/bin/python -m pytest -q \
  tests/test_cloud_voxel_style.py \
  tests/test_earth_guide_render.py \
  tests/test_cli_cloud_and_overlay_args.py \
  tests/test_sky_disc_viewport_and_ground_tint.py
```
