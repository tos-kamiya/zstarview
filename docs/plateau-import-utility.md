# PLATEAU ZIP 取込ユーティリティ

最終更新: 2026-03-14

## 1. この文書の位置づけ

この文書は、`zstarview-import-plateau-zip` の利用方法をまとめたものである。  
`zstarview` 本体の機能仕様は `docs/specification.md`、内部設計は `docs/design.md` を参照する。

## 2. 目的

`zstarview-import-plateau-zip` は、利用者がローカルにダウンロードした PLATEAU CityGML ZIP から、
urban outline overlay 用の derived building tiles を生成して、実行中の `zstarview` パッケージ配下へ書き込むための補助 CLI である。

## 3. 入力と出力

### 3.1 形式

```bash
zstarview-import-plateau-zip ZIP_PATH [options]
```

### 3.2 入力

- `ZIP_PATH`
  - PLATEAU CityGML ZIP ファイル 1 個
- `--workers N`
  - `udx/bldg/*.gml` の変換を並列化する worker 数
  - 実データは大きいことが多いため、通常は指定してよい
- `--min-building-height-m METERS`
  - derived tile 生成時に保持する最低建物高さ
  - 既定値は `40`
  - `0` を指定した場合は高さによる除外を行わない

### 3.3 出力

- 実行中の `zstarview` パッケージ配下 `zstarview/data/plateau_derived/<city>/bldg`
  - derived tile JSON
  - `tile_index.json`

## 4. 挙動

- ZIP を一時展開し、その中の `udx/bldg` を検出して変換する
- 1 回の呼び出しで処理する ZIP は 1 個だけとする
- import 時に `--min-building-height-m` で残した建物を、runtime 側の都市アウトライン描画では高さ再フィルタせず利用する

## 5. 使用例

通常は `--workers` を付けて使う。

```bash
zstarview-import-plateau-zip \
  ~/Downloads/12100_chiba-shi_city_2024_citygml_1_op.zip \
  --workers 8
```

低層建物も残したい場合:

```bash
zstarview-import-plateau-zip \
  ~/Downloads/32201_matsue-shi_city_2024_citygml_1_op.zip \
  --workers 8 \
  --min-building-height-m 0
```

## 6. エラー条件

次の場合はエラー終了する。

- ZIP が存在しない
- ZIP 名から自治体コードを推定できない
- ZIP 展開結果から `udx/bldg` を検出できない
