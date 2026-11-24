# Handoff: PLATEAU キャッシュの現行の問題点と改修方針

作成日: 2026-11-25。対象: `zstarview-download-plateau-buildings` と PLATEAU 派生キャッシュの
読込。LOD2 壁面の利用とは独立した改修として扱う（壁面の利用は、データ量の見積りが
ないため一時的に見送った。`plateau-lod2-wall-cache-design-2026-11-25.md` を参照）。

主なコード:
- `src/zstarview/data/plateau_buildings.py`: 準備CLI（`_prepare_city_code`、
  `_complete_cache_metadata`、`_cache_matches_catalog`、`_list_plateau_caches`、
  `build_arg_parser`）
- `src/zstarview/data/plateau_building_cache.py`: `is_valid_plateau_cache`、
  `find_plateau_building_dirs`
- `src/zstarview/data/building_source.py`: PLATEAU/Overture の選択

## 問題点（2026-11-25 時点の実装で確認）

### 1. 読込対象のディレクトリを名前で絞っていない（重複読込）

`find_plateau_building_dirs` と `_list_plateau_caches` は、root 直下の全ディレクトリを
`is_valid_plateau_cache`（メタデータのスキーマ・`status`・索引の有無）に掛けるだけで、
ディレクトリ名を見ない。次のものが有効なキャッシュとして同時に読まれ得る。

- `.outdated-*`／`.replaced-*`: 中身が現行スキーマ3なら有効。手元環境の
  `.outdated-*` が今読まれていないのは、スキーマ1・2で不一致になっているためにすぎない
  （例: `13103_2025.outdated-20260715142530` はスキーマ1、`...145855` はスキーマ2）。
  今後スキーマ3のキャッシュを更新すると、退避された旧版が新版と重複して読まれる。
- 同じ自治体の別年度（例: `13103_2024` と `13103_2025`）。
- `--input-zip` を `--year` なしで使ったときに作られる `<code>_latest`（問題点4）。

重複すると、同じ建物の屋根線・壁が二重に描かれ、候補数上限も無駄に消費する。

### 2. カタログ照合で退避済みのメタデータを拾う（不要な再ダウンロード）

`_complete_cache_metadata` は `output_root.glob(f"{city_code}_*")` の最初の完成済み
メタデータを返す。この glob は `<code>_<year>.outdated-*` にも一致し、順序は
ファイルシステム依存である。手元環境で実行した結果:

```text
13101 -> derived_tile_schema_version 1   (.outdated-* を返した)
13103 -> 3
13110 -> 1                                (.outdated-* を返した)
32201 -> 3
```

スキーマ1のメタデータは `_cache_matches_catalog` で必ず不一致になるため、キャッシュが
最新でも「更新が必要」と判定され、再ダウンロード（数GB）が走る。

### 3. カタログ更新時、新版の完成前に旧版を退避する

`_prepare_city_code` のカタログ経路は、既存の完成済みキャッシュを staging の生成
**前**に `<dir>.outdated-<時刻>` へ移す。取得・変換・容量不足・中断で失敗すると、
正式な名前のディレクトリがなくなる。

現状は問題点1のおかげで `.outdated-*` がそのまま読まれ、表面化しにくい。問題点1だけを
直すと、失敗時にその自治体の PLATEAU データが消えたように見える。**問題点1と3は同時に
直す必要がある。**

一方、`--overwrite` 経路は staging 完成後に旧版を `.replaced-*` へ移しており、順序は
正しい。

### 4. `--year` が機能しておらず、`<code>_latest` を作る

- ヘルプは「PLATEAU preparation year or latest」だが、カタログ URL
  （`.../citygml/{city_code}?types=bldg`）は年度を含まず、常に最新版を返す。
  年度を選ぶ機能はない。
- カタログ経路では、カタログに `year` がない場合の代替値にのみ使われる。手元の
  24自治体はすべてカタログの年度で作られていた。
- `--input-zip` では `--year` の値がそのままディレクトリ名になる。既定値 `latest` の
  ままだと `<code>_latest` ができ、問題点1により `<code>_2025` 等と重複する。
- 年度入りの URL 雛形 `DEFAULT_API_URL`・`build_download_url`・`--api-url` は、
  テスト（`tests/data/test_plateau_buildings.py`）以外から呼ばれていない。
- README・`docs/` に `--year` の記載はない。

### 5. `--input-zip` 周りの不整合

- 既存の完成済みキャッシュがあると、`--overwrite` なしでは `FileExistsError` で止まる。
  `--overwrite` を付けると旧版は `.replaced-*` に残る（問題点1の対象）。
- `--keep-zip` はダウンロード時にしか ZIP をコピーしない。`--input-zip` で
  キャッシュ内の `source-citygml.zip` を指定して再変換すると、新版に ZIP が残らない。
- 元 ZIP を識別する情報がない（`cache_meta.json` は `archive_url` と
  `archive_size_bytes` のみ。`--input-zip` ではどちらも `None`）。保存済み ZIP が
  既存キャッシュと同じ元データかを判定できない。

### 6. 退避ディレクトリが溜まり続ける

`.outdated-*`／`.replaced-*`／`.incomplete-*` は削除されない。手元環境では東京23区と
松江市のほぼ全自治体に `.outdated-*` が1〜2個ずつ残っており、それぞれ元キャッシュと
同程度の容量を使う（港区の `bldg/` は約176 MB）。削除の手段や案内もない。

### 参考: 屋根面の抽出が LOD を区別しない

`_roof_surfaces` は `RoofSurface` 以下の全 `posList` を走査し、LOD2/LOD3・ポリゴン・
穴を区別しない。港区2025年の `53393599` タイルでは、11棟が LOD2 と LOD3 の両方の
境界面を持つ。キャッシュ管理とは別の問題であり、詳細屋根の40 m閾値の撤廃
（決定済み）と合わせて影響を評価する。

## 決定済みの改修方針

利用者と合意済み。CLI の非互換変更は 2.2.0 で出す。

- **自治体ごとに1世代**: キャッシュを `<5桁の自治体コード>/` に置き、名前に年度を
  入れない。新形式どうしの重複は構造的に起きない（問題点1・4）。
- **元 ZIP のチェックサム**: `cache_meta.json` に `source_zip_sha256` と
  `source_zip_size_bytes` を記録する。ダウンロード中のストリーム、または ZIP 指定時の
  読込時に計算する。カタログ由来の年度等は情報として記録するが、ディレクトリ名や
  更新判定には使わない（問題点5）。
- **差し替えと削除**: 新しいデータを準備したら、年度は比べずに置き換える。
  1. `.<code>-<unique-id>/` の staging で生成し、読み戻して検証する。失敗時は staging を
     消し、旧版をそのまま残す。
  2. 既存の `<code>/` を `.` 始まりの名前に退避し、staging を `<code>/` に移す。
  3. 退避した旧版と、同じ自治体の旧形式（`<code>_<year>/` とその `.outdated-*`／
     `.replaced-*`／`.incomplete-*`）を削除し、削除したパスを表示する。他の自治体には
     触れない。
  4. 2と3の間で中断して残った `.` 始まりのディレクトリは読まず、次回の準備時に消す。

  `.outdated-*`／`.replaced-*` へのバックアップは廃止する（問題点3・6）。
- **旧形式の読込**: 新形式 `<code>/` がない自治体に限り、旧形式
  `<code>_<4桁年度>/` のうち年度が最大のもの1つを読み、他はログに出す。それ以外
  （`.` 始まり、`.outdated-*`、`.replaced-*`、`.incomplete-*`、`_latest` 等）は読まない
  （問題点1）。カタログ照合で既存キャッシュを探す処理も同じ規則で選ぶ（問題点2）。
- **準備CLIの作り直し**: `--year`、`--api-url`、`--overwrite` は引き継がない。
  `DEFAULT_API_URL` と `build_download_url` は削除する。
  - カタログ経由: カタログが既存の記録と一致すればダウンロードせずに終了する。
    それ以外は推定ダウンロード量を表示し、確認のうえ取得・置換する。
  - ZIP 指定: チェックサムが既存と一致すれば、そのまま再変換する。不一致や既存なしの
    場合は、置き換える旨を表示して確認する。
  - 確認は対話の y/N と、非対話用の確認省略オプションに統一する。
  - 元 ZIP を残すオプションは、ダウンロード・ZIP 指定のどちらでも ZIP を staging に
    コピーする。
  - 一覧と削除の操作を持たせる。一覧では、選ばれなかった旧形式や `.` 始まりの残骸を
    「読込対象外・削除可能」として示す。
- 元 ZIP のハッシュを使わず年度を ZIP 内のメタデータから読み取る案は、命名の揺れと
  別年度の副メタデータの同居（つくば市2023年）があり、規則が複雑になるため不採用。
  調査結果は `dev-notes/session-2026-11-25.md` を参照。

## 未決定（仕様で決めること）

1. コマンド名とサブコマンド（準備・一覧・削除）、確認省略オプションの名前。
   現行の補助オプション（`--temp-dir`、`--catalog-url`、`--max-download-bytes`、
   `--min-building-height-m`、`--jsonl`）の取捨。複数自治体の一括準備の可否。
2. 準備成功時の旧形式・残骸の削除を、準備の確認1回にまとめるか、別に確認するか。
3. 旧形式が残っている自治体を起動時に案内するか（頻度と表示先）。

## 検証事項

1. 新形式と旧形式の併存、旧形式の複数年度、`.outdated-*`、`.replaced-*`、
   `.incomplete-*`、`_latest`、`.` 始まりの残骸がある状態で、自治体ごとに1つだけ
   読み込まれること（`find_plateau_building_dirs` と一覧の両方）。
2. カタログ照合が、退避済みディレクトリのメタデータを拾わないこと（問題点2の回帰試験）。
3. カタログ経由の新規準備・置換・スキップ。置換時の旧形式と残骸の削除。
4. ZIP 指定で、チェックサム一致時の再変換、不一致・既存なし時の確認、元 ZIP 保持の
   有無、キャッシュ内の `source-citygml.zip` からの再変換。
5. 取得・変換の失敗、中断（置換手順の各段階）、ディスク不足、置換中の読込で、旧版が
   失われないこと。Windows でのディレクトリ移動。
6. 公開時に `release-notes.md`（2.2.0: 準備CLIの非互換と移行手順）、
   `docs/specification.md`（準備CLI）、`docs/design.md`（キャッシュ構成）を更新する。
