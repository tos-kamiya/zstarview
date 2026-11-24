# zstarview 設計書

最終更新: 2026-10-10

この文書は、`zstarview` の内部設計の入口である。
`docs/design/` 以下に、責務ごとに分割した詳細文書を置く。

## 実行環境診断

通常のGUI起動時は、起動ログへASCII JSON形式のruntime diagnosticsを1回記録
する。診断項目はapp version/revision、PID・親PID・session ID、Python実装・版・
実行ファイル、GIL状態、OS・アーキテクチャ、主要依存ライブラリのdistribution
版である。NumPyやAstropy等の重いモジュールを診断のためだけにimportせず、
`importlib.metadata`から版を取得する。

同じ情報は `zstarview-diagnose-runtime` で手動採取できる。通常出力は端末向け
key-value、`--json`は報告書向けJSONとし、どちらも非ASCII文字を直接出力しない。
Git管理外の配布物ではrevisionを`unknown`とし、`ZSTARVIEW_APP_REVISION`で運用時に
明示できる。診断情報は原因特定や無クラッシュの証明ではなく、実行条件を比較
可能にするための記録として扱う。

## GUIの循環GCと寿命管理

`gui.viewer.main()` は `GuiGarbageCollector` の実行範囲でGUIを起動する。
GUI動作中はPythonの自動循環GCを停止し、参照カウントによる解放は維持する。
アプリケーションを親とする1秒間隔のQTimerで、GUI所有スレッドから若い世代を
回収する。割り当て圧が元のGCしきい値に達した場合を対象とし、全世代も30秒間隔で
回収する。開始時にGCが無効だった場合は手動回収も追加しない。

回収時は共有 `native_work_lock` を非ブロッキングで取得し、科学計算が実行中なら
そのtickを見送る。full回収の期限は、実際に回収できた後で更新する。
終了・起動失敗・ダイアログ取消時には共有executorを停止してworkerを待ち、
その後にtimerを停止して元のGC有効状態を復元する。起動処理と資源回収を分けるため、
公開 `main()` が寿命を管理し、従来の起動処理は内部 `_main()` で行う。

これはworkerの割り当てをきっかけにQt関連の循環参照が別スレッドで回収される
経路への対策であり、特定のSIGSEGV原因の確定ではない。第三者コードの明示的な
`gc.collect()` や、GC以前に発生したヒープ破壊を隔離するものでもない。
根拠となる寿命制約は [QtのThreads and QObjects](https://doc.qt.io/qt-6/threads-qobject.html)、
GC制御の動作は [Pythonのgc API](https://docs.python.org/3/library/gc.html) を参照する。

## プロセス分離の移行基盤

GUIクラッシュの根本原因は未確定のため、既存のGUI内計算経路をこの段階では
置き換えない。移行の第一段階として、`zstarview.processes` にQt非依存の
要求・結果schemaと、イベントループから `poll()` するプロセス監督を置く。

監督は1つのプロセスを実行中に保ち、追加要求は最新の1件だけを待機させる。
要求はセッション、worker世代、要求ID、layer/view世代、入力revisionを含む
JSONで作業ディレクトリへatomic renameする。結果manifestは同じ識別情報と
statusを持ち、監督はprotocol version、世代、成果物の作業領域内包含、宣言
サイズを検証する。worker異常終了、期限超過、壊れたmanifest、古い結果は
通常の成功結果と区別する。

`TropicalCycloneController.update()` はこの基盤を使い、取得・キャッシュ読み込み・
snapshot/polygonの辞書化を `python -m zstarview.tropical_cyclones.worker` へ移す。
GUIはmanifest検証後に `payload.json` を読み、snapshot collectionを所有する。
cyclone controllerは共有GUI executorへ処理を委譲しない。worker失敗時にGUI内
経路へ自動フォールバックしない。

成果物はGUIが読み込み後に解放し、解放されない場合もsession単位の件数・容量
上限で古いjobディレクトリを回収する。workerの異常終了とtimeoutには、cyclone
controllerで最大2回、短いbackoff付きの再起動を行う。通常の取得失敗やmanifest
検証失敗は自動再試行せず、UIへ失敗を返す。既存のcloud workerが生成するpickleは、
この新しい汎用IPCの成果物として扱わない。

GMNメテオ軌跡のデータ取得、キャッシュ、座標固定処理は
[meteor-trails.md](design/meteor-trails.md)に記載する。

`zstarview` には次の 3 つのアプリケーション入口がある。

- `zstarview`
  - CLI 引数で初期状態を与えて GUI を起動する
- `zstarview-gui`
  - 起動前ダイアログを先に開いてから GUI を起動する
- `zstarview-export-image`
  - 起動間隔の抑制は `cli/export_image_cooldown.py` がキャッシュ直下に開始時刻と
    起動ごとの識別子をJSONとして保存する。60秒以内の記録では起動を拒否し、処理終了時に
    識別子が変わっていれば重複の可能性を警告する。JSONの置換は一時ファイルを使うが、
    読み取りから書き込みまでは排他しないため、同時起動の防止や確実な検出は保証しない。
    デバッグ用の回避オプションは記録の読み書き自体を行わない。
  - GUI を起動せず、1 枚の画像を headless で書き出す

補助的に `zstarview-debug` という console script 版の GUI ランチャーもある。
これは `zstarview` と同じ `main()` を呼ぶが、Windows では terminal を伴って起動しやすく、起動時ログを見たい診断用途に向く。

補助的に `zstarview-install-overturemaps-exe-cli` という staging 用の console script も置いてよい。
これはダウンロード済みの `overturemaps` Windows 実行ファイルを `CACHE_PATH/overturemaps.exe` にコピーするだけの単純な入口として扱う。
`overturemaps` 本体の探索やダウンロードは行わず、ファイル配置だけを責務にする。

これら 3 つは、地点解決、時刻解釈、描画、キャッシュ、外部データ取得の核心を共有する。
差分は「どの入口から始まるか」「対話 UI を持つか」「1 枚の画像で終わるか」にある。
GUIの一時的な表示モードである`Inverted City`はこの共有核心には含めず、通常の都市データと投影結果をGUI側で表示する直前に変換する。したがって、設定保存、共有scene data、`zstarview-export-image`には反映しない。

通常ビューアのdisplay tone curveは、レイヤー固有のvisibility boostとは分離した
最終表示補償として扱う。完成済みpresent frameとvolatile UI overlayの間へ置き、
補正が有効な場合だけ追加の表示用サーフェイスを生成する。キャリブレーションUI、
LUT、表示用キャッシュの詳細は[rendering-pipeline.md](design/rendering-pipeline.md)と
[gui-screen-update-and-cache.md](design/gui-screen-update-and-cache.md)に記載する。

Landscape viewの注釈フェードはGUIの100ms timerでカーソル位置と最終移動時刻を監視し、
`RenderHudState.landscape_annotation_opacity`として描画pipelineへ渡す。カーソルが
Landscape viewのウィンドウ外、ウィンドウ境界上、画面端、または約0.35秒間静止した時はopacityを
毎秒0.2ずつ下げ、約5秒で非表示にする。カーソルが動き始めるかウィンドウ内へ戻った時は1.0へ即時復帰させる。画面端判定には
カーソル位置上の画面geometryを使い、ウィンドウ内判定にはhostのenter/leave eventを使う。
ウィンドウ外ではグローバルカーソル座標を参照しない。Wayland等で退出後の座標が
更新されない場合にもフェードを進める。モード切替時にはopacityを即時初期化する。
メニューが表示されている間はカーソルがメニューへ移っても注釈を表示状態に保ち、
閉じた後に通常のフェード条件へ戻す。
ガイド、アステリズム、
惑星・月を含むラベル候補とHUD注釈の描画だけにQPainter opacityを適用し、
星空や地形などのシーンレイヤーはフェード対象に含めない。opacityはframe cache keyにも
含め、補間中の描画結果が古い注釈濃度で再利用されないようにする。

### 夜間光データの配布

1.35.0 では、夜間光の元データを EOG の 2025 年次 VIIRS Nighttime Lights
VNL v2.2 `average_masked` とする。元配布物をリリース時の前処理でGeoTIFF
に変換し、PyPIパッケージには含めず、GitHub Release `night-lights-2025`の
リリースアセットとして配布する。実行時の固定入口は次のmanifest URLとする。

`https://github.com/tos-kamiya/zstarview/releases/download/night-lights-2025/manifest.json`

実行時はmanifestを取得してタイル一覧とSHA-256を確定し、観測地点から最大
夜間光距離までのレイが通過しうるタイルだけを取得する。ダウンロードした
ファイルは一時ファイルへ書き込み、SHA-256、GeoTIFFのCRS・解像度・寸法・
範囲を検証した後にキャッシュへatomic renameする。GitHub Pagesやアプリ作者の
Webサイトを実行時の必須依存にはしない。

manifest・各タイル・NOTICEは同じReleaseのAssetsとして管理する。manifestには
データセット識別子、データ年、各Assetの相対ファイル名・URL・SHA-256、元データ
URL・元ファイルのSHA-256、変換スクリプトのGitコミットを記録する。

配布アセットには、データ年・製品名・変換方法・SHA-256 checksum を記録した
manifest または README を添付する。EOG VNL は CC BY 4.0 の条件に従い、
EOG の帰属、GeoTIFF への変換を行った旨、関連する出典論文を配布物と
アプリのクレジットに記載する。

通常の`zstarview`では、都市アウトラインの線と面塗りに淡い青白色`(207, 229, 255)`を使い、道路ランプの暖色と分離したLED照明風の都市光として表現する。

### OvertureMaps exe staging

`import_overture_buildings.py` は `overturemaps` CLI を呼び出して Overture building data を取得する。現行の呼び出し経路は `--overturemaps-bin` を受け取り、`shutil.which(...)` で解決した実行ファイルか、明示されたパスをそのまま使う。

`zstarview-install-overturemaps-exe-cli` は、この経路に対して手元の Windows 版 exe を用意するための補助である。
目的は、GitHub Releases から取得した `v1.0.1` 以上の `*-windows-x86_64.exe` のような資産を、固定名 `overturemaps.exe` として `CACHE_PATH` 配下に staging することにある。
コピー後の実行パスを `--overturemaps-bin` に渡せば、既存の import パイプラインや GUI 側の呼び出しと接続できる。
この helper はリリース資産名の suffix を解釈しないため、バージョン番号やアーキテクチャ名は destination には残らない。
将来もし自動探索を足すなら、`CACHE_PATH/overturemaps.exe` を first-class に見る lookup を別途設計してもよいが、この helper 自体には含めない。

実行ファイルの解決順は、呼び出し側が `--overturemaps-bin` を明示した場合を最優先にし、次に `CACHE_PATH/overturemaps.exe` の staging 済みファイルを見て、最後に `PATH` 上の `overturemaps` を探索する形が自然である。
この順序なら、利用者が明示したパスを壊さず、Arm64 Windows の回避策として staging を使え、従来のインストール形態もそのまま残せる。

### PLATEAU building preparation

都市アウトラインとreplayの建物統合では、標高サンプリングや候補選択の前に、IDと形状・高さ等の全フィールドが一致する建物を先着順で重複排除する。同一IDでも内容が異なるレコードは保持する。欠損IDや旧インポーターの生成ID（`bldg-` / `feature-` と数字）はタイル内でしか一意でないため重複排除せず、親建物のIDとの照合にも使わない。建物部分による親の抑制は重複排除後に行う。

PLATEAUは、通常起動時の暗黙のダウンロードソースではなく、利用者が明示的に準備する任意の高密度ソースとして扱う。
準備入口は `zstarview-download-plateau-buildings` とし、APIからCityGMLを取得し、既存のデータソース非依存なderived building tile形式へ変換して `CACHE_PATH/plateau_buildings` 配下へ保存する。

CityGMLの取得・展開・変換はこの準備CLIプロセスの責務であり、通常のGUI workerや描画スレッドに持ち込まない。
CLIは大容量データを扱うため、取得前に対象データセット、年度、推定サイズ、必要空き容量を確認し、利用者の明示確認後に処理を開始する。カタログからZIPと展開対象ファイルのサイズを取得できる場合は、`shutil.disk_usage()`で選択した一時ファイルシステムの空き容量を確認し、ZIP・展開後ファイル・安全余裕の合計が収まらなければダウンロードを開始しない。ユーザー単位のディスククォータは対象外とし、事前確認後の競合による`ENOSPC`は実行中にも処理する。
ZIPと展開ディレクトリは作業ディレクトリに置き、変換が正常完了した場合だけ完成済みキャッシュへ反映する。失敗時は部分キャッシュを有効化せず、可能なら一時データを削除する。

PLATEAU派生キャッシュの完了判定はディレクトリの存在だけに依存しない。`cache_meta.json` のメタデータスキーマ、派生タイルスキーマ、データソース、変換完了状態と `bldg/tile_index.json`、参照先タイルを検証し、観測地点と検索半径を含むキャッシュだけを選択する。

準備CLIの`--list`は、同じキャッシュ検証を使う読み取り専用の全件一覧経路である。通常出力は自治体コード・データ年度・保存パスの1行形式とし、必要な絞り込みは標準的なテキスト検索ツールへ委ねる。`--jsonl`では`cache_meta.json`の詳細に保存パスを加え、1キャッシュ1行で出力する。カタログAPIには接続せず、無効なキャッシュは一覧から除外する。
取得済みの元CityGMLを残すことは既定にせず、検査や再変換が必要な場合だけ明示オプションで保持する。
準備CLIを再実行した場合は、PLATEAUカタログの整備年度、登録年度、仕様、建物ファイル数、建物ファイルサイズ合計を既存の`cache_meta.json`と照合する。全項目が一致するキャッシュは再利用し、いずれかが異なる場合は利用者の確認後に再ダウンロードして置き換える。通常の`zstarview`起動時にはこのカタログ照合を行わない。

PLATEAU建物形状はLOD0/LOD1を基本のフォールバックとし、利用可能なLOD2については`RoofSurface`だけを追加利用する段階方式とする。LOD2の壁面全体やLOD3/LOD4は、空を背景とした建物シルエット表示に対してデータ量と処理量が過大になるため対象外とする。PLATEAUのCityGML APIが提供するLOD別地物数を、対象データの詳細度確認に利用してよい（[CityGML API](https://docs.plateauview.mlit.go.jp/api/rest/operations/datacatalogcitygmlconditions/)）。

派生タイルの建物レコードは、従来の外周リング・高さ・範囲を後方互換の基本情報として保持する。LOD1では、壁面を多数の線として直接描画せず、単純立体の屋根面を観測地点から見たシルエットへ変換し、既存の屋根ライン表示と同じ描画経路へ渡す。LOD2では、CityGMLの`RoofSurface`ポリゴンと各頂点の標高を保存し、各屋根面の外周を観測地点基準のAlt/Azへ投影する。同一建物内で標高差が小さい面群は、投影後に共有辺を相殺して外周リングへ連結し、完全包含される内側の面や穴の境界は描画対象から除外する。標高差の許容値は観測地点を0mとし、距離に比例して増加させ、0kmで-1m、1kmで0m、4kmで3m、4km以上では3mに固定する。許容値が負の場合は面の結合を行わず、個別のRoofSurfaceをそのまま描画する。壁面内部や底面は描画しない。表示密度を抑えるため、LOD2屋根面の描画は当面40m以上の建物に限定し、それ以外はLOD1/LOD0へフォールバックする。

LOD1の初期描画は、箱の上面を屋根ポリゴンとして扱う。上面の各頂点を観測地点基準の東・北・上座標へ変換し、Alt/Az座標上へ投影したポリゴン境界を、既存の建物外周ラインと同じ`UrbanOutlinePolyline`経路へ渡す。壁面内部や底面は描画対象にしない。方位角の0度/360度境界、視野境界、検索半径境界の分割、エッジのサンプリング、細い輪郭の簡略化は既存のOverture Maps用処理を再利用する。これにより、LOD1の立体情報を使いながら画面上の線の密度を増やさない。

形状情報を追加する派生タイルでは、`DERIVED_TILE_SCHEMA_VERSION`を3へ更新し、`geometry_lod`、`geometry_mode`、LOD別建物数、`roof_surfaces`、詳細サーフェス数をメタデータまたはタイルレコードへ保存する。スキーマ1/2のキャッシュは新しい形状方式に対応済みとはみなさず、準備CLIの再実行時に再生成する。キャッシュ更新判定では、既存の整備年度・登録年度・仕様・ファイル数・ファイルサイズの照合を引き続き適用する。

形状方式の実験対象は、東京タワーを含む港区の自治体コード`13103`とする。`13103`のCityGMLカタログでLOD2の提供状況を確認し、`RoofSurface`が存在する建物について派生タイルを生成する。東京タワー周辺では、Overture Mapsの従来表示、PLATEAUのLOD0、LOD1、LOD2屋根面シルエットを比較し、輪郭の再現性、画面の混雑、変換時間、キャッシュ容量、起動後の描画時間を評価する。

### Building source selection

建物ソースの選択はGUI側のcontrollerが行い、各ソースの取得・変換処理とは分離する。
選択順は、観測地点に対応する有効なPLATEAU派生キャッシュ、通常のOverture派生キャッシュ、最後にOvertureのオンデマンド取得である。PLATEAUが選択された場合は、近距離建物と遠距離skyscraper補助のOverture取得を行わず、PLATEAU派生タイルだけを共通の建物読込・輪郭計算経路へ渡す。
PLATEAU対象地域かどうかだけではPLATEAU取得を開始せず、PLATEAU準備済みキャッシュの有無を優先条件とする。

PLATEAUキャッシュが有効な地点では、Overtureデータを無条件に重ねない。これにより同一建物の二重表示を避け、PLATEAUを選択した利用者のデータソース意図を保持する。
PLATEAUキャッシュがない場合は、対象地域内であっても既存のOverture経路を使い、起動時の大容量ダウンロードを避ける。

観測地点の地面標高は地点解決時に近傍DEMから求め、都市アウトラインと地形horizonは
その値を観測地点の共通基準として使う。DEMの範囲切り出しは元データの画素格子に揃え、
読み込み半径によって同一地点の補間結果が変わらないようにする。建物位置の標高は
引き続きそれぞれの代表点でDEMをサンプリングし、DEMが得られない場合は観測地点の標高を使う。

PLATEAUとOverture Mapsの通常近距離レイヤーでは、既存の派生輪郭と高さから観測地点まで5 m以上、設定した半径以内の近似壁を生成する。
壁の半径は既定0.2 km、0で無効とし、都市アウトラインの取得半径を上限とする。
観測地点を含む建物は除外する。地面標高と下端・上端高さから壁面を作り、2 mを基本間隔とし、
大きな面では分割数を抑えた小領域を観測地点基準のAlt/Azへ変換する。既存の地点基準結果とともに
キャッシュし、画面描画時に視野へ投影する。建物輪郭までの距離により5〜10 mで不透明度を増やし、
10 mから設定半径まで減らす。半径が10 m以下の場合は、5 mからのフェードインと半径への
フェードアウトを組み合わせる。不透明度ごとのパスを低い順に画面用画像へ上書きし、最大の不透明度を適用する。
分割前に観測地点から各壁面の四隅へ光線をたどり、同等以上の不透明度を持つ単一の別壁面の
内側を通って完全に隠れると判定できる場合だけ、その壁面を除外する。方位角の範囲で候補を絞り、
境界や部分遮蔽など判定が曖昧な場合は残す。除外結果は観測地点基準の結果とともに再利用する。
壁の画像は描画先の幅・高さをそれぞれ半分にして作成し、滑らかに拡大して合成する。
屋根の輪郭線は従来の解像度で描画する。
壁の視野判定と投影は四角形をまとめて計算し、観測地点基準の輪郭リスト・視線・視野・画面形状・
反転状態が同じ間は投影済みパスを再利用する。表示色や不透明度は合成時に適用する。
遠距離補助の経路では壁を生成しない。

`--use-building-top`はGUI構築前に観測高度を確定する必要があるため、通常の非同期urban controllerとは別の起動前処理を維持する。ただし建物キャッシュの選択と鮮度条件は通常経路に揃える。選択順は、準備済みPLATEAU、同じ観測地点を中心とする有効な通常2.5km Overture `building` / `building_part`キャッシュ、有効な0.15km専用キャッシュ、最後に0.15kmの不足feature typeだけの取得とする。Overtureキャッシュはメタデータの中心と半径が要求0.15km範囲を包含し、30日TTL内で、既知の最新Overtureリリースと一致することを確認する。これにより通常urban outlineを表示済みの地点では起動前の追加取得を避け、初回だけ小範囲取得へフォールバックする。

### 外部 API の識別

外部 HTTP API へのリクエストは、`zstarview/<app-version> (+service)` 形式の識別可能な `User-Agent` を持つ。  
共通基底は現行版 `zstarview/1.32.11` で、各クライアントは短い suffix を足して区別する。
将来の版更新時は、`zstarview/<current-version>` の基底部分だけを差し替えればよい。

- `build_user_agent("water-overlay")` -> `zstarview/1.32.11 (+water-overlay)`
- `build_user_agent("nominatim")` -> `zstarview/1.32.11 (+nominatim)`
- `build_user_agent("night-lights")` -> `zstarview/1.32.11 (+night-lights)`
- `build_user_agent("overture-release")` -> `zstarview/1.32.11 (+overture-release)`
- `build_user_agent("geosatellite")` -> `zstarview/1.32.11 (+geosatellite)`
- `build_user_agent("tropical-cyclone")` -> `zstarview/1.32.11 (+tropical-cyclone)`
- `build_user_agent("ip-api")` -> `zstarview/1.32.11 (+ip-api)`
- `build_user_agent("opensky")` -> `zstarview/1.32.11 (+opensky)`
- `build_user_agent("gmn-meteors")` -> `zstarview/1.32.11 (+gmn-meteors)`
- `build_user_agent("satellites-celestrak")` -> `zstarview/1.32.11 (+satellites-celestrak)`
- `build_user_agent("satellites-horizons")` -> `zstarview/1.32.11 (+satellites-horizons)`
- `build_user_agent("satellites-wheretheiss")` -> `zstarview/1.32.11 (+satellites-wheretheiss)`
- `build_user_agent("copernicus-dem")` -> `zstarview/1.32.11 (+copernicus-dem)`
- `build_user_agent("s3")` -> `zstarview/1.32.11 (+s3)`
- `build_user_agent("skyfield-loader")` -> `zstarview/1.32.11 (+skyfield-loader)`
- `build_user_agent("solar-hover")` -> `zstarview/1.32.11 (+solar-hover)`

この方針は、サービス運営側のトラフィック識別を助けつつ、障害調査でどの経路が使われたかを追いやすくする。  
仕様上の公開一覧は `docs/specification.md` に置き、実装の増減があってもそこへ反映する。

### 24時間キャッシュの履歴表示

雲・航空機・人工天体の履歴キャッシュを、設定可能な期間保持する。保持期間は起動前設定ダイアログから1〜720時間で指定し、既定値は24時間とする。過去表示の対象時刻がこの期間内なら、天体計算は選択時刻を使い、これらのデータレイヤーは通常表示と同じく現在時刻に基づく取得・キャッシュ更新を続け、最新データを描画する。Geo-satelliteは元画像PNGを撮影時刻・画像種別ごとに保存し、ライブ経路は現在利用可能な最新画像を使う。過去の天体表示と最新のデータレイヤーを同時に表示するため、データレイヤーは過去状態の再現を意味しない。ライブ表示のfresh/stale条件は保存期間から独立させる。

`--place` 起動時は、入力queryと解決後の座標を含む起動地点metadataを時刻付きcache directoryへ記録する。これはsession IDを使わず、地点の来歴とデータscopeを確認するための情報とする。各レイヤーのsource時刻・取得時刻・observer/bbox scopeを区別する。過去時刻の天体表示中もデータ取得は停止せず、現在値に基づく通常の更新を続ける。

設定値は `config.json` へ保存し、起動metadataおよび各履歴キャッシュのcleanupと、通常ビューアの過去表示中にライブレイヤーを有効にする期間が同じ値を参照する。人工衛星要素の表示鮮度やライブデータのstale fallbackには適用しない。保存先、JSON schema、欠損時の表示仕様は[キャッシュ履歴と過去時刻表示](design/cache-history-replay.md)に定める。GUIでは `--datetime` や負の `--hours` で天体表示の対象時刻を指定する。

### オフラインモード

`--offline` は `zstarview.offline` に集約する。

- 状態はプロセス環境変数 `ZSTARVIEW_OFFLINE=1` で持つ。`enable_offline_mode()` が設定し、`Popen` で起動するワーカー（雲ソース、雲 voxel、熱帯低気圧）にも継承される。`zstarview/__init__.py` が import 時に `install_from_env()` を呼ぶので、`python -m zstarview...` で起動した子プロセスも同じガードを持つ。
- 安全網として `socket.getaddrinfo` と `socket.socket.connect` / `connect_ex` を差し替え、ループバック以外への接続を `OfflineError` で失敗させる。`OfflineError` は `OSError` の派生なので、`urllib` では `URLError` として届き、既存の「取得失敗時は古いキャッシュを使う」経路（DEM、道路の夜間光、水面、Nominatim キャッシュなど）をそのまま通る。Skyfield や astropy の内部ダウンロードもここで止まる。
- astropy は `iers.conf.auto_download = False`、`auto_max_age = None`、`iers_degraded_accuracy = "warn"`、`data.conf.allow_internet = False` にする。同梱 IERS-A 予測値の鮮度切れで変換が例外になるのを避けるためで、UT1 の誤差は描画上無視できる。
- ガードだけに頼ると取得失敗の警告と再試行が続くため、主要な入口では通信の手前で判断する。
  - `apply_offline_overrides()`：短期間データのレイヤー（航空機、人工衛星、予報降水量、雲、熱帯低気圧、GMN メテオ）の opacity を 0 にし、`--moon-style image` を `sphere` にする。Open-Meteo 確認ダイアログより前に適用する。
  - 地点解決：`_resolve_auto_location()` はオフラインでは即座に失敗し、省略時の順序（auto → last → Tokyo）で前回地点へ進む。`search_place_candidates()` はレート制御スロットと通信を飛ばし、キャッシュだけを返す。
  - `resolve_search_targets()` は手元カタログで見つからなければ人工衛星と JPL の callback を呼ばない。
  - 月と太陽のホバー画像コントローラーは取得を開始しない。
  - Overture：リリース確認は失敗をメタデータへ記録せず、最後に確認したリリースを返す。`overturemaps` サブプロセスは起動前に `OfflineError` を送出する。起動時の `overturemaps` 存在確認も行わない。
  - 夜間光：sky worker はオフラインでタイルを得られなかった地点を記録し、同じ地点では再計算しない。
- export-image はオフライン時に `allow_partial_data` を真にし、キャッシュ不在のレイヤーを省略して保存する。
- キャッシュ不在は `OfflineCacheMissError`（`OfflineError` の派生）で表す。`is_offline_failure()` は例外の `__cause__` / `__context__` / `URLError.reason` をたどってオフライン起因かを判定し、`failure_log_level()` がオフライン起因の失敗を INFO に下げる。
  - 地形（`fetch_copernicus_dem`）はオフラインではダウンロードを試さず、TTL 切れのタイルも使い、不足タイル数を INFO 1 行にまとめる。タイルが 1 枚もない場合はオンライン時と同じ例外文言にして、既存の平坦地面フォールバックを通す。GUI のステータスは `not cached`。
  - 道路の夜間光は 10 km / 5 km のスコープキャッシュを TTL に関係なく使い、なければ `OfflineCacheMissError`。半径を縮めた再試行は行わない。
  - 水面は Overpass に問い合わせず、キャッシュがなければ海域マスクだけで描く。走査半径は観測者の高さで変わるため、完全一致のスコープがなければ `load_water_overlay_cache_for_location()` で同じ地点の最新キャッシュを使う（GUI と export-image で共通）。
  - `-B` の建物上端解決（`_resolve_building_top_height_m`）は、オフラインでは TTL 切れの Overture キャッシュも使い、0.15km の取得は行わない。キャッシュがなければ INFO を出して地表基準にする。
  - 都市アウトラインは派生データがあれば TTL に関係なく使い、なければ `OfflineCacheMissError`。超高層タイルは存在するものだけ使う。
  - export-image の `_await_background_task_result()` はオフライン起因の失敗を `layer_failures` に入れず、INFO で「skipped offline」と記録する。
- replay の `LocationDataPolicy` は地点解決だけを対象とする仕組みで、レイヤー全体を止めるオフラインモードとは独立に残す。

### Nominatim 地名検索とオフラインキャッシュ

Nominatim 地名検索は、CLI の `--place` から使う地点解決経路と、起動前設定ダイアログの候補選択経路で共通の検索サービスを使う。
通信、レート制御、候補の正規化、キャッシュ保存、通信不能時のフォールバックをこの共通層へ集約し、各入口はトップ候補の採用または候補選択だけを担当する。
既存の `last_city` は最後に選択した地点を復元する設定データであり、再利用可能な候補一覧を保持する地名検索キャッシュとは分離する。

#### オンライン優先フロー

利用者が検索を確定した場合、キャッシュを先に読まず Nominatim へ問い合わせる。
1 件以上の有効な候補を取得した場合だけ、候補一覧と取得日時をキャッシュへ保存する。
0 件結果、HTTP エラー、応答解析失敗、通信失敗では、既存の正常なキャッシュを更新しない。

DNS 解決失敗、接続失敗、タイムアウトを通信不能として扱い、その場合だけ一致するキャッシュを検索する。
`HTTPError` は `URLError` より先に分類し、HTTP エラー、レート制限、応答解析失敗、オンラインでの 0 件結果を通信不能フォールバックへ含めない。
一致するエントリーを使用した場合は、キャッシュ由来であることと取得日時を呼び出し側へ返し、UI または起動時状態表示で利用者に示す。
通信不能かつキャッシュ不在の場合は、空の `LocationResolveError` へ変換せず、ネットワーク接続が必要で利用可能なキャッシュがないことを示す利用者向けメッセージを保持する。

公開 Nominatim への要求は利用者操作によってのみ開始し、入力途中の自動補完には使用しない。
`filelock` を使った共有ロックと最終要求開始日時のマーカーをキャッシュルートに置き、複数スレッドおよび複数プロセスを通じて要求開始間隔を 1 秒以上にする。
共有ロックを保持した状態で直前の要求開始日時を確認し、必要な待機後、今回の開始日時を記録してから通信を開始する。

#### キャッシュキー

キャッシュキーは次のフィールドから作る。

- キースキーマバージョン
- Nominatim provider の識別子
- 正規化した検索文字列
- 正規化した国コード。未指定は `null`
- 正規化した応答言語
- Nominatim へ要求する最大候補数

検索文字列は Unicode NFKC で正規化し、前後の空白を除去し、連続する Unicode 空白を単一の ASCII space にする。
国コードは前後空白を除去して ASCII lowercase にし、言語は前後空白を除去する。
通常のアプリ入口では検索結果数を 5 件とする。既存の検索 API は候補数を変更できるため、異なる件数の結果を混同しないよう最大候補数もキャッシュキーへ含める。
異なる生入力が同一の完全な正規化キーになった場合は同じ検索とみなし、オンライン検索成功時に同じエントリーを更新する。
`zstarview-export-replay` の `-p` も通常GUIと同じonline-first検索経路を使う。
過去の起動metadataや検索cacheから独自に座標を先取りせず、成功時のキャッシュ上書きと
通信障害時のフォールバックを共有する。replayのローカルデータpolicyは標高・建物読込と
検索cache rootの指定だけに適用し、検索順序や候補のランキングを変更しない。

キーオブジェクトを `sort_keys=True`、`separators=(",", ":")`、`ensure_ascii=False` の JSON として直列化し、UTF-8 bytes の SHA-256 hex digest を計算する。
検索文字列をファイル名へ直接含めず、digest だけを使うことで、パス区切り、予約名、制御文字、Unicode、長い入力を含む検索も各対応ファイルシステムで安全に扱う。

#### JSON バケット

保存先は次の構造とする。

```text
CACHE_PATH/nominatim/schema-1/
  <sha256>.json
```

SHA-256 はエントリーの一意識別子ではなく、JSON バケットを選ぶための値として扱う。
各 JSON ファイルは次のように複数エントリーを保持できる。

```json
{
  "schema_version": 1,
  "hash_algorithm": "sha256",
  "entries": [
    {
      "key": {
        "key_schema_version": 1,
        "provider": "nominatim.openstreetmap.org",
        "query": "松江駅",
        "countrycode": "jp",
        "language": "ja",
        "limit": 5
      },
      "original_query": "松江駅",
      "fetched_at_utc": "2026-08-02T07:15:00Z",
      "results": [
        {
          "name": "松江駅",
          "lat": 35.4641778,
          "lon": 133.0628539,
          "category": "railway",
          "type": "station",
          "importance": 0.7
        }
      ]
    }
  ]
}
```

読み込み時はファイル名の digest だけを信用せず、各エントリーの完全な `key` を要求キーと比較する。
各 `key` から digest を再計算し、バケットのファイル名と一致することも検証する。
完全なキーが一致すればそのエントリーを使用し、digest が同じでも完全なキーが異なれば SHA-256 collision として別エントリーのまま保持する。

オンライン検索成功時は、完全なキーが一致するエントリーだけを置き換える。
一致するエントリーがなければ同じバケットの `entries` に追加し、他のエントリーを保持する。
バケット単位の `filelock` を取得して read-modify-write を行い、同じディレクトリに一意な一時ファイルを書いてから atomic replace する。
JSON が破損している、スキーマが未知である、または既存エントリーの検証に失敗した場合は、そのファイルを空バケットとして上書きせず、警告を記録して保存を中止する。

取得日時は UTC の timezone-aware datetime として生成し、JSON には `YYYY-MM-DDTHH:MM:SSZ` 形式で保存する。
キャッシュ利用時の表示では利用者のローカル日時へ変換し、タイムゾーンを識別できる表記にする。
この用途は少数の完全一致検索であり、部分一致、集計、関係検索を必要としないため、SQLite などの DB は導入しない。

- 夜間光は GeoTIFF 由来の base layer と、DEM と距離だけで作る edge glow layer に分ける。base layer は `night_light_opacity`、edge glow layer は `ridge_glow_opacity` で別々に調整する。
- night lightのbase layerは、屋根ポリゴン面塗りの実効係数から昼間用の下限 `0.05` を差し引き、0未満を0へ丸めた係数を指定opacityへ乗算する。このため太陽高度 `-9` 度以下かつ深夜活動係数が1のときは指定opacityの95%、`-4` 度以上では0となる。太陽が最低高度を通過した上昇側では屋根面塗りと同じ活動係数を含める。edge/ridge glowは自然の大気光として人工照明係数から除外する。
- Diffuse skyの指定opacityは最大値として扱い、night lightと同じ`-9`度から`-4`度の太陽高度係数を乗算する。night light、road night lights、または建物屋根面が有効かつ対応データを持つ場合は、下方南中前の人工照明活動量1.0でDiffuse skyを40%へ抑え、下方南中後に活動量が0.45へ下がるのと反比例して100%まで滑らかに戻す。これらの人工照明が描画されない場合はこの減光を外す。ridge glowは人工照明として数えない。
- Gaia/AKARI共通のDiffuse skyには上記の人工照明減光を適用する。これは静的な全天データへの視覚的補正であり、データ自体の時間変化ではない。
- Diffuse sky と空色ディスクは縮小描画サーフェイスでは地形地平線下まで描き、両者と雲を合成した後に、フル解像度の画素逆投影・方位補間で地形地平線下をクリップする。境界には表示1.5画素相当のsoftnessを持たせ、地形プロファイルがない場合は従来どおり全域を描画する。
- 夜間のsky discのアンビエント成分はDiffuse skyから独立させ、従来どおり`-9`度から`-18`度の区間で暗くし、`-18`度以下で最低値へ固定する。これは人工照明とは独立した、完全な黒を避けるための空色成分である。
- base layer は副稜線レイヤー配列の順序をそのまま使う。`0` 番は最初の副稜線であり、主稜線は入力しない。
- edge glow は夜間光の色を借りず、`GlowMask` の固定 tint で描画する。ridge glow はまず粗い全域モデルを作り、その上に主稜線近傍だけ高解像度の補助モデルを重ねる二段構成としてよい。細密モデルの外側は外挿で埋めず、粗いモデルを背景値として使ってよい。
- `night_light` と `ridge glow` は別の強度成分として扱うが、どちらも最終的には同じ `GlowMask` 系へ折りたたんで描画する。night light は従来どおり直接描画し、ridge glow は粗い基盤モデルと局所的な補助モデルを合成してから描画時にマスクする。粗密の切り替えは硬く切らず、境界では滑らかにつないでよい。
- 建物アウトラインの閉じた屋根ポリゴンは、輪郭線とは別の面塗り成分として扱う。従来の基準アルファを1.2倍し、その最終強度へnight lightと同じ太陽高度係数（`-4` 度で0、`-9` 度で1、その間はsmoothstep補間）を単純に乗算する。太陽が最低高度を通過した上昇側だけは、Road Lightsの連続光と同じ活動係数をさらに乗算する。両係数の積には0.05の下限を設け、昼間も夜間最大光量の5%に相当する面色を残す。輪郭線は地物可読性のため、面塗りとは独立した強度で残す。
- ただし glow の可視化は cloud のダウンロード完了と結びつけず、sky snapshot 完了時点で base sky/glow を先に再描画してよい。cloud overlay は別トリガーで後から重ねてよい。
- `atmosphere.py` の現行RGB近似では、入射太陽光を `[0.45, 0.625, 1.0]` とする表示校正を除き、散乱積算後の独立したチャンネルバランス補正を適用しない。代表波長から求めた散乱・消散係数、オゾン吸収、薄明多重散乱近似を表示色の入力とし、多波長化・CIE XYZ変換へ移行する際にはこれらの3成分近似を置き換える。
- `zstarview-export-image` は単発処理なので、GUI の splash warm-up と同じ順序に固定しなくてよい。短命な aircraft / satellites を先に開始し、cloud は独立した取得経路として早めに走らせ、取得後に (alt, az) グリッドへの取り込みを行い、その後に terrain / urban / water / night light をまとめて揃える設計を許容してよい。
- 通常GUIの雲取得間隔は `cloud_refresh.CloudRefreshSchedule` で決める。
  `CloudImageState` / `GeoSatelliteState` が独立した判定履歴を持ち、成功時の
  `window_update_cloud.py` callbackで、観測時刻と取得処理の完了時刻を渡す。
  通常600秒、30分以上古い新たな観測結果なら300秒、同じ配信元・観測時刻の再取得なら600秒とする。
  source不明・未来時刻は前倒しせず、取得失敗時は従来の600秒再試行を維持する。
  cloud stateをclearしたときは判定履歴も破棄する。単発exportやreplayには適用しない。

### 観測者中心の雲 alt/az グリッド

cloud source worker は、取得した B13/C13 の輝度温度から観測者中心の `(alt, az)` グリッド（`CloudAltAzGrid`）を作る。内部シェルは地球半径に `3 km`〜`11 km` を1 km間隔で加えた9層で、低層・中層・高層の3グループへ集約した `shell_amounts` と被覆率を持つ。

- 雲の画面描画は voxel 方式だけで行い、このグリッドは描画に使わない。GUI と画像出力は、観測地点、被覆率、メタデータを得るためにこのグリッドを参照する。`scripts/prototype_cloud_cell_shading.py` の `--cloud-model altaz` も、このグリッドの `physical_shell_amounts` を使う。
- 全体雲量からシェル重みを決める処理は、9内部シェルを3グループへ分けて使う。グループ間は `(0.0, 1.0, 0.0)` から `(0.2, 0.6, 0.2)` へsmoothstepで補間し、各グループ内は3高度へ均等に `1/3` ずつ分配する。
- 以前の shell 描画方式では、B16/C16 の companion band を取得して、B13-B16 差でグループ間の配分を補正していた。shell 方式の廃止に合わせて B16/C16 の取得と補正を削除した。B16 を含む旧形式の雲履歴キャッシュは、`CloudSourceData.__setstate__` が `auxiliary_bands` を読み捨てて、そのまま読み込める。

### 雲 voxel 描画方式の統合

voxel RGBA は画像半径を `content_fov_deg` に対応させた正方形として生成する。画面のジオメトリ半径は `edge_fov_deg` に対応するため、合成時の画像半径を `geometry.radius * content_fov_deg / edge_fov_deg` として、ジオメトリ中心に配置する。生成と合成の角度尺度を一致させ、縦長ウィンドウでも content FOV 内の空に画像境界による未描画帯が生じないようにする。GUI の画像半径は従来の 256 px（一辺 513 px）を維持する。

雲は voxel 方式だけで描画する（以前の `--cloud-mode voxel|shell` と shell 描画は削除した）。3次元 voxel 内の光の散乱と透過を計算した後、`cloud_voxel_style.apply_cloud_voxel_marks()` で透明マークをRGBA画像に適用する。マークの種類は `CloudVoxelMarkStyle(mode, value)` で表し、`--cloud-stripe holes[,RADIUS]|cutout[,COUNT]` から作る。既定の `holes` は、alt/az に固定した同心円配置の点（`CLOUD_VOXEL_HOLE_COUNT` 個。天頂の1点と、高度を等間隔に分けた輪からなり、最下段の輪は地平線から半ステップ上に置く。各輪の点数は周長に比例させ、最大剰余法で合計を合わせる。k 番目の輪（0始まり）の並びの始まりを `k * 黄金角` だけ回し（その輪の点間隔で剰余をとる）、同心円のままひまわり状の渦の筋をつくる。`cloud_hole_altaz_points()` で一度だけ計算する）を `project_altaz_to_cloud_image()` で雲画像の画素に順投影し、各点の周りの小さな窓だけアルファ値を下げる。この順投影は worker が使う `render.ground_mask.inverse_project_disc()`（正距方位図法、画像半径＝`content_fov_deg`）の逆変換であり、テストで往復一致を確認する。半径は513x513基準の画素数で、画像サイズに比例させる。`cutout` は従来の513x513基準の疎な右下がり45度の線（既定16本・1px）である。この後処理はGUIと画像出力・replayで共通化し、voxel worker protocolには含めない。穴の位置は視点（`view_center`、`content_fov_deg`）に依存するため、replay では各チャンクの画像取り込み後に viewer の視点を渡して適用する。voxel 方式には追加の雲量しきい値を適用しない。`--cloud-voxel-threshold` は廃止し、指定時は CLI 入力エラーとする。

実装は次の構成とする。

親側の `cloud_voxel_worker` transport module は、起動・入力のpickle化・結果画像の
検証と読み込みを担当する。voxel renderer、JIT関数の登録、地表光の取得・サンプリング、
逆投影のimportは `_worker_main()` 内に置き、親がtransportをimportする際には実行しない。
この境界は雲描画専用処理に対するもので、night lightなどの既存Numba経路は親に残る。
GUIからは共有native lockを渡し、親の入力pickle化・NumPy画像読み込み・昼夜係数・
マーク後処理だけをロック内で行う。subprocess起動・待機中にはそのロックを保持しない。

1. **入力と互換性**: `--cloud-stripe holes[,RADIUS]|cutout[,COUNT]` を parser で `(mode, value, 1.0)` に解釈し、`CloudVoxelMarkStyle.from_cli()` で `SkyWindowRuntimeOptions.cloud_voxel_mark_style`、GUI初期化、`zstarview-export-image`、replay まで渡す。保存済みプロファイルに残る shell 用の値（`cloud_mode`、`cloud_missing_tint_opacity`、`holes`・`cutout` 以外の `cloud_stripe`）は `launch_profile` の読み込み時に捨てる。
2. **データ保持**: cloud source workerが返す `CloudSourceData` を描画入力まで保持する。データ取得と衛星画素から密度ボリュームを作る処理は、GUIプロセスやpaint経路で行わない。
3. **別プロセス計算**: 一回実行worker subprocessへpickle化したsourceと、地点・視点・描画サイズ・太陽高度/方位・月の高度/方位/離角・zstarview本体が決める昼光混合率と環境光基準RGBを渡す。protocol version、request id、source key、出力寸法を結果manifestで検証する。workerが視点別のvoxel cloud RGBAを計算し、親は完了後に結果を取り込む。GUI側は既存の非同期controllerからworkerを起動し、要求世代/source identityが最新の結果だけ採用する。新しいカメラ要求は実行中workerの後に最新の1件を処理し、古い完了結果は表示しない。exportは同じworkerを同期的に呼び、layer deadlineを適用する。
4. **レイヤー順**: 空色と雲を合成した後に星をすべて描く。voxel cloud RGBAはworkerで計算し、cloud opacityはpaint時の合成へ一度適用する。雲は星を隠さない。
5. **失敗と資源管理**: timeout、worker異常終了、source不在、投影範囲外は既存cloud layerの失敗/partial-data扱いへ流し、GUIでは直前のsourceに対する新しい要求だけを継続する。workerの一時入力・出力は専用一時ディレクトリに置き、正常終了・失敗・取消後に回収する。native Astropy/Skyfield workはこのworkerで増やさず、雲numpy/numba処理と地表光の取得・サンプリングを分離する。

voxel workerは native CloudDisc sourceでは `shade_native_voxels` を使い、Geo-satellite sourceでは同じ `_render` 光線追跡・照明・透過処理を共有する `shade_geo_satellite_voxels` を使う。native sourceは高度ごとの衛星画素 footprint を近似するが、Geo-satellite sourceは投影済み雲量画像の水平グリッドを一度だけ観測者近傍へ線形化し、1〜9 kmの全層で同一のセル基底と原点を使う。Geo画像の雲量を既存の低・中・高層重みから各3層へ均等配分する。画像は物理的な輝度温度ではなく、層配分も実測雲高度や光学的厚さではない。worker境界のprotocolにはversion、request id、source key、結果寸法を持たせ、壊れた出力や古い要求を通常結果として受け入れない。

night light の有効条件は terrain horizon の生成結果の有無に合わせる。terrain horizon がまだない間は夜間光の alpha grid を作らず、terrain horizon が用意できた時点で 1 回だけ alpha grid を生成して保持する。以後は同じ terrain 条件ではその grid を使い回し、terrain horizon が再計算されたときだけ night light 側も再生成する。

Geo-satellite voxelでは、観測点から水平約200 kmの描画対象へ切り出した全画素を `_refine_geo_cloud_amount` で縦横3分割する。各小セル中心の雲量は有効な原画素の双線形補間で得る。欠損画素に属する小セルは空とし、範囲外の参照座標は端に固定する。レイ追跡は全域で小セル幅を使い、視線と太陽光の光学的厚さには雲量と実際のkm通過距離を掛ける。雲量を分割数で割らない。内部比較用に `refinement_factor=1` を指定できる。結果metadataは倍率、対象窓、元と細分化後のセル基底を記録する。

### Open-Meteo モデル予報降水雨線

降水雨線は Open-Meteo Weather Forecast API の `minutely_15` に含まれる
`precipitation`、`rain`、`showers` を入力とする。現在時刻の前後にある15分区間のうち、
区間中心が現在時刻に最も近いものを選ぶ。これは複数の気象モデルを地域に応じて選択・接続した予報値で
あり、実況観測データとして命名または表示しない。API 応答の対象時刻、15分スロット、
取得時刻、単位、欠測を降水値と分けて保持する。ネットワーク取得、JSON検証、DEM
sampling、Alt/Az投影、および描画 primitive の準備は worker thread で行う。UI/render
thread は準備済みの48地点の雨線群を QPainter で描画するだけとする。

provider は `open-meteo-noncommercial` に固定し、`off` は opacity `0` で表す。
Commercial API、API key、provider選択はこの仕様の対象外とする。Free API endpoint は
`https://api.open-meteo.com` とし、利用者が非商用条件、利用上限、帰属条件を明示的に
確認した後だけ通信する。

同意状態は通常設定と同じ共通 `config.json` に
`open_meteo_noncommercial_terms_version` として保存する。現在の要求versionとの完全一致を
有効な同意とし、versionを上げた場合は再確認を要求する。この値は秘密情報ではない。
`zstarview-gui` の設定ダイアログは opacity を保存するだけで、同意状態を変更せず通信も
行わない。

対話的起動では、effective opacity が正かつ有効な同意がない場合、viewer本体、worker、
network requestを生成する前に確認ダイアログを表示する。同意時はversionを原子的に保存して
起動を続行する。キャンセル時は保存済みopacityを維持したまま起動を中止する。
`zstarview-gui` が起動元の場合は、既存の起動失敗通知経路で
`Open-Meteo Free API terms were not accepted. zstarview was not started.` と同等の
メッセージを最後に表示し、launcherも終了する。降水だけを無効化してviewerを続行する
fallbackは設けない。

HelpメニューのOpen-Meteo利用条件actionは常時有効とする。降水表示toggleはaction自体を
常時表示し、保存済みopacityが `0` ならdisabled、正ならenabledとする。起動完了時には
同意済みであることが保証されるため、toggle処理はnetwork同意を扱わず、runtimeの
enabled状態だけを変更する。Offにしても保存済みopacityは変更しない。

Helpメニューのライセンス一覧は、実行中versionと最新READMEのずれを避けるため、
`zstarview.data/licenses_and_data_sources.md`をpackage dataとして同梱し、`QTextBrowser`で
表示する。文書中の外部URLだけをdesktop browserへ渡す。本文はread-onlyだが選択・コピーを
許可し、`Copy All`は`toPlainText()`をclipboardへ書き込む。

`zstarview-export-image` は同じ precipitation opacity option と共通設定を読むが、同意を
得る対話UIや同意を書き込むCLI optionは持たない。opacityが正なら、argument検証後かつ
network、renderer、出力一時ファイルの準備前に同意versionを検査する。無効ならASCIIの
案内をstderrへ出して非0終了し、出力を残さない。有効なら周辺48地点を取得し、通常viewerと
同じ投影・雨線style・`Forecast: Open-Meteo` 帰属を単発画像へ描画する。opacityが `0` なら
同意状態を読まず、降水処理を開始しない。

表示サンプルは観測地点中心の環状領域へ黄金角フィロタキシスで配置する。既定値は
`min_distance_km = 0`、
`max_distance_km = 32`、`sample_count = 48`、黄金角
`137.507764°` とする。サンプル番号 `n = 0..N-1` の方位角と距離は次式で求める。

```text
azimuth_n = (azimuth_offset + n * 137.507764°) mod 360°
t_n = (n + 0.5) / N
distance_n = sqrt(min_distance_km^2
                  + (max_distance_km^2 - min_distance_km^2) * t_n)
```

`sqrt` による半径を使い、環状領域内の点密度を面積に対してほぼ一様にする。
`azimuth_offset` は北基準の固定値 `0°` とし、時刻やデータ更新によって配置を回転
させない。観測地点の移動時は同じローカル配置を新しい地点へ連続的に移し、格子線
への近さを条件に全点を半セル移動するような不連続な切り替えは行わない。周辺配置には
中心点を生成せず、`0 km` から地理的に投影する雨線を置く。最内側の点も中心点とは一致
しないため、黄金角によるひまわり状の配置を維持する。

周辺48地点を1回の複数座標リクエストへまとめ、
`cell_selection=nearest` を指定
して、既定の land cell 選択による表示点の意図しない移動を避ける。応答は要求順との
対応、座標数、時刻、単位を検証する。各値は Open-Meteo が選択した気象モデル格子から
補間・downscaleした地点予報であり、表示点を元モデルの格子点とは扱わない。欠測または
降水強度が `0.1 mm/h` 未満の表示点は描画しない。

Open-Meteo の `minutely_15.precipitation` は各時刻までの直前15分間の積算量 `amount_mm`
として扱い、次式（`interval_seconds = 900`）で1時間当たりの強度へ正規化する。選択区間は
時刻から15分を引いた範囲として保持する。15分スロット、未知の単位、
または必須 interval の欠落は応答全体の検証失敗とする。個別地点の欠測はその地点だけ
を描画対象外とし、`0 mm` と区別する。

```text
r_mm_h = amount_mm * 3600 / interval_seconds
```

雨線の基準位置は表示点の DEM 地表高とし、既存の地理オーバーレイと同じ地球曲率、屈折、
観測者高度、Alt/Az 変換を使う。DEM が利用できない場合は基準楕円体面または既存の
地表フォールバックを使い、降雨取得まで失敗扱いにしない。雨線の上端は物理高度として
投影せず、基準位置から距離に応じた表示角 `height_deg` だけ仰角を増やす。

```text
t = clamp(distance_km / 32, 0, 1)
height_deg = 7.0 + (1.6 - 7.0) * t
opacity_factor = 1.0 + (0.35 - 1.0) * t
streak_count = clamp(ceil(log2(1 + r_mm_h / 2)), 1, 6)
```

各表示点では、DEM地表高（DEM利用不可時は標高0m）の投影点を下端中央とする画面空間の
正方形タイルを上方へ作り、タイルを基準に
「／」方向（右上から左下へ向かう）の平行な実線を描く。タイルの一辺は既存の距離ベースの
表示スケールから得る表示高さを使い、距離が遠いほど小さくする。上下方向はこのタイルの
高さをクリップ範囲として維持し、左右方向は正方形の境界でクリップしない。これにより各線分の
上下端の高さをそろえつつ、左右へ自然にはみ出させる。これは物理的な雨柱の投影ではなく、方向・距離・
降水強度を読み取るための模式記号である。

雨線群のローカルな「上」は画面全体の上に固定せず、各表示点を同じ方位のまま少し高い高度へ
移動した点を投影して求める。その方向が画面上のローカル天頂方向になるため、天頂を見上げる
視野では雨線が画面中心方向を向き、通常の視野では従来の斜め方向に近い外観を保つ。

降水強度は線の本数で表す。線の中心オフセットは、線間隔を `s` としたとき、奇数 `n` では
`(-floor(n/2), ..., 0, ..., floor(n/2)) * s`、偶数 `n` では
`(-(n-1)/2, ..., -0.5, 0.5, ..., (n-1)/2) * s` とする。`s` はタイル一辺の `0.25` 倍を
基準とし、viewer の line-width scale を乗算する。これにより奇数本では中央線が
現れ、偶数本では中央が空く。例えば1本は中央、2本は中央を挟む左右、3本は中央とその
左右に配置する。

雨線は単一の明るい青色の実線 pen で描き、現在の雨滴表現にある点線の重ね描きは行わない。線端は
角型とする。表示線幅はタイル一辺の `0.08` 倍を基準とし、600px基準表示で
`0.45` pxを下限としてサブピクセル化による消失を防ぐ。タイル寸法と線幅の両方を
距離に応じて縮小し、近距離ほど大きく太く、遠距離ほど小さく細くする。線 alpha は従来の2倍にする
（最大255へクランプ）。降水強度は線の本数、距離はタイルの大きさと
不透明度で表す。最終 alpha は
レイヤー全体の opacity に `opacity_factor` を乗算する。線幅と線間隔には viewer の
line-width scale を乗算し、距離補正後の alpha と既存の色体系を使う。表示角は物理的な
降水高度ではない。

柱が同じ画面方位へ密集する場合は、画面空間または小さな方位ビン内で弱い柱から間引いて
よいが、強い降水地点を弱い地点で置き換えない。欠測と `0 mm/h` は全処理段階で区別する。

`minutely_15` は現在のWeather Forecast APIの同じ `/v1/forecast` endpointから取得し、
別のAPIサービスやCommercial APIへ切り替えない。日本などネイティブな15分モデルの対象外
地域では、Open-Meteoが時間予報から補間した値になることをUIおよびライセンス説明で明示する。
成功応答は process 内 memory cache だけに保持し、disk cache は作らない。同一 scope
の freshness TTL は10分とし、timezone-aware UTC の `fetched_at_utc` で判定する。
`time.monotonic()` は request timeout や worker shutdown deadline にだけ使用する。
cache key は schema version、固定 provider mode、順序付きの丸め済み48座標、要求変数、単位、
`cell_selection` から作り、API key を含めない。同一 key の同時 request は1本へまとめる。
観測地点または provider mode が変われば別 scope として即時取得する。更新失敗時は以前の
雨線へ fallback せず、現在の雨線を消して `Precipitation: unavailable` とする。

更新時はフィロタキシスの表示位置を固定し、雨線の本数だけを新しい対象時刻へ更新する。
これにより、気象モデルの格子配列、同心円、または更新ごとの配置回転が視覚上の偽の
移動として読まれることを避ける。単一時刻の柱配置から降雨の接近・遠ざかりを推定
しない。レイヤーが有効な間は対象時刻と `Forecast: Open-Meteo` に相当するコンパクトな
statusを表示する。intervalは応答検証と内部データに保持するがstatusには表示しない。
レイヤーが無効な間は画面上に帰属を表示しない。完全な帰属、CC BY 4.0 link、
およびモデル予報値を表示用の柱へ変換した旨をREADMEとライセンス表示に記載する。

### Water polygon and boundary overlay

水面の取得経路は、Overpass API の `natural=water`、`waterway=riverbank`、および対応する multipolygon を `WaterPolygonFootprint` として保持する。各フットプリントは外周リング・内周リング・水面種別・タグを持つ。`waterway=river` などの中心線だけの要素は、閉じた水域ポリゴンへ復元できない場合、水面フットプリントとして扱わない。

取得後は観測地点中心のローカル座標へ変換し、距離に応じたグリッドでリングを簡略化する。簡略化済みのフットプリントをキャッシュおよび後続処理の共通入力とし、取得元の高密度な頂点列を描画経路へ渡さない。これにより、水面線を描画しても元データの頂点数に比例した描画負荷を避ける。

水面表示は同じ簡略化済みフットプリントから境界折れ線を作る。また、レイ走査の水面サンプルから円周線を作る。サンプル点そのものは画面に描画しない。

- 水面サンプル: 観測地点からのレイを走査し、ポリゴン内部に入ったサンプルを `WaterOverlayPoint` に変換する。これは円周線の頂点データであり、点マーカーとしては描画しない。
- 境界折れ線: 外周リング（必要に応じて内周リング）を観測地点基準の Alt/Az へ投影し、表示可能な隣接頂点を `WaterOverlayPolyline` 相当のフラグメントへまとめる。ポリゴンを面として塗りつぶすのではなく、低密度の輪郭線として描画する。

境界折れ線の生成では、まずローカル XY 上で境界線用の最大10kmの円にリングをクリップする。水面サンプルは従来どおり、観測高度に応じた水面走査範囲を使う。次に、距離円との交点を補間したリングを Alt/Az へ変換し、前方半球と表示 FOV の判定を行う。方位角の 0 度/360 度境界、視野境界、地平線下の表示限界をまたぐ部分は分割し、画面外の線分を描画しない。

境界線の太さと透明度は距離に応じて調整し、遠距離の川・湖が輪郭として読めることと、都市アウトラインや地形線が過密にならないことを両立する。

海面マスクの距離帯は、それぞれの下限を除外し、上限までのサンプルだけを生成する。近距離の海面を異なる解像度の帯域から重複取得しない。

水面の円周線は、走査サンプルから生成する。`water_surface_rings.build_water_surface_ring_polylines()` は `scan_distance_m`、水面種別、水域IDでサンプルをグループ化し、各グループを方位順に並べる。距離は投影後の斜距離 `distance_km` ではなく、走査時の地表距離を使う。半径の選別や新しい水域判定・標高推定は行わず、既存サンプルの仰角・方位角・高さの投影結果をそのまま線の頂点として利用する。

隣接する走査方位の点だけを結び、欠落した方位サンプル、水面種別または水域IDの変化で線を分割する。方位間隔が走査間隔の1.5倍を超える点を結ばない。0度/360度の境界では隣接点を短い別フラグメントとして結び、全周水面の場合にも途切れない線を作る。孤立点と `scan_distance_m` を持たない点は円周線にしない。点の重複は同一方位ごとに除去する。

GUIでは `WaterOverlayController` の背景処理でサンプルから円周線を組み立て、既存の水域境界線と併せてスコープキャッシュへ格納する。円周線の組み立ては純Pythonのグループ化・ソートであり、新しいネイティブ演算は行わない。内陸水域の高さも元のサンプルの投影を引き継ぐ。画像出力でもサンプル取得後に同じ生成関数を使う。水面線の更新はフレームキャッシュキーに反映し、既存の水面線レンダラーで視野クリップと地形遮蔽を適用する。線幅は水面サンプルの走査時の地表距離を優先し、境界線など走査距離がない場合は投影距離を使う。500mを基準に距離の `-1/4` 乗で変化するカーブを従来の80%に縮小し、境界線と円周線に同じ距離依存の線幅を適用する。円周線は観測地点に近いほど、画面Y座標へ最大振幅1.25pxの上下変位を加える。変位は方位の24度周期と8度周期の正弦波を重ねる。最短線の振幅には観測地点からの距離による `exp(-distance_km / 0.5km)` を適用し、最短線からさらに遠い点には追加で `exp(-Δd / 0.2km)` を乗じる。位相は水域ID・種別・走査距離から決定的に作るため、再描画でも形が固定される。fast-mode と水域境界線には適用しない。アルファ値はレイヤーの指定不透明度を上限とし、最短水面線距離を `d_min` として `opacity × (1 - 0.8 × 2^(-(distance_km - d_min) / 2))` により変化させる。最短点で上限の20%となり、その距離から2km進むごとに残りの差が半分になる。地形遮蔽によるアルファ低下はこの距離カーブに乗算する。円周線は全頂点に `scan_distance_m` があることから境界線と区別し、円周線と境界線は同じ距離なら同じアルファ値で描く。

ビューポート操作中の fast-mode は円周線だけを描画し、ポリゴン境界線を省略する。これにより、サンプル由来の水面形状を維持しつつ、操作中の線密度と描画量を抑える。

#### 海岸線ベクタタイル配布・キャッシュ

海岸線の近距離境界を描画するため、OSM Water Polygonsの海岸ポリゴンから生成した静的ベクタタイルを、既存の `water_overlay` キャッシュとは分離して保持する。アプリ本体には同梱せず、GitHub Releasesのデータリリースから利用者が明示的に取得する。

既存の32×16全球グリッドをトップレベルのアドレスとして使い、通常タイルと、容量が大きいトップレベルタイルだけを4×4に分割した子タイルを混在させる。親タイルを子タイルへ分割した場合、親の重いデータファイルは保存せず、manifestだけで子タイルの存在を示す。海岸線境界の読み込み・クリップ範囲は観測地点から最大10kmとする。川・湖などの内陸水域の境界も最大10kmとするが、水面サンプルの取得範囲は既存の高度依存設定を維持する。

キャッシュの最終構造は次のとおりとする。

```text
CACHE_PATH/coastline/
  osm-water-polygons/20260725/schema-1/
    manifest.json
    READY
    grid-32x16/
      y02/x17/tile.geojson.gz
      y12/x09/children.json
      y12/x09/children/q00.geojson.gz
      y12/x09/children/q11.geojson.gz
      y12/x09/children/q22.0
```

配布アセットはx列単位とする。1つのアセットには、その経度列の全緯度行を格納する。

```text
coastline-manifest-20260725.json
coastline-grid-x00-20260725.zip
coastline-grid-x01-20260725.zip
...
coastline-grid-x31-20260725.zip
```

緯度・経度のタイル境界と配布アセットの境界は独立している。アーカイブは配布上の単位であり、展開後のキャッシュは必ず `grid-32x16/yNN/xNN/` の構造へ戻す。x列単位にすることで、経度範囲から必要なアセットを直接選択できる。

`manifest.json`には、データソースと更新日、座標系、32×16グリッド、子タイル分割係数、各トップレベルタイルの `single`/`4x4` 配置、アセット名、サイズ、SHA-256を保存する。経度列の選択に使えるよう、アセットごとのx列番号も記録する。子タイルの境界で線が途切れないよう、生成時に境界上の線分を隣接タイルへ少し重複させるか、クリップ端点を連結可能な形で保存する。

ダウンロードは専用CLI `zstarview-download-coastline` だけが行う。通常のGUI起動時に自動ダウンロードしてはならない。CLIは `--lon-min` と `--lon-max` で指定された経度範囲に重なるx列を選択し、グリッド境界まで範囲を拡張して必要なアセットを取得する。全世界を取得する場合は `--all` を使う。`--all` は全海岸線列に加えて25m全球水面マスクも取得する。既定のデータ版はGitHub Releaseタグ `coastline-data-20260725` であり、`--cache-dir` でキャッシュの基底ディレクトリを変更できる。海岸線と水面マスクのダウンロード入口は分散させない。

ダウンロードした各アセットは一時ファイルへストリーム保存し、サイズとSHA-256を検証する。展開時は絶対パス、`..`、シンボリックリンクなどの危険なエントリを拒否する。選択した全x列の展開が完了した後にだけ `READY` を作成する。更新前には既存の `READY` を外し、各x列は一時展開後に旧ディレクトリを一時退避して置き換える。置換に失敗した場合は旧ディレクトリを復元する。これにより、既存の非空キャッシュにも再実行でき、未取得・不完全・ハッシュ不一致の更新を有効なキャッシュとして扱わない。リーダーも `READY` のないキャッシュを読み込まない。

GitHub Releaseのデータ版は、アプリのリリースとは独立したタグで管理する。アプリはデータ版を明示的に指定して取得し、未取得または取得失敗時は現在の水面サンプル経路へフォールバックする。古いデータ版は新しい版の検証完了後も直ちには削除せず、利用中でない版だけを後のキャッシュ整理対象とする。

海岸線ベクタタイルは海岸境界専用であり、川・湖・池などの内陸水域の取得経路を置き換えない。GeoJSON試作では1タイルが未圧縮252MB、4×4分割後も最大約73MBとなったため、座標量子化・形状簡略化・gzipまたはバイナリベクタ形式は、配布データ生成時に別途評価する。

### Road night lights

道路夜間光（`Road Night Lights`）は、都市の交通量や街路灯の実測データを追加する機能ではなく、OSMの道路中心線から生成する補助的な視覚レイヤーとする。道路の`way`に付いた`highway`タグを分類に使い、初期対象は`motorway`、`trunk`、`primary`、`secondary`、`tertiary`に限定する。`motorway_link`などの`*_link`は初期対象から除外する。

Overpassへの要求は道路種別ごとに分けず、観測地点から10kmの範囲を正規表現フィルタで一括取得する。応答のwayとノード列は、取得後に`highway`値ごとのグループへ分ける。要求半径、対象タグ集合、データスキーマのバージョン、観測地点を正規化したキャッシュキーを使い、同じ範囲の再表示で外部APIへ再要求しない。道路キャッシュのTTLは取得日時から30日とし、取得日時はキャッシュメタデータの`fetched_at_utc`で判定する。

道路取得の更新フローは次のとおりとする。

1. 30日以内の10kmキャッシュがあれば、それを使用する。
2. 10kmキャッシュがない、または期限切れの場合は、10km範囲を取得する。
3. 10km取得がHTTP 504または通信タイムアウトで失敗した場合だけ、同じ地点の5km範囲を1回取得する。5kmのキャッシュキーは半径を含むため、10kmキャッシュとは分離する。
4. 10kmまたは5kmの取得が成功した場合は、成功した範囲の結果を保存し、描画へ渡す。
5. 更新に失敗した場合は、期限切れであっても既存の10kmまたは5kmキャッシュを stale fallback として使用してよい。利用可能なstaleキャッシュもない場合は、道路レイヤーを省略して処理を継続する。

空の正常レスポンス（`elements: []`）は取得成功として扱い、HTTPエラーや通信タイムアウトのフォールバック条件とは区別する。時刻変更では道路を再取得せず、キャッシュ済みジオメトリを再投影し、描画強度とランプ表示状態だけを更新する。

取得したwayは観測地点中心のローカルXYへ変換し、0.5kmの内側円と10kmの外側円による環状範囲へクリップする。wayの線分が距離円を横切る場合は交点を補間し、表示対象部分だけを残す。10kmの要求範囲は水域境界線の表示範囲と揃え、0.5km未満を除外することで、観測地点近傍の道路光が画面を過度に占有することを防ぐ。

クリップ前に、海岸線ポリラインと同じ距離依存グリッドで道路wayの頂点を簡略化する。遠距離ほど大きなグリッドを使い、連続する同一グリッド点を除去する。連続ストロークはこの簡略化済みwayを使い、街路灯風の光点は線長に沿った別のサンプリング結果として扱う。

道路座標をAlt/Azへ投影する前に、連続ストロークの頂点とランプ光点の各座標でDEM地表高をサンプリングする。投影に渡す目標高は、車のライトを抽象化した連続ストロークで`DEM地表高 + 0.5m`、道路ランプの光点で`DEM地表高 + 8m`とする。連続ストロークの0.5mは車灯の位置の抽象化と地表への埋没防止を兼ね、ランプ光点の8mは主要道路の一般的な灯具高を抽象化する。両者の高さを分けて保持し、楕円体高の固定値`0m`を道路座標の目標高として使用しない。

描画は、道路中心線に沿った約6000Kの白色ストロークと、中心線を約240m間隔でサンプリングした約3000Kの小さな光点の二層とする。ストロークは`road_light_opacity`を基準光量として、既存の`night lights`と同じ太陽高度フェードを適用する。さらに現在の太陽高度・方位と観測緯度から太陽赤緯と下方南中時の最低高度を求め、太陽が最低高度を通過した上昇側だけ、最低高度から`-9`度まで活動係数を1.00から0.45へsmoothstepで減衰させる。光点は太陽高度`0`度で消灯し、`0`度から`-4`度にかけてフェードし、`-4`度以下でストローク基準光量の約0.8倍として描画する。上昇側の活動係数は光点へ適用しない。光点間隔を従来の約120mから倍増させることで、オレンジ色のランプ光の密度を半分に抑える。いずれも実際の街路灯位置や交通量を意味するデータとして扱わない。

ストロークと光点には、観測地点からの距離に応じた共通の滑らかな減衰を加える。0.5km付近では係数1.0、10km付近では係数0.2を下限とし、遠距離の道路光が地平線付近を過度に占有しないようにする。

道路候補は `way` 単位で、既定値 `5000` の上限を設ける。環状範囲などの幾何学的な絞り込み後、各 `way` の観測地点からの最短距離だけで近い順に並べ、上位候補だけを道路形状の簡略化、DEMサンプリング、Alt/Az投影へ渡す。角長、道路の頂点数、`highway` 種別、取得順は選別基準にしない。1本の採用 `way` から生成される連続ストロークとランプ光点は同じ候補に属する。固定容量の候補ヒープと安定した識別子による決定的なタイブレークを使う。`--road-light-max-candidates 0` は道路夜間光を無効化する。

候補上限は取得要求やキャッシュの保存件数には適用しない。キャッシュには取得可能な道路ジオメトリを保持し、観測地点や上限値の変更時に同じキャッシュから再選別できるようにする。視野中心は候補選別に使わないため、視野変更だけでは道路候補を再選別しない。

道路データの取得・解析はGUIスレッドで行わず、既存の静的補助レイヤー準備と同じバックグラウンド経路へ置く。道路のジオメトリはキャッシュ済みの取得結果から再投影できるよう、外部座標と正規化済みタグを保持し、時刻変更では再取得せず描画強度と光点の状態だけを更新する。

Road Lightsの表示基準値は`road_light_opacity`として、既存のVIIRS Night Lights用`night_light_opacity`とは分離する。起動時のCLI値を`SkyWindowUserOptions`から`SkyWindow`へ渡し、`0`の場合は`_road_light_toggle_supported`を無効にしてメニューからの再有効化を禁止する。メニューで一時的にOFFにした場合は、起動時に指定したopacityを復元する。OFF中は`start_background_road_night_lights_update()`が即時に戻るため、道路データ取得を開始しない。

Road Lightsの通常のステータス行は、都市アウトラインと同じくデータソースや利用状態だけを示し、候補数や上限適用後の採用数を表示しない。候補数などの診断情報は必要な場合にログへ出力する。

#### 25m水面マスクの任意配布・圧縮キャッシュ

125m、250m、500mの水面マスクはアプリに同梱する既存経路として維持する。海岸付近の表示位置を改善するための25m水面マスクは、アプリ本体には同梱せず、GitHub Releaseから利用者が明示的に取得する任意データとする。キャッシュが有効な場合は0〜250mの水面走査に25mを使い、250m〜2kmは125mへ戻す。未取得、取得失敗、破損時には0〜2kmを従来の125m経路で処理する。

25mデータは、OSM Water Polygonsのデータ日付とスキーマを含むバージョン付きキャッシュに保存する。

```text
CACHE_PATH/water/
  osm-water-polygons/20260725/schema-1/resolution-25m/
    manifest.json
    READY
    water-tiles-25m-global-20260725.zip
```

現行の配布アセットは既存の32×16全球グリッド全体を一つのZIPに格納する。ZIPには、`tile_yNN_xNN.tif`、水面が全くないタイルの`.0`、全面水面のタイルの`.1`を格納する。manifestにはタイル範囲、サイズ、SHA-256、データ日付、CRS、解像度を記録する。将来、必要な経度列だけを取得できるよう、同じmanifest形式を保ったx列単位の分割へ移行してよい。

25m ZIPは展開済みTIFFのディレクトリを通常キャッシュとして作らない。アプリはZIP内の必要なTIFFエントリだけをメモリ上で読み込んで水面判定を行う。全タイルを画素配列へ一括展開せず、必要なタイルのGeoTIFFをメモリ上で開き、走査終了時に解放する。これにより、約137MBの配布ZIPを保持しつつ、展開後の約1GB規模のディスク占有を避ける。

25m水面マスクのダウンロードは `zstarview-download-coastline --water-25m` で有効化する。現行アセットは全球ZIPなので、この指定では全球ZIPを取得する。`--all` は全海岸線列に加えて25m全球ZIPも取得し、`--all --water-25m` も同じ結果となる。ダウンロード中は既存の海岸線アセットと同じ進捗表示、サイズ検証、SHA-256検証、一時ファイル、`READY` マーカーの手順を使う。ZIPは展開せず、検証済みの完成ファイルをバージョン付きキャッシュへ配置する。

25m帯では最短サンプル距離を50mとし、125m以上の帯域では従来の最短125mを維持する。25m経路は海岸付近の水面サンプル位置の改善を目的とし、川・湖の境界ベクタデータや内陸水域のOverpass取得経路を置き換えない。25mラスタとベクタ海岸線の位置が一致しない場合もあるため、25mを導入しただけで海岸線の幾何精度が保証されるとは扱わない。

## Earth guide の描画

`render/earth_guide.py` は、地平線下の画素ごとに地球を透視した先の地表点を求め、
陸・海の色と日照による明るさで塗った RGBA 画像を作る。続けて、夜側の陸の
輪郭線を白い細線で重ねる。合成処理
（`gui/composite.py` の `_overlay_earth_guide`）がこれを空のディスクへ重ねる。

処理の流れは次のとおり。

1. `ground_mask.inverse_project_disc` で表示の1/2の解像度の画素を alt/az へ
   逆投影する。`fast_mode` では塗り自体を省略する。
2. 方位ごとの可視上限（地形地平線と、観測高度による伏角＋1度の低い方）より
   下の画素だけを残し、上限から1度かけてフェードさせる。
3. alt/az を観測者の east/north/up 基底で地球固定座標の方向ベクトルにし、
   観測者位置から単位球と交わる遠い側の交点（`far_surface_points`）を求める。
   観測高度が高く、地平線直下の視線が地球に当たらない場合は無効とする。
4. 太陽の alt/az も同じ基底で地球固定座標へ変換する。太陽光は平行光とみなせる
   ので、地表点 `p` での太陽高度の正弦は `p·s` になる。観測地点では `p = up`
   となり、自分の空の太陽高度と一致する。`earth_guide_sunlight_level` がこれを
   昼・薄明・夜の明るさへ写す。
5. 地表点の経緯度で陸地率マスクを双線形補間し、陸色と海色を混ぜる。
6. 観測地点からの地表距離（`acos(p·up)`）が500km以内の点は消す
   （`earth_guide_near_fade`）。塗りは500kmで段差として切り、輪郭線は800kmまでに
   滑らかに現す。塗りをフェードさせると、明るさで表している日照と区別できないため。近い陸地は地平線直下に圧縮されて
   筋状になるため、角度ではなく距離で切る。

夜間光は同梱の0.1度格子PNGを一度読み込み、上記の地表点の経緯度で双線形補間する。値は2025年VIIRS年次輝度の0.05度ごとの平均から、各2×2画素の平均75%と最大25%を混ぜ、`log1p(value) / log1p(20)` を0〜255へ量子化した表示用強度である。地表点の太陽高度-8度から-14度への夜の重み、塗りと同じ距離・地平線マスク、地球ガイドopacityを使い、近景のNight Lightsと共通の色で塗りに合成する。近景のNight Lightsとは別のデータと表示条件を使う。PNGの読み込みに失敗した場合は夜間光だけを省略する。高速表示では塗り自体を省略するため夜間光も描かない。

都市光の不透明度には1.75倍の表示係数を掛け、1.0を上限とする。

夜側の輪郭線は、陸地リングの辺を0.5度以下に細分した線分（`load_coastline_segments`、
LOD ごとにキャッシュ）から作る。日付変更線と南極点に沿った辺はデータの継ぎ目なので
除く。各端点で太陽高度から夜の強さ（-8度で0、-14度で1）を、上と同じ地平線・
距離・視野の条件から可視度を求め、両端の小さい方を線分の強さとする。画面上で極端に
長い線分は捨てる。線分は強さを8段階に量子化し、段階ごとに `drawLines` でまとめて
描く。塗りの薄明は-18度まで続き、輪郭線は-8度から現れ始めるので、-8度から
-18度で塗りと輪郭線が重なり、昼から夜へ表示が途切れずにつながる。輪郭線の開始高度は
`EARTH_GUIDE_NIGHT_OUTLINE_START_DEG`、塗りの薄明の終わりは
`EARTH_GUIDE_TWILIGHT_LIMIT_DEG` で別々に持つ。

`fast_mode` では塗りを作らず、`compute_night_outline_segments(night_only=False)` で
LOD 0 の輪郭線を昼夜の重みなしに描く。塗りの逆投影と陸地マスク参照がなくなるため、
2400×1400 の表示でも1回約20〜30msで済む。

陸地率マスクは LOD 3 の陸地リングを Pillow で経緯度平面へ塗りつぶして作る
（0.25度格子を2倍で描いて平均し、海岸を滑らかにする）。同梱リングは日付変更線で
分割済みで、南極は南極点の辺まで閉じているため、単純な多角形塗りで正しい形になる。
マスクは初回に一度だけ作り、プロセス内でキャッシュする。LOD 0〜2 のデータは
マスクには使わないが、データ保守と後方互換のため同梱したままにする。

太陽の alt/az は従来から合成キャッシュのキーに含まれているため、時刻の変化で
日照が変わると合成画像が作り直される。大きなウィンドウ（2400×1400、半径1100px）
での1回の描画は約0.3秒で、旧実装の輪郭線描画と同程度である。

## SkyWindow ミックスインの共有状態と型チェック

`SkyWindow` は `SkyWindowCoreMixin`（`gui/window.py`）が入力・メニュー操作・各種更新のミックスイン
（`window_input.py`、`window_actions.py`、`window_update_*.py`）を合成して構成する。
ミックスイン間で読み書きする状態（`viewer_data`、`state`、レイヤーの表示フラグと不透明度、
`_geo_satellite_enabled` など）の実体は `SkyWindowCoreMixin.__init__` で初期化する。

- `gui/window_shared_state.py` の `SkyWindowSharedState` は、これらの共有状態の型だけを宣言する
  クラスで、上記ミックスインが継承する。値は持たず、実行時の挙動は変えない。ミックスイン単体でも
  mypy が属性の型を解決できるようにするための契約である。
- 共有状態を新たにミックスインから参照・更新する場合は、`__init__` での初期化に加えて
  `SkyWindowSharedState` に型を宣言する。
- メニューの action 属性（`_action_*`）は `QAction | None` と宣言する。`__init__` で `None` を入れた後に
  `_build_window_menu()` が作成し、利用できないレイヤーの action は `None` のまま残るため、
  参照側は `None` を確認する。
- 描画キャッシュは `_render_cached_frame_image` が属性名の文字列（`_present_frame_cache_key` など）で
  キーと画像を保存する。キーを計算するメソッドはこの属性名と重ならない名前
  （`_render_frame_cache_key`、`_compute_present_frame_cache_key`）にする。同名にすると
  インスタンス属性がメソッドを隠す。
- ミックスインが別のミックスインや `SkyWindowCoreMixin` のメソッドを呼ぶ場合は、`gui/window_shared_methods.py` の
  `SkyWindowSharedMethods`（`SkyWindowSharedState` の基底）にシグネチャを宣言する。宣言は `TYPE_CHECKING` の
  ブロック内に置き、実行時には存在させない。`SkyWindowClientWidget` は描画ミックスインを合成し、
  未定義の属性を `__getattr__` で所有ウィンドウに転送するため、実行時にスタブがあると転送が働かなくなる。
  実体の定義とシグネチャが食い違うと mypy が報告する。
- 描画ミックスイン（`SkyWindowRenderCacheMixin`、`SkyWindowRenderMixin`）は `SkyWindowCoreMixin` と
  `SkyWindowClientWidget` の両方に合成される。実行時には `SkyWindowSharedState` を継承せず、型チェック時だけ
  `SkyWindowSharedState` と `QWidget` を基底として見せる（`window_render_cache.py` の `_RenderMixinBase`）。
  クライアントウィジェットでだけ使う `_owner` は `SkyWindowRenderMixin` に注釈だけで宣言する。
- `SkyWindowInputMixin` と `SkyWindowActionsMixin` は QWidget 派生のウィンドウにのみ合成されるため、
  型チェック時だけ `QWidget` を基底として見せ、`super().eventFilter()` などの協調呼び出しや、
  ウィンドウを親にする `QMenu`/`QAction` の生成を解決する。
- テストでミックスインを継承したスタブに偽のコントローラなどを入れる場合は、スタブのクラスで
  その属性を `Any` と宣言する（例: `tests/test_scheduler_idle.py` の `_SchedulerProbe`）。
- ミックスイン内のメソッド呼び出しは `self` 経由で行い、`SkyWindow.<method>(self, ...)` の
  ようなクラス経由の呼び出しはしない。テストで `SimpleNamespace` などのスタブを `self` に
  渡す場合は、必要な実メソッドをスタブに束縛する。

mypy の設定（`pyproject.toml`）では、src にモジュール単位の除外はない。`tests.*` に限り `attr-defined`、`call-overload`、`arg-type`、`method-assign` を無効化する。
テストはダックタイピングのテストダブルを具象型の引数に渡したり、インスタンスのメソッドを
差し替えたりするため、これらの検査は意図したパターンを報告するだけになる。

## 文書構成

- [cache-sequence-export.md](design/cache-sequence-export.md)
  - 履歴キャッシュの時刻範囲照会（実装済み）と毎分PNG出力・HUD・再開の設計案
- [overview.md](design/overview.md)
  - 設計方針、全体アーキテクチャ、3 つのアプリの共通前提
- [runtime.md](design/runtime.md)
  - スレッドモデル、GUI 状態更新、処理フロー、エラー処理、キャッシュ方針
- [gui-screen-update-and-cache.md](design/gui-screen-update-and-cache.md)
  - GUI の更新トリガー、再描画の流れ、フレーム/合成/状態キャッシュの整理
- [data-model.md](design/data-model.md)
  - 主要データ構造、状態オブジェクト、アプリ間で共有する scene/state の境界
- [rendering-pipeline.md](design/rendering-pipeline.md)
  - 描画パイプライン、オーバーレイ合成、ラベル、外部依存
- [legacy-archive.md](design/legacy-archive.md)
  - 分割前の単一ファイル版の記録。参照用の履歴として残す

利用者向けの機能説明は [specification.md](specification.md) を参照する。

## 起動時のレイヤー分担

- スプラッシュ表示中は、観測地点と観測高度が確定したあとに変化しない静的レイヤーを先行して warm-up してよい。
- この静的レイヤーには、DEM、地形地平線、夜間光と ridge glow、水面、都市アウトラインを含めてよい。
- これらの準静的レイヤーは、既存キャッシュの再利用、欠損時のダウンロード、初回描画用キャッシュの作成までを含めてよい。
- 雲、Geo-satellite 雲、航空機、人工衛星、台風・サイクロンのような動的レイヤーは、通常 GUI が見えてから遅延させてよい。
- 初期データロード完了前のウィンドウリサイズは、初期ロードの世代やキャッシュを破棄せず、最新サイズでの再計算要求を保留する。初期ロード完了後に保留要求を一度だけ処理し、起動中のリサイズ連続通知でスプラッシュ表示が停滞しないようにする。
