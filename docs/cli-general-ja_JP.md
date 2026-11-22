# 一般 CLI オプション

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `-h`, `--help` | ヘルプメッセージを表示して終了します。 | |
| `--window-geometry restore\|X,Y,W,H` | 初期ウィンドウ位置と大きさを指定します。`restore` で前回終了時の位置/サイズを復元し、`X,Y,W,H` で整数値を直接指定できます。Wayland ではウィンドウ位置の復元は利用できません（サイズ復元は有効です）。 | |
| `--window-frame {frameless,window}` | ウィンドウ装飾モードを選びます。`frameless` は従来の枠なし表示、`window` は OS 標準のタイトルバーと枠を使います。 | `frameless` |
| `--observation-info auto\|top\|bottom\|off` | 起動時の観測情報ブロックの表示モードを指定します。 | `bottom` |
| `--initial-display-mode MODE` | GUIの初期表示モードを選び、Spaceによる表示状態循環の開始位置を決めます。`normal`（別名`default`）、`inverted-city`（`inverted`）、`simple-no-labels`（`simple`）、`simple-labels`（`labels`）、`landscape`を指定できます。Urban Outlineが利用できない場合、`Inverted City`は`normal`へフォールバックします。このオプションは`zstarview`専用です。 | `normal` |
| `--include-direction-grid` | `zstarview-export-image` 専用です。出力画像に方向グリッドを含めます。 | |
| `-t`, `--theme {night,transparent,transparent-10..90}` | ウィンドウ背景と空色ディスクの不透明度のテーマを指定します。`night` は既定の暗いテーマです。`transparent` は `transparent-40` の別名で、`transparent-10` から `transparent-90` までは 10 刻みの透明度プリセットです。 | `night` |
| `--visibility-boost MULTIPLIER` | 薄い補助レイヤーの見やすさを持ち上げる倍率です。`1.0` より大きい値で、地形地平線・地球ガイド・都市アウトライン・空/雲ディスク・人工衛星・航空機・アステリズムなどの不透明度を底上げします。 | `1.0` |
| `--overlay-font-size POINTS` | ウィンドウ内に描画するラベルとHUD文字だけの基本フォントサイズを指定します。ステータス行・メニュー・ダイアログ・標準Qtウィジェットには影響しません。 | `11` |
| `--offline` | ネットワークに接続せずに起動します。長期間キャッシュするデータ（地形、都市アウトライン、夜間光、道路の夜間光、水面）はキャッシュがあれば使い、なければ描きません。航空機、人工衛星、予報降水量、衛星雲画像、熱帯低気圧、メテオ軌跡、月と太陽のホバー画像は使いません。場所を省略すると前回の場所、なければ `Tokyo` を使い、`--place` と検索はキャッシュと手元のデータだけを使います。`zstarview-export-image` と起動ダイアログ（`Offline`）でも指定できます。 | |
| `--clear-long-lived-cache` | トラブルシュート用オプションです。起動前に長寿命の DEM / 都市アウトラインキャッシュを削除します。3 日以内に再度使うと起動を拒否し、再実行可能日時を表示します。 | |

#### 脚注

※3 Geo-satellite が実験扱いなのは、配信元の画像が生の雲データではなく人間向けに加工された表示用プロダクトだからです。海岸線・国境・マスク領域などに穴や欠損ができることがあり、inpainting などの補完や代替処理が必要になるため、標準の GOES/Himawari 経路ほど堅牢ではありません。
