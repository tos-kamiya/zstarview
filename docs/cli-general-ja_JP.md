# 一般 CLI オプション

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `-h`, `--help` | ヘルプメッセージを表示して終了します。 | |
| `--window-geometry restore\|X,Y,W,H` | 初期ウィンドウ位置と大きさを指定します。`restore` で前回終了時の位置/サイズを復元し、`X,Y,W,H` で整数値を直接指定できます。Wayland ではウィンドウ位置の復元は利用できません（サイズ復元は有効です）。 | |
| `--window-frame {frameless,window}` | ウィンドウ装飾モードを選びます。`frameless` は従来の枠なし表示、`window` は OS 標準のタイトルバーと枠を使います。 | `frameless` |
| `--observation-info auto\|top\|bottom\|off` | 起動時の観測情報ブロックの表示モードを指定します。 | `bottom` |
| `--include-direction-grid` | `zstarview-export-image` 専用です。出力画像に方向グリッドを含めます。 | |
| `-t`, `--theme {night,day,white,black,transparent,transparent-10..90}` | 背景と星の見え方のテーマを指定します。`transparent` は `transparent-40` の別名で、`transparent-10` から `transparent-90` までは 10 刻みの透明度プリセットです。 | `night` |
| `--visibility-boost MULTIPLIER` | 薄い補助レイヤーの見やすさを持ち上げる倍率です。`1.0` より大きい値で、地形地平線・地球ガイド・都市アウトライン・空/雲ディスク・人工衛星・航空機・地面 tint などの不透明度を底上げします。 | `1.0` |
| `--clear-long-lived-cache` | トラブルシュート用オプションです。起動前に長寿命の DEM / 都市アウトラインキャッシュを削除します。3 日以内に再度使うと起動を拒否し、再実行可能日時を表示します。 | |
