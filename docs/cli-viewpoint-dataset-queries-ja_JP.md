# 観測地点名参照ツール

GUI を起動せずに、同梱タワー/展望地点データと山頂ビューポイントデータを参照できます。

| オプション | 説明 | デフォルト |
| :--- | :--- | :--- |
| `-h`, `--help` | ヘルプメッセージを表示して終了します。 | |
| `--list-viewpoints {t,m}` | 同梱タワー (`t`) または山 (`m`) の主表示名を出力して終了します。各行は `t/NAME` または `m/NAME` 形式で、利用可能な場合は ASCII 代替名を優先します。 | |
| `--list-viewpoint-names {t,m}` | 同梱タワー (`t`) または山 (`m`) の名前を、多言語名と ASCII 代替名込みで一覧出力して終了します。各行は `t/NAME` または `m/NAME` 形式です。 | |
| `--show-viewpoint-json NAME` | 指定名で同梱ビューポイントを解決し、利用可能な場合は `ascii_name` を含む JSON メタデータを出力して終了します。`t/` または `m/` を付けると対象 kind を明示できます。 | |

```bash
zstarview --list-viewpoints t
zstarview --list-viewpoint-names t
zstarview --show-viewpoint-json "t/Tokyo Skytree"
zstarview --list-viewpoints m
zstarview --show-viewpoint-json "m/Mount Fuji"
```

これらのオプションは相互排他で、`location` 引数や時刻・描画オプションとは併用できません。
`--list-viewpoints` では、利用可能な場合は ASCII 代替名を優先表示します。
`--list-viewpoint-names` では、元の綴りと ASCII 代替綴りの両方を含みます。
prefix なしの `--show-viewpoint-json` で tower と mountain の両方に完全一致した場合は、`t/...` / `m/...` 候補を列挙して曖昧一致エラーにします。
