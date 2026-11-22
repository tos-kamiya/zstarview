# mypy の除外設定の段階的解除の引き継ぎ

作成日: 2026-11-21。状態: 除外設定の解除は完了。src にモジュール単位の mypy 除外は残っていない。
PySide6 の stub と実行時の食い違いは調査済み（下の「PySide6 の stub と実行時の食い違い」）。
`QImage.save` は解決済み。`setCurrentItem(None)` の `type: ignore` だけが残っている。

## 目的と範囲

`pyproject.toml` の mypy 設定に残る除外（`ignore_errors` と `attr-defined` などの無効化）を段階的に外し、
SkyWindow まわりも型チェックの対象にする。各段階で `mypy src/zstarview tests` をエラー 0 件に保ち、
段階ごとに独立したコミットにする。実行時の挙動は変えない。挙動に関わるバグを見つけた場合は、
回帰テストと一緒に別コミットで直す（前回もこの方針で3件のバグを修正した）。

範囲外: `tests.*` で無効にしている `arg-type`、`method-assign`、`attr-defined`、`call-overload`
（テストダブルの意図したパターンのため。理由は `docs/design.md` の「SkyWindow ミックスインの共有状態と型チェック」）。

## 前提（2026-11-21 時点）

- `9589f506 refactor: resolve mypy errors and warnings` で、現行設定ではエラー・注記とも 0 件。
  経緯と判断理由は `dev-notes/session-2026-11-21.md`「mypy のエラー・警告の解消」。
- 現在の除外設定:
  - 全体: `disable_error_code = ["import-untyped", "attr-defined", "call-overload"]`
  - `zstarview.gui.*`、`zstarview.render.*`: `attr-defined` を無効化（全体と重複）
  - `zstarview.render.draw`、`zstarview.gui.window`、`zstarview.gui.sky_worker`: `ignore_errors = true`
- `gui/window_shared_state.py` の `SkyWindowSharedState` は、ミックスイン間の共有状態のうち、
  前回 `has-type` が出た約25属性だけを宣言している。

## 実測値（設定を一時コピーして計測。`pyproject.toml` は未変更）

| シナリオ | 変更内容 | エラー数 |
|---|---|---|
| A | `ignore_errors` の3モジュールを外す | 44件（2ファイル） |
| B | A + src 全体で `attr-defined` を有効化 | 891件（32ファイル） |
| C | B + src 全体で `call-overload` を有効化 | 915件（40ファイル） |

### A の内訳（44件）

- `render/draw.py`: 0件。すぐに除外を外せる。
- `gui/sky_worker.py`: 3件。`image_size` の `tuple[int, ...]`、`DeepSkyTable` と `dict[str, ndarray]` の不一致。
- `gui/window.py`: 41件。
  - 22件: `self._action_* = ... | None` の代入。`SkyWindowActionsMixin` の推論では `QAction` だが、
    `window.py` では `None` を入れている。メニューの action 属性の型をどちらかに宣言して揃える。
  - 2件: `window.py:677` で、`_present_frame_cache_key` のインスタンス属性が
    `SkyWindowRenderCacheMixin` の同名メソッドを上書きしている。このため src 側は
    `SkyWindowRenderMixin._present_frame_cache_key(self, ...)` とクラス経由で呼んでいる。
    属性かメソッドの改名が必要（テストも追随が必要）。
  - 7件: `X | None` を `X` 型の変数に代入（1963、2009、2062-2064 行など）、`tuple[float, ...]`（1127 行）。
  - 3件: `has-type`（`_client_geometry_sync_done`、`_landscape_annotation_opacity`、`_landscape_fade_was_active`）。
  - `window.py:208`: 既存の `type: ignore[arg-type]` では `type-var` と `return-value` を覆えていない。

### B の内訳（`attr-defined` 847件）

- 大半はミックスインの `self` 上の属性: `SkyWindowRenderMixin` 133、`SkyWindowActionsMixin` 110、
  `SkyWindowRenderCacheMixin` 64、`SkyWindowSkyUpdatesMixin` 61、`SkyWindowCloudUpdatesMixin` 61、
  `SkyWindowInputMixin` 57、`SkyWindowUpdatesMixin` 32 など。
  ミックスインから参照される未宣言の属性名は197種類。
- ファイル別の上位: `window_update_overlays.py` 155、`window_render.py` 139、`window_actions.py` 110、
  `window_render_cache.py` 66、`window_update_sky.py` 61、`window_update_cloud.py` 61、`window_input.py` 57。
- `object` に対する属性参照 74件。主に `cli/export_image_inputs.py`（69件）で、`args: object` を
  受け取って属性を読んでいる。`argparse.Namespace` か専用の型にすれば大半が消える見込み。
- そのほか `gui/composite.py` 15、`gui/viewer.py` 13、`gui/cloud_render.py` 8 など少数。

### C の追加分（`call-overload` 24件）

16ファイルに散在。`window_actions.py` に7件、他は1-2件ずつ。

## 進捗（2026-11-21 後半）

- 完了: `ignore_errors` は全廃。`attr-defined` と `call-overload` も src 全体で有効（ミックスインを含む）。
  ミックスインは選択肢 1 で対応し、`gui/window_shared_state.py`（属性）と `gui/window_shared_methods.py`
  （`TYPE_CHECKING` 限定のメソッドスタブ）に宣言した。詳細と判断理由は `dev-notes/session-2026-11-21.md`、
  設計上の規則は `docs/design.md`「SkyWindow ミックスインの共有状態と型チェック」。
- 残る除外は `tests.*` の override（`attr-defined`、`call-overload`、`arg-type`、`method-assign`）だけ（範囲外）。
- 以下の「方針の選択肢」「提案する進め方」「実測値」は作業前の記録として残している。

## PySide6 の stub と実行時の食い違い（調査済み、2026-11-21）

### 調査方法

PySide6 6.10.0、6.10.3、6.11.0、6.11.3、6.12.0（2026-11-21 時点で PyPI の最新）を、
`QT_QPA_PLATFORM=offscreen uv run --no-project --python 3.12 --with "PySide6==X" python -I check.py` で切り替えて確認した。
stub（`QtGui.pyi`、`QtWidgets.pyi`）は 6.10.0 と 6.12.0 を比べた。

### 不整合 1: `QImage.save(fileName, format)` の format の型

- stub（調べた全バージョンで同じ）:
  `def save(self, fileName: str, /, format: bytes | bytearray | memoryview | None = ..., quality: int = ...) -> bool`
- 実行時（調べた全バージョンで同じ）:

  | 呼び出し | 結果 |
  |---|---|
  | `im.save("t.png", "PNG")`（str） | True。stub では型エラー（`call-overload`） |
  | `im.save("t.png", b"PNG")`（bytes） | `ValueError: 'PySide6.QtGui.QImage.save' called with wrong argument values`。stub どおりなのに拒否される |
  | `im.save("t.png")`（format を省略） | True（拡張子で判定） |
  | `im.save("t.png.tmp")`（format を省略） | False（拡張子から判定できない） |

- 同じ QtGui の `QPixmap.save` の stub は `format: str | None` で、実行時の挙動と合っている。
  食い違っているのは `QImage.save` の stub だけと見られる。
- 型チェックと実行時の両方を満たす `QImage.save(path, format)` の書き方はない。
- 対応（`08aac648`）: `QImageWriter(path, b"PNG").write(image)` に置き換えた。
  - `QImageWriter` は stub も実行時も format を bytes で受け取り、両者が一致している。
  - str を渡すと `QImageWriter(path, "PNG")` は実行時に `TypeError` になり、stub でも型エラーになるので、取り違えは mypy で検出できる。
  - 拡張子に関係なく format を指定できるので、`.tmp` の一時ファイルにも使える。
  - 共通ヘルパーは `render/qt_image.save_qimage_png()`。上流の stub が直っても、この書き方を変える必要はない。
- 対象外: device を受け取るオーバーロードの `QImage.save(QBuffer, "PNG")`（テストで使用）。
  型エラーにならず、実行時にも動く。

### 不整合 2: `QListWidget.setCurrentItem(None)`

- stub（全バージョンで同じ）: `def setCurrentItem(self, item: QListWidgetItem, /) -> None`（`None` を受け付けない）。
- 実行時（全バージョンで同じ）: `None` を渡すと選択解除として動く（`currentRow() == -1`、選択数 0）。
- 同じ効果を持つ代替として `setCurrentRow(-1)` がある。結果も、発行されるシグナルとその順序も同じだった
  （`itemSelectionChanged` → `currentItemChanged(None, 前の項目)` → `currentRowChanged(-1)`）。
- 対応: 利用者の判断で、当面は `type: ignore[call-overload]  # stubs omit None` のまま残す。
  対象は `gui/place_search_dialog.py`（1か所）と `gui/famous_star_search_dialog.py`（2か所）。
  外すときは `setCurrentRow(-1)` に置き換える。

### 上流の状況

どちらの不整合についても、そのものを扱った issue は見つからなかった（2026-11-21 時点）。類似の報告:
- PYSIDE-1968: `QPixmap.save` に bytes の format を渡すと同じ `ValueError` が出る件。6.3.0 で修正済み。`QImage` は対象外。
- PYSIDE-2235: `QWidget.setParent(None)` を stub が拒否していた件。`Optional` を追加して修正された。
  `setCurrentItem` も同じ種類の問題。

### 調査前の記録

PySide6 6.11.1 で、次の2点が stub と実行時で食い違っていた。現在はどちらも
`type: ignore[call-overload]` と理由のコメントで対応している。

- `QImage.save(path, format)`: stub は `format: bytes | bytearray | memoryview | None` だが、
  6.11.1 は `b"PNG"` を `ValueError`（wrong argument values）で拒否し、`"PNG"`（str）だけ受け付けた。
  対象: `src/zstarview/replay/sequence.py`（テストなし）、`scripts/prototype_cloud_cell_shading.py`。
- `QListWidget.setCurrentItem(None)`: stub は `None` を受け付けないが、6.11.1 では選択解除として動作する。
  対象: `gui/place_search_dialog.py`、`gui/famous_star_search_dialog.py`。

`pyproject.toml` の対応範囲は `PySide6>=6.10,<7` のため、次を確認する。

- 6.10 系（下限）と最新の 6.x で、`QImage.save` が str と bytes のどちらを受け付けるか。
  str を拒否するバージョンがあれば、バージョンに依存しない書き方（例: format 引数を省き拡張子で判定させる、
  `QImageWriter` を使う）を検討する。`replay/sequence.py` は一時ファイル名が `.tmp` で終わるため、
  拡張子での判定には変更が必要。
- 同じく `setCurrentItem(None)` の挙動。
- stub 側の修正状況（PySide6 の issue tracker）。修正されたら `type: ignore` を外す。
- 確認用の最小コード（offscreen で実行）:

  ```python
  from PySide6.QtGui import QImage
  from PySide6.QtWidgets import QApplication, QListWidget
  app = QApplication([])
  im = QImage(4, 4, QImage.Format.Format_RGBA8888)
  print(im.save("/tmp/t.png", "PNG"))
  try:
      print(im.save("/tmp/t.png", b"PNG"))
  except ValueError as exc:
      print("bytes rejected:", exc)
  lw = QListWidget(); lw.addItem("a"); lw.setCurrentRow(0); lw.setCurrentItem(None); print(lw.currentRow())
  ```

  `QT_QPA_PLATFORM=offscreen uv run --with "PySide6==6.10.*" python check.py` のように
  バージョンを切り替えて実行する。
- テストがない `replay/sequence.py` の PNG 保存に回帰テストを追加することも検討する。

## 方針の選択肢（次のセッションで最初に決める）

ミックスインの `attr-defined`（800件超）の扱いが最大の論点。

1. **`SkyWindowSharedState` の拡張**（前回採用した方式の延長）: 197属性とミックスイン間で呼ぶメソッドを
   宣言する。メソッドは `TYPE_CHECKING` 下の宣言か Protocol で表す。実行時の挙動は変わらない。
   属性が多いので、分類（状態・コントローラ・Qt widget・タイマー・action など）してファイルを分ける余地がある。
2. **ホスト Protocol と `self` 型注釈**: `class SkyWindowHost(Protocol)` を定義し、ミックスインのメソッドに
   `self: SkyWindowHost` を付ける。契約は明確になるが、注釈の追加が全メソッドに及ぶ。
3. **ミックスインの分解**: 状態やコントローラを別オブジェクトにして委譲する。設計としては最も良いが、
   大規模な挙動変更のリスクがあり、mypy 対応の範囲を超える。前回は見送った。

前回の相談では 1 を選んだ（2 と 3 は規模が大きいため見送り）。続ける場合も 1 が自然だが、
宣言が200行規模になる点を利用者に確認してから進める。

## 提案する進め方

各段階でエラー 0 件・全テスト成功を確認してコミットする。

1. `render/draw.py` の `ignore_errors` を外す（0件のため設定変更のみ）。
2. `gui/sky_worker.py` の `ignore_errors` を外し、3件を直す。
3. `gui/window.py` の `ignore_errors` を外し、41件を直す。`_present_frame_cache_key` の名前衝突と
   `_action_*` の型の統一は設計判断を含むので、方針を利用者に説明してから進める。
4. `cli/export_image_inputs.py` の `args: object` を型付けする（`attr-defined` 69件の解消に向けた準備）。
5. src 全体の `attr-defined` 有効化。ミックスインの扱い（上の選択肢）を決めてから、
   ファイル単位またはミックスイン単位で進める。途中段階では、対象外のモジュールだけ
   overrides で無効化を残し、徐々に範囲を広げる。
6. `call-overload` の有効化（24件）。

## 計測方法（再現用）

`pyproject.toml` の `[tool.mypy]` と overrides を一時ファイルにコピーし、除外を外したうえで実行する。

```bash
.venv/bin/python -m mypy --config-file <tmp>.toml --cache-dir <tmp-cache> src/zstarview tests
```

`tests.*` の override（`attr-defined`、`call-overload`、`arg-type`、`method-assign` の無効化）は残して計測した。

## 注意点

- AGENTS.md: 横断的な変更は、範囲とリスクを説明して確認を取ってから進める。
  `dev-notes/session-YYYY-MM-DD.md` に判断理由を残す。
- 前回、ミックスインの呼び出しを `self` 経由に変えたところ、`SimpleNamespace` のスタブを使うテストが
  失敗した。テスト側は `_bind_window_methods`（`tests/terrain/test_terrain_horizon_ui.py`）や
  `_WindowStub`（`tests/_window_render_support.py`）の委譲メソッドで対応している。
- `ruff check --select I --fix` をディレクトリ単位で実行すると、無関係なファイルの import 順も
  変わる。対象ファイルだけに実行する。
