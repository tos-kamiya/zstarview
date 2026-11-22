# Atlas 削除の引き継ぎ

作成日: 2026-11-21。状態: 実施済み（2026-11-21、詳細は session-2026-11-21.md「Atlas の削除」）。利用者の決定は「Atlas 自体を削除する」。
雲の shell 方式削除（`a08c4472`）の際に範囲を調査し、雲と独立した大きな変更なので別作業に分けた。

## 目的と範囲

白背景の別入口アプリ Atlas（`zstarview.gui.atlas`、表示 ID `instrument`）と、Atlas だけが使う
分岐・設定・テーマ・文書を削除する。通常ビューア（`scenic`）と `zstarview-export-image`、replay の
見た目と挙動は変えない。

前提:
- Atlas は配布対象外。`pyproject.toml` の `[project.scripts]` にエントリはなく、
  `docs/specification.md` 冒頭にも「製品配布対象から一時的に除外」と書かれている。
- Atlas の雲描画（`render_atlas_cloud_layer` / `_draw_instrument_cloud_layer`）は
  `a08c4472` で削除済み。現状の Atlas は雲を描かない。

## 削除対象（2026-11-21 時点の調査）

### 入口と専用モジュール
- `src/zstarview/gui/atlas.py`（`apply_atlas_profile`、`_ATLAS_DEFAULTS`、`main`）
- `src/zstarview/render/atlas_pipeline.py`（`InstrumentSkyPresentation` ほか）
- `src/zstarview/render/instrument_background.py`
  - `draw_instrument_background` は Atlas 専用。
  - `draw_instrument_time_of_day_marker` は通常ビューア（`zstarview_pipeline.py`）も使うので残す。
    通常経路は常に `tint_rgba` を渡すため、`tint_rgba=None` の場合の `atlas_background_tint_rgba`
    フォールバックは不要になる。モジュール名・定数名（`ATLAS_TIME_OF_DAY_MARKER_*`）の改名も検討する。
- `render/background.py` の `atlas_background_tint_rgba` と `ATLAS_*_TINT_RGBA`

### 表示 ID と分岐（`presentation_id == "instrument"`）
- `render/pipeline.py`: `_is_instrument_presentation` と、その分岐（base scene の振り分け、
  流星の芯の色、簡易表示のラベル、航空機・惑星への `instrument_presentation` 引数、星ラベルの色）
- `render/solar_system.py`: `instrument_presentation` 引数と `_ATLAS_PLANET_OUTLINE_*`
- `render/aircraft.py`: `instrument_presentation` 引数
- `render/meteors.py`: `METEOR_ATLAS_CORE_COLOR`
- `gui/window_render.py`、`gui/window.py`、`gui/window_updates.py`、`gui/window_update_sky.py`、
  `gui/window_render_cache.py`: `presentation_id` による分岐とキャッシュキー
- `presentation_id` の保持: `gui/window_inputs.py`、`gui/window_display_settings.py`、
  `gui/window_shared_state.py`、`render/render_types.py`、`gui/viewer.py`
  - 「scenic 以外がない」状態になるので、フィールドごと消すか残すかを決める
    （消すとテストのダブルにも波及する）。

### Atlas だけが設定する値
- `star_data_policy`（`positional_static` / `scenic_view_scoped`）: 保持されるだけで挙動に使われていない。
  `window_inputs.py`、`window.py`、`viewer.py`
- `light_background_star_outline`: Atlas の既定値だけが True にする。
  `render/stars.py` の `light_background_outline` 経路（`draw_stars_*` の引数と分岐、1017 行付近）、
  `render/pipeline.py`、`gui/window_render.py`（星サーフェイスのキャッシュキーを含む）、
  `window_inputs.py`、`window_display_settings.py`、`window.py`、`viewer.py`、`render_types.py`
  - 星描画は 1.62.14 で最適化したばかりなので、削除後に星描画テストと見た目を必ず確認する。

### CLI と起動
- `cli/args.py`: `allow_atlas_theme`、`_parse_theme` の Atlas 分岐、`include_scenic_arguments`、
  `vmag_limit_max`、`default_overrides`（Atlas 以外に呼び出し元がないか確認する）
- `gui/viewer.py` の `main` / 内部関数: `apply_app_profile`、`parser_default_overrides`、
  `include_scenic_arguments`、`vmag_limit_max`
- `gui/startup_dialog.py`: `include_twinkle_options`（Atlas だけが False で呼ぶ）

### テーマ
- `paths.py`: `ATLAS_THEME_PRESET`（`atlas-white`）、`ATLAS_GUIDE_STYLE`、`ATLAS_OVERLAY_STYLES`、
  `PALETTE_ATLAS_DIRECTION_GUIDE_RGB`、`THEME_STYLES_BY_PRESET` の Atlas エントリ
- `OverlayStyles.cloud` と `CloudLayerStyle`: shell 雲削除の後は Atlas の雲描画専用だった。
  現在は参照する描画コードがないので、このとき一緒に削除できる。

### テスト
Atlas を直接扱うもの: `tests/test_atlas_profile.py`、`tests/test_atlas_pipeline_layers.py`。
Atlas の分岐やテーマを含むもの（該当テストだけを削る）: `test_window_render_context.py`、
`test_window_render_hud.py`、`test_theme_window_background.py`、`test_earth_guide_render.py`、
`test_color_palette.py`、`test_startup_log.py`、`test_planet_marker_solar.py`、
`test_planet_marker_artificial.py`、`test_deep_sky_objects_render.py`、`meteors/test_render.py`、
`test_sky_disc_viewport_and_ground_tint.py`、`test_cli_arg_helpers.py`、
`gui/test_startup_dialog_helpers.py`（`test_startup_dialog_hides_twinkle_for_atlas`）。

### 文書
- 削除: `docs/atlas-specification.md`、`docs/design/atlas.md`
- 参照の整理: `docs/specification.md`（8 行目の注記、275、580、593、598、601、724 行付近）、
  `docs/design.md`（133、665、733 行付近）、`docs/design/rendering-pipeline.md`（74〜82、405 行付近）、
  `docs/design/meteor-trails.md`（30、69 行）、`docs/cli-celestial.md` / `-ja_JP.md`（asterism opacity の Atlas 言及）
- `docs/implementation-archive.md` は更新しない（AGENTS.md）。
- Atlas は配布対象外なので、release notes は原則不要。CLI の受理値が変わる場合だけ判断する。

## 進め方の提案

1. 入口・専用モジュール・Atlas テストを削除し、`presentation_id` の分岐を通常経路に固定する。
2. `light_background_star_outline` と `star_data_policy` を削除する（星描画テストと見た目を確認）。
3. CLI・viewer・起動ダイアログの Atlas 用フックを削除する。
4. テーマ（`atlas-white`、`CloudLayerStyle`）と文書を整理する。

各段階で次を通す:
`QT_QPA_PLATFORM=offscreen uv run -p .venv/bin/python pytest`、
`uv run -p .venv/bin/python mypy src/zstarview tests`、`ruff check src/zstarview tests`
（既存の E402 9件は対象外）。

## 関連して残っている作業

- Geo-satellite パイプライン（`geosatellite/pipeline.py`）が、shell 時代の `disc_gray` 投影を
  まだ計算している。結果はどこからも使われていない。`project_gray_image_to_disc` と
  `geosatellite/cache.py` の投影キャッシュも含めて削除できる。
- README の screenshot10 は廃止した halftone 雲の画像なので、撮り直しが必要。
- 雲描画の並列化（`_render` の `prange` 化）は利用者の判断で見送った（2026-11-21）。
  試作の測定値は `dev-notes/session-2026-11-21.md`「雲レイヤー生成時間の調査」にある。
