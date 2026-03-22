# uv Workflow

`zstarview` を `uv` 環境で開発・実行するための手順です。

## 前提

- `uv` がインストール済み
- リポジトリルート（`zstarview/`）でコマンド実行

## 1. 環境構築

```bash
uv venv --python 3.12
uv pip install -p .venv/bin/python -e ".[dev]"
```

補足:
- 既存環境を作り直す場合は `uv venv --clear --python 3.12` を使用。
- Python 3.10+ が対象です。
- `.[dev]` を付けると開発用依存（例: `pytest`）も同時に導入されます。

## 2. アプリ実行

```bash
uv run -p .venv/bin/python zstarview
```

モジュール実行:

```bash
uv run -p .venv/bin/python -m zstarview.gui.viewer
```

例:

```bash
uv run -p .venv/bin/python zstarview Tokyo
uv run -p .venv/bin/python zstarview \"35.68;139.76\" --datetime \"2026-02-26 21:00 JST\"
```

## 3. テストと検証

テスト:

```bash
uv run -p .venv/bin/python pytest
```

カバレッジ:

```bash
uv run -p .venv/bin/python coverage run -m pytest
uv run -p .venv/bin/python coverage report
```

型チェック:

```bash
uv run -p .venv/bin/python mypy --install-types --non-interactive src/zstarview tests
```

簡易コンパイル確認:

```bash
uv run -p .venv/bin/python -m compileall src/zstarview
```

## 4. ビルド

`build` が未導入なら追加:

```bash
uv pip install -p .venv/bin/python build
```

配布物作成:

```bash
uv run -p .venv/bin/python -m build
```

## 5. トラブルシューティング

- `pytest` で `botocore` が見つからない場合:
  - 依存が未インストールの可能性があるため、再度 `uv pip install -p .venv/bin/python -e ".[dev]"` を実行。
- Linux/X11 で Qt 起動に失敗する場合:
  - `README.md` の X11 セクション（`libxcb-cursor0`）を参照。

## 6. 関連ドキュメント

- 星カタログ再生成: `docs/developer/star-catalog-generation.md`
