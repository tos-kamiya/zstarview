# CPython Test Matrix

`zstarview` を CPython 3.10 / 3.11 / 3.12 / 3.13 / 3.14 で確認するための手順です。

`uv` を使い、各 Python バージョンごとに独立した仮想環境を作ってテストします。

自動実行用スクリプト:

```bash
scripts/run_python_matrix.sh
```

このスクリプトは `.venv-3.10` 〜 `.venv-3.14` を作成/再利用し、各 venv の `python` を直接使って既定で `pytest` と `compileall` を順に実行します。

`mypy` は既定では実行しません。必要なときだけ代表バージョン（3.10 と 3.14）または全バージョンで追加実行します。

## 前提

- `uv` がインストール済み
- リポジトリルート（`zstarview/`）でコマンド実行
- 対象は CPython:
  - `3.10`
  - `3.11`
  - `3.12`
  - `3.13`
  - `3.14`

## 1. 対象 Python を用意する

`uv` 管理の Python をまだ入れていない場合:

```bash
uv python install 3.10 3.11 3.12 3.14
```

確認:

```bash
uv python list
```

## 2. バージョンごとの仮想環境を作る

各バージョンごとに専用 venv を作ります。

```bash
uv venv .venv-3.10 --python 3.10
uv venv .venv-3.11 --python 3.11
uv venv .venv-3.12 --python 3.12
uv venv .venv-3.14 --python 3.14
```

作り直す場合:

```bash
uv venv --clear .venv-3.10 --python 3.10
uv venv --clear .venv-3.11 --python 3.11
uv venv --clear .venv-3.12 --python 3.12
uv venv --clear .venv-3.14 --python 3.14
```

## 3. 各環境に開発依存を入れる

```bash
uv pip install -p .venv-3.10/bin/python -e ".[dev]"
uv pip install -p .venv-3.11/bin/python -e ".[dev]"
uv pip install -p .venv-3.12/bin/python -e ".[dev]"
uv pip install -p .venv-3.14/bin/python -e ".[dev]"
```

## 4. テストを回す

最小構成の確認:

```bash
uv run -p .venv-3.10/bin/python pytest
uv run -p .venv-3.11/bin/python pytest
uv run -p .venv-3.12/bin/python pytest
uv run -p .venv-3.14/bin/python pytest
```

必要に応じて型チェックも回します。

代表バージョン（推奨）:

```bash
uv run -p .venv-3.10/bin/python mypy --install-types --non-interactive src/zstarview tests
uv run -p .venv-3.14/bin/python mypy --install-types --non-interactive src/zstarview tests
```

全バージョンで回したい場合:

```bash
uv run -p .venv-3.10/bin/python mypy --install-types --non-interactive src/zstarview tests
uv run -p .venv-3.11/bin/python mypy --install-types --non-interactive src/zstarview tests
uv run -p .venv-3.12/bin/python mypy --install-types --non-interactive src/zstarview tests
uv run -p .venv-3.14/bin/python mypy --install-types --non-interactive src/zstarview tests
```

簡易コンパイル確認:

```bash
uv run -p .venv-3.10/bin/python -m compileall src/zstarview
uv run -p .venv-3.11/bin/python -m compileall src/zstarview
uv run -p .venv-3.12/bin/python -m compileall src/zstarview
uv run -p .venv-3.14/bin/python -m compileall src/zstarview
```

## 5. 推奨チェック順

普段の軽い確認:

1. 変更を 1 つの基準環境で確認する
2. 仕上げ前に 3.10〜3.14 の `pytest` を一通り回す
3. 必要に応じて `compileall` を回す
4. 型を確認したい変更では、代表バージョン（3.10 と 3.14）で `mypy` を回す

リリース候補前の確認:

1. 4 バージョンすべてで `pytest`
2. 4 バージョンすべてで `compileall`
3. 代表バージョン（通常は 3.10 と最新 3.14）で `mypy`
4. 代表バージョンでアプリ起動確認

## 6. 一括実行例

シェルでまとめて回す例:

```bash
for py in 3.10 3.11 3.12 3.14; do
  uv run -p ".venv-$py/bin/python" pytest
done
```

`mypy` も同様:

```bash
for py in 3.10 3.14; do
  uv run -p ".venv-$py/bin/python" mypy --install-types --non-interactive src/zstarview tests
done
```

スクリプト利用例:

```bash
scripts/run_python_matrix.sh
scripts/run_python_matrix.sh --with-mypy
scripts/run_python_matrix.sh --mypy-only
scripts/run_python_matrix.sh --mypy-all
scripts/run_python_matrix.sh --skip-install
```

## 7. 注意点

- `uv.lock` を更新した後は、各環境で `uv pip install -p ... -e ".[dev]"` を再実行して依存を揃える。
- `PySide6`, `rasterio`, `satpy` などは環境差分の影響を受けやすいので、リリース前は最低でも 1 回は新しい venv で確認する。
- 失敗が 1 バージョンだけに出る場合は、まずその環境の `pip install` をやり直してから再確認する。

## 8. 関連ドキュメント

- `docs/developer/uv-workflow.md`
- `docs/specification.md`
- `docs/design.md`
