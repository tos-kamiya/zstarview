#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PY_VERSIONS=("3.10" "3.11" "3.12" "3.13")
MYPY_VERSIONS=("3.10" "3.13")
RUN_PYTEST=1
RUN_MYPY=0
RUN_COMPILEALL=1
SKIP_INSTALL=0
MYPY_ALL=0

usage() {
  cat <<'EOF'
Usage:
  scripts/run_python_matrix.sh [options]

Options:
  --pytest-only     Run only pytest across 3.10/3.11/3.12/3.13
  --with-mypy       Also run mypy on representative versions (3.10 and 3.13)
  --mypy-only       Run only mypy on representative versions (3.10 and 3.13)
  --mypy-all        Run mypy across 3.10/3.11/3.12/3.13
  --compileall-only Run only compileall across 3.10/3.11/3.12/3.13
  --skip-install    Reuse existing venvs without reinstalling .[dev]
  -h, --help        Show this help

Default:
  Create/update .venv-3.10 .. .venv-3.13, install .[dev], then run
  pytest and compileall for each version.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pytest-only)
      RUN_PYTEST=1
      RUN_MYPY=0
      RUN_COMPILEALL=0
      ;;
    --mypy-only)
      RUN_PYTEST=0
      RUN_MYPY=1
      RUN_COMPILEALL=0
      ;;
    --with-mypy)
      RUN_MYPY=1
      ;;
    --mypy-all)
      RUN_MYPY=1
      MYPY_ALL=1
      ;;
    --compileall-only)
      RUN_PYTEST=0
      RUN_MYPY=0
      RUN_COMPILEALL=1
      ;;
    --skip-install)
      SKIP_INSTALL=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
  shift
done

echo "==> Ensuring CPython interpreters are available"
uv python install "${PY_VERSIONS[@]}"

for py in "${PY_VERSIONS[@]}"; do
  venv=".venv-$py"
  pybin="$venv/bin/python"

  echo
  echo "==> Preparing environment for CPython $py"
  if [[ ! -x "$pybin" ]]; then
    uv venv "$venv" --python "$py"
  fi

  if [[ "$SKIP_INSTALL" -eq 0 ]]; then
    uv pip install -p "$pybin" -e ".[dev]"
  fi

  if [[ "$RUN_PYTEST" -eq 1 ]]; then
    echo "==> [$py] pytest"
    "$pybin" -m pytest
  fi

  if [[ "$RUN_MYPY" -eq 1 ]]; then
    run_mypy_for_this_version=0
    if [[ "$MYPY_ALL" -eq 1 ]]; then
      run_mypy_for_this_version=1
    else
      for mypy_py in "${MYPY_VERSIONS[@]}"; do
        if [[ "$py" == "$mypy_py" ]]; then
          run_mypy_for_this_version=1
          break
        fi
      done
    fi

    if [[ "$run_mypy_for_this_version" -eq 1 ]]; then
      echo "==> [$py] mypy"
      "$pybin" -m mypy --install-types --non-interactive src/zstarview tests
    fi
  fi

  if [[ "$RUN_COMPILEALL" -eq 1 ]]; then
    echo "==> [$py] compileall"
    "$pybin" -m compileall src/zstarview
  fi
done

echo
echo "Matrix run completed."
