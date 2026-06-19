# Contributing

Thank you for helping improve zstarview!

- Start here: read `AGENTS.md` for the full contributor guide (structure, commands, style, testing, and PR norms).

## Quick Start
- Create venv + editable install:
  - `python -m venv .venv && source .venv/bin/activate`
  - `pip install -U pip && pip install -e .`
- Run locally: `zstarview [options] [city]` or `python -m zstarview.zstarview`.

## Before You Commit
- Type check: `mypy --install-types --non-interactive src/zstarview tests`.
- Tests + coverage: `coverage run -m pytest && coverage report`.
- Keep GUI code responsive; move heavy work off the UI thread.
- Follow naming and layout in `AGENTS.md` (e.g., tests in `tests/test_*.py`).

## Commits & Pull Requests
- Commits: use Conventional Commits (e.g., `feat:`, `fix:`, `docs:`, `chore:`).
- PRs: include a clear description, reference issues (e.g., `Closes #123`), provide repro/validation steps, and screenshots/short clips for UI changes. Keep scope focused.

For additional build and packaging steps, see `AGENTS.md` (type checking, build artifacts, and GNOME launcher helpers).
