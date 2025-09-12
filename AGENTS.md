# Repository Guidelines

## Project Structure & Module Organization
- `src/zstarview/`: Main package. CLI entry in `zstarview.py`; desktop helper in `make_desktop_file.py`.
- `src/zstarview/data/`: Runtime assets (catalogs, fonts, icons).
- `tests/`: Add `test_*.py` here.
- `docs/`: Project documentation and images.
- Root: `pyproject.toml`, `uv.lock`, `README.md`.

## Build, Test, and Development Commands
- Create venv + editable install:
  - `python -m venv .venv && source .venv/bin/activate`
  - `pip install -U pip && pip install -e .`
- Run the app locally: `zstarview [options] [city]` or `python -m zstarview.zstarview`.
- Type check: `mypy --install-types --non-interactive src/zstarview tests`.
- Tests + coverage: `coverage run -m pytest && coverage report`.
- Build wheel/sdist: `python -m build` (requires `pip install build`).
- GNOME launcher: `zstarview-make-desktop-file [--write]`.

## Coding Style & Naming Conventions
- Python 3.10+. Follow PEP 8 with 4‑space indentation; add type hints where practical.
- Names: modules/files `snake_case.py`; functions/variables `snake_case`; classes `CamelCase`; constants `UPPER_SNAKE`.
- Keep the GUI responsive: avoid heavy work on the UI thread; move I/O/compute off the GUI thread. Place assets under `src/zstarview/data/`.

## Testing Guidelines
- Framework: pytest. Focus on pure logic (coordinate transforms, projections, phase angle math).
- Location & names: put tests in `tests/` as `test_*.py`.
- Rendering: skip or mock PyQt5 UI code. Run `pytest` or use the coverage command above.

## Commit & Pull Request Guidelines
- Commits: use Conventional Commits (e.g., `feat:`, `fix:`, `docs:`, `chore:`).
- Pull requests: include a clear description and rationale; reference issues (e.g., `Closes #123`); add repro/validation steps and screenshots or short clips for UI changes; keep scope small and focused.

## Security & Configuration Tips
- App config: `~/.config/zstarview/config.json` (auto‑managed). Do not commit local configs.
- Cache: ephemeris and data cached via `appdirs`; safe to delete when troubleshooting.

## Session Logging (Raw‑ish Notes in `dev-notes/`)
- Goal: Leave a human‑readable, raw‑ish session record under `dev-notes/session-YYYY-MM-DD.md` (not a full transcript, but close).
- When starting work:
  - Create `dev-notes/session-YYYY-MM-DD.md` and add a short header (date, scope).
- While working (for each meaningful action):
  - Append a small block capturing what was done with a transcript feel:
    - Time (HH:MM)
    - Ran: the exact command (in backticks)
    - Output (excerpt) as a fenced block, keep raw lines; truncate if too long
    - Optional: Why / Result / Next (1–2 lines each)
  - Example snippet to append:
    
    ```
    - 10:22 Ran: `git log --oneline -n 3`
      Output:
      ```text
      abcd123 feat: add X
      efgh456 fix: Y
      ....
      ```
      Why: Inspect recent commits
      Result: Verified last change
    ```

- Group related commands under a short heading when it helps scanability.
- Redact sensitive info (tokens, secrets, private paths) if they appear in output.
- Prefer brevity for rationale; keep the output raw‑ish.
- At the end, add a short summary (commits, tags, pushes, follow‑ups).
- Optional: If you also capture a full terminal log (`script`, `tee`, etc.), reference the file at the top of the markdown and paste only key excerpts.

Notes:
- Historical session logs may exist under `docs/` (e.g., `docs/session-YYYY-MM-DD.md`). New logs should go to `dev-notes/`.
