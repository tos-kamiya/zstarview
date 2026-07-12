# Repository Guidelines

## Project Structure & Module Organization
- `src/zstarview/`: Main package. GUI entry in `gui/viewer.py`; CLI tools, including the desktop helper, are under `cli/`.
- `src/zstarview/data/`: Runtime assets (catalogs, fonts, icons).
- `tests/`: Add `test_*.py` here.
- `docs/`: Project documentation and images.
- Root: `pyproject.toml`, `uv.lock`, `README.md`.

## Documentation Roles
- `docs/specification.md`: User-facing functional specification. Keep it more detailed than `README.md`, but focused on externally visible behavior, supported inputs, UI behavior, constraints, and failure behavior.
- `docs/design.md`: Internal design reference. Document architecture, module responsibilities, processing flows, thread model, and key data structures here.
- `release-notes.md`: User-facing public release history. Keep this focused on shipped behavior changes.
- `dev-notes/session-YYYY-MM-DD.md`: Internal development history. Use this for dated decision logs, rationale, command transcripts, TODO, and INPROGRESS notes.
- `docs/implementation-archive.md`: Historical archive of older implementation notes. Do not add new entries here or use it as a source of truth for current behavior.
- When updating docs, preserve this separation of concerns instead of mixing user-visible behavior, internal design, public release history, internal session history, and archive material in the same file.

## Build, Test, and Development Commands
- Create or update the development environment:
  - `uv venv --python 3.12`
  - `uv pip install -p .venv/bin/python -e ".[dev]"`
- Run the app locally: `uv run -p .venv/bin/python zstarview [options] [city]` or
  `uv run -p .venv/bin/python -m zstarview.gui.viewer`.
- Type check: `uv run -p .venv/bin/python mypy --install-types --non-interactive src/zstarview tests`.
- Tests + coverage: `uv run -p .venv/bin/python coverage run -m pytest && uv run -p .venv/bin/python coverage report`.
- Build wheel/sdist: `uv run -p .venv/bin/python -m build` (install `build` with `uv pip install` if needed).
- GNOME launcher: `uv run -p .venv/bin/python zstarview-make-desktop-file [--write]`.
- Import cleanup: use `uv run -p .venv/bin/python ruff check --select I --fix src/zstarview tests` for
  import sorting, then run `uv run -p .venv/bin/python ruff check src/zstarview tests` to catch any
  residual import-related errors such as `E402` in path-hacked scripts.

## Coding Style & Naming Conventions
- Python 3.10+. Follow PEP 8 with 4‑space indentation; add type hints where practical.
- Names: modules/files `snake_case.py`; functions/variables `snake_case`; classes `CamelCase`; constants `UPPER_SNAKE`.
- Keep the GUI responsive: avoid heavy work on the UI thread; move I/O/compute off the GUI thread. Place assets under `src/zstarview/data/`.
- Cross-platform text safety:
  - Any text that can be written directly to a terminal, console, log, CLI help, exception message, or subprocess stdout/stderr should be ASCII-only unless there is a documented reason not to.
  - If code needs to recognize non-ASCII text, prefer Unicode escape sequences like `"\u2019"` instead of embedding the character literally in source.
  - UI-only strings may use non-ASCII when required and validated on supported platforms.

## Testing Guidelines
- Framework: pytest. Focus on pure logic (coordinate transforms, projections, phase angle math).
- Location & names: put tests in `tests/` as `test_*.py`.
- Rendering: skip or mock PySide6 UI code. Run `pytest` or use the coverage command above.

## Commit & Pull Request Guidelines
- Commits: use Conventional Commits (e.g., `feat:`, `fix:`, `docs:`, `chore:`).
- Pull requests: include a clear description and rationale; reference issues (e.g., `Closes #123`); add repro/validation steps and screenshots or short clips for UI changes; keep scope small and focused.
- When bumping the package version, add or update the corresponding entry in `release-notes.md` in the same change unless the bump is purely internal and intentionally undocumented.
- When drafting release notes, review the change set back to the previous version-bump commit so the entry covers the full release window, not only the immediate parent commit.

### Commit Message Format
- Subject: a single line using Conventional Commits (imperative, concise). Example:
  - `feat: improve robust PDF text matching`
- Body: optional, a few short lines separated by a blank line from the subject. Use bullets for specifics; wrap at ~72 chars.
  - Example:
    - `- Normalize dashes/quotes; collapse line-break hyphenation`
    - `- Add auto phrase-split fallback (order, gap, ratio)`
    - `- Add tests; help epilog`
- Tips:
  - No trailing period in the subject; keep it ≤ 50 chars when possible.
  - Use present tense, active voice ("add", "fix", "update").
  - Reference issues in the body (e.g., `Closes #123`).
  - One logical change per commit; split unrelated changes.

## Security & Configuration Tips
- App config: `~/.config/zstarview/config.json` (auto‑managed). Do not commit local configs.
- Cache: ephemeris and data cached via `appdirs`; safe to delete when troubleshooting.

## Session & Design Log (`dev-notes/`)

- File: `dev-notes/session-YYYY-MM-DD.md`
- Primary Goal: Capture the why behind key decisions made during a session. Design discussions and their outcomes are the main focus.
- Secondary Goal: Supplement decisions with the how by logging relevant command-line transcripts for context, debugging, or reproducibility.

### Logging Format

Entries should use one of the following formats, choosing the one that best fits the context.

#### 1. Design Decision (Default)

- Structure:
  ```markdown
  - Topic: {Brief description of the feature or problem}
    - Decision: A clear and concise summary of the final decision.
    - Rationale: The core reasons why this decision was made. Explain benefits and trade-offs.
    - Alternatives Considered (Optional): Briefly mention other options and why they were not chosen.
  ```

#### 2. Command Transcript (When Necessary)

- Use when specific commands, their sequence, and their output are critical (e.g., troubleshooting, verifying a change).
- Structure:
  ```markdown
  - {Time HH:MM} Ran: `{command}`
    Output:
    ```text
    {A concise excerpt of the output}
    ```
    Result: {Brief one-line explanation of the outcome}
  ```

### General Guidance
- Prioritize logging decisions over raw transcripts; group related command transcripts under a single Topic/Decision when possible.
- Keep outputs concise and redact sensitive info (tokens, secrets, private paths) if they appear.
- At the start of work: create the session file with a short header (date, scope).
- At the end: add a short summary (commits, tags, pushes, follow-ups).
- Optional: If you also capture a full terminal log (`script`, `tee`, etc.), reference the file at the top and paste only key excerpts.

## Decision Safeguards
- For risky, destructive, security-sensitive, cross-cutting, or materially ambiguous changes, state the issue, risk, and proposed scope before proceeding. Ask for confirmation when the action would create external side effects, destroy data, or materially expand the requested scope.
- For ordinary in-scope local changes, use the repository context and reasonable decision criteria to proceed, then validate the result. Do not require a special approval exchange for every uncertainty.
- Trigger examples:
  - Ambiguity in requirements or library choice
  - Irreversible/destructive actions (history rewrites, large deletes)
  - Large refactors without tests or clear rollback
  - Security-sensitive or dependency changes without validation
  - Cross‑cutting API changes impacting many files
  - Conflicts with testing conventions or this AGENTS.md
- Process:
  - Summarize the issue and risk, and propose a safer scope or incremental alternative.
  - Request explicit confirmation when the risk or scope requires it.
  - Prefer incremental changes with minimal tests and/or feature flags.
  - Log the decision point and commands in `dev-notes/session-YYYY-MM-DD.md`.
  - Example warning:

    ```text
    ⚠️ Confirmation needed — Risky or ambiguous change detected
    - Issue: Potentially broad or unclear change could cause breakage.
    - Risk: High scope + unclear intent; potential breakage across components.
    - Proposal: Clarify target behavior, add a minimal test, then proceed incrementally.
    Please confirm the proposed scope or provide adjusted constraints.
    ```
  - Scope: Applies to the entire repository. More-nested AGENTS.md files may refine but not weaken these safeguards.

Notes:
- Historical session logs may exist under `docs/` (e.g., `docs/session-YYYY-MM-DD.md`). New logs should go to `dev-notes/`.
