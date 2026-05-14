# Search Objects at Startup

| Option | Description | Default |
| :----- | :---------- | :------ |
| `--search QUERY` | Resolve a named object at startup. Uses the same local-first matching as GUI `Search Objects...`; known artificial satellites use the app-side current position first and fail instead of falling back to JPL if their current position cannot be obtained, while CLI/export-image JPL lookup may include major bodies and small bodies; bare `QUERY` searches by label or ID (e.g. `Ceres`, `2000001`). | |
| `--search label=QUERY` | Search labels only (e.g. `label=Ceres`). | |
| `--search id=QUERY` | Search IDs only (e.g. `id=2000001`). | |
| `--search-keep-marker` | Keep the selected target as marker plus label. | |
| `--list` | `zstarview-export-image` only. List candidates and exit without rendering. | |

If `-A` or `-Z` is also given, that axis stays fixed and the search result fills the other axis. The search result altitude is clamped to `-5°`.
