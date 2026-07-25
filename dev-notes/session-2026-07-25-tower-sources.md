# Session 2026-07-25

- Scope: Expand tower viewpoints with field-level provenance and overseas observation viewpoints.

- Topic: Multiple data sources for tower viewpoints
  - Decision: Keep Wikidata as the primary source for facility identity, coordinates, names, and structural height. Add an item-level `sources` array so observation-deck names and heights can be attributed to official facility websites.
  - Rationale: A tower's structural height and the observer's deck height are different facts and are often maintained by different sources. Field-level provenance keeps the distinction explicit without replacing the existing Wikidata-based pipeline.

- Topic: Initial overseas expansion
  - Decision: Add Burj Khalifa, Taipei 101, Sky Tower Auckland, Empire State Building, Lotte World Tower, and Petronas Towers with curated observation heights and source URLs.
  - Validation: Tower data tests passed (10 tests); all six names resolved through the runtime viewpoint resolver; `git diff --check` passed.

- Summary: Updated the data-generation schema, bundled tower dataset, resolver metadata keys, developer documentation, and tests. Existing unrelated working-tree changes were preserved.
