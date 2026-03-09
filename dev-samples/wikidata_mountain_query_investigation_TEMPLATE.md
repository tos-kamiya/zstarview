# Wikidata mountain query investigation (TEMPLATE)

## Goal

Find representative mountains suitable for summit-view use in zstarview.

Important: the goal is not to build a complete mountain gazetteer.  
The goal is to assemble a curated list of mountain viewpoints that are:

- regionally representative,
- commonly recognized,
- stable enough for name-based startup resolution,
- and usable as summit-view locations.

## Files

- Query memo: `dev-samples/wikidata_mountain_query_investigation_YYYY-MM-DD.md`
- Raw WDQS result(s): `dev-samples/wikidata_mountain_query_raw_result_YYYY-MM-DD.json`
- Curated seed JSON: `dev-samples/wikidata_mountain_candidates_YYYY-MM-DD.json`
- Normalizer script: `dev-samples/build_curated_mountain_viewpoints.py`
- Normalized output: `dev-samples/wikidata_mountain_viewpoints_YYYY-MM-DD.json`

## Discovery strategy

Use WDQS and/or Wikidata item pages only for discovery and enrichment.
Do not rely on automatic ranking by elevation alone.

Selection should be curated by:

- regional representativeness,
- general recognizability,
- stable naming,
- and coordinate quality.

## Notes from investigation

- Record which query patterns worked.
- Record which patterns produced too much noise.
- Record any classes or properties that were useful.
- Record why specific mountains were kept or excluded.

## Practical takeaway

Keep the split:

- Wikidata / WDQS: discovery and enrichment
- curated seed JSON: editorial selection
- local script: normalization into zstarview data shape

This keeps the final dataset reviewable and avoids overfitting the final list to
fragile or overly broad SPARQL.

## Candidate seed contract

The curated seed JSON should contain a top-level `candidates` list.

Each candidate may contain:

- `qid`
- `name`
- `names`
- `labels`
- `latitude_deg`
- `longitude_deg`
- `elevation_m`
- `region_hint`
- `country_codes`
- `wikipedia_titles`

Seed values override Wikidata-derived values when both are present.

## Query notes

Paste WDQS queries or links here.

### Query 1: Broad mountain discovery

Use this to get a broad list of mountains with coordinates and elevation.

```sparql
SELECT ?item ?itemLabel ?coord ?elevation ?country ?countryLabel WHERE {
  ?item wdt:P31 ?instance .
  ?instance wdt:P279* wd:Q8502 .   # mountain

  ?item wdt:P625 ?coord .
  OPTIONAL { ?item wdt:P2044 ?elevation . }
  OPTIONAL { ?item wdt:P17 ?country . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "ja,en". }
}
ORDER BY DESC(?elevation) ?itemLabel
```

### Query 2: Highest point of each country

Use this as a practical starting point for regionally representative candidates.

```sparql
SELECT ?country ?countryLabel ?item ?itemLabel ?coord ?elevation WHERE {
  ?country wdt:P31 wd:Q3624078 .   # sovereign state
  ?country wdt:P610 ?item .        # highest point

  ?item wdt:P31 ?instance .
  ?instance wdt:P279* wd:Q8502 .   # mountain

  ?item wdt:P625 ?coord .
  OPTIONAL { ?item wdt:P2044 ?elevation . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "ja,en". }
}
ORDER BY ?countryLabel
```

### Query 3: Representative peaks by mountain range

Use this to complement the country-highest-point list with famous peaks from
major ranges.

```sparql
SELECT ?range ?rangeLabel ?item ?itemLabel ?coord ?elevation WHERE {
  VALUES ?range {
    wd:Q513   # Himalayas
    wd:Q5463  # Andes
    wd:Q5462  # Alps
  }

  ?item wdt:P4552 ?range .         # mountain range
  ?item wdt:P31 ?instance .
  ?instance wdt:P279* wd:Q8502 .   # mountain
  ?item wdt:P625 ?coord .
  OPTIONAL { ?item wdt:P2044 ?elevation . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "ja,en". }
}
ORDER BY ?rangeLabel DESC(?elevation)
```

These queries are intended for manual use in `https://query.wikidata.org/`.
Treat them as discovery aids; the final candidate list should still be curated
manually before normalization.
