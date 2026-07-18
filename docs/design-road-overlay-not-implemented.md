# Road Overlay Design Proposal (Not Implemented)

Status: not implemented. The design below is retained for comparison and must
not be interpreted as an active runtime architecture.

## Candidate pipeline

1. Download Overture Transportation `segment` features for a small bbox around
   the observer.
2. Select `subtype=road` segments and retain only explicit `width_rules`.
3. Expand each width-rule interval to a strip on both sides of the centerline.
4. Clip the strip to 50m ground distance from the observer.
5. Project the strip into the existing observer-relative coordinate system.
6. Omit acquisition, processing, and drawing during fast-mode and simplified
   display.

The retrieval bbox could be about 100m while the drawing radius remains 50m,
so segments crossing the boundary are not lost. A future cache key would need
to include the observer, radius, Overture release, and width-processing schema
version.

## Feasibility observation

The Tokyo Station probe retrieved 108 road segments in a 50m bbox and 266 road
segments in a 100m bbox. Neither sample contained a segment with
`width_rules`. The GeoJSON outputs were approximately 789KB and 1.0MB,
respectively, so data volume was not the limiting factor; width-data coverage
was.

## Decision

Do not add a runtime road controller, road cache, or renderer based on this
proposal yet. First compare a source that provides road edges or another
reliable width representation, and only then decide whether this design should
be revived.
