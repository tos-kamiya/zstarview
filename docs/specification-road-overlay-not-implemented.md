# Road Overlay Proposal (Not Implemented)

Status: not implemented. This file is a preserved investigation note, not a
current user-facing feature specification.

## Proposed behavior

- Treat roads as an optional reference layer, not as a complete map or
  navigation dataset.
- Limit the drawing range to 50m ground distance from the observer.
- Retrieve a slightly larger bbox, such as 100m, and clip segments to the 50m
  drawing range.
- Use Overture Transportation road segments as the candidate source.
- Draw only segments that contain explicit `width_rules` data.
- Do not infer a width from road class, subclass, or a typical road width.
- If there are no usable width rules, draw no road layer.
- Omit the layer entirely during fast-mode or simplified display.

## Investigation result

The probe successfully downloaded Overture segment data around Tokyo Station:

| Bbox radius | Features | Road segments | Segments with `width_rules` | GeoJSON size |
| --- | ---: | ---: | ---: | ---: |
| 50m | 156 | 108 | 0 | 789,439 bytes |
| 100m | 314 | 266 | 0 | 1,005,166 bytes |

The result shows that the data volume is modest, but no road width data was
present in either Tokyo Station sample. This makes an explicit-width-only road
overlay unlikely to be useful without another data source or a broader survey.

## Current decision

Do not implement this road overlay in the current application. Keep the probe
and this proposal for future comparison with Japanese road-edge or 3D city
model sources.
