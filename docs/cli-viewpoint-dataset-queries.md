# Viewpoint dataset queries

You can inspect the bundled tower/viewpoint and mountain/viewpoint datasets without launching the GUI.

| Option | Description | Default |
| :----- | :---------- | :------ |
| `-h`, `--help` | Show this help message and exit. | |
| `--list-viewpoints {t,m}` | Print bundled tower (`t`) or mountain (`m`) primary names and exit. Output lines use `t/NAME` or `m/NAME`; ASCII fallback names are preferred when available. | |
| `--list-viewpoint-names {t,m}` | Print bundled tower (`t`) or mountain (`m`) names including localized and ASCII-fallback names, then exit. Output lines use `t/NAME` or `m/NAME`. | |
| `--show-viewpoint-json NAME` | Resolve a bundled viewpoint and print its JSON metadata, including `ascii_name` when available, then exit. Prefix `NAME` with `t/` or `m/` to force tower-only or mountain-only resolution. | |

```bash
zstarview --list-viewpoints t
zstarview --list-viewpoint-names t
zstarview --show-viewpoint-json "t/Tokyo Skytree"
zstarview --list-viewpoints m
zstarview --show-viewpoint-json "m/Mount Fuji"
```

These options are mutually exclusive, do not accept the `location` argument, and cannot be combined with time/rendering options.
`--list-viewpoints` prefers ASCII fallback names when available.
`--list-viewpoint-names` includes both the original spellings and ASCII fallback spellings.
`--show-viewpoint-json` reports an ambiguity error with prefixed candidates if an unprefixed name matches both a tower and a mountain exactly.
