# Figure font

This directory used to bundle a commercial Helvetica font family, which
cannot legally be redistributed in a public repository. It has been
replaced with **Liberation Sans**, a free, metric-compatible sans-serif
(SIL Open Font License 1.1, see `OFL.txt`), so all plotting notebooks and
scripts that load a font from here produce figures with the same layout
and spacing as before.

Code that references the font family by name (`rcParams['font.family']`
in the notebooks, `base_family` in `workflow/scripts/00_R_presets.R`) has
been updated from `"Helvetica"` to `"Liberation Sans"` accordingly.

If you need pixel-identical Helvetica rendering (e.g. to exactly match
the typeset journal PDF), install a licensed copy of Helvetica locally
and point the relevant `font_path`/`font_import()` calls at it instead —
just don't commit the font files themselves.
