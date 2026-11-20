# Research atlas

[Open the atlas](https://gravitas-epistemic-atlas.f-firmenich.chatgpt.site). The hosted publication requires the owner account. The repository standalone HTML works offline.

Six questions, a pausable tour, draggable map, filters and downloadable evidence. The summary remains readable without JavaScript.

Print-map typography uses sentence-case headings and increased spacing between headings and cards.

DOM tests: `npm install --no-save jsdom` · `node tools/test_epistemic_dom.cjs`.
Browser tests (Playwright + Chromium): `python tools/test_epistemic_graph.py`.

## Workflow

The default view is a navigable canvas: move individual nodes, select to illuminate connections and double-click to focus. Use the minimap to navigate. The previous atlas, summary and tour remain available.

Print diagrams use Open Sans, bundled with its license in `docs/fonts`. To include TikZ with XeLaTeX, load `docs/graph_fonts.tex` after fontspec from the project root.
