# Atlas de investigación

[ Abrir el atlas ](https://gravitas-epistemic-atlas.f-firmenich.chatgpt.site). La publicación alojada requiere la cuenta propietaria. El HTML autónomo del repositorio funciona sin conexión.

Seis preguntas, recorrido con pausa, mapa arrastrable, filtros y evidencia descargable. El resumen puede leerse sin JavaScript.

La tipografía del mapa impreso usa Open Sans, títulos en mayúsculas y separación ampliada entre encabezados y cajas.

DOM tests: `npm install --no-save jsdom` · `node tools/test_epistemic_dom.cjs`.
Browser tests (Playwright + Chromium): `python tools/test_epistemic_graph.py`.

## Workflow

La entrada principal es un lienzo navegable: mové cada nodo, seleccioná para iluminar conexiones y usá doble clic para enfocar. El minimapa permite recorrer el conjunto. Se conservan la vista anterior, el resumen y el tour.

Los gráficos impresos usan Open Sans, incluida con su licencia en `docs/fonts`. Para incluir el TikZ con XeLaTeX, cargá `docs/graph_fonts.tex` después de fontspec, desde la raíz del proyecto.
