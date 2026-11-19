from pathlib import Path
import json
from streamlit.testing.v1 import AppTest
R=Path(__file__).resolve().parents[1]
rows=[]
for lang in ['Español','English']:
 app=AppTest.from_file(str(R/'ui/app.py'),default_timeout=50).run()
 app.sidebar.selectbox[0].set_value(lang).run()
 options=list(app.sidebar.radio[0].options)
 for name in options:
  app.sidebar.radio[0].set_value(name).run()
  errors=[str(e.message) for e in app.exception]
  rows.append({'language':lang,'page':name,'errors':errors})
  assert not errors,(lang,name,errors)
 app.sidebar.radio[0].set_value(options[9]).run();assert not app.exception
 choices=[x for x in app.selectbox if x.label in ['Experimento','Experiment']]
 choices[0].set_value('Emulador condicional' if lang=='Español' else 'Conditional emulator').run();assert not app.exception
 rows.append({'language':lang,'action':'Topology conditional emulator controls','errors':[]})
 selectors=[x for x in app.selectbox if x.label in ['ID de pasada','Pass ID']]
 expected=sorted({r['pass_id'] for r in json.loads((R/'research/fractalgo/results/cassini.json').read_text())['rows']})
 assert list(map(int,selectors[0].options))==expected
 selectors[0].set_value(20).run();assert not app.exception
 assert len(app.dataframe[0].value)==2
 rows.append({'language':lang,'action':'Real pass IDs and final pass 20','errors':[]})
 # Exercise the window switch, including the retained adverse long-window result.
 app.sidebar.radio[0].set_value(options[8]).run()
 selectors=[x for x in app.selectbox if 'trigger' in x.label]
 selectors[0].set_value(1000).run();assert not app.exception
 rows.append({'language':lang,'action':'Fermi 1000 s selector','errors':[]})
(R/'records/ui_tests_v1_11.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False)+'\n')
print('PASS',len(rows),'UI checks')
