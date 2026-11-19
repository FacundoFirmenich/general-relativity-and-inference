from pathlib import Path
import json
from audit import figures as audit_figures
from cassini import figures as cassini_figures
R=Path(__file__).resolve().parents[2]
if __name__=='__main__':
    audit_figures(json.loads((R/'research/revision_v17/results/audit.json').read_text()))
    cassini_figures(json.loads((R/'research/revision_v17/results/cassini.json').read_text()))
    print('Four bilingual revision figures regenerated from saved records.')
