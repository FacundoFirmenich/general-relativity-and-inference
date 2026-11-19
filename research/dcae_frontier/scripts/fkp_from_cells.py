#!/usr/bin/env python3
"""Load ONLY source-backed DESI/SDSS cell catalog; never invent missing cells."""
import argparse,json,pathlib,sys
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fkp import load_cells,cumulative_fkp
p=argparse.ArgumentParser(description='FKP-derived DCAE for observed binned galaxy survey cells')
p.add_argument('--input',type=pathlib.Path,required=True,help='CSV with measured cells, nbar, power, radial bins')
p.add_argument('--output',type=pathlib.Path,required=True)
p.add_argument('--survey-protocol',required=True,help='Permanent link to audited survey selection and nbar pipeline')
p.add_argument('--power-protocol',required=True,help='Permanent link/version for adopted P(k,z)')
a=p.parse_args()
if not a.input.exists():
    raise SystemExit('NOT_ESTIMABLE: authentic cell input not accessible; no synthetic data substituted')
cells,digest=load_cells(a.input)
result={'status':'ANALYTICAL_FKP_DERIVED_NOT_CAUSAL_VALIDATION',
 'input_sha256':digest,'survey_protocol':a.survey_protocol,'power_protocol':a.power_protocol,
 'units':{'cell_volume':'cMpc^3','nbar':'1/cMpc^3','P_k':'cMpc^3','k':'1/cMpc'},
 'bins':cumulative_fkp(cells),'independent_external_validation':False}
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(a.output)
