"""Reproduce the screenshot notebook on the preserved official snapshot, without network or plot display."""
from pathlib import Path
import json,io,contextlib,gzip,warnings,hashlib
import numpy as np,pandas as pd,requests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad
R=Path(__file__).resolve().parents[1]
p=R/'originals/Hunter_1B_1BkRZtxyfU9_jwEuejU5e0bOlKG0_L.ipynb'
nb=json.loads(p.read_text());source=''.join(nb['cells'][0]['source'])
raw=(R/'raw/heasarc_fermilgrb.tdat.gz').read_bytes()
class Snapshot:
 content=raw
 def raise_for_status(self):pass
old_get=requests.get;old_show=plt.show
requests.get=lambda *args,**kwargs:Snapshot();plt.show=lambda:None
out=io.StringIO();scope={}
try:
 with contextlib.redirect_stdout(out),warnings.catch_warnings(record=True) as warn:
  warnings.simplefilter('always');exec(compile(source,str(p),'exec'),scope)
finally:requests.get=old_get;plt.show=old_show;plt.close('all')
(R/'results/hunter_exact_replay_stdout.txt').write_text(out.getvalue())
correct=pd.read_csv(R/'results/NASA_HEASARC_FERMILGRB.csv')
original=scope['df'];oldselected=scope['clean_data']
first=original.iloc[0]
result={'original_id':'1B_1BkRZtxyfU9_jwEuejU5e0bOlKG0_L','snapshot_sha256':hashlib.sha256(raw).hexdigest(),'source_executed_unchanged':True,'network_replaced_with_verified_snapshot':True,'figure_display_suppressed':True,'warning_messages':sorted(set(str(w.message) for w in warn)),'warning_count':len(warn),'parsed_rows_original':len(original),'parsed_rows_correct':len(correct),'original_dataframe_index_first':str(original.index[0]),'original_first_NAME':str(first['NAME']),'correct_first_name':str(correct.iloc[0]['name']),'original_first_REDSHIFT':float(first['REDSHIFT']),'correct_first_redshift':float(correct.iloc[0]['redshift']),'original_first_TL0':float(first['TL0']),'correct_first_tl0':float(correct.iloc[0]['tl0']),'selected_rows_original':len(oldselected),'X_minmax_original':[float(scope['X'].min()),float(scope['X'].max())],'Y_minmax_original':[float(scope['Y'].min()),float(scope['Y'].max())],'original_slope':float(scope['m_fit']),'original_slope_se':float(scope['m_err']),'fit_rows_original':int(np.sum(scope['Y']<200)),'interpretation':'Exact archived parser and fit replay on current verified catalogue; descriptive audit, not physical inference.'}
# No silent parameter fallback. Conditional OLS is displayed only as a parser/cut sensitivity diagnostic.
valid=correct[(correct.redshift>0)&(correct.tl0>0)].copy()
om=float(scope['cosmo'].Om0);ode=float(scope['cosmo'].Ode0)
valid['K']=valid.redshift.apply(lambda z:quad(lambda u:(1+u)/np.sqrt(om*(1+u)**3+ode),0,z)[0])
for label,frame in [('all_positive_onsets',valid),('archived_cut_tl0_lt_200',valid[valid.tl0<200])]:
 design=np.column_stack([frame.K,np.ones(len(frame))]);coef=np.linalg.lstsq(design,frame.tl0,rcond=None)[0]
 result[label]={'n':len(frame),'slope_s':float(coef[0]),'intercept_s':float(coef[1]),'scope':'Unweighted descriptive regression only; not a LIV estimator or bound.'}
aligned=original.loc[correct['name']]
result['column_shift_checks']={
 'original_REDSHIFT_equals_catalog_luminosity_distance':bool(np.allclose(pd.to_numeric(aligned.REDSHIFT,errors='coerce'),correct.luminosity_distance,equal_nan=True)),
 'original_TL0_equals_catalog_TL100':bool(np.allclose(pd.to_numeric(aligned.TL0,errors='coerce'),correct.tl100,equal_nan=True)),
 'valid_redshift_range':[float(valid.redshift.min()),float(valid.redshift.max())],
 'original_redshift_range':[float(oldselected.REDSHIFT.min()),float(oldselected.REDSHIFT.max())]
}
def finite_json(x):
 if isinstance(x,dict):return {k:finite_json(v) for k,v in x.items()}
 if isinstance(x,list):return [finite_json(v) for v in x]
 if isinstance(x,float) and not np.isfinite(x):return None
 return x
result=finite_json(result)
(R/'results/hunter_exact_replay.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(result,indent=2))
