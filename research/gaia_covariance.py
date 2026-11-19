"""Recompute Gaia FPR transit covariance contrasts from ESA CCD records."""
from pathlib import Path
import csv,json,hashlib,datetime
import numpy as np
ROOT=Path(__file__).resolve().parent
def main():
    output=[];summaries=[];checks=[]
    for number,name in [(21,'Lutetia'),(1917,'Cuyo'),(1566,'Icarus'),(161989,'Cacus')]:
        file=ROOT/'raw'/f'gaia_{number}.csv'
        rows=list(csv.DictReader(file.open()))
        transits=sorted(set(r['transit_id'] for r in rows))
        rejected=0;single=0;local=[]
        for tid in transits:
            group=[r for r in rows if r['transit_id']==tid and r['is_rejected'].strip().lower() in ['false','0','f']]
            if not group:rejected+=1;continue
            if len(group)==1:single+=1;continue
            angle=np.radians([float(r['position_angle_scan']) for r in group])
            H=np.c_[np.sin(angle),np.cos(angle)]
            def matrices(kind):
                a=np.array([float(r['ra_error_'+kind]) for r in group])
                b=np.array([float(r['dec_error_'+kind]) for r in group])
                rho=np.array([float(r['ra_dec_correlation_'+kind]) for r in group])
                M=np.empty((len(group),2,2));M[:,0,0]=a*a;M[:,1,1]=b*b;M[:,0,1]=M[:,1,0]=rho*a*b
                return M
            R=matrices('random');S=matrices('systematic')
            d=np.einsum('ni,nij,nj->n',H,R,H)
            s2=np.einsum('ni,nij,nj->n',H,S,H)
            assert np.all(d>0) and np.all(s2>=0)
            s=np.sqrt(s2);one=np.ones(len(group))
            shared=np.diag(d)+np.outer(s,s)
            independent=np.diag(d+s2)
            two_dimensional=np.diag(d)+H@S[0]@H.T
            uncertainty=lambda C:float(1/np.sqrt(one@np.linalg.solve(C,one)))
            sh,ind,two=map(uncertainty,[shared,independent,two_dimensional])
            woodbury=np.sum(1/d)-np.sum(s/d)**2/(1+np.sum(s*s/d))
            checks.append(abs(woodbury/(one@np.linalg.solve(shared,one))-1))
            record=dict(number_mp=number,name=name,transit_id=tid,n_ccd=len(group),
                        sigma_shared_mas=sh,sigma_independent_mas=ind,
                        ratio=sh/ind,two_dimensional_relative_change=two/sh-1)
            output.append(record);local.append(record)
        ratios=np.array([r['ratio'] for r in local])
        summaries.append(dict(number_mp=number,name=name,raw_ccd=len(rows),total_transits=len(transits),
                              rejected_transits=rejected,single_accepted_ccd_transits=single,compared_transits=len(local),
                              ratio_min=float(ratios.min()),ratio_median=float(np.median(ratios)),ratio_max=float(ratios.max()),
                              two_dimensional_max_relative_change=max(abs(r['two_dimensional_relative_change']) for r in local),
                              input_sha256=hashlib.sha256(file.read_bytes()).hexdigest()))
    assert max(checks)<1e-10
    result=dict(executed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                scope='Uncertainty of a common along-scan transit displacement, not a Yukawa-amplitude fit',
                selection='ESA is_rejected=false, at least two accepted CCDs per transit',
                max_dense_woodbury_relative_difference=max(checks),objects=summaries)
    (ROOT/'results/gaia_covariance.json').write_text(json.dumps(result,indent=2)+'\n')
    with (ROOT/'results/gaia_covariance_transits.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
    radar=[]
    for number in [1566,161989]:
        f=ROOT/'raw'/f'radar_{number}.json';j=json.loads(f.read_text())
        rr=[dict(zip(j['fields'],r)) for r in j['data']]
        selected=[r for r in rr if r['epoch'].startswith('2015-' if number==1566 else '2022-')]
        radar.append(dict(number_mp=number,signature=j['signature'],all_delay_records=len(rr),
                          selected_records=len(selected),dates=sorted(set(r['epoch'][:10] for r in selected)),
                          stations=sorted(set(r['rcvr'] for r in selected)),records=selected,
                          input_sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    (ROOT/'results/radar_inventory.json').write_text(json.dumps(radar,indent=2)+'\n')
    print(json.dumps(result,indent=2));print(json.dumps(radar,indent=2))
if __name__=='__main__':main()
