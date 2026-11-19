"""Operator, force, and independent numerical checks."""
from pathlib import Path
import json,hashlib
import numpy as np
import h5py
from scipy.linalg import cholesky,solve_triangular,svd
from .physics import validate as physics_validate
R=Path(__file__).resolve().parents[1]

def main():
    out=R/'results';out.mkdir(exist_ok=True)
    f=next((R/'inputs').glob('likelihood_*.h5'))
    assert hashlib.sha256(f.read_bytes()).hexdigest()=='db72d5fb076ad732c48a29032ee86c85f68faef14c4b0ee055f1ffad44ecf50d'
    with h5py.File(f) as h:C=h['covariance/value'][...];W=h['window/value'][...]
    L=cholesky(C,lower=True);A=solve_triangular(L,W,lower=True)
    U,s,Vh=svd(A,full_matrices=True)
    threshold=max(A.shape)*np.finfo(float).eps*s[0];rank=int(np.sum(s>threshold))
    witness=Vh[rank];null_res=np.linalg.norm(A@witness)/(np.linalg.norm(A)*np.linalg.norm(witness))
    assert null_res<1e-13
    # Visible and null-space witnesses in the same input Euclidean norm.
    response=[float(np.linalg.norm(A@Vh[0])**2),float(np.linalg.norm(A@witness)**2)]
    desi=dict(source='DESI DR1 BGS public likelihood, z=0.1--0.4, theta=0.05',
        raw_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),shape=list(W.shape),
        numerical_rank=rank,nullity=W.shape[1]-rank,svd_threshold=float(threshold),
        whitened_singular_max=float(s[0]),whitened_singular_min=float(s[-1]),
        relative_null_residual=float(null_res),equal_input_norm_squared_responses=response,
        scope='Linear observation operator; singular values depend on stated input units and discretization. No cosmological fit or f_R0 bound is derived.')
    (out/'desi_operator.json').write_text(json.dumps(desi,indent=2)+'\n')
    np.savez_compressed(out/'desi_operator_arrays.npz',singular=s,visible=Vh[0],null=witness)
    phy=physics_validate();(out/'physics_validation.json').write_text(json.dumps(phy,indent=2)+'\n')
    print(json.dumps(dict(physics=phy,desi=desi),indent=2))

if __name__=='__main__':main()
