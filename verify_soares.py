"""Independent monoclinic audit. Reads local data; no network or external services."""
import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import platform
import numpy as np

ROOT = Path(__file__).resolve().parent
PERMS = [(p, (-1) ** sum(p[i] > p[j] for i in range(4) for j in range(i+1, 4)))
         for p in itertools.permutations(range(4))]

def read_rows(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        sample = f.read(2048); f.seek(0)
        rows = list(csv.DictReader(f, delimiter=';' if ';' in sample.splitlines()[0] else ','))
    if 'H' in rows[0]:
        return np.array([[float(r[k]) for k in ('H','K','L','2THETA')] for r in rows])
    return np.array([[float(r[k]) for k in ('h','k','l','two_theta_deg')] for r in rows])

def metric(rows, wavelength):
    if not np.all(np.isfinite(rows)) or not np.all(rows[:,:3] == np.round(rows[:,:3])):
        raise ValueError('Require finite input and integer indices')
    if np.any(np.all(rows[:,:3] == 0, axis=1)) or np.any((rows[:,3] <= 0) | (rows[:,3] >= 180)):
        raise ValueError('Invalid reflection or angle')
    h,k,l = rows[:,:3].astype(np.int64).T
    x = np.column_stack((h*h, k*k, l*l, 2*h*l))
    y = (2*np.sin(np.deg2rad(rows[:,3]/2))/wavelength)**2
    return x,y

def exact_det(m):
    # Input indices here are small: int64 products are bounded safely.
    if np.max(np.abs(m), initial=0) > 10000:
        raise ValueError('Index matrix exceeds exact int64 determinant bound')
    d = np.zeros(len(m), dtype=np.int64)
    for p,s in PERMS:
        d += s * m[:,0,p[0]] * m[:,1,p[1]] * m[:,2,p[2]] * m[:,3,p[3]]
    return d

def audit(rows, wavelength=1.5418, batch_size=50000):
    x,y = metric(rows, wavelength)
    combos = itertools.combinations(range(len(rows)),4)
    counts = dict(attempted=0, singular=0, nonphysical=0, accepted=0)
    center = None
    sums = np.zeros(6); squares = np.zeros(6)
    max_residual = 0.0
    while True:
        ix = np.array(list(itertools.islice(combos,batch_size)), dtype=np.int64)
        if not len(ix): break
        counts['attempted'] += len(ix)
        m = x[ix]; rhs = y[ix]
        good = exact_det(m) != 0
        counts['singular'] += int(np.sum(~good))
        m = m[good].astype(float); rhs = rhs[good]
        if not len(m): continue
        q = np.linalg.solve(m, rhs[...,None])[...,0]
        A,B,C,D = q.T
        good = np.all(np.isfinite(q),axis=1) & (A>0) & (B>0) & (C>0) & (A*C-D*D>0)
        counts['nonphysical'] += int(np.sum(~good))
        A,B,C,D = q[good].T
        if not len(A): continue
        w = -D/np.sqrt(A*C)
        if np.any(np.abs(w)>=1): raise AssertionError('Angle outside physical domain')
        a = 1/np.sqrt(A*(1-w*w)); b = 1/np.sqrt(B); c = 1/np.sqrt(C*(1-w*w))
        beta = np.rad2deg(np.arccos(w)); vol = a*b*c*np.sqrt(1-w*w)
        vals = np.column_stack((a,b,c,beta,vol,w))
        if center is None: center = vals[0].copy()
        shifted = vals-center
        sums += shifted.sum(axis=0); squares += (shifted*shifted).sum(axis=0)
        counts['accepted'] += len(vals)
        residual = np.abs(np.einsum('nij,nj->ni',m[good],q[good])-rhs[good])
        max_residual = max(max_residual,float(residual.max(initial=0)))
    n = counts['accepted']
    if not n: raise ValueError('No admissible solutions')
    mean = center+sums/n
    sd = np.sqrt(np.maximum(0,squares/n-(sums/n)**2))
    return {'counts':counts,'statistics':{k:{'mean':float(v),'population_sd':float(s)}
            for k,v,s in zip(('a_A','b_A','c_A','beta_deg','volume_A3','cos_beta'),mean,sd)},
            'beta_conventions':{'acos_mean_cos_deg':float(np.rad2deg(np.arccos(mean[5]))),
               'degrees_sd_cos':float(np.rad2deg(sd[5])),
               'asin_sd_cos_deg':float(np.rad2deg(np.arcsin(sd[5])))},
            'max_linear_equation_residual_A_minus2':max_residual}

def verify(sr, full):
    assert sr['counts'] == dict(attempted=35,singular=5,nonphysical=0,accepted=30)
    assert abs(sr['statistics']['beta_deg']['population_sd']-0.038406)<0.0000005
    assert abs(sr['beta_conventions']['asin_sd_cos_deg']-0.035719)<0.0000005
    if full is not None:
        assert full['counts'] == dict(attempted=9381251,singular=486135,nonphysical=180,accepted=8894936)
        for k,expected in [('a_A',12.328400),('b_A',7.143199),('c_A',10.882509),('volume_A3',891.137456)]:
            assert abs(full['statistics'][k]['mean']-expected)<0.0000005, k
        for k,expected in [('a_A',0.139674),('b_A',0.038774),('c_A',0.161765),('volume_A3',14.748165)]:
            assert abs(full['statistics'][k]['population_sd']-expected)<0.0000005, k
        assert abs(full['statistics']['beta_deg']['mean']-111.582198)<0.0000005
        assert abs(full['statistics']['beta_deg']['population_sd']-0.341156)<0.0000005
        assert abs(full['beta_conventions']['asin_sd_cos_deg']-0.307476)<0.0000005
    return 'PASS'

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--full',action='store_true',help='Enumerate all 9,381,251 four-row combinations')
    p.add_argument('--output',type=Path,help='Write JSON to this path; otherwise print JSON only')
    args = p.parse_args()
    small = read_rows(ROOT/'data'/'SrSiO3_email.csv')
    sr = audit(small)
    full = audit(read_rows(ROOT/'data'/'SrSiO3_powderpattern_table_collCode#130895.csv')) if args.full else None
    out = {'runtime':{'python':platform.python_version(),'numpy':np.__version__},
           'wavelength_A':1.5418,'sd_ddof':0,'seven_reflections':sr,'full_124_reflections':full,
           'verification':verify(sr,full)}
    text = json.dumps(out,indent=2,allow_nan=False)
    if args.output: args.output.write_text(text+'\n',encoding='utf-8')
    print(text)

if __name__ == '__main__': main()
