"""Trial reference solver for indexed monoclinic and hexagonal spacings.

Uses ordinary reciprocal-metric crystallography. No proprietary methods.
This is an independent reference implementation, not a reproduction of a
historical published mean after changing the accepted subset population.
Run: python lattice_metric_trial.py --self-test
"""
import argparse,csv,itertools,json,math,pathlib,platform
import numpy as np

def determinant_integer(matrix):
    """Exact Bareiss determinant for a small square integer matrix."""
    a=[list(map(int,r)) for r in matrix];n=len(a);sign=1;previous=1
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        value=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*value-a[i][k]*a[k][j]
                assert numerator%previous==0
                a[i][j]=numerator//previous
            a[i][k]=0
        previous=value
    return sign*a[-1][-1]

def validate_rows(rows,wavelength):
    if not math.isfinite(wavelength) or wavelength<=0:raise ValueError('Wavelength must be positive and finite.')
    normalized=[]
    for row in rows:
        if len(row)!=4:raise ValueError('Expected h,k,l,two_theta_deg.')
        h,k,l,t=map(float,row)
        if not all(math.isfinite(v) for v in (h,k,l,t)):raise ValueError('Nonfinite input.')
        if not all(v.is_integer() for v in (h,k,l)):raise ValueError('Miller indices must be integers.')
        if h==k==l==0:raise ValueError('000 is not a diffraction reflection.')
        if not 0<t<180:raise ValueError('Require 0 < two_theta_deg < 180.')
        normalized.append((int(h),int(k),int(l),t))
    return normalized

def solve_subset(rows,system,wavelength=1.5418):
    rows=validate_rows(rows,wavelength)
    count={'monoclinic':4,'hexagonal':2}.get(system)
    if count is None:raise ValueError('Supported systems: monoclinic, hexagonal.')
    if len(rows)!=count:raise ValueError(f'Expected {count} rows.')
    ds=[float(wavelength/(2.0*math.sin(math.radians(t/2.0)))) for h,k,l,t in rows]
    rhs=np.array([1/d**2 for d in ds])
    if system=='monoclinic':
        integer_matrix=[[h*h,k*k,l*l,2*h*l] for h,k,l,t in rows]
        scaled_rhs=rhs
    else:
        # Multiply the standard hexagonal equations by 3 for an integer matrix.
        integer_matrix=[[4*(h*h+h*k+k*k),3*l*l] for h,k,l,t in rows]
        scaled_rhs=3*rhs
    det=determinant_integer(integer_matrix)
    result={'rows_h_k_l_two_theta':rows,'determinant_integer':det,'d_A':ds}
    if det==0:return dict(result,status='singular_index_system')
    matrix=np.array(integer_matrix,dtype=float)
    result['condition_number']=float(np.linalg.cond(matrix))
    try:q=np.linalg.solve(matrix,scaled_rhs)
    except np.linalg.LinAlgError:return dict(result,status='numerical_solve_failure')
    result['coefficients_A_minus2']=q.tolist()
    if not np.all(np.isfinite(q)):return dict(result,status='nonfinite_coefficients')
    if system=='monoclinic':
        A,B,C,D=map(float,q)
        if A<=0 or B<=0 or C<=0 or A*C-D*D<=0:
            return dict(result,status='nonphysical_metric')
        w=-D/math.sqrt(A*C)
        if not -1<w<1:return dict(result,status='nonphysical_angle')
        a=1/math.sqrt(A*(1-w*w));b=1/math.sqrt(B);c=1/math.sqrt(C*(1-w*w))
        values={'a_A':a,'b_A':b,'c_A':c,'cos_beta':w,'beta_deg':math.degrees(math.acos(w))}
        predicted=[(h*h/a**2+k*k*(1-w*w)/b**2+l*l/c**2-2*h*l*w/(a*c))/(1-w*w) for h,k,l,t in rows]
    else:
        if not np.all(q>0):return dict(result,status='nonphysical_metric')
        a,c=1/math.sqrt(float(q[0])),1/math.sqrt(float(q[1]))
        values={'a_A':a,'c_A':c}
        predicted=[(4/3)*(h*h+h*k+k*k)/a**2+l*l/c**2 for h,k,l,t in rows]
    residuals=np.array(predicted)-rhs
    max_residual=float(max(abs(residuals)))
    result.update(values=values,residuals_original_equations_A_minus2=residuals.tolist(),max_residual=max_residual)
    # Explicit numerical validation; no peak selection or outlier filtering.
    if not np.all(abs(residuals)<=1e-12+1e-10*abs(rhs)):
        return dict(result,status='original_equation_residual_failure')
    return dict(result,status='accepted')

def summarize(solutions,system):
    accepted=[r for r in solutions if r['status']=='accepted']
    if not accepted:return {'accepted':0}
    keys=('a_A','b_A','c_A','beta_deg') if system=='monoclinic' else ('a_A','c_A')
    stats={k:{'mean':float(np.mean([r['values'][k] for r in accepted])),
              'population_sd':float(np.std([r['values'][k] for r in accepted],ddof=0))} for k in keys}
    out={'accepted':len(accepted),'statistics':stats,'maximum_original_equation_residual':max(r['max_residual'] for r in accepted)}
    if system=='monoclinic':
        ws=np.array([r['values']['cos_beta'] for r in accepted])
        out['legacy_beta_convention']={'degrees_acos_mean_cos':math.degrees(math.acos(float(ws.mean()))),
            'degrees_asin_population_sd_cos':math.degrees(math.asin(float(ws.std()))),
            'warning':'The second quantity is not the SD of the subset beta angles.'}
    return out

def analyze(rows,system,wavelength=1.5418):
    rows=validate_rows(rows,wavelength)
    n={'monoclinic':4,'hexagonal':2}[system]
    if len(rows)<n:raise ValueError('Insufficient input rows.')
    results=[]
    for number,ix in enumerate(itertools.combinations(range(len(rows)),n),1):
        r=solve_subset([rows[i] for i in ix],system,wavelength)
        r.update(system_number=number,input_rows_1based=[i+1 for i in ix]);results.append(r)
    summary=summarize(results,system)
    summary.update(input_count=len(rows),attempted=len(results),status_counts={s:sum(r['status']==s for r in results) for s in sorted({r['status'] for r in results})})
    return {'system':system,'wavelength_A':wavelength,'summary':summary,'subsets':results}

def load_csv(path):
    with pathlib.Path(path).open(newline='',encoding='utf-8-sig') as f:
        return [[float(r[k]) for k in ('h','k','l','two_theta_deg')] for r in csv.DictReader(f)]

def synthetic_rows(rows,system,a,b,c,beta,wavelength=1.5418):
    out=[];w=math.cos(math.radians(beta))
    for h,k,l,t in rows:
        r=(h*h/a**2+k*k*(1-w*w)/b**2+l*l/c**2-2*h*l*w/(a*c))/(1-w*w) if system=='monoclinic' else (4/3)*(h*h+h*k+k*k)/a**2+l*l/c**2
        out.append((h,k,l,2*math.degrees(math.asin(wavelength*math.sqrt(r)/2))))
    return out

def self_test(sr,zn):
    checks=[]
    def check(name,condition):
        if not condition:raise AssertionError(name)
        checks.append(name)
    # Independent exact determinant oracle for row-pivot and singular cases.
    rng=np.random.default_rng(20261029)
    for size in (2,4):
        for _ in range(50):
            matrix=rng.integers(-3,4,size=(size,size)).tolist()
            direct=0
            for perm in itertools.permutations(range(size)):
                inversions=sum(perm[i]>perm[j] for i in range(size) for j in range(i+1,size))
                direct+=(-1 if inversions%2 else 1)*math.prod(matrix[i][perm[i]] for i in range(size))
            check('exact_determinant_against_independent_permutation_expansion',determinant_integer(matrix)==direct)
    for system,rows,expected in [('monoclinic',sr,30),('hexagonal',zn,42)]:
        result=analyze(rows,system);check(system+'_accepted_count',result['summary']['accepted']==expected)
        valid=[r for r in result['subsets'] if r['status']=='accepted']
        for r in valid:
            base=r['values']
            for perm in itertools.permutations(r['rows_h_k_l_two_theta']):
                alt=solve_subset(perm,system)
                check(system+'_subset_row_permutation',alt['status']=='accepted' and all(abs(alt['values'][k]-v)<1e-9 for k,v in base.items()))
            reversed_indices=[(-h,-k,-l,t) for h,k,l,t in r['rows_h_k_l_two_theta']]
            alt=solve_subset(reversed_indices,system)
            check(system+'_whole_index_reversal',all(abs(alt['values'][k]-v)<1e-9 for k,v in base.items()))
        for beta in ([70,90,111.58,125] if system=='monoclinic' else [90]):
            truth={'a_A':12.328,'c_A':10.882} if system=='monoclinic' else {'a_A':3.2511,'c_A':5.208}
            if system=='monoclinic':truth.update(b_A=7.143,beta_deg=beta)
            gen=synthetic_rows(rows,system,truth['a_A'],truth.get('b_A',1),truth['c_A'],beta)
            rr=analyze(gen,system)
            for r in rr['subsets']:
                if r['status']=='accepted':check(system+'_synthetic_recovery',all(abs(r['values'][k]-v)<1e-9 for k,v in truth.items()))
        try:validate_rows([(0,0,0,30)],1.5418)
        except ValueError:check('reject_000',True)
        else:check('reject_000',False)
    # Deliberately impossible spacing relationships yield a nonphysical metric.
    bad=[(1,0,0,10),(0,1,0,10),(0,0,1,10),(1,0,1,120)]
    check('reject_nonphysical_metric',solve_subset(bad,'monoclinic')['status']=='nonphysical_metric')
    # A false old 3x3 gate must not reject full-rank system 6.
    rr=analyze(sr,'monoclinic');s6=rr['subsets'][5]
    check('recover_full_rank_case_rejected_by_old_3x3_gate',s6['status']=='accepted' and s6['determinant_integer']!=0)
    check('retain_genuine_singular_cases',rr['summary']['status_counts']['singular_index_system']==5)
    return {'status':'PASS','checks':len(checks),'distinct_checks':sorted(set(checks))}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--self-test',action='store_true');parser.add_argument('--output',default='trial_results.json');args=parser.parse_args()
    root=pathlib.Path(__file__).resolve().parent
    sr=load_csv(root/'SrSiO3_email.csv');zn=load_csv(root/'ZnO_email.csv')
    results={'scope':'Reference implementation of standard crystallographic equations; unchanged supplied selections.','runtime':{'python':platform.python_version(),'numpy':np.__version__},'SrSiO3':analyze(sr,'monoclinic'),'ZnO':analyze(zn,'hexagonal')}
    if args.self_test:results['self_test']=self_test(sr,zn)
    path=pathlib.Path(args.output);path.write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps({'SrSiO3':results['SrSiO3']['summary'],'ZnO':results['ZnO']['summary'],'self_test':results.get('self_test'),'output':str(path)},indent=2))
