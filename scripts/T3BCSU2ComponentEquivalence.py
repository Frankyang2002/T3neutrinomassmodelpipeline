#!/usr/bin/env python3
"""Exact component-level residual constraints for Workstream C (B(-1), C(-1)).

Reads the full-residual export ZIP, not the untrusted source Lagrangian as Python.
No physical Wilson coefficients or counterterms are computed. All cancellations
are conditional on registered native tensors satisfying the explicitly audited
conditions, and on Bose symmetry of the source tensor projections.
"""
import argparse
import re,collections, itertools, fractions, pathlib, hashlib, json, sys, zipfile

# Tiny parser of Matchete InputForm (read only) 
TOK = re.compile(r'\s*([A-Za-z$][A-Za-z0-9$]*|\d+|[\[\]\{\}\(\),+\-*/^])')
class ParseError(ValueError): pass
class Parser:
  def __init__(self,source):
    matches = list(TOK.finditer(source)); self.ts=[m.group(1) for m in matches]; self.i=0
    if ''.join(m.group(0) for m in matches).strip() != source.strip():
      at=sum(len(m.group(0)) for m in matches); raise ParseError('unknown text at '+repr(source[at:at+80]))
  def pop(self,p=None):
    if self.i>=len(self.ts):raise ParseError('end of input')
    x=self.ts[self.i]; self.i+=1
    if p is not None and x!=p:raise ParseError('expected '+p+', got '+x)
    return x
  def peek(self):return self.ts[self.i] if self.i<len(self.ts) else None
  def parse(self):
    node=self.sum()
    if self.peek() is not None:raise ParseError('trailing: '+repr(self.ts[self.i:self.i+10]))
    return node
  def sum(self):
    first=self.product(); seq=[first]
    while self.peek() in ('+','-'):
      op=self.pop(); n=self.product(); seq.append(n if op=='+' else ('neg',n))
    return first if len(seq)==1 else ('add',tuple(seq))
  def product(self):
    first=self.power(); seq=[first]
    while self.peek() in ('*','/'):
      op=self.pop(); n=self.power();seq.append(n if op=='*' else ('inv',n))
    return first if len(seq)==1 else ('mul',tuple(seq))
  def power(self):
    node=self.unary()
    if self.peek()=='^':self.pop();node=('pow',node,self.power())
    return node
  def unary(self):
    if self.peek()=='-':self.pop();return ('neg',self.unary())
    if self.peek()=='+':self.pop();return self.unary()
    return self.atom()
  def atom(self):
    t=self.pop()
    if t=='(':
      x=self.sum();self.pop(')');return x
    if t=='{':
      xs=[]
      if self.peek()!='}':
        xs.append(self.sum())
        while self.peek()==',':self.pop();xs.append(self.sum())
      self.pop('}');return ('list',tuple(xs))
    if t.isdigit():return ('num',int(t))
    if not re.fullmatch('[A-Za-z$][A-Za-z0-9$]*',t):raise ParseError('bad atom '+t)
    if self.peek()=='[':
      self.pop();xs=[]
      if self.peek()!=']':
        xs.append(self.sum())
        while self.peek()==',':self.pop();xs.append(self.sum())
      self.pop(']');return ('call',t,tuple(xs))
    return ('sym',t)

def iscall(x,name=None):return isinstance(x,tuple) and x[0]=='call' and (name is None or x[1]==name)
def sym(x):
  if x[0]=='sym':return x[1]
  raise ParseError('not sym '+repr(x))
def call(x,name):
  if not iscall(x,name):raise ParseError('not '+name+' '+repr(x))
  return x[2]
def lst(x):
  if x[0]!='list':raise ParseError('not list '+repr(x))
  return x[1]
def name_bar(x):
  if iscall(x,'Bar'):
    a=call(x,'Bar')
    if len(a)!=1:raise ParseError('Bar')
    return (sym(a[0]),True)
  return sym(x),False

def flatten_sum(x):
  if x[0]=='add':
    return sum((flatten_sum(y) for y in x[1]),[])
  if x[0]=='neg':return [('neg',y) for y in flatten_sum(x[1])]
  return [x]

def get_factors(x):
  sgn=fractions.Fraction(1); factors=[]
  def walk(node,modifier=1):
    nonlocal sgn
    k=node[0]
    if k=='neg':sgn*=-1;walk(node[1],modifier)
    elif k=='mul':
      for q in node[1]:walk(q,modifier)
    elif k=='inv':walk(node[1],-modifier)
    elif k=='num':sgn*=fractions.Fraction(node[1])**modifier
    elif k=='pow':
      b,e=node[1:]
      if e[0]!='num':raise ParseError('unexpected power')
      if b[0]=='num':sgn*=fractions.Fraction(b[1])**(modifier*e[1])
      else:
        for _ in range(e[1]):walk(b,modifier)
    else:
      if modifier==-1:raise ParseError('unexpected inverse '+repr(node))
      factors.append(node)
  walk(x)
  return sgn,factors

def index(x):
  if iscall(x,'Bar'):x=call(x,'Bar')[0]
  k,typ=call(x,'Index');return (sym(k),repr(typ))

def labels(factors):
  cgs=[];fields=[];couplings=[]
  for f in factors:
    barred=False
    if iscall(f,'Bar'):
      f=call(f,'Bar')[0];barred=True
    if iscall(f,'CG'):
      n,slots=call(f,'CG');cgname,barname=name_bar(n)
      if barred: raise ParseError('Bar CG outer?')
      cgs.append((cgname,barname,tuple((index(s),iscall(s,'Bar')) for s in lst(slots))))
    elif iscall(f,'Field'):
      n,kind,indices,extra=call(f,'Field');fields.append((sym(n),barred,index(lst(indices)[0])))
    elif iscall(f,'Coupling'):
      n,indices,order=call(f,'Coupling'); couplings.append((sym(n),barred))
    else:raise ParseError('unexpected factor '+repr(f)[:300])
  return cgs,fields,couplings

def eps(a,b):return (a==1 and b==2)-(a==2 and b==1)
def inv_tensor_value(slots,matching,values):
  score=1
  for i,j in matching:
    s,t=slots[i],slots[j]
    x,y=values[i],values[j]
    score *= eps(x,y) if s==t else (x==y)
  return score
MATCHINGS=(((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))
# For identical bosons, invariant CG tensors must be symmetric in the indicated slots.
SYM={
  'T3Scalar1SelfCGInv1':((0,2),(1,3)),
  'T3MENativeS2SelfCGInv1':((0,2),(1,3)),
  'T3MENativeExtraQ000112I1CG':((2,3),),
  'T3MENativeExtraQ000202I1CG':((0,1),(2,3)),
  'T3MENativeExtraQ001201I1CG':((1,2),),
  'T3MENativeExtraQ020020I1CG':((0,1),(2,3)),
  'T3MENativeExtraQ020200I1CG':((0,1),(2,3)),
}

def project_components(v,swaps):
  if not swaps:return v
  accum={tuple(k):fractions.Fraction(x) for k,x in v.items()}
  for a,b in swaps:
    new={}
    for idx in itertools.product((1,2),repeat=4):
      idxswap=list(idx);idxswap[a],idxswap[b]=idxswap[b],idxswap[a]
      new[idx]=(accum[idx]+accum[tuple(idxswap)])/2
    accum=new
  return accum

def independent_basis(slots,name):
  import sympy as sp
  assert len(slots)==4
  candidates=[]
  for matching in MATCHINGS:
    raw={seq:inv_tensor_value(slots,matching,seq) for seq in itertools.product((1,2),repeat=4)}
    projected=project_components(raw,SYM.get(name,()))
    candidates.append(projected)
  basis=[];vectors=[]
  for c in candidates:
    vec=[c[k] for k in sorted(c)]
    if sp.Matrix(vectors+[vec]).rank()>len(vectors):
      basis.append(c);vectors.append(vec)
  return basis

# Evaluate symbolic polynomial for one candidate tensor per CG label. 
# Monomial = sorted coupling factors, sorted scalar (bar/label/component) factors; coefficient rational.
def polynomial(terms,tensor=None,mode='all'):
  res=collections.defaultdict(fractions.Fraction)
  for term in terms:
    coef,factors=get_factors(term)
    cg,fs,cs=labels(factors)
    if len(cg)>1:raise ValueError('multiple CGs')
    if mode=='base' and cg:continue
    if mode=='cg' and not cg:continue
    dummies=sorted(set(([v[0] for v in cg[0][2]] if cg else [])+[f[2] for f in fs]))
    for ass in itertools.product((1,2),repeat=len(dummies)):
      mapping=dict(zip(dummies,ass))
      cgval=1
      if cg:
        if tensor is None:raise ValueError('missing CG tensor')
        cgindices=tuple(mapping[q[0]] for q in cg[0][2])
        cgval=tensor[cgindices]
      if not cgval:continue
      fkey=tuple(sorted((field,barred,mapping[idx]) for field,barred,idx in fs))
      ckey=tuple(sorted(cs))
      res[(ckey,fkey)]+=coef*cgval
  return {k:v for k,v in res.items() if v}

class Evidence:
  def __init__(self, archive):
    self.archive=zipfile.ZipFile(archive)
    self.names=set(self.archive.namelist())
  def text(self, name):
    if name not in self.names:raise ValueError('Missing archive member: '+name)
    return self.archive.read(name).decode('utf-8')
  def sha(self,name):
    if name not in self.names:raise ValueError('Missing archive member: '+name)
    return hashlib.sha256(self.archive.read(name)).hexdigest()

def audit_model(evidence,model):
  pre=f'{model}_full_eft1_residual'
  manifest=json.loads(evidence.text(pre+'_full_residual_manifest.json'))
  summary=[]; groups_by_cg={}
  for group in manifest['coupling_groups']:
    text=evidence.text(group['file'])
    terms=flatten_sum(Parser(text).parse())
    if len(terms)!=group['additive_terms']:raise ValueError('Term mismatch '+str((len(terms),group)))
    cgs=[]
    for t in terms:
      cg,_,_=labels(get_factors(t)[1]);cgs+=cg
    labels_set=set(cg[0] for cg in cgs)
    if len(labels_set)!=1:raise ValueError('multiple CG labels '+str(labels_set))
    label=next(iter(labels_set))
    representative=next(x for x in cgs if not x[1])
    slots=tuple(x[1] for x in representative[2])
    basis=independent_basis(slots,label)
    # Construct a component tensor satisfying the independently solved
    # conditions.  This checks every ORIGINAL term, including the CG-free
    # delta contractions, against all tensor-bearing terms exactly.
    if label=='T3MENativeCrossCGInv1':coeffs=(1,-1)
    elif label=='T3MENativeCrossCGInv2':coeffs=(1,1)
    else:coeffs=(1,)
    if len(coeffs)!=len(basis):raise ValueError('Unexpected source CG rank')
    exact_tensor={key:sum(c*array[key] for c,array in zip(coeffs,basis))
      for key in basis[0]}
    if polynomial(terms,exact_tensor):
      raise ValueError('Exact conditional component cancellation failed: '+group['structural_sector'])
    base=polynomial(terms,mode='base')
    projection=[polynomial(terms,tensor,mode='cg') for tensor in basis]
    keys=set(base).union(*(set(q) for q in projection))
    rows=[([q.get(k,0) for q in projection],-base.get(k,0),k) for k in keys]
    summary.append(dict(sector=group['structural_sector'],raw_terms=len(terms),cg=label,rank=len(basis),n_equations=len(rows),n_base_monomials=len(base),n_cg_monomials_per_basis=[len(p) for p in projection]))
    if label not in groups_by_cg:groups_by_cg[label]=dict(rows=[],rank=len(basis),basis=basis)
    if groups_by_cg[label]['rank']!=len(basis):raise ValueError('basis rank changed')
    groups_by_cg[label]['rows']+=rows
  import sympy as sp
  results={}
  for label,data in groups_by_cg.items():
    rows=data['rows'];rank=data['rank'];A=sp.Matrix([[sp.Rational(x) for x in row[0]] for row in rows]);b=sp.Matrix([sp.Rational(row[1]) for row in rows]);
    rA=A.rank();rAb=A.row_join(b).rank()
    if rAb>rA:status='INCONSISTENT'; sol=None
    else:
      sol_set=sp.linsolve((A,b));sol=next(iter(sol_set));status='UNIQUE' if rA==rank else 'NONUNIQUE'
    failrows=[]
    if status=='INCONSISTENT':
      # Find minimal incremental set of inconsistent row constraints to illustrate contradiction
      sel=[]
      for ix,row in enumerate(rows):
        trial=sel+[ix];ta=sp.Matrix([[sp.Rational(z) for z in rows[i][0]] for i in trial]);tb=sp.Matrix([sp.Rational(rows[i][1]) for i in trial]);
        if ta.row_join(tb).rank()>ta.rank():failrows=trial;break
        if ta.rank()>len(sel):sel.append(ix)
      if failrows:
        failrows=[dict(coefficients=[str(k) for k in rows[i][0]],rhs=str(rows[i][1]),monomial=str(rows[i][2])[:700]) for i in failrows]
    results[label]=dict(status=status,n_equations=len(rows),rank=rank,matrix_rank=rA,augmented_rank=rAb,required_coefficients=[str(s) for s in sol] if sol is not None else None,contradiction=failrows)
  if sum(q['raw_terms'] for q in summary)!=137 or len(summary)!=25:
    raise ValueError('Expected exactly 137 additive terms across 25 sectors')
  return summary,results

def self_test():
  for s in ['-(1*2/4)+X[]^2','Foo[Bar[Index[d$$1,T3BSMd2]],{2}]']:
    try:Parser(s).parse()
    except Exception as e: raise AssertionError(('parser',s,str(e))) from e
  slots=(False,True,True,False)
  b=independent_basis(slots,'T3MENativeCrossCGInv1')
  assert len(b)==2
  assert eps(1,2)==1 and eps(2,1)==-1
  assert len(independent_basis((False,False,False,False),'T3MENativeExtraQ020020I1CG'))==1
  print('PASS: expression parser, SU(2) invariant basis rank and bosonic symmetry projection')

def check_native_files(path,model):
  if not path:return None
  p=pathlib.Path(path)
  matches=sorted(p.glob(f'{model}_*native_cross_cg_check.json'))
  if not matches:return {'available':False,'reason':'Native CG check JSON not found'}
  if len(matches)!=1:raise ValueError('Ambiguous native CG checks '+str(matches))
  record=json.loads(matches[0].read_text(encoding='utf-8'))
  if record.get('format')!='t3-bc-native-cross-cg-exact-component-check-v1':
    raise ValueError('Unexpected native CG check format '+str(matches[0]))
  expected={'B_m1':'B_-1','C_m1':'C_-1'}[model]
  if record.get('model')!=expected:raise ValueError('Incorrect model in native CG check')
  arrpath=matches[0].with_name(matches[0].name.replace('_native_cross_cg_check.json','_native_cross_cg_arrays.inputform.txt'))
  if not arrpath.is_file() or arrpath.stat().st_size<10:
    raise ValueError('Missing exact native CG tensor arrays: '+str(arrpath))
  # Independently crosscheck the solved basis coefficients rather than
  # accepting the Wolfram script's boolean alone.
  a1=record.get('CrossInv1_coeff_A_inputform')
  b1=record.get('CrossInv1_coeff_B_inputform')
  a2=record.get('CrossInv2_coeff_A_inputform')
  b2=record.get('CrossInv2_coeff_B_inputform')
  coefficients_consistent=(a1=='1' and b1=='-1' and isinstance(a2,str)
     and a2==b2 and a2 not in ('0','0.','0.0'))
  native_truth=(record.get('native_cross_conditions_met') is True and
     coefficients_consistent and record.get('CrossInv1_in_singlet_span') is True and
     record.get('CrossInv2_in_singlet_span') is True and
     record.get('CrossInv1_nonzero') is True and record.get('CrossInv2_nonzero') is True)
  return {'available':True,
          'file':str(matches[0]),
          'native_arrays_file':str(arrpath),
          'native_arrays_sha256':hashlib.sha256(arrpath.read_bytes()).hexdigest(),
          'independently_checked_coefficients_consistent':coefficients_consistent,
          'CrossInv1_exact_required_A_minus_B':record.get('CrossInv1_exact_required_A_minus_B') is True,
          'CrossInv2_required_equal_A_B_coefficients':record.get('CrossInv2_required_equal_A_B_coefficients') is True,
          'native_cross_conditions_met':native_truth,
          'native_coefficients':{k:v for k,v in record.items() if k.endswith('_inputform')}}

def run(args):
  ev=Evidence(args.archive)
  audit_summary=json.loads(ev.text('summary.json'))
  if audit_summary.get('format')!='t3-bc-complete-eft1-residual-audit-v1':
    raise ValueError('Invalid attached full residual summary format')
  if len(audit_summary.get('models',[]))!=2:raise ValueError('Expected B/C residual summary')
  out={'format':'t3-bc-su2-residual-component-constraints-v1',
       'source_zip':str(pathlib.Path(args.archive).resolve()),
       'status':'ConditionalComponentCancellationNativeTensorEvidencePending',
       'assumptions':['SU(2) doublet fundamental/anti-fundamental epsilon and delta conventions',
          'Four-boson tensors symmetrized on identical source-field legs',
          'Conjugate CG component tensors use the conjugate source values',
          'All four component-field variables and all spurion/coupling variables are symbolic and independent'],
       'physical_matching_certified':False,
       'renormalized_finite_C5_computed':False,
       'models':[]}
  for model in ('B_m1','C_m1'):
    observed=next((j for j in audit_summary['models'] if j['model']=={'B_m1':'B_-1','C_m1':'C_-1'}[model]),None)
    if observed is None:raise ValueError('Missing benchmark model in summary')
    fpre=model+'_full_eft1_residual_'
    for key,filename in [('saved_eomeft','saved_eomeft_full.inputform.txt'),
       ('reconstructed_candidate','reconstructed_candidate_full.inputform.txt'),
       ('symbolic_residual','full_symbolic_residual.inputform.txt')]:
      got=ev.sha(fpre+filename);want=observed['full_expression_sha256'][key]
      if got!=want:raise ValueError(f'{model} hash mismatch for {key}: {got} != {want}')
    sectors,solved=audit_model(ev,model)
    family_keys=set(solved)
    if family_keys != set(SYM)|{'T3MENativeCrossCGInv1','T3MENativeCrossCGInv2'}:
      raise ValueError('Unexpected native tensor family inventory: '+str(family_keys))
    # Exact source-independent solutions for the mixed-tensor families.
    inv1=solved['T3MENativeCrossCGInv1']
    inv2=solved['T3MENativeCrossCGInv2']
    if inv1['status']!='UNIQUE' or inv1['required_coefficients']!=['1','-1']:
      raise ValueError('CrossInv1 basis requirement changed: '+str(inv1))
    if inv2['status']!='NONUNIQUE' or inv2['matrix_rank']!=1 or not any('tau0' in str(x) for x in inv2['required_coefficients']):
      raise ValueError('CrossInv2 necessary subspace changed: '+str(inv2))
    # Seven families have zero component polynomial for EVERY invariant allowed
    # by the source Bose symmetries; their normalization is irrelevant.
    other={k:v for k,v in solved.items() if k not in ("T3MENativeCrossCGInv1","T3MENativeCrossCGInv2")}
    if any(v['matrix_rank']!=0 or v['augmented_rank']!=0 for v in other.values()):
      raise ValueError('Non-cross residual not identically zero')
    native=check_native_files(args.native_checks_directory,model)
    out['models'].append({'model':observed['model'],
      'residual_sha256':observed['full_expression_sha256']['symbolic_residual'],
      'raw_terms':sum(s['raw_terms'] for s in sectors),
      'sector_count':len(sectors),'tensor_family_count':len(solved),
      'constructive_all_25_sectors_zero_under_required_tensors':True,
      'noncross_bose_symmetric_families_cancelled':len(other),
      'exact_mixed_tensor_condition_1':'CrossInv1 = A - B',
      'exact_mixed_tensor_condition_2':'CrossInv2 = gamma*(A + B), gamma nonzero',
      'basis_A':'delta(i,j)*delta(k,l)',
      'basis_B':'delta(i,k)*delta(j,l)',
      'linear_constraints':solved,'sector_inventory':sectors,
      'native_tensor_test':native,
      'component_equivalence_proven':(native is not None and native.get('native_cross_conditions_met') is True),
      'matchete_lagrangian_check_original_passed':False,
      'physical_matching_certified':False})
  if out['models'][0]['residual_sha256']!=out['models'][1]['residual_sha256']:
    raise ValueError('The two benchmark residual expressions differ')
  if all(m['component_equivalence_proven'] for m in out['models']):
    out['status']='ExactComponentEquivalenceUnderNativeCGAndBoseSymmetries'
  elif args.native_checks_directory:
    out['status']='NativeCGEvidenceIncompleteOrConditionsFailed'
  pathlib.Path(args.output).parent.mkdir(parents=True,exist_ok=True)
  pathlib.Path(args.output).write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
  for item in out['models']:
    print(item['model']+': 137 terms/25 sectors; 7 Bose tensor families cancel; CrossInv1=A-B, CrossInv2=gamma(A+B); native='+('MET' if item['component_equivalence_proven'] else 'PENDING/FAILED'))
  print('PASS: hashes and conditional SU(2) component identities; report='+str(args.output))

if __name__=='__main__':
  parser=argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--archive',help='T3BC_FullEFT1Residual_Outputs.zip')
  parser.add_argument('--output',help='Path to JSON summary')
  parser.add_argument('--native-checks-directory',help='Optional new native CG check outputs')
  parser.add_argument('--self-test',action='store_true')
  args=parser.parse_args()
  try:
    if args.self_test:self_test()
    else:
      if not args.archive or not args.output:parser.error('--archive and --output are required')
      run(args)
  except (ValueError,ParseError,KeyError,AssertionError) as e:
    print('FAIL (not certified): '+str(e),file=sys.stderr);sys.exit(1)

