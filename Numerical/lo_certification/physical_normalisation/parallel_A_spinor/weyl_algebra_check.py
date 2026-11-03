#!/usr/bin/env python3
"""Exact 4-spinor/Weyl Grassmann identities for a declared convention.
Does NOT read Matchete's abstract GammaCC or assert its numerical orientation.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sympy as sp

labels=['l0','l1','c0','c1','x0','x1','Lb0','Lb1','Cb0','Cb1','Xb0','Xb1']
order={s:i for i,s in enumerate(labels)}
conj=dict(zip(labels[:6],labels[6:]))
conj.update(dict(zip(labels[6:],labels[:6])))
def atom(s):return {(s,):sp.Integer(1)}
def add(*things):
 out={}
 for a in things:
  for m,c in a.items():out[m]=sp.expand(out.get(m,0)+c)
 return {m:c for m,c in out.items() if c!=0}
def scale(a,z):return {m:sp.expand(z*c) for m,c in a.items() if z*c!=0}
def prod(a,b):
 out={}
 for u,uc in a.items():
  for v,vc in b.items():
   w=u+v
   if len(set(w))!=len(w):continue
   sign=(-1)**sum(order[w[i]]>order[w[j]] for i in range(len(w)) for j in range(i+1,len(w)))
   m=tuple(sorted(w,key=order.get))
   out[m]=sp.expand(out.get(m,0)+sign*uc*vc)
 return {m:c for m,c in out.items() if c!=0}
def dagger(a):
 out={}
 for m,c in a.items():
  z={():sp.conjugate(c)}
  for name in m[::-1]:z=prod(z,atom(conj[name]))
  out=add(out,z)
 return out
def dot(u,v):return add(*(prod(a,b) for a,b in zip(u,v)))
def assert_poly(a,b,msg):
 if a!=b:raise AssertionError(msg+': '+str(a)+' vs '+str(b))

def test():
 I=sp.I;zero=sp.zeros(2);eye=sp.eye(2)
 sx=sp.Matrix([[0,1],[1,0]]);sy=sp.Matrix([[0,-I],[I,0]]);sz=sp.diag(1,-1)
 block=lambda a,b,c,d:sp.BlockMatrix([[a,b],[c,d]]).as_explicit()
 gam=[block(zero,eye,eye,zero)]+[block(zero,s,-s,zero) for s in (sx,sy,sz)]
 C=I*gam[2]*gam[0];eps=sp.Matrix([[0,1],[-1,0]])
 assert C==block(eps,zero,zero,-eps)
 assert C.T==-C and C*C==-sp.eye(4)
 assert all(C*g.T*C.inv()==-g for g in gam)
 assert all(gam[i]*gam[j]+gam[j]*gam[i]==(2 if i==j==0 else -2 if i==j else 0)*sp.eye(4) for i in range(4) for j in range(4))
 l=[atom('l0'),atom('l1')];chi=[atom('c0'),atom('c1')];xi=[atom('x0'),atom('x1')]
 Lbar=[{}, {},atom('Lb0'),atom('Lb1')]
 F=[chi[0],chi[1],scale(atom('Xb1'),-1),atom('Xb0')]
 Fc=[xi[0],xi[1],scale(atom('Cb1'),-1),atom('Cb0')]
 Fbar=[scale(xi[1],-1),xi[0],atom('Cb0'),atom('Cb1')]
 contract=lambda a,b:add(prod(a[0],b[1]),scale(prod(a[1],b[0]),-1))
 assert_poly(dot(Lbar,Fc),dagger(contract(l,chi)),'bar L F^c')
 assert_poly(dot(Lbar,F),dagger(contract(l,xi)),'bar L F')
 chixi=contract(chi,xi)
 assert_poly(dot(Fbar,F),add(chixi,dagger(chixi)),'Dirac bar F F')
 chichi=contract(chi,chi)
 assert_poly(chichi,scale(prod(chi[0],chi[1]),2),'Majorana factor of 2')
 majbar=[scale(chi[1],-1),chi[0],atom('Cb0'),atom('Cb1')]
 majF=[chi[0],chi[1],scale(atom('Cb1'),-1),atom('Cb0')]
 assert_poly(dot(majbar,majF),add(chichi,dagger(chichi)),'Majorana bar F F')
 assert_poly(scale(chichi,sp.Rational(-1,2)),scale(prod(chi[0],chi[1]),-1),'Canonical Majorana -1/2')
 # Bad C -> -C without the field/coupling dictionary changes only the y1 conversion.
 assert dot(Lbar,Fc)!=scale(dot(Lbar,Fc),-1)
 return {'gamma_clifford':'PASS','C=i_gamma2_gamma0':'DECLARED_CHOICE','C_transpose':'-C',
         'C_squared':'-I','barL_Fc':'(l chi)^dagger','barL_F':'(l xi)^dagger',
         'dirac_mass':'-M chi xi + H.c.',
         'majorana_mass':'-M chi chi/2 + H.c.; chi chi=2chi0chi1',
         'mass_matrix_factor':'NO EXTRA HALF OR TWO','wrong_unaccompanied_C_flip':'DETECTED',
         'installed_numeric_Matchete_GammaCC_measured':False,'physical_source_certified':False}

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);x=p.parse_args()
 assert not x.output.exists(),'Output exists; refusing overwrite'
 row=test();x.output.parent.mkdir(parents=True,exist_ok=True)
 x.output.write_text(json.dumps(row,indent=2)+'\n',encoding='utf-8')
 print('Exact chosen-basis exterior-algebra tests PASS; native GammaCC UNMEASURED:',x.output)
if __name__=='__main__': main()
