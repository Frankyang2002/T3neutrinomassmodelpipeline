"""Read-only source-evidence collection for A(-4)/D(-2) fermion matching.

Intended path: scripts/CollectT3FermionMatchingEvidence.py
Does not evaluate or claim renormalised Wilson matching. No source changes.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

MODELS = {'A_-4': 'A_m4', 'D_-2': 'D_m2'}
SOURCES = {
 'A_-4': 'output/general_observable/fermion_first_all16_v1/A_m4/data/eft1_after_F_wilson_seed.json',
 'D_-2': 'output/general_observable/fermion_first_self_fix_v2/D_m2/data/eft1_after_F_wilson_seed.json',
}
C5 = 'output/general_matchete_matching/{tag}/c5_coefficient.txt'
SEARCH_DIRS = (
 'output/parallel_t3/full_general_batch2',
 'output/parallel_t3/thesis_full_general_joint_scalar_batch9_v1',
 'output/parallel_t3/thesis_hard_general_hhss_batch10_v1',
 'output/parallel_t3/full_general_batch1',
 'output/parallel_t3',
 'output/general_matchete_matching',
)
# Large generated arrays and unrelated campaign data are not necessary for source inspection.
SUFFIXES = {'.json','.txt','.wl','.wls','.m','.log','.md','.csv'}
MAX_FILE = 12 * 1024 * 1024
MAX_TOTAL = 65 * 1024 * 1024
MAX_FILES = 250


def sha256(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda: f.read(1024*1024),b''): h.update(block)
 return h.hexdigest()


def source_ids(path):
 if not path.is_file():return None
 try:
  data=json.loads(path.read_text(encoding='utf-8-sig'))
  # Seed file's operator list is under source-specific field; raw manifest is preserved in archive.
  counts={k:len(v) for k,v in data.items() if isinstance(v,list) and ('operator' in k.lower() or 'term' in k.lower())}
  return {'top_level_fields':list(data)[:40], 'list_counts':counts}
 except (ValueError,UnicodeError):return {'status':'not_json'}


def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--project',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True)
 args=ap.parse_args()
 root=args.project.resolve(); dest=args.output.resolve()
 if not root.is_dir():ap.error('project directory missing')
 if dest.exists():ap.error('output already exists; refusing overwrite')
 if root in dest.parents:ap.error('output must be outside project to prevent self-collection')
 chosen={}; missing=[]; omitted=[]
 def add(rel,required=False):
  p=root/rel
  if not p.is_file():
   if required:missing.append(rel)
   return
  if p.suffix.lower() not in SUFFIXES: return
  size=p.stat().st_size
  if size>MAX_FILE:
   omitted.append({'path':rel,'reason':'file size cap','bytes':size});return
  chosen[rel]=p
 for model,tag in MODELS.items():
  add(SOURCES[model],True)
  add(C5.format(tag=tag),True)
  for rel in [f'output/general_matchete_matching/{tag}/general_matching.log',
              f'output/parallel_t3/full_general_batch2/thesis_full_general_batch2_laurent_v2/{model}/manifest.json']:
   add(rel)
 # Include only model-linked source ledgers, manifests, audit reports and compact provenance.
 for base in SEARCH_DIRS:
  folder=root/base
  if not folder.is_dir():continue
  for p in folder.rglob('*'):
   if not p.is_file():continue
   rel=p.relative_to(root).as_posix(); low=rel.lower()
   if rel in chosen or p.suffix.lower() not in SUFFIXES:continue
   if not any(t in low for t in ('a_-4','d_-2','a_m4','d_m2')):continue
   if not any(t in low for t in ('ledger','manifest','laurent','hard','matching','source','operator','audit','coefficient','matchete','raw','input')):continue
   add(rel)
 ordered=sorted(chosen)
 total=0; included=[]
 for rel in ordered:
  size=chosen[rel].stat().st_size
  if len(included)>=MAX_FILES or total+size>MAX_TOTAL:
   omitted.append({'path':rel,'reason':'total/count cap','bytes':size});continue
  included.append(rel);total+=size
 report={
  'status':'read_only_source_collection_not_renormalised_matching',
  'project':str(root),'included_bytes':total,
  'files':[{'path':rel,'bytes':chosen[rel].stat().st_size,'sha256':sha256(chosen[rel])} for rel in included],
  'missing_required':missing,'omitted':omitted,
  'tree_seed_summaries':{model:source_ids(root/rel) for model,rel in SOURCES.items()},
  'known_unresolved':['GammaCC-to-Weyl sign and symmetry factor','UV-hard minus EFT1-soft finite counterterm subtraction','complete operator basis and scalar normalisation'],
 }
 dest.parent.mkdir(parents=True,exist_ok=True)
 with ZipFile(dest,'x',compression=ZIP_DEFLATED,compresslevel=6) as z:
  for rel in included:z.write(chosen[rel],arcname=rel)
  z.writestr('EVIDENCE_MANIFEST.json',json.dumps(report,indent=2,allow_nan=False)+'\n')
 print('Created:',dest)
 print('Source files:',len(included),'uncompressed bytes:',total)
 print('Missing required:',missing)
 print('Omitted due to caps:',len(omitted))
 if missing or omitted:print('WARNING: evidence set incomplete; inspect EVIDENCE_MANIFEST.json')
 return 0

if __name__=='__main__':raise SystemExit(main())
