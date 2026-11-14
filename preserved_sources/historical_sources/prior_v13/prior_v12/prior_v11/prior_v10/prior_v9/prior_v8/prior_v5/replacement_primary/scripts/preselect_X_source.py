from pathlib import Path
import csv,gzip,json,re
N=Path(__file__).resolve().parents[1];X=N.parent/'X_source_reanalysis_20260910';O=N/'source_X_candidates';O.mkdir(exist_ok=True)
paths=[X/'private_provider/acaf_threshold.chrX.pvar']+[N/'source_X_counts'/x for x in ['cases.acount','controls.acount','cases.vmiss','controls.vmiss']]
def rows(path):
 f=path.open();return csv.DictReader((s for s in f if not s.startswith('##')),delimiter='\t')
rds=[rows(p) for p in paths];n=0;sites=0;kept=0;unsupported=0
with gzip.open(O/'count_candidates.tsv.gz','wt') as out:
 w=csv.writer(out,delimiter='\t');w.writerow(['source_index','allele_index','CHROM','POS','REF','ALT','case_AC','case_AN','control_AC','control_AN','case_callrate_preliminary','control_callrate_preliminary'])
 for idx,parts in enumerate(zip(*rds,strict=True)):
  pv,ca,co,cm,om=parts;n+=1
  assert ca['REF']==co['REF']==pv['REF'] and ca['ALT']==co['ALT']==pv['ALT']
  assert ca['ID']==co['ID']==cm['ID']==om['ID']==pv['ID']
  alts=pv['ALT'].split(',');ac=list(map(float,ca['ALT_CTS'].split(',')));bc=list(map(float,co['ALT_CTS'].split(',')));assert len(alts)==len(ac)==len(bc)
  an=float(ca['OBS_CT']);bn=float(co['OBS_CT']);cr=1-float(cm['MISSING_CT'])/float(cm['OBS_CT']);br=1-float(om['MISSING_CT'])/float(om['OBS_CT'])
  if cr<.98 or br<.98:continue
  any_kept=False
  for ai,(alt,a,b) in enumerate(zip(alts,ac,bc),1):
   if not re.fullmatch('[ACGT]+',pv['REF']) or not re.fullmatch('[ACGT]+',alt):unsupported+=1;continue
   if min(a,an-a)<20 or min(b,bn-b)<20 or min(a+b,an+bn-a-b)<100:continue
   w.writerow([idx,ai,pv['#CHROM'],pv['POS'],pv['REF'],alt,a,an,b,bn,cr,br]);kept+=1;any_kept=True
  sites+=any_kept
  if n%100000==0:print('SOURCE_X_COUNT_SCREEN',n,kept,flush=True)
assert n==2153825
(O/'COUNT_SCREEN_COMPLETE.json').write_text(json.dumps({'source_sites':n,'candidate_alleles':kept,'candidate_sites':sites,'unsupported_alleles_at_called_sites':unsupported,'status':'preliminary count-only candidates; not final QC; exact masked-genotype counts, biological MAF, site filters and female-control HWE still required'},indent=2))
print('SOURCE_X_COUNT_SCREEN_COMPLETE',n,kept,flush=True)
