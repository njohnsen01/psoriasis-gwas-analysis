#!/usr/bin/env python3
"""Reconstruct psoriasis follow-up matrices; never use saved model input pickles.
X: unaffected provider multiallelic PGEN, with whole-genotype male heterozygosity masking.
Autosomes: original raw PLINK exports (historical fixed lead panel).
"""
from pathlib import Path
import argparse,hashlib,importlib.metadata as im,json,sys
import numpy as np,pandas as pd
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--audit-root',type=Path,required=True)
p.add_argument('--cohort-dir',type=Path,required=True)
p.add_argument('--labs-dir',type=Path,required=True)
p.add_argument('--old-raw-root',type=Path,required=True)
p.add_argument('--new-raw-root',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
X=a.audit_root/'X_source_reanalysis_20260910';rr=a.audit_root/'relatedness_sensitivity_20260910'
sys.path.insert(0,str(X/'python_dependencies'));import pgenlib
inputs={};proof=[]
def record(path):
 path=Path(path).resolve()
 if str(path) not in inputs:
  h=hashlib.sha256()
  with path.open('rb') as f:
   for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
  inputs[str(path)]={'bytes':path.stat().st_size,'sha256':h.hexdigest()}
 return path
def read(path,**kw):return pd.read_csv(record(path),sep='\t',**kw)
ph=read(a.cohort_dir/'psoriasis_pheno.tsv',dtype={'IID':str,'FID':str}).set_index('IID')
cv=read(a.cohort_dir/'psoriasis_covariates.tsv',dtype={'IID':str,'FID':str}).set_index('IID')
assert ph.index.is_unique and cv.index.is_unique and set(ph.index)==set(cv.index)
assert np.array_equal(ph.FID,cv.FID.reindex(ph.index))
d=ph.join(cv.drop(columns='FID'),validate='one_to_one')
d['SEX']=d.sex.map({1:1,0:2});assert d.SEX.notna().all()
metrics=read(X/'private_provider_qc/genomics_metrics_Dec142023_1859_02_tz0000.tsv',dtype={'research_id':str}).set_index('research_id')
assert metrics.index.is_unique and d.index.isin(metrics.index).all()
d=d.join(metrics,validate='one_to_one')
coverage=read(X/'private_provider/coverage_failures.tsv',dtype=str)
flags=read(X/'private_provider_qc/flagged_samples.tsv',dtype=str)
anc=read(a.audit_root/'phenotype_validation_20260910/private_full_ancestry.tsv',dtype={'research_id':str}).set_index('research_id')
assert anc.index.is_unique and d.index.isin(anc.index).all()
d=d.join(anc[['ancestry_pred']],validate='one_to_one')
d['coverage_eligible']=~d.index.isin(coverage.research_id)
d['canonical_match']=((d.SEX==1)&(d.dragen_sex_ploidy=='XY'))|((d.SEX==2)&(d.dragen_sex_ploidy=='XX'))
d['provider_qc_flag']=d.index.isin(flags.s)
d['strict_eligible']=d.coverage_eligible & d.canonical_match & ~d.provider_qc_flag
removed=read(rr/'samples_relatedness_flagged_samples.tsv',dtype=str).iloc[:,0]
d['unrelated']=~d.index.isin(removed)
pairs=read(rr/'samples_relatedness.tsv',dtype={'i.s':str,'j.s':str})
kept=d.index[d.strict_eligible & d.unrelated]
assert not (pairs['i.s'].isin(kept)&pairs['j.s'].isin(kept)).any()
M=read(rr/'variant_map.tsv');assert len(M)==25 and M.variant_id.is_unique
ps=read(X/'private_provider/acaf_threshold.chrX.psam',dtype=str)
sids=pd.Index(ps['#IID']);assert sids.is_unique and d.index.isin(sids).all()
ix=sids.get_indexer(d.index)
pvp=record(X/'private_provider/acaf_threshold.chrX.pvar');nskip=0
with pvp.open() as f:
 for line in f:
  if not line.startswith('##'):break
  nskip+=1
pv=pd.read_csv(pvp,sep='\t',skiprows=nskip,dtype=str)
offsets=np.r_[np.uint64(0),(pv.ALT.str.count(',').to_numpy(np.uint64)+2).cumsum()].astype(np.uintp)
pgen=record(X/'private_provider/acaf_threshold.chrX.pgen')
G=np.empty(len(ps),dtype=np.int8);A=np.empty(len(ps)*2,dtype=np.int32)
rsmap={'chrX:5998724:A:G':'rs183018930','chrX:22917619:C:T':'rs5970952','chrX:67736589:T:C':'rs5964607','chrX:67842269:G:T':'rs5919457'}
with pgenlib.PgenReader(str(pgen).encode(),raw_sample_ct=len(ps),allele_idx_offsets=offsets) as reader:
 assert reader.get_variant_ct()==len(pv)
 for v in M[M.variant_id.str.startswith('chrX:')].itertuples():
  chrom,pos,ref,ea=v.variant_id.split(':');assert ea==v.effect_allele
  candidates=pv.index[(pv['#CHROM']==chrom)&(pv.POS==pos)&(pv.REF==ref)]
  found=[(int(i),pv.at[i,'ALT'].split(',').index(ea)+1) for i in candidates if ea in pv.at[i,'ALT'].split(',')]
  assert len(found)==1,(v.variant_id,found)
  j,ai=found[0];reader.read(j,G,allele_idx=ai);reader.read_alleles(j,A)
  alleles=A.reshape(-1,2);expected=(alleles==ai).sum(1).astype(np.int8);expected[(alleles<0).any(1)]=-9
  assert np.array_equal(G,expected),'PGEN target-allele projection disagrees with GT'
  g=G[ix].astype(float);het=(alleles[ix,0]!=alleles[ix,1])&(alleles[ix]>=0).all(1)
  g[(g<0)|(het & d.SEX.eq(1).to_numpy())]=np.nan
  assert set(g[d.SEX.eq(1)&np.isfinite(g)])<=set([0.,2.])
  d[rsmap[v.variant_id]]=g
  proof.append({'variant':v.variant_id,'source':'provider_chrX_PGEN','variant_index':j,'allele_index':ai,'counted_allele':ea,'whole_GT_projection_verified':True,'male_heterozygotes_masked':True})
print('X_RECONSTRUCTED',flush=True)
d=d.loc[d.strict_eligible].copy();d.to_pickle(a.out/'X_inputs.private.pkl')
base=read(a.labs_dir/'psoriasis_biomarker_analysis_base_no_genotypes.tsv',dtype={'IID':str,'FID':str,'person_id':str}).set_index('IID')
assert base.index.is_unique and set(base.index)==set(ph.index)
bio=base.loc[base.index.isin(d.index)].copy()
rawcache={}
for v in M.itertuples():
 if v.variant_id.startswith('chrX:'):
  bio[v.aligned_column]=d[rsmap[v.variant_id]].reindex(bio.index);continue
 root=a.new_raw_root if v.source=='new_targeted_extraction' else a.old_raw_root
 matches=[]
 for path in sorted(root.glob('*.raw')):
  cols=path.open().readline().split()
  for col in cols:
   if col.startswith(v.variant_id+'_'):matches.append((path,col))
 assert len(matches)==1,(v.variant_id,matches)
 path,col=matches[0]
 if path not in rawcache:
  rawcache[path]=read(path,dtype={'IID':str,'FID':str}).set_index('IID')
 raw=rawcache[path];assert raw.index.is_unique and bio.index.isin(raw.index).all()
 ea=col[len(v.variant_id)+1:];ref=v.variant_id.split(':')[2]
 assert ea in [v.effect_allele,ref]
 g=pd.to_numeric(raw[col],errors='coerce').reindex(bio.index)
 assert g.dropna().between(0,2).all()
 bio[v.aligned_column]=g if ea==v.effect_allele else 2-g
 proof.append({'variant':v.variant_id,'source':str(path),'raw_counted_allele':ea,'effect_allele':v.effect_allele,'transformation':'identity' if ea==v.effect_allele else '2-minus-raw'})
bio['unrelated']=d.unrelated.reindex(bio.index);bio.to_pickle(a.out/'biomarker_inputs.private.pkl')
M.to_csv(a.out/'variant_map.tsv',sep='\t',index=False)
(a.out/'input_manifest.private.json').write_text(json.dumps(inputs,indent=2))
(a.out/'allele_provenance.json').write_text(json.dumps(proof,indent=2))
(a.out/'COMPLETE.json').write_text(json.dumps({'status':'reconstructed','phenotype':'original definition with post-hoc coverage/canonical-ploidy/provider-flag exclusions','cohort_n':len(d),'cases':int(d.psoriasis_case.sum()),'pruned_n':int(d.unrelated.sum()),'pruned_cases':int(d.loc[d.unrelated,'psoriasis_case'].sum()),'genotype_source':'provider multiallelic X PGEN and historical raw autosomal PLINK exports; not a fresh source-wide GWAS','pgenlib':im.version('pgenlib'),'numpy':np.__version__,'pandas':pd.__version__,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
print('FOLLOWUP_INPUTS_RECONSTRUCTED',len(d),flush=True)
