from pathlib import Path
import datetime,hashlib,importlib.util,itertools,json,os,shutil,subprocess,sys,traceback
import numpy as np,pandas as pd,openpyxl,zstandard
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';O=I/'conditional_current_v1'
sys.path.insert(0,str(N.parent/'X_source_reanalysis_20260910/python_dependencies'));import pgenlib
spec=importlib.util.spec_from_file_location('primary_auto',N/'scripts/run_autosomes.py');auto=importlib.util.module_from_spec(spec);spec.loader.exec_module(auto)
sha=auto.sha
def save(p,d):
 p.write_text(json.dumps(d,indent=2,default=str))
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def state(stage,**kw):
 save(O/'STATUS.private.json',dict(stage=stage,updated_utc=stamp(),pid=os.getpid(),**kw));print(stage,kw,flush=True)
def run(cmd,p):
 save(p.with_suffix('.command.json'),cmd)
 with p.open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=auto.ENV)
 if r.returncode:raise RuntimeError('Command failed '+str(r.returncode)+' '+str(p))
def tab(w,name):
 it=w[name].iter_rows(values_only=True);h=next(it);return pd.DataFrame(list(it),columns=h)
def group(r):
 c=int(r.CHROM);p=int(r.GENPOS)
 if c==5 and 131000000<=p<=133000000:return 'chr5_131_133Mb'
 if c==5 and 159000000<=p<=160000000:return 'chr5_159Mb'
 if c==18 and 54000000<=p<=55000000:return 'chr18_54Mb'
 if c==19 and 10000000<=p<=11000000:return 'chr19_10Mb'
 return None
def readpv(p):return pd.read_csv(p,sep='\t',dtype={'ID':str,'REF':str,'ALT':str})
def readps(p):
 d=pd.read_csv(p,sep='\t',dtype=str);d=d.rename(columns={'#FID':'FID','#IID':'IID'});assert d.IID.is_unique;return d
def pgvec(prefix,want,subset=None):
 pv=readpv(Path(str(prefix)+'.pvar'));ps=readps(Path(str(prefix)+'.psam'));kw={} if subset is None else {'sample_subset':np.array(subset,dtype=np.uint32)}
 reader=pgenlib.PgenReader(str(prefix).encode()+b'.pgen',**kw);n=len(ps) if subset is None else len(subset);out={}
 try:
  for v in want:
   j=int(np.flatnonzero(pv.ID.eq(v))[0]);x=np.empty(n,np.float64);reader.read_dosages(j,x);assert np.all((x<0)|((x>=0)&(x<=2)));out[v]=x
 finally:reader.close()
 return ps,out
def main():
 assert not O.exists(),'Preserve existing conditional directory; review before a new attempt';O.mkdir()
 state('initializing')
 protocol={'created_utc':stamp(),'status':'post-hoc current-primary conditional sensitivity, not independent replication','selection':'Every current QC-passing reporting-group lead in four explicitly recorded multi-lead autosomal regions; coordinate grouping is for this test plan, not evidence of signal independence. MHC and TNIP1 single reporting leads not part of this within-region comparison. Historical X targets not promoted through this analysis.','genotypes':'Frozen original PGEN chunks; source generation, CRC32C, size and SHA256 checked; PLINK2 make-pgen and pmerge preserve source dosage, with exact participant-aligned dosage comparison. No new QC relaxation.','model':'Original primary REGENIE 4.1.2 Step2 binary model, original keep, phenotype, numeric covariates and completed Step1 LOCO predictions. Firth approximate correction at pThresh .01; Firth does not itself model relatedness. Existing Step1 whole-genome predictions supply the relatedness-aware framework.','baseline_gate':{'N_alleles_positions':'exact source equality','BETA_SE':{'rtol':2e-5,'atol':1e-7},'LOG10P':{'rtol':2e-5,'atol':1e-5},'stop_if_mismatch':True},'families':{'strongest_lead_conditioned':{'attempted_tests':9,'bonferroni_threshold':.05/9},'mutually_conditioned':{'attempted_tests':13,'bonferroni_threshold':.05/13}},'interpretation':'Both families are exploratory post-selection sensitivities. Also display conventional genome-wide threshold. Neither family alone establishes causal variant, fine-mapping, or independent replication. Convergence/collinearity failures remain unavailable, not null results.','LD':'Pair-complete ALT-dosage Pearson r squared, current primary provider-pruned controls in EUR/AFR/AMR; N >=500, both MAC >=20, both call rates >=.98, nonzero variance. Pruning reflects available released kinship pairs, not a guarantee of absence of all relatives. All exact counts remain private pending coordinated disclosure review.','documentation':['https://rgcgithub.github.io/regenie/options/','https://www.cog-genomics.org/plink/2.0/data']}
 save(O/'PROTOCOL.private.json',protocol)
 versions={'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__}
 for tool in ['regenie','plink2']:
  r=subprocess.run([tool,'--version'],capture_output=True,text=True);versions[tool]=(r.stdout+r.stderr).strip()
 assert '4.1.2' in versions['regenie'];save(O/'runtime_versions.json',versions)
 inputs={str(p):sha(p) for p in [Path(__file__),N/'scripts/run_autosomes.py',N/'inputs/keep.tsv',N/'inputs/phenotype.tsv',N/'inputs/covariates.tsv',N/'step1/model_pred.list',N/'step1/model_1.loco',N/'step1/VALIDATED.json',I/'supplement_staging_v6/S1_reconciled_tables_HELD.xlsx']}
 valid=json.loads((N/'step1/VALIDATED.json').read_text());assert valid['status']=='validated';assert inputs[str(N/'step1/model_1.loco')]==valid['loco_sha256']
 assert (N/'step1/model_pred.list').read_text().split()==['psoriasis_case',str(N/'step1/model_1.loco')]
 keep=pd.read_csv(N/'inputs/keep.tsv',sep='\t',header=None,names=['FID','IID'],dtype=str);assert len(keep)==valid['n'] and keep.IID.is_unique
 save(O/'INPUTS.private.json',inputs)
 w=openpyxl.load_workbook(I/'supplement_staging_v6/S1_reconciled_tables_HELD.xlsx',read_only=True,data_only=True);lead=tab(w,'Autosomal_reporting_groups');lead['region']=lead.apply(group,axis=1);lead=lead[lead.region.notna()].copy();assert len(lead)==13 and lead.ID.is_unique
 lead.sort_values(['CHROM','GENPOS'],inplace=True);lead.to_csv(O/'planned_leads.private.tsv',sep='\t',index=False);ids=set(lead.ID);found=set();chunks=[];baseline=[];sources=[]
 for p in sorted((N/'results/autosomes').glob('*/passed_variants.txt')):
  if not p.parent.name.startswith(('chr5_','chr18_','chr19_')):continue
  hit=ids.intersection(p.read_text().splitlines())
  if not hit:continue
  assert not (hit&found);found|=hit;chunks.append((p.parent,sorted(hit)))
  c=json.loads((p.parent/'COMPLETE.json').read_text());assert c['status']=='validated';assoc=p.parent/'assoc_psoriasis_case.regenie';assert sha(assoc)==c['results_sha256']
  q=pd.read_csv(p.parent/'variant_QC.private.tsv.gz',sep='\t').set_index('ID');assert q.loc[sorted(hit),'PASS'].all();q.loc[sorted(hit)].to_csv(O/(p.parent.name+'_target_QC.private.tsv'),sep='\t')
  a=pd.read_csv(assoc,sep=r'\s+');a=a[a.ID.isin(hit)].copy();assert len(a)==len(hit);a['source']=str(assoc);a['source_sha256']=sha(assoc);baseline.append(a)
 assert found==ids and len(chunks)==5
 baseline=pd.concat(baseline,ignore_index=True).set_index('ID');baseline.to_csv(O/'frozen_primary_targets.private.tsv',sep='\t')
 assert shutil.disk_usage(O).free>128*1024**3
 prefixes=[];source_dosages={}
 for source,hit in chunks:
  state('preparing_source',chunk=source.name)
  d=O/'sources'/source.name;d.mkdir(parents=True);rec=json.loads((source/'source_receipts.json').read_text());receipt={}
  for ext in ['pgen','pvar','psam']:
   path=d/('input.'+('pvar.zst' if ext=='pvar' else ext));receipt[ext]=auto.download(rec[ext]['uri'],path)
   for k in ['generation','crc32c','bytes','sha256']:assert str(receipt[ext][k])==str(rec[ext][k]),(source.name,ext,k)
  with (d/'input.pvar.zst').open('rb') as f,(d/'input.pvar').open('wb') as g:zstandard.ZstdDecompressor().copy_stream(f,g)
  save(d/'verified_source_receipts.private.json',receipt);sources.append({'chunk':source.name,'receipts':receipt,'original_receipt_sha256':sha(source/'source_receipts.json')})
  pv=readpv(d/'input.pvar').set_index('ID');assert pv.index.is_unique
  for v in hit:
   assert pv.loc[v,'REF']==str(baseline.loc[v,'ALLELE0']) and pv.loc[v,'ALT']==str(baseline.loc[v,'ALLELE1']) and int(pv.loc[v,'POS'])==int(baseline.loc[v,'GENPOS'])
  (d/'selected.txt').write_text('\n'.join(hit)+'\n');prefix=d/'slice'
  run(['plink2','--pfile',str(d/'input'),'--keep',str(N/'inputs/keep.tsv'),'--extract',str(d/'selected.txt'),'--make-pgen','--sort-vars','--threads','4','--memory','4096','--out',str(prefix)],d/'slice.console.log')
  sp=readps(d/'input.psam');subset=np.flatnonzero(sp.IID.isin(set(keep.IID))).astype(np.uint32);assert len(subset)==len(keep)
  _,raw=pgvec(d/'input',hit,subset);sliceps,sliced=pgvec(prefix,hit);order=pd.Index(sp.iloc[subset].IID).get_indexer(sliceps.IID);assert (order>=0).all();assert set(sliceps.IID)==set(keep.IID)
  for v in hit:assert np.array_equal(raw[v][order],sliced[v]),('source_to_slice_dosage',v);source_dosages[v]=pd.Series(sliced[v],index=sliceps.IID)
  prefixes.append(prefix)
 save(O/'SOURCE_RECEIPTS.private.json',sources);(O/'merge_list.txt').write_text('\n'.join(map(str,prefixes))+'\n')
 state('merging_selected_dosages')
 run(['plink2','--pmerge-list',str(O/'merge_list.txt'),'--make-pgen','--threads','4','--memory','4096','--out',str(O/'current_leads')],O/'merge.console.log')
 ps,G=pgvec(O/'current_leads',ids);assert len(ps)==len(keep) and set(ps.IID)==set(keep.IID);pv=readpv(O/'current_leads.pvar');assert set(pv.ID)==ids and pv.ID.is_unique
 for v in ids:assert np.array_equal(G[v],source_dosages[v].loc[ps.IID].to_numpy()),('slice_to_merge_dosage',v)
 save(O/'GENOTYPE_CONCORDANCE.private.json',dict(status='exact_source_slice_merge_dosage_equality',variants=len(ids),samples=len(ps),outputs={str(p):sha(p) for p in [O/'current_leads.pgen',O/'current_leads.pvar',O/'current_leads.psam']}))
 (O/'all_targets.txt').write_text('\n'.join(lead.ID)+'\n')
 template=json.loads((chunks[0][0]/'assoc.console.command.json').read_text())
 def fit(label,targets,condition):
  f=O/'fits'/label;f.mkdir(parents=True);(f/'targets.txt').write_text('\n'.join(targets)+'\n');cmd=template.copy()
  for flag,value in [('--pgen',str(O/'current_leads')),('--extract',str(f/'targets.txt')),('--out',str(f/'assoc'))]:cmd[cmd.index(flag)+1]=value
  if condition:(f/'condition.txt').write_text('\n'.join(condition)+'\n');cmd+=['--condition-list',str(f/'condition.txt')]
  state('fitting',model=label);run(cmd,f/'console.log');a=pd.read_csv(f/'assoc_psoriasis_case.regenie',sep=r'\s+');expected=set(targets)-set(condition);assert set(a.ID)==expected and a.ID.is_unique
  log=(f/'assoc.log').read_text();assert 'number of individuals used in analysis = '+str(valid['n']) in log;assert str(valid['cases'])+' cases and '+str(valid['controls'])+' controls' in log;assert 'LOCO predictions' in log
  a['analysis']=label;a['conditioning_IDs']=';'.join(condition);a['output']=str(f/'assoc_psoriasis_case.regenie');a['output_sha256']=sha(f/'assoc_psoriasis_case.regenie');a['log_sha256']=sha(f/'assoc.log');return a
 b=fit('unconditioned',lead.ID.tolist(),[]).set_index('ID');audit=[]
 for v in lead.ID:
  for k in ['CHROM','GENPOS','ALLELE0','ALLELE1','N','TEST']:assert str(b.loc[v,k])==str(baseline.loc[v,k]),('source_identity',v,k)
  for k in ['BETA','SE','LOG10P']:
   x=float(b.loc[v,k]);y=float(baseline.loc[v,k]);ok=bool(np.isclose(x,y,rtol=2e-5,atol=1e-5 if k=='LOG10P' else 1e-7));audit.append(dict(ID=v,field=k,primary=y,reproduced=x,absolute_difference=abs(x-y),passed=ok))
 pd.DataFrame(audit).to_csv(O/'unconditioned_concordance.private.tsv',sep='\t',index=False);assert all(x['passed'] for x in audit),'Unconditioned mismatch: conditional fits not launched'
 save(O/'BASELINE_VERIFIED.private.json',dict(status='passed',checks=len(audit),variants=len(ids),protocol_sha256=sha(O/'PROTOCOL.private.json')));state('baseline_verified')
 results=[];failures=[]
 for region,leads in lead.groupby('region',sort=True):
  groupids=leads.sort_values('LOG10P',ascending=False).ID.tolist();strongest=groupids[0]
  plans=[('strongest_lead_conditioned',region+'_strongest',groupids,[strongest])]+[('mutually_conditioned',region+'_mutual_'+str(j+1),groupids,[x for x in groupids if x!=v]) for j,v in enumerate(groupids)]
  for family,label,targets,condition in plans:
   try:
    a=fit(label,targets,condition);a['region']=region;a['family']=family;a['family_attempted']=protocol['families'][family]['attempted_tests'];a['family_threshold']=protocol['families'][family]['bonferroni_threshold'];results.append(a)
   except Exception as e:
    failures.append(dict(region=region,family=family,analysis=label,targets=';'.join(sorted(set(targets)-set(condition))),condition=';'.join(condition),reason=str(e)));save(O/'fit_failures.private.json',failures)
 a=pd.concat(results,ignore_index=True) if results else pd.DataFrame();assert len(a)+sum(len(x['targets'].split(';')) for x in failures)==22
 if len(a):
  for k in ['BETA','SE','LOG10P']:a[k]=pd.to_numeric(a[k],errors='coerce')
  a['estimable']=np.isfinite(a[['BETA','SE','LOG10P']]).all(axis=1)&a.SE.gt(0)&a.LOG10P.ge(0);a['OR']=np.exp(a.BETA);a['CI95_Wald_low']=np.exp(a.BETA-1.959963984540054*a.SE);a['CI95_Wald_high']=np.exp(a.BETA+1.959963984540054*a.SE);a['P_REGENIE']=10.**(-a.LOG10P);a['family_significant']=a.estimable&a.P_REGENIE.lt(a.family_threshold);a['genome_wide_significant']=a.estimable&a.P_REGENIE.lt(5e-8);a['effect_unit']='per ALT-coded dosage unit; autosomal 0/1/2';a.to_csv(O/'conditional_estimates.private.tsv',sep='\t',index=False)
 save(O/'fit_failures.private.json',failures)
 state('computing_control_LD')
 ldpath=Path('/home/jupyter/psoriasis_reproducibility_20260910/reconstructed_ancestry_v2/ancestry_autosomal_inputs.private.pkl');D=pd.read_pickle(ldpath);D.index=D.index.astype(str);assert D.index.is_unique;D=D.loc[D.index.intersection(ps.IID)].copy();assert D.unrelated.notna().all();anc=D.ancestry_pred.astype(str).str.upper();eligible=D.unrelated.astype(bool)&D.psoriasis_case.eq(0)&anc.isin(['EUR','AFR','AMR']);LD=[]
 for region,L in lead.groupby('region'):
  for u,v in itertools.combinations(L.ID,2):
   for ancestry in ['EUR','AFR','AMR']:
    groupids=D.index[eligible&anc.eq(ancestry)];ix=pd.Index(ps.IID).get_indexer(groupids);assert (ix>=0).all();x=G[u][ix];y=G[v][ix];ox=x>=0;oy=y>=0;pair=ox&oy;n=int(pair.sum());cx=float(ox.mean()) if len(x) else np.nan;cy=float(oy.mean()) if len(y) else np.nan;sx=float(x[pair].sum());sy=float(y[pair].sum());mx=min(sx,2*n-sx);my=min(sy,2*n-sy);reasons=[]
    if n<500:reasons.append('pair_N_below_500')
    if mx<20 or my<20:reasons.append('pair_MAC_below_20')
    if not np.isfinite(cx) or not np.isfinite(cy) or min(cx,cy)<.98:reasons.append('callrate_below_0.98')
    if n<2 or np.var(x[pair])<=0 or np.var(y[pair])<=0:reasons.append('zero_or_unavailable_variance')
    r2=float(np.corrcoef(x[pair],y[pair])[0,1]**2) if not reasons else np.nan
    LD.append(dict(region=region,ID1=u,ID2=v,ancestry=ancestry,reference='current primary provider-pruned controls; combined sex',reference_N=len(ix),pair_N=n,MAC1=mx,MAC2=my,callrate1=cx,callrate2=cy,r2=r2,available=not reasons,reason=';'.join(reasons),dosage='ALT dosage read from source-verified PGEN; autosomal 0/1/2',ancestry_source=str(ldpath),ancestry_source_sha256=sha(ldpath)))
 pd.DataFrame(LD).to_csv(O/'ancestry_control_LD.private.tsv',sep='\t',index=False);assert len(LD)==60
 for p,h in inputs.items():assert sha(Path(p))==h,('input_changed_during_run',p)
 complete=dict(status='complete_secure_conditional_sensitivity',updated_utc=stamp(),baseline_variants=len(ids),conditional_attempted=22,conditional_output_rows=len(a),fit_failures=failures,LD_pairs_by_ancestry=len(LD),outputs={p.name:sha(p) for p in O.glob('*.private.tsv')},interpretation='Post-hoc discovery-overlapping conditional sensitivity; no independent replication or automatic independence claim. Private counts require disclosure review.')
 save(O/'COMPLETE.private.json',complete);state('complete',conditional_output_rows=len(a),failed_models=len(failures),LD_rows=len(LD))
if __name__=='__main__':
 try:main()
 except Exception as e:
  if O.exists():save(O/'FAILED.private.json',dict(error=str(e),traceback=traceback.format_exc(),updated_utc=stamp()))
  traceback.print_exc();raise
