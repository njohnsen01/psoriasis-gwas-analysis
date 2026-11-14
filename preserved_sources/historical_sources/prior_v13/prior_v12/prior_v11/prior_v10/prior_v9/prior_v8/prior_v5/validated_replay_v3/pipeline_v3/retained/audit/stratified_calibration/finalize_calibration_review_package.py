from pathlib import Path
import json,time,hashlib,numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime,timezone
R=Path(__file__).resolve().parent;D=R/'diagnostics';S=json.loads((R/'strata.json').read_text())
while True:
 state=json.loads((R/'pipeline_status.json').read_text())
 if state['status']=='complete':break
 if state['status']=='failed':
  (D/'review_package_failure.txt').write_text(state.get('error','Pipeline failed'));raise SystemExit(1)
 time.sleep(30)
cal=json.loads((D/'genome_calibration.json').read_text());table=pd.DataFrame(cal['groups'])
receipts=[json.loads(p.read_text()) for p in (R/'receipts').glob('*.json')]
assert len(receipts)==1149 and len({r['stem'] for r in receipts})==1149
qc=[]
for name in S:
 rs=[r['strata'][name] for r in receipts]
 n_valid=sum(r['n_valid'] for r in rs);n_rows=sum(r['n_rows'] for r in rs)
 assert n_valid==int(table[(table.stratum==name)&(table.group=='all')].n_tests.iloc[0])
 qc.append({'stratum':name,'cohort_N':S[name]['n'],'cases':S[name]['cases'],'controls':S[name]['controls'],'genotype_chunks':len(rs),'source_variants':sum(r['n_source_variants'] for r in receipts),'association_rows':n_rows,'valid_tests':n_valid,'invalid_tests':sum(r['n_invalid'] for r in rs),'source_variants_without_row':sum(r['n_source_without_row'] for r in rs),'failed_firth_tests':sum(r['failed_firth_tests'] for r in rs),'firth_tests':sum(r['firth_tests'] for r in rs),'ignored_low_MAC':sum(r['ignored_low_MAC'] for r in rs)})
qct=pd.DataFrame(qc);qct.to_csv(D/'Stratum_Test_Accounting.tsv',sep='\t',index=False)
assert all(qct.source_variants==13255796)
assert np.allclose(table.lambda_GC,table.lambda_from_CHISQ,rtol=2e-5,atol=2e-5)
table=table.merge(qct[['stratum','cohort_N','cases','controls']],on='stratum',validate='many_to_one')
table.to_csv(D/'Supplementary_Table_S10_Full_Output_Calibration.tsv',sep='\t',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7})
for stem,names in [('Supplementary_Figure_S4_Sex_QQ',['female','male']),('Supplementary_Figure_S5_Ancestry_QQ',['EUR','AFR','AMR'])]:
 fig,axes=plt.subplots(len(names),2,figsize=(7.08,3.1*len(names)),squeeze=False)
 for i,name in enumerate(names):
  q=pd.read_csv(D/(name+'_genome_exact_QQ.tsv'),sep='\t')
  for j,key in enumerate(['autosomes','X']):
   ax=axes[i,j];series=[(key,'#4C78A8' if key=='autosomes' else '#C94F70','All tests')]
   if key=='autosomes':series.append(('autosomes_no_MHC','#14866D','MHC excluded'))
   for group,color,label in series:
    v=q[q.group==group];t=table[(table.stratum==name)&(table.group==group)].iloc[0]
    assert len(v)>0 and v['rank'].max()==t.n_tests and np.allclose(v.expected,-np.log10((v['rank']-.5)/t.n_tests),atol=2e-14,rtol=0)
    ax.scatter(v.expected,v.observed,s=2.2,c=color,alpha=.65,linewidths=0,rasterized=True,label=label+rf', $\lambda_{{GC}}={t.lambda_GC:.3f}$')
   t=table[(table.stratum==name)&(table.group==key)].iloc[0];v=q[q.group==key];xm=float(v.expected.max())+.2;ym=max(float(v.observed.max())*1.08,xm)
   ax.plot([0,xm],[0,xm],color='#7C8795',lw=.8,ls=(0,(4,3)));ax.set(xlim=(0,xm),ylim=(0,ym),xlabel=r'Expected $-\log_{10}(P)$',ylabel=r'Observed $-\log_{10}(P)$')
   label=name.capitalize() if name in ['female','male'] else name
   ax.set_title(f'{label} — '+('autosomes' if key=='autosomes' else 'chromosome X')+'\n'+f'{int(t.n_tests):,} valid tests',fontsize=8.5,pad=6)
   ax.text(-.17,1.10,chr(97+2*i+j),transform=ax.transAxes,fontweight='bold',fontsize=10)
   ax.legend(loc='upper left',frameon=False,fontsize=6.8,handletextpad=.2);ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',lw=.4,color='#D9DDE3');ax.set_axisbelow(True)
 fig.subplots_adjust(left=.10,right=.98,bottom=.07,top=.94,hspace=.58,wspace=.36)
 fig.savefig(D/(stem+'.png'),dpi=300,facecolor='white');fig.savefig(D/(stem+'.pdf'),facecolor='white');plt.close(fig)
lines=['# Complete sex and ancestry calibration results','', 'Computed from complete, separately fitted REGENIE outputs. QQ ranks and lambda use unit weight per valid test; points are thinned only after full ranks are computed. These diagnostics do not constitute independent replication.','', '| Group | Cohort N | Cases | Controls | Autosomal tests | Autosomal lambda | X tests | X lambda | Non-PAR X lambda | PAR lambda |','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
for name in S:
 t=table[table.stratum==name].set_index('group');c=S[name]
 lines.append(f'| {name} | {c["n"]:,} | {c["cases"]:,} | {c["controls"]:,} | {int(t.loc["autosomes","n_tests"]):,} | {t.loc["autosomes","lambda_GC"]:.3f} | {int(t.loc["X","n_tests"]):,} | {t.loc["X","lambda_GC"]:.3f} | {t.loc["X_nonPAR","lambda_GC"]:.3f} | {t.loc["X_PAR","lambda_GC"]:.3f} |')
lines+=['','Methods draft: We independently refitted REGENIE Step 1 in female, male, EUR, AFR, and AMR groups using the frozen pruned autosomal marker panel, then tested the complete primary Step 2 variant inputs with each group\'s own LOCO predictions. Female and male models adjusted for age, age squared, and principal components 1–10; constant sex and age-by-sex terms were omitted. Ancestry-group models retained the primary covariates. Approximate Firth testing was applied at nominal P < 0.01, with minimum minor-allele count 5. Chromosome X testing used hg38 pseudoautosomal boundaries. Inflation and QQ ranks were derived from all valid association P values with unit weights, separately for autosomes and X, with MHC exclusion and PAR/non-PAR summaries. Variant counts and failed-test accounting are provided in Supplementary Table S10 and Stratum_Test_Accounting.tsv.','', 'Figure captions: Supplementary Figure S4 presents female and male autosomal and X QQ plots; Supplementary Figure S5 presents EUR, AFR, and AMR plots. The tested variant count and full-output lambda are shown in each panel. Green overlays exclude chromosome 6:25–34 Mb from the autosomal analysis. The dashed diagonal denotes the uniform null. Exact full-output ranks precede display thinning; no sampling weights are used.','', 'Delivery remains pending an approved export process, visual inspection of the new stratum figures, and insertion into the manuscript/supplement. Do not resolve the PI comment solely because these files exist.']
(D/'Calibration_Results_and_Methods_for_Review.md').write_text('\n'.join(lines))
manifest={'completed_utc':datetime.now(timezone.utc).isoformat(),'all_1149_chunks_accounted_for':True,'all_valid_test_totals_match_diagnostics':True,'lambda_P_and_CHISQ_numerically_agree':True,'visual_review_required':True,'export_performed':False,'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in D.iterdir() if p.is_file() and (p.name.startswith('Supplementary_') or p.name in ['Stratum_Test_Accounting.tsv','Calibration_Results_and_Methods_for_Review.md'])}}
(D/'review_package_ready.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
