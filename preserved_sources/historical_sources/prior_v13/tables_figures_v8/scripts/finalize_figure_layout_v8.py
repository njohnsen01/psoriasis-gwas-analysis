from update_supplement_figures_v8 import *

def rebuild_selected():
 # Changes below affect layout only; reuse the already verified association display data.
 r=json.loads((O/'MANHATTAN_SOURCE.private.json').read_text());d=pd.read_csv(O/'Manhattan_display.private.tsv.gz',sep='\t');figure1(d,r['counts'])
 reuse_regions();conditional_figures();sensitivity_figures()
 # Complete pooled MHC-excluded ranks were already generated from the full frozen outputs.
 qp=N/'submission_assets_20260911/pooled_autosomes_MHC_excluded_autosomal_QQ_ranks.private.tsv.gz'
 m=json.loads((N/'submission_assets_20260911/MANIFEST.json').read_text())
 assert sha(qp) in str(m),'Pooled MHC-excluded rank checksum is not in the original manifest'
 q=pd.read_csv(qp,sep='\t');c=D['Calibration'];c=c[(c.group=='pooled')&(c.domain=='autosomes_MHC_excluded')];assert len(c)==1;c=c.iloc[0]
 assert int(q.N_valid_tests.iloc[0])==int(c.valid_tests)
 cp=O/'pooled_MHC_excluded_calibration.private.json';savejson(cp,{'lambda_GC_from_full_output_P':c.lambda_GC,'N_valid_tests':int(c.valid_tests),'source':'Calibration and source-manifest-verified full-rank table','interval':'chr6:25000000-34000000 inclusive; GRCh38','weights':'unit'})
 fig,ax=plt.subplots(figsize=(87.5/25.4,3.3));fig.subplots_adjust(left=.19,right=.96,bottom=.18,top=.84);qq(ax,qp,cp,'Pooled autosomes excluding MHC','a')
 export_fig(fig,'Supplementary_Figure_S4f','Pooled MHC-excluded complete-output calibration','The pooled autosomal QQ plot excludes chromosome 6 positions 25,000,000–34,000,000 inclusive on GRCh38. Every valid retained test contributes unit weight. Expected quantiles and observed values use complete-output ranks, and display thinning does not rerank points or alter λGC. Shading denotes a pointwise 95% beta order-statistic envelope.','Manifest-verified pooled MHC-excluded full-output ranks; Calibration')
 savejson(O/'LAYOUT_REVISION_COMPLETE.private.json',{'figures':figures,'sources':sources,'code_sha256':sha(Path(__file__))})

if __name__=='__main__':
 try:rebuild_selected()
 except Exception as e:savejson(O/'LAYOUT_REVISION_FAILED.private.json',{'error':str(e),'traceback':traceback.format_exc()});raise

