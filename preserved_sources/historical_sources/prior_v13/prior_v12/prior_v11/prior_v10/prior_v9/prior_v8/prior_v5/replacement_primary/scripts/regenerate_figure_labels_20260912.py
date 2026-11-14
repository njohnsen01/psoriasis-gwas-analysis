from pathlib import Path
import sys,importlib.util,json,hashlib,traceback
N=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(N/'scripts'))
spec=importlib.util.spec_from_file_location('label_assets',N/'scripts/build_submission_assets_labels_20260912.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
def main():
 for group in a.GROUPS:
  d=a.load(group,'autosomes');no=d.loc[~(d.CHROM.eq(6)&d.GENPOS.between(25000000,34000000))]
  fig=a.plt.figure(figsize=(7.2,7.1));gs=fig.add_gridspec(2,2,height_ratios=[1,1.15])
  a.manhattan(fig.add_subplot(gs[0,:]),d,'A  '+(group if group in ['EUR','AFR','AMR'] else group.capitalize())+' autosomal association')
  a.qq(fig.add_subplot(gs[1,0]),d,'B  All autosomal tests');a.qq(fig.add_subplot(gs[1,1]),no,'C  MHC excluded')
  fig.subplots_adjust(left=.1,right=.98,bottom=.09,top=.94,hspace=.42,wspace=.35);a.export_figure(fig,'Autosomal_'+group)
  print('LABEL_FIGURE_COMPLETE','Autosomal_'+group,flush=True)
  del d,no
  a.chromosome_x(group)
 files=[{'name':p.name,'bytes':p.stat().st_size,'sha256':a.sha(p)} for p in sorted(a.O.iterdir()) if p.is_file() and p.name not in ['MANIFEST.json','COMPLETE.json','FAILED.json']]
 a.save(a.O/'MANIFEST.json',{'files':files,'source_snapshot_sha256':a.sha(N/'freeze_review_20260912/SNAPSHOT_INDEX.private.json'),'scope':'Figure labels only: Tests replaces ambiguous N, ancestry codes capitalized. Same source association outputs, ranks, curves and thresholds. Prior regional plot unchanged.','release':'HELD'})
 a.save(a.O/'COMPLETE.json',{'status':'complete','figure_sets':12,'visual_QA':'pending','script_sha256':a.sha(Path(__file__)),'module_sha256':a.sha(N/'scripts/build_submission_assets_labels_20260912.py')})
 print('FIGURE_LABEL_UPDATE_COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:a.save(a.O/'FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});raise
