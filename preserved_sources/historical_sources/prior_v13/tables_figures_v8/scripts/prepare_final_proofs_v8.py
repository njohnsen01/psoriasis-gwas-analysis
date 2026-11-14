from update_supplement_figures_v8 import *
from update_main_tables_v8 import render_table
supplementary_calibration(only='female')
savejson(O/'FEMALE_LABEL_FINAL.private.json',{'figures':figures,'sources':sources,'code_sha256':sha(Path(__file__))})
r=json.loads((O/'TABLES_COMPLETE.private.json').read_text())
for t in r['tables']:
 t['data']=pd.DataFrame(t['data']);render_table(t)
proofs=[]
for p in sorted(F.glob('*.png')):
 with Image.open(p) as im0:
  im1=im0.convert('RGB');ratio=min(1,1450/im1.height,2500/im1.width)
  if ratio<1:im1=im1.resize((round(im1.width*ratio),round(im1.height*ratio)),Image.Resampling.LANCZOS)
  dest=P/('final_'+p.name);im1.save(dest)
  proofs.append({'source':str(p),'source_sha256':sha(p),'proof':str(dest),'proof_sha256':sha(dest),'scale':ratio})
savejson(P/'FINAL_PROOFS.private.json',proofs)
print('FINAL_PROOFS',len(proofs),flush=True)
