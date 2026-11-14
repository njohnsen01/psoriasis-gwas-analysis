from pathlib import Path
import json,hashlib,shutil
from PIL import Image
O=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912/tables_figures_v8')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
names=['FIGURE1_COMPLETE','FORESTS_COMPLETE','SUPPLEMENT_FIGURES_COMPLETE','LAYOUT_REVISION_COMPLETE','FEMALE_LABEL_FINAL','FIGURE3_FINAL','CONDITIONAL_LAYOUT_FINAL']
figs={}
for n in names:
 for r in json.loads((O/(n+'.private.json')).read_text())['figures']:figs[r['id']]=r
changed=[];history=O/'build_history/rgba_tiff_originals';history.mkdir(parents=True,exist_ok=True)
for r in figs.values():
 p=Path(r['files']['tiff']['file'])
 with Image.open(p) as im:
  if im.mode=='RGB':continue
  assert im.mode=='RGBA' and im.getchannel('A').getextrema()==(255,255),(p,im.mode)
  oldhash=sha(p);dpi=im.info['dpi'];pixels=im.convert('RGB');expected=hashlib.sha256(pixels.tobytes()).hexdigest()
 shutil.copy2(p,history/p.name);pixels.save(p,compression='tiff_lzw',dpi=dpi)
 with Image.open(p) as check:assert check.mode=='RGB' and hashlib.sha256(check.tobytes()).hexdigest()==expected
 r['files']['tiff']['sha256']=sha(p);r['tiff_normalization']={'original_sha256':oldhash,'visible_RGB_pixels':'identical; fully opaque alpha removed','original_retained':str(history/p.name)};changed.append(r)
(O/'TIFF_NORMALIZATION_FINAL.private.json').write_text(json.dumps({'figures':changed,'code_sha256':sha(Path(__file__))},indent=2))
p=O.parent/'scripts/audit_tables_figures_v8.py';s=p.read_text();s=s.replace("'FIGURE3_FINAL','CONDITIONAL_LAYOUT_FINAL']","'FIGURE3_FINAL','CONDITIONAL_LAYOUT_FINAL','TIFF_NORMALIZATION_FINAL']");p.write_text(s)
print('TIFF_RGB_NORMALIZED',len(changed))

