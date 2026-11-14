"""Non-sensitive operational status only; no participant records or association findings."""
from pathlib import Path
import datetime,json,re
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013')
print('\033[2J\033[H',end='')
print('PSORIASIS X-CHROMOSOME ANALYSES')
print(datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC'))
print()
print('Full-X scan: %d / 128 sections finished' % len(list((R/'full_x').glob('chunk_*/DONE.json'))))
print('Whole-scan male coding audit: %d / 128 sections finished' % len(list((R/'full_x_ploidy_audit').glob('chunk_*.private.json'))))
print('Rare-variant annotation: %d / 1189 sections finished' % len(list((R/'rare_vat_v2').glob('part_*.json'))))
for title,marker in [('Annotation preparation','rare_vat_v2/ANNOTATION_COMPLETE.json'),('Rare genotype QC/export','rare_genotypes/GENOTYPES_READY.json'),('Rare gene tests','rare_gene_tests/TESTS_COMPLETE.json'),('Rare aggregate review','rare_gene_tests/PRIMARY_REVIEW.json'),('Full-X lead verification','full_x_followup/FOLLOWUP_COMPLETE.json'),('Full-X aggregate review','aggregate_review/FULL_X_REVIEW.json')]:
 print(title+': '+('checkpoint present; inspect before interpretation' if (R/marker).exists() else 'pending'))
for stem in ['run_full_x.production','audit_full_x_ploidy','chain_rare_preparation','chain_rare_association','chain_rare_review','chain_full_x_followup','chain_full_x_review']:
 p=R/(stem+'.pid')
 if p.exists():
  pid=p.read_text().strip();stat=Path('/proc')/pid/'stat'
  state=stat.read_text().split(')')[-1].strip().split()[0] if stat.exists() else 'ended'
  print(stem+': '+('running' if state not in ['ended','Z'] else 'ended'))
p=R/'prepare_rare_genotypes.console.log'
if p.exists():
 text=p.read_text(errors='replace')
 stages=re.findall(r'\[Stage [^\r\n]+',text)
 if stages:print('Genotype processing '+stages[-1])
 print('Genotype log error marker: '+str('Traceback (most recent call last)' in text))
p=R/'full_x_ploidy_audit/GLOBAL_PLOIDY_AUDIT.json'
if p.exists():print('Male coding audit result: '+json.loads(p.read_text())['status'])
print()
print('Results remain provisional until validation and scientific review finish.')
print('Protected participant data remain inside Verily.')
