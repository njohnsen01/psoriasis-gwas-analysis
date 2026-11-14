import re
from collections import defaultdict
EXCLUDE=4|256|512|1024|2048
TARGETS=[(5995921,'G','A'),(5996494,'GCT','G'),(5998724,'A','G')]
def cigar_map(start,cigar,seq,qual):
 bases,insertions,deleted={},{},set();r,q=start,0
 for ns,op in re.findall(r'(\d+)([MIDNSHP=X])',cigar):
  n=int(ns)
  if op in 'M=X':
   for k in range(n):bases[r+k]=(seq[q+k].upper(),ord(qual[q+k])-33 if qual!='*' else -1,q+k)
   r+=n;q+=n
  elif op=='I':
   insertions[r-1]=[(seq[q+k].upper(),ord(qual[q+k])-33 if qual!='*' else -1,q+k) for k in range(n)];q+=n
  elif op=='D':deleted.update(range(r,r+n));r+=n
  elif op=='N':r+=n
  elif op=='S':q+=n
 return bases,insertions,deleted
def span_haplotype(bases,insertions,deleted,lo,hi):
 if lo not in bases or hi not in bases:return None
 vals=[]
 for pos in range(lo,hi+1):
  if pos in bases:vals.append(bases[pos])
  elif pos not in deleted:return None
  vals.extend(insertions.get(pos,[]))
 if not vals:return None
 return ''.join(x[0] for x in vals),min(x[1] for x in vals),vals[0][2]
def observe(fields,target,reference_window=None):
 flag,mq,start=int(fields[1]),int(fields[4]),int(fields[3])
 if flag&EXCLUDE or fields[2]!='chrX' or fields[9]=='*':return None
 pos,ref,alt=target;bases,ins,deleted=cigar_map(start,fields[5],fields[9],fields[10])
 if len(ref)==len(alt)==1:
  if pos not in bases:return None
  allele,bq,cycle=bases[pos];call='REF' if allele==ref else 'ALT' if allele==alt else 'OTHER'
 else:
  lo,hi,ref_hap,alt_hap=reference_window;v=span_haplotype(bases,ins,deleted,lo,hi)
  if v is None:return None
  allele,bq,cycle=v;call='REF' if allele==ref_hap else 'ALT' if allele==alt_hap else 'OTHER'
 tags={t.split(':',2)[0]:t.split(':',2)[2] for t in fields[11:] if len(t.split(':',2))==3}
 return dict(template=tags.get('RG','')+'|'+fields[0],call=call,MQ=mq,BQ=bq,reverse=bool(flag&16),cycle=(len(fields[9])-cycle if flag&16 else cycle+1),read_length=len(fields[9]),seq=fields[9],soft_clipped=sum(int(n) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',fields[5]) if op=='S'),mate_chrY=fields[6]=='chrY',alternative_chrY=('chrY,' in tags.get('SA','') or 'chrY,' in tags.get('XA','')))
def unique_high_quality(rows):
 by_template=defaultdict(list)
 for r in rows:
  if 30<=r['MQ']<255 and r['BQ']>=20:by_template[r['template']].append(r)
 chosen,conflicts=[],0
 for vals in by_template.values():
  top=max(v['BQ'] for v in vals);best=[v for v in vals if v['BQ']==top]
  if len({v['call'] for v in best})!=1:conflicts+=1;continue
  chosen.append(sorted(best,key=lambda v:(-v['MQ'],v['cycle']))[0])
 return chosen,conflicts
def self_test():
 def sam(cigar,seq,start=100,flag=0,mq=60):return ['read',str(flag),'chrX',str(start),str(mq),cigar,'=',str(start),'0',seq,'I'*len(seq)]
 assert observe(sam('5M','AAGAA'),(102,'G','A'))['call']=='REF'
 assert observe(sam('5M','AAAAA'),(102,'G','A'))['call']=='ALT'
 assert observe(sam('5M','AAAAA',flag=1024),(102,'G','A')) is None
 ref='AAAAAGCTCTAAAAA';alt=ref[:6]+ref[8:];window=(100,114,ref,alt)
 for cigar in ['6M2D7M','7M2D6M','8M2D5M']:assert observe(sam(cigar,alt),(105,'GCT','G'),window)['call']=='ALT'
 assert observe(sam('15M',ref),(105,'GCT','G'),window)['call']=='REF'
 assert span_haplotype(*cigar_map(100,'6M2N7M',alt,'I'*13),100,114) is None
 a=observe(sam('5M','AAAAA'),(102,'G','A'));b=observe(sam('5M','AAGAA'),(102,'G','A'))
 assert unique_high_quality([a,b])==([],1)
 assert len(unique_high_quality([a,a])[0])==1
 assert not unique_high_quality([observe(sam('5M','AAAAA',mq=255),(102,'G','A'))])[0]
 return True
if __name__=='__main__':
 assert self_test();print('CIGAR_INDEL_OVERLAP_AND_QUALITY_TESTS_PASSED')

