xqcmasks={};xqcsummary=[];xqchwe={}
for an in ['EUR','AFR','AMR']:
 hx=pd.read_csv(xpanelqc/(an+'.female_control_hwe.hardy.x'),sep=r'\s+').set_index('ID');hp=pd.read_csv(xpanelqc/(an+'.female_control_hwe.hardy'),sep=r'\s+').set_index('ID')
 h=pd.concat([hx.FEMALE_ONLY_P,hp.P]);assert len(h)==376171 and h.index.is_unique
 xqchwe[an]=h
xqchwe['pooled']=pd.concat(xqchwe,axis=1).min(axis=1,skipna=False)
for an in ['pooled','EUR','AFR','AMR']:
 pair=[]
 for state in ['cases','controls']:
  d=pd.read_csv(xpanelqc/(an+'.'+state+'.gcount'),sep=r'\s+').set_index('ID');m=pd.read_csv(xpanelqc/(an+'.'+state+'.vmiss'),sep=r'\s+').set_index('ID');assert len(d)==376171 and d.index.equals(m.index) and d.index.is_unique
  ac=d.HET_REF_ALT_CTS+2*d.TWO_ALT_GENO_CTS+d.HAP_ALT_CTS
  ancount=2*(d.HOM_REF_CT+d.HET_REF_ALT_CTS+d.TWO_ALT_GENO_CTS)+d.HAP_REF_CT+d.HAP_ALT_CTS
  assert ((ac>=0)&(ac<=ancount)).all();assert (d.MISSING_CT==m.MISSING_AND_HETHAP_CT).all();assert (d.OBS_CT+d.MISSING_CT==m.OBS_CT).all()
  pair.append(pd.DataFrame({'AC':ac,'AN':ancount,'MAC':np.minimum(ac,ancount-ac),'CALL_RATE':1-m.MISSING_AND_HETHAP_CT/m.OBS_CT}))
  keepids=set(pd.read_csv(xpanelqc/(an+'.'+state+'.keep.tsv'),sep='\t',header=None,dtype=str)[1]);sm=xps.IID.isin(keepids).to_numpy()
  for vi,(_,row) in enumerate(xpvar.iterrows()):
   g=np.asarray(xgt[vi]);valid=sm&(g>=0)&~((xsex==1)&(g==1));manual_ac=float(np.where(xsex[valid]==1,g[valid]/2,g[valid]).sum());manual_an=int(np.where(xsex[valid]==1,1,2).sum());sid=row['ID'];assert int(ac.loc[sid])==manual_ac and int(ancount.loc[sid])==manual_an,(an,state,sid)
 a,c=pair;ac=a.AC+c.AC;nn=a.AN+c.AN;mac=np.minimum(ac,nn-ac);maf=mac/nn
 mask=(a.MAC>=20)&(c.MAC>=20)&(mac>=100)&(maf>=.01)&(a.CALL_RATE>=.98)&(c.CALL_RATE>=.98)
 strict=mask&xqchwe[an].reindex(mask.index).ge(1e-15)
 xqcmasks[an]=pd.DataFrame({'count_call_pass':mask,'plus_two_sided_HWE_pass':strict,'biological_MAF':maf});xqcmasks[an].to_csv(xpanelqc/(an+'.variant_QC_mask.private.tsv'),sep='\t')
 xqcsummary.append({'stratum':an,'original_panel':len(mask),'count_call_pass':int(mask.sum()),'plus_two_sided_HWE_pass':int(strict.sum()),'four_leads':[{**{'ID':sid},**xqcmasks[an].loc[sid,['count_call_pass','plus_two_sided_HWE_pass']].astype(bool).to_dict()} for sid in xpvar.ID]})
(xpanelqc/'variant_QC_summary.json').write_text(json.dumps(xqcsummary,indent=2));print(json.dumps(xqcsummary,indent=2));print('Four-lead direct biological allele counts verified in all 8 case/control groups')