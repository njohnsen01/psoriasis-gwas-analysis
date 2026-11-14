import os
import pandas as pd, numpy as np, hashlib, json
from scipy.stats import chi2
from pathlib import Path
full_path=Path(f"/home/jupyter/psoriasis_regenie_loco_step2/attempt_{os.environ['PSO_PRIMARY_RUN_ID']}/final_aggregate/psoriasis_corrected_regenie_all_results.tsv.gz")
audit_dir=Path("/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09")
parts=[]
for chunk in pd.read_csv(full_path,sep="	",usecols=["CHROM","POS","LOG10P","P","CHISQ","N","TEST","EXTRA"],dtype={"CHROM":str},chunksize=500000):
    parts.append(chunk)
frame=pd.concat(parts,ignore_index=True)
valid=np.isfinite(frame.LOG10P)&(frame.LOG10P>=0)
mask_mhc=frame.CHROM.eq("6")&frame.POS.between(25000000,34000000)
masks={"all":valid,"all_no_MHC":valid&~mask_mhc,"autosomes":valid&frame.CHROM.isin([str(i) for i in range(1,23)]),"autosomes_no_MHC":valid&frame.CHROM.isin([str(i) for i in range(1,23)])&~mask_mhc,"X":valid&frame.CHROM.isin(["X","23"])}
results=[];points=[]
for label,mask in masks.items():
    a=np.sort(frame.loc[mask,"LOG10P"].to_numpy(float))[::-1];n=len(a)
    centers=a[[(n-1)//2,n//2]]
    lam=float(np.mean(chi2.isf(10.0**(-centers),1))/chi2.ppf(.5,1))
    stat=frame.loc[mask,"CHISQ"].to_numpy(float)
    results.append({"group":label,"n":n,"n_GWS":int((a>-np.log10(5e-8)).sum()),"lambda_GC":lam,"lambda_from_CHISQ":float(np.nanmedian(stat)/chi2.ppf(.5,1)),"n_p_equal_1":int((a==0).sum()),"max_log10p":float(a[0]),"min_variant_N":int(frame.loc[mask,"N"].min()),"max_variant_N":int(frame.loc[mask,"N"].max())})
    ranks=np.unique(np.r_[np.arange(1,min(n,2000)+1),np.rint(np.geomspace(min(n,2001),n,2500)).astype(int)])
    for rank in ranks:points.append({"group":label,"rank":int(rank),"n":n,"expected":float(-np.log10((rank-.5)/n)),"observed":float(a[rank-1])})
    np.save(audit_dir/(label+"_sorted_logp.npy"),a)
summary={"source":str(full_path),"n_records":len(frame),"n_valid_LOG10P":int(valid.sum()),"n_invalid_LOG10P":int((~valid).sum()),"test_types":frame.TEST.value_counts(dropna=False).to_dict(),"n_invalid_P":int((~np.isfinite(frame.P)|(frame.P<0)|(frame.P>1)).sum()),"groups":results,"method":"all valid association LOG10P values; unit weight per test; QQ display thinning only after exact full-data ranks"}
h=hashlib.sha256()
with full_path.open("rb") as f:
    for block in iter(lambda:f.read(8*1024*1024),b""):h.update(block)
summary["sha256"]=h.hexdigest()
pd.DataFrame(points).to_csv(audit_dir/"exact_qq_points.tsv",sep="	",index=False)
(audit_dir/"full_calibration_summary.json").write_text(json.dumps(summary,indent=2))
