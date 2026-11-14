from pathlib import Path
import os,json,pandas as pd,numpy as np
from google.cloud import bigquery
root=Path("/home/jupyter")
audit_dir=root/"psoriasis_manuscript_cohort_calibration_2026-09-09"
frozen=root/f"psoriasis_regenie_loco_indexed_range_step1/attempt_{os.environ['PSO_PRIMARY_RUN_ID']}/frozen_inputs"
ph=pd.read_csv(frozen/"psoriasis_step1_pheno.tsv",sep="	",dtype={"IID":str,"FID":str})
co=pd.read_csv(frozen/"psoriasis_step1_covariates.tsv",sep="	",dtype={"IID":str,"FID":str})
anc=pd.read_csv(root/"psoriasis_phase2_reproducibility/secure_participant_level/ancestry_preds.tsv",sep="	",usecols=["research_id","ancestry_pred"],dtype={"research_id":str})
base=ph.merge(co,on=["FID","IID"],validate="one_to_one").merge(anc,left_on="IID",right_on="research_id",how="left",validate="one_to_one")
assert len(base)==392155 and base.psoriasis_case.sum()==5805
ages=base.groupby("psoriasis_case").age.agg(n="count",median="median",q1=lambda x:x.quantile(.25),q3=lambda x:x.quantile(.75),mean="mean",sd="std")
ages.to_csv(audit_dir/"age_summary.tsv",sep="	")
client=bigquery.Client(project=os.environ["GOOGLE_PROJECT"])
cdr=os.environ["WORKSPACE_CDR"]
dem=client.query(f"SELECT person_id,dob,has_ehr_data FROM `{cdr}.cb_search_person`").to_dataframe();dem["person_id"]=dem.person_id.astype(str)
ehr=client.query((audit_dir/"ehr_observation_span.sql").read_text()).to_dataframe();ehr["person_id"]=ehr.person_id.astype(str)
clinical=base.merge(dem,left_on="IID",right_on="person_id",how="left",validate="one_to_one").merge(ehr,on="person_id",how="left",validate="one_to_one")
clinical["span_years"]=(pd.to_datetime(clinical.last_ehr_date)-pd.to_datetime(clinical.first_ehr_date)).dt.days/365.25
clinical.loc[clinical.has_ehr_data!=1,"span_years"]=np.nan
summary=clinical.groupby("psoriasis_case").span_years.agg(n_available="count",median="median",q1=lambda x:x.quantile(.25),q3=lambda x:x.quantile(.75),mean="mean",sd="std",zero_days=lambda x:(x==0).sum())
summary["n_total"]=clinical.groupby("psoriasis_case").size();summary["n_unavailable"]=summary.n_total-summary.n_available
summary.to_csv(audit_dir/"ehr_span_summary.tsv",sep="	")
# Never export participant-level frames or unsuppressed small-cell ancestry counts.
