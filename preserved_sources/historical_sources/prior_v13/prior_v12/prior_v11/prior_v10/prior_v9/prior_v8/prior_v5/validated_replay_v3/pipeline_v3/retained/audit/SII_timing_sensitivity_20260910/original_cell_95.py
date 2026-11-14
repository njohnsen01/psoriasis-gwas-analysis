# ============================================================
# Step 3: Extract biomarker labs from measurement table
# ============================================================

import os
import pandas as pd
import numpy as np

try:
    import pandas_gbq
    read_gbq = pandas_gbq.read_gbq
except ImportError:
    read_gbq = pd.read_gbq

GOOGLE_PROJECT = os.environ.get("GOOGLE_PROJECT")
WORKSPACE_CDR = os.environ.get("WORKSPACE_CDR", "wb-silky-artichoke-2408.C2024Q3R8")

print("GOOGLE_PROJECT =", GOOGLE_PROJECT)
print("WORKSPACE_CDR =", WORKSPACE_CDR)

cohort_ids = (
    analysis_base["person_id"]
    .astype(str)
    .drop_duplicates()
    .tolist()
)

print("[INFO] Cohort IDs:", len(cohort_ids))

# Candidate lab concept IDs from prior AoU project + additional lipid/CBC IDs
MEASUREMENT_CONCEPT_IDS = [
    3001123, 3002030, 3003215, 3004327, 3007461, 3008342,
    3013498, 3013650, 3017501, 3017732, 3018010, 3019198,
    3024386, 3024929, 3027368, 3031586, 3037234, 3037511,
    3038058, 3043111, 37028517, 37031129, 37037425, 37040154,
    37045722, 37048421, 37049637, 37056087, 37060904, 37062698,
    40772531, 40779159, 40785788, 43055366
]

CHUNK_SIZE = 10000
all_lab_chunks = []

for i in range(0, len(cohort_ids), CHUNK_SIZE):
    chunk = cohort_ids[i:i + CHUNK_SIZE]
    ids_sql = ",".join([f"'{x}'" for x in chunk])
    concept_ids_sql = ",".join(map(str, MEASUREMENT_CONCEPT_IDS))

    labs_sql = f"""
    WITH cohort AS (
      SELECT person_id
      FROM UNNEST([{ids_sql}]) AS person_id
    )
    SELECT
      CAST(m.person_id AS STRING) AS person_id,
      m.measurement_concept_id,
      c.concept_name AS lab_name,
      c.concept_code AS lab_code,
      c.vocabulary_id AS lab_vocabulary,
      m.measurement_datetime,
      m.value_as_number,
      m.value_source_value,
      m.unit_concept_id,
      u.concept_name AS unit_name,
      m.unit_source_value,
      m.measurement_source_value,
      m.measurement_source_concept_id,
      sc.concept_name AS source_lab_name,
      sc.concept_code AS source_lab_code,
      sc.vocabulary_id AS source_vocabulary
    FROM `{WORKSPACE_CDR}.measurement` m
    JOIN cohort
      ON CAST(m.person_id AS STRING) = cohort.person_id
    LEFT JOIN `{WORKSPACE_CDR}.concept` c
      ON m.measurement_concept_id = c.concept_id
    LEFT JOIN `{WORKSPACE_CDR}.concept` u
      ON m.unit_concept_id = u.concept_id
    LEFT JOIN `{WORKSPACE_CDR}.concept` sc
      ON m.measurement_source_concept_id = sc.concept_id
    WHERE m.value_as_number IS NOT NULL
      AND (
        m.measurement_concept_id IN ({concept_ids_sql})
        OR LOWER(c.concept_name) LIKE '%cholesterol in ldl%'
        OR LOWER(c.concept_name) LIKE '%ldl cholesterol%'
        OR LOWER(c.concept_name) LIKE '%cholesterol in hdl%'
        OR LOWER(c.concept_name) LIKE '%hdl cholesterol%'
        OR LOWER(c.concept_name) LIKE '%triglyceride%'
        OR LOWER(c.concept_name) LIKE '%platelet%'
        OR LOWER(c.concept_name) LIKE '%neutrophil%'
        OR LOWER(c.concept_name) LIKE '%lymphocyte%'
        OR LOWER(sc.concept_name) LIKE '%ldl%'
        OR LOWER(sc.concept_name) LIKE '%hdl%'
        OR LOWER(sc.concept_name) LIKE '%triglyceride%'
        OR LOWER(sc.concept_name) LIKE '%platelet%'
        OR LOWER(sc.concept_name) LIKE '%neutrophil%'
        OR LOWER(sc.concept_name) LIKE '%lymphocyte%'
      )
    """

    lab_chunk = read_gbq(
        labs_sql,
        project_id=GOOGLE_PROJECT,
        dialect="standard",
        use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
        progress_bar_type=None
    )

    all_lab_chunks.append(lab_chunk)
    print(f"[INFO] Finished chunk {i//CHUNK_SIZE + 1}: rows={lab_chunk.shape[0]}")

labs_raw = pd.concat(all_lab_chunks, ignore_index=True)

labs_raw.to_csv(
    "psoriasis_biomarker_labs_raw.tsv",
    sep="\t",
    index=False
)

print("[INFO] Raw lab rows:", labs_raw.shape)
display(labs_raw.head())