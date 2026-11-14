# ============================================================
# Block 6 — Recreate raw laboratory checkpoint if absent
# ============================================================

from google.cloud import bigquery

REBUILD_RAW_LABS = not RAW_LABS.exists()

MEASUREMENT_CONCEPT_IDS = [
    3001123, 3002030, 3003215, 3004327, 3007461, 3008342,
    3013498, 3013650, 3017501, 3017732, 3018010, 3019198,
    3024386, 3024929, 3027368, 3031586, 3037234, 3037511,
    3038058, 3043111, 37028517, 37031129, 37037425, 37040154,
    37045722, 37048421, 37049637, 37056087, 37060904, 37062698,
    40772531, 40779159, 40785788, 43055366,
]

if not REBUILD_RAW_LABS:
    print("[SKIP] Retained raw laboratory checkpoint found:", RAW_LABS)
else:
    client = bigquery.Client(project=GOOGLE_PROJECT)
    cohort_ids = cohort["person_id"].drop_duplicates().tolist()
    concept_ids_sql = ",".join(map(str, MEASUREMENT_CONCEPT_IDS))

    query = f"""
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
    LEFT JOIN `{WORKSPACE_CDR}.concept` c
      ON m.measurement_concept_id = c.concept_id
    LEFT JOIN `{WORKSPACE_CDR}.concept` u
      ON m.unit_concept_id = u.concept_id
    LEFT JOIN `{WORKSPACE_CDR}.concept` sc
      ON m.measurement_source_concept_id = sc.concept_id
    WHERE CAST(m.person_id AS STRING) IN UNNEST(@cohort_ids)
      AND m.value_as_number IS NOT NULL
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

    chunk_size = 10_000
    first_write = True
    total_rows = 0

    for start in range(0, len(cohort_ids), chunk_size):
        chunk_ids = cohort_ids[start:start + chunk_size]
        config = bigquery.QueryJobConfig(
            query_parameters=[
                bigquery.ArrayQueryParameter(
                    "cohort_ids",
                    "STRING",
                    chunk_ids,
                )
            ]
        )
        chunk = client.query(query, job_config=config).to_dataframe()
        chunk.to_csv(
            RAW_LABS,
            sep="\t",
            index=False,
            mode="w" if first_write else "a",
            header=first_write,
        )
        first_write = False
        total_rows += len(chunk)
        print(
            f"Chunk {start // chunk_size + 1}: "
            f"{len(chunk):,} rows; cumulative {total_rows:,}"
        )

    print("Raw laboratory rows written:", total_rows)
