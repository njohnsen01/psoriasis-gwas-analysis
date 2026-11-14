
    WITH cohort AS (
      SELECT person_id
      FROM UNNEST(@cohort_ids) AS person_id
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
    FROM `wb-silky-artichoke-2408.C2024Q3R8.measurement` m
    JOIN cohort
      ON CAST(m.person_id AS STRING) = cohort.person_id
    LEFT JOIN `wb-silky-artichoke-2408.C2024Q3R8.concept` c
      ON m.measurement_concept_id = c.concept_id
    LEFT JOIN `wb-silky-artichoke-2408.C2024Q3R8.concept` u
      ON m.unit_concept_id = u.concept_id
    LEFT JOIN `wb-silky-artichoke-2408.C2024Q3R8.concept` sc
      ON m.measurement_source_concept_id = sc.concept_id
    WHERE m.value_as_number IS NOT NULL
      AND (
        m.measurement_concept_id IN (3001123,3002030,3003215,3004327,3007461,3008342,3013498,3013650,3017501,3017732,3018010,3019198,3024386,3024929,3027368,3031586,3037234,3037511,3038058,3043111,37028517,37031129,37037425,37040154,37045722,37048421,37049637,37056087,37060904,37062698,40772531,40779159,40785788,43055366)
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
    