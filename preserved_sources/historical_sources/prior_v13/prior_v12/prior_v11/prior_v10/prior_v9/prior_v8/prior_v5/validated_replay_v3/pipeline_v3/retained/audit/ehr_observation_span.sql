WITH dates AS (
SELECT person_id, MIN(condition_start_date) lo, MAX(condition_start_date) hi FROM `wb-silky-artichoke-2408.C2024Q3R8.condition_occurrence` GROUP BY person_id
UNION ALL SELECT person_id, MIN(procedure_date), MAX(procedure_date) FROM `wb-silky-artichoke-2408.C2024Q3R8.procedure_occurrence` GROUP BY person_id
UNION ALL SELECT person_id, MIN(drug_exposure_start_date), MAX(drug_exposure_start_date) FROM `wb-silky-artichoke-2408.C2024Q3R8.drug_exposure` GROUP BY person_id
UNION ALL SELECT person_id, MIN(visit_start_date), MAX(visit_start_date) FROM `wb-silky-artichoke-2408.C2024Q3R8.visit_occurrence` WHERE visit_type_concept_id != 44818519 GROUP BY person_id
) SELECT person_id, MIN(lo) first_ehr_date, MAX(hi) last_ehr_date FROM dates GROUP BY person_id