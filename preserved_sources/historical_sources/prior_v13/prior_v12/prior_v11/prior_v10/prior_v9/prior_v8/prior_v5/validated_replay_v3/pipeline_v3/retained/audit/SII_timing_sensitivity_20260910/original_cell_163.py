# ============================================================
# Block 8 — Harmonize units and calculate participant medians
# ============================================================

labs = pd.read_csv(
    CLASSIFIED_LABS,
    sep="\t",
    dtype={"person_id": str, "measurement_concept_id": str},
    low_memory=False,
)

labs["value_as_number"] = pd.to_numeric(
    labs["value_as_number"],
    errors="coerce",
)
labs = labs.dropna(subset=["value_as_number"]).copy()
labs = labs[labs["value_as_number"] > 0].copy()

labs["unit_all_lower"] = (
    labs["unit_name"].astype(str).str.lower()
    + " | "
    + labs["unit_source_value"].astype(str).str.lower()
)
labs["value_std"] = labs["value_as_number"].copy()

ldl_hdl_mmol = (
    labs["biomarker"].isin(["LDL", "HDL"])
    & labs["unit_all_lower"].str.contains("mmol", na=False)
)
tg_mmol = (
    labs["biomarker"].eq("TG")
    & labs["unit_all_lower"].str.contains("mmol", na=False)
)

labs.loc[ldl_hdl_mmol, "value_std"] = (
    labs.loc[ldl_hdl_mmol, "value_as_number"] * 38.67
)
labs.loc[tg_mmol, "value_std"] = (
    labs.loc[tg_mmol, "value_as_number"] * 88.57
)

cbc = labs["biomarker"].isin(["ANC", "ALC", "platelet_count"])
cell_per_ul = labs["unit_all_lower"].str.contains(
    r"/ul|/µl|/mcl|per microliter|cells/ul|cells per microliter",
    case=False,
    na=False,
    regex=True,
)
already_thousands = labs["unit_all_lower"].str.contains(
    r"10\*3|10\^3|thousand|k/ul|10\*9/l|billion per liter",
    case=False,
    na=False,
    regex=True,
)
convert_cells_to_thousands = cbc & cell_per_ul & ~already_thousands
labs.loc[convert_cells_to_thousands, "value_std"] = (
    labs.loc[convert_cells_to_thousands, "value_as_number"] / 1000
)

plausibility_ranges = {
    "LDL": (1, 500),
    "HDL": (1, 250),
    "TG": (1, 2000),
    "platelet_count": (1, 2000),
    "ANC": (0.01, 100),
    "ALC": (0.01, 100),
}

clean_parts = []
for biomarker, (lower, upper) in plausibility_ranges.items():
    part = labs[labs["biomarker"].eq(biomarker)].copy()
    part = part[part["value_std"].between(lower, upper)]
    clean_parts.append(part)

labs_clean = pd.concat(clean_parts, ignore_index=True)

labs_person = (
    labs_clean
    .groupby(["person_id", "biomarker"])["value_std"]
    .median()
    .reset_index()
    .pivot(index="person_id", columns="biomarker", values="value_std")
    .reset_index()
)
labs_person.columns.name = None

print("Clean laboratory rows:", len(labs_clean))
print("Participants with at least one retained biomarker:", len(labs_person))
print("Nonmissing source components:")
print(labs_person.notna().sum().to_dict())
