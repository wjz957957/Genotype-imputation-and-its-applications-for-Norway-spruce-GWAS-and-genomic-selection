import pandas as pd

xlsx = "F1146F1147-Pedigree-V2.xlsx"
ids_file = "target1297.ids"

df = pd.read_excel(xlsx, sheet_name="F1146F1147")
target_ids = pd.read_csv(ids_file, header=None, names=["Sheet_ID"])

print("Excel rows:", len(df))
print("Target samples:", len(target_ids))

merged = target_ids.merge(
    df,
    on="Sheet_ID",
    how="left"
)

print("Matched rows:", merged["Trial"].notna().sum())
print("Unmatched:", merged["Trial"].isna().sum())

print("\nGenetic groups:")
print(merged["Groups"].value_counts(dropna=False))

print("\nGroup_FFF:")
print(merged["Group_FFF"].value_counts(dropna=False))

print("\nWith both coordinates:")
print(merged[["Lat_O", "Long_O"]].notna().all(axis=1).sum())

print("\nTrial:")
print(merged["Trial"].value_counts(dropna=False))

merged.to_csv(
    "target1297_subgroup_metadata.tsv",
    sep="\t",
    index=False
)

unmatched = merged.loc[
    merged["Trial"].isna(),
    ["Sheet_ID"]
]

unmatched.to_csv(
    "target1297_unmatched.ids",
    sep="\t",
    index=False,
    header=False
)

print("\nSaved: target1297_subgroup_metadata.tsv")
print("Saved: target1297_unmatched.ids")
