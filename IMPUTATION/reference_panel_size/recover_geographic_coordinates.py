#!/usr/bin/env python3

import pandas as pd
import numpy as np

xlsx = "F1146F1147-Pedigree-V2.xlsx"
target_file = "target1297_subgroup_metadata.tsv"

# Read data
ped = pd.read_excel(xlsx, sheet_name="F1146F1147")
target = pd.read_csv(target_file, sep="\t")

print("Pedigree rows:", len(ped))
print("Target rows:", len(target))

# Normalize IDs to strings
for col in ["Sheet_ID", "Mum", "RAPiD_Genomics_Sample_Code"]:
    if col in ped.columns:
        ped[col] = ped[col].astype("string").str.strip()

for col in ["Sheet_ID", "Mum", "RAPiD_Genomics_Sample_Code"]:
    if col in target.columns:
        target[col] = target[col].astype("string").str.strip()

# Existing coordinates in target
has_direct = target["Lat_O"].notna() & target["Long_O"].notna()

print("\nDirect coordinates:")
print("With direct coordinates:", has_direct.sum())
print("Without direct coordinates:", (~has_direct).sum())

# Build possible parent coordinate lookup tables

# 1. Sheet_ID -> coordinate
sheet_lookup = (
    ped.loc[
        ped["Lat_O"].notna() & ped["Long_O"].notna(),
        ["Sheet_ID", "Lat_O", "Long_O", "Groups", "Group_FFF"]
    ]
    .dropna(subset=["Sheet_ID"])
    .drop_duplicates(subset=["Sheet_ID"])
    .set_index("Sheet_ID")
)

# 2. RAPiD sample code -> coordinate
rapid_lookup = (
    ped.loc[
        ped["Lat_O"].notna() & ped["Long_O"].notna(),
        ["RAPiD_Genomics_Sample_Code", "Lat_O", "Long_O", "Groups", "Group_FFF"]
    ]
    .dropna(subset=["RAPiD_Genomics_Sample_Code"])
    .drop_duplicates(subset=["RAPiD_Genomics_Sample_Code"])
    .set_index("RAPiD_Genomics_Sample_Code")
)

# Start with original coordinates
target["Geo_Lat"] = target["Lat_O"]
target["Geo_Long"] = target["Long_O"]
target["Geo_source"] = np.where(has_direct, "direct", "missing")

# Try recovery via Mum -> Sheet_ID
recovered_sheet = 0

for i, row in target.loc[~has_direct].iterrows():
    mum = row.get("Mum")

    if pd.notna(mum) and mum in sheet_lookup.index:
        target.at[i, "Geo_Lat"] = sheet_lookup.at[mum, "Lat_O"]
        target.at[i, "Geo_Long"] = sheet_lookup.at[mum, "Long_O"]
        target.at[i, "Geo_source"] = "Mum_to_Sheet_ID"
        recovered_sheet += 1

# For still missing, try Mum -> RAPiD code
still_missing = target["Geo_Lat"].isna() | target["Geo_Long"].isna()

recovered_rapid = 0

for i, row in target.loc[still_missing].iterrows():
    mum = row.get("Mum")

    if pd.notna(mum) and mum in rapid_lookup.index:
        target.at[i, "Geo_Lat"] = rapid_lookup.at[mum, "Lat_O"]
        target.at[i, "Geo_Long"] = rapid_lookup.at[mum, "Long_O"]
        target.at[i, "Geo_source"] = "Mum_to_RAPiD"
        recovered_rapid += 1

final_has = target["Geo_Lat"].notna() & target["Geo_Long"].notna()

print("\nRecovered through mother mapping:")
print("Mum -> Sheet_ID:", recovered_sheet)
print("Mum -> RAPiD:", recovered_rapid)

print("\nFinal geographic coverage:")
print("With coordinates:", final_has.sum())
print("Without coordinates:", (~final_has).sum())
print("Coverage:", round(final_has.mean(), 4))

print("\nCoordinate source:")
print(target["Geo_source"].value_counts(dropna=False))

# Coordinate ranges
geo = target.loc[final_has].copy()

print("\nLatitude summary:")
print(geo["Geo_Lat"].describe())

print("\nLongitude summary:")
print(geo["Geo_Long"].describe())

print("\nGenetic groups among samples with geographic coordinates:")
print(geo["Groups"].fillna("Unknown").value_counts())

# Save files
target.to_csv(
    "target1297_geographic_metadata.tsv",
    sep="\t",
    index=False
)

geo[
    ["Sheet_ID", "Groups", "Mum", "Geo_Lat", "Geo_Long", "Geo_source"]
].to_csv(
    "target_with_coordinates.tsv",
    sep="\t",
    index=False
)

print("\nSaved:")
print("target1297_geographic_metadata.tsv")
print("target_with_coordinates.tsv")
