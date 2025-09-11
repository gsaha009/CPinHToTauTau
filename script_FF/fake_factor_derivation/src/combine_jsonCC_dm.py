'''
Author : @oponcet
Date : 29-01-2025
Script to combine multiple closure_correction JSON files into one,
grouped by njets (0,1,2), without dm categorization.
'''
import json
import glob
import os

year = "2023"  # Change this to the desired year

# Find input files
input_files = glob.glob(
    f"/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/script_FF/fake_factor_derivation/outputs/{year}/closure_correction/*/*.json"
)

merged_structure = {
    "schema_version": 2,
    "description": "Closure corrections for the httcp analysis (no dm categorization)",
    "corrections": [
        {
            "name": "closure_corrections_fit",
            "description": "Fit closure corrections merged across all input files",
            "version": 1,
            "inputs": [
                {
                    "name": "met_var_qcd_h1",
                    "type": "real",
                    "description": "Transverse momentum of the tau"
                },
                {
                    "name": "njets",
                    "type": "int",
                    "description": "Number of jets in the event (0, 1, 2)"
                },
                {
                    "name": "syst",
                    "type": "string",
                    "description": "Systematic variations: 'nom', 'up', 'down'"
                }
            ],
            "output": {
                "name": "closure_correction",
                "type": "real",
                "description": "Closure correction to apply to data-MC"
            },
            "data": {
                "nodetype": "category",
                "input": "njets",
                "content": []
            }
        }
    ]
}

# Dictionary keyed by njets (0,1,2), each holding syst content
njets_content = {}

for file in input_files:
    with open(file, "r") as f:
        data = json.load(f)

    # Extract njets blocks
    njets_blocks = data["corrections"][0]["data"]["content"]

    for nj in njets_blocks:
        nj_key = nj["key"]       # 0, 1, 2
        syst_block = nj["value"] # syst category

        if nj_key not in njets_content:
            # Initialize syst container for this njets
            njets_content[nj_key] = {
                "nodetype": "category",
                "input": "syst",
                "content": []
            }

        # Merge syst entries into this njets
        njets_content[nj_key]["content"].extend(syst_block["content"])

# # Rebuild content list sorted by njets key
merged_structure["corrections"][0]["data"]["content"] = [
    {"key": nj, "value": value} for nj, value in sorted(njets_content.items())
]

# Save the merged JSON
output_file = f"/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/script_FF/fake_factor_derivation/outputs/{year}/closure_correction/json/closure_correction_{year}_merged.json"
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as f:
    json.dump(merged_structure, f, indent=4)

print(f"Merged JSON saved to {output_file}")
