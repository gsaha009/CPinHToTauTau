'''
Author : @oponcet
Date : 29-01-2025
Script to combine the json files into one.
'''
import json
import glob


year = "2023_postBPix"  # Change this to the desired year
# Define the files and the structure
input_files = glob.glob(f"/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/script_FF/fake_factor_derivation/outputs/{year}/closure_correction/*/*.json")  # Replace with the path to your files
merged_structure = {
    "schema_version": 2,
    "description": "Closure corrections for the httcp analysis",
    "corrections": [
        {
            "name": "closure_corrections_fit",
            "description": "Fit closure corrections for all decay modes and jet categories",
            "version": 1,
            "inputs": [
                {
                    "name": "met_var_qcd_h1",
                    "type": "real",
                    "description": "Transverse momentum of the tau"
                },
                {
                    "name": "dm",
                    "type": "int",
                    "description": "Reconstructed tau decay mode of leading tau: pi_1, rho_1, a1dm11_1, a1dm10_1, a1dm2_1"
                },
                {
                    "name": "njets",
                    "type": "int",
                    "description": "Number of jets in the event (has_0j, has_1j, has_2j)"
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
                "input": "dm",
                "content": []
            }
        }
    ]
}

# Process each file and merge into the main structure
dm_content = {}

for file in input_files:
    with open(file, "r") as f:
        data = json.load(f)

    correction_data = data["corrections"][0]["data"]
    dm = correction_data["content"][0]["key"]  # Get the decay mode
    njets_content = correction_data["content"][0]["value"]["content"]  # Get jet-specific data

    # Add the decay mode and its associated jet categories to the main content
    if dm not in dm_content:
        dm_content[dm] = {"nodetype": "category", "input": "njets", "content": njets_content}
    else:
        # If the decay mode exists, extend its `njets` content
        dm_content[dm]["content"].extend(njets_content)

# Add all decay modes to the merged structure
merged_structure["corrections"][0]["data"]["content"] = [
    {"key": dm, "value": value} for dm, value in dm_content.items()
]


# Save the merged JSON
output_file = f"/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/script_FF/fake_factor_derivation/outputs/{year}/closure_correction/json/closure_correction_{year}.json"
# guranty that the output directory exists
import os
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as f:
    json.dump(merged_structure, f, indent=4)

print(f"Merged JSON saved to {output_file}")
