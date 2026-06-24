print("Starting Spills Filter Algorithm")

import uproot
import numpy as np
import pandas as pd
import sys
import os
sys.path.append(os.path.abspath("/mnt/netapp2/Store_uni/home/usc/ie/dcr/software/hk"))
import functions_spills
import glob
import os
import argparse
import pickle

# Create parser
parser = argparse.ArgumentParser(description="Partition to analyse")

# Adding optional arguments 
parser.add_argument(
    "--partition",
    help="Specify the partition number (if not specified, all partitions will be used)",
    default="all"
)

# Parsing arguments
args = parser.parse_args()
partition = args.partition

# Other arguments for the analysis
run_number = "2384"  # Run number
drun = "2386"
nhits_threshold = 300  # Threshold for nHits
nhits_window = 5000  # Window for nHits
death_time = 6000  # Death time for nHits

version = 'v0_5'

# Paths for files #############################################################################################
root_dir = f"/data/halmazan/WCTE/data/{run_number}/{version}/"
root_file_path = f"{root_dir}WCTE_offline_R{run_number}S0P{partition}.root"
output_path = f"/scratch/halmazan/WCTE/files/filtered_files/"

source_pos = []
sources_pos = {

    '2386': [0.0, 729.725, 0.0],
    '2387': [0.0, 119.725, 0.0],
    '2388': [0.0, -160.275, 0.0],
    '2389': [-476.4, 729.725, 291.7],
    '2390': [-476.4, 119.725, 291.7],

    #'2386': [0.0, +30.5, 0.0],
    #'2387': [0.0, -30.5, 0.0],
    #'2388': [0.0, -58.5, 0.0],
    #'2389': [-47.64, +30.5, 29.17],
    #'2390': [-47.64, -30.5, 29.17], 
}
if run_number in sources_pos.keys():
    source_pos = sources_pos[run_number]
else:
    source_pos = sources_pos[drun]

# Showing results
if partition == "all":
    print("Analysing all partitions.")

    print(f"Loading Run {run_number}...")
    root_files = sorted(glob.glob(os.path.join(root_dir, "*.root")))
    root_files = sorted(root_files, key=lambda file_path: int(file_path.split("P")[-1].split(".")[0]))
    root_files = root_files[:-1]

    print(f"Found {len(root_files)} ROOT files.")

    times_branch_sorted, times_branch_sorted_TOF, charge_branch_sorted, mpmt_id_branch_sorted, pmt_id_branch_sorted, event_number_branch, _ = functions_spills.multiple_partition(root_files, source_pos)

    print("Runs loaded.")
    N_events = max(event_number_branch) + 1

print(f"Total number of events in run {run_number}: {N_events}")

# Filter spills using nHits threshold ###################################################################################
print(f"Applying filter for Run {run_number}...")

(
    times_branch_modified_TOF,
    times_branch_modified,
    charge_branch_modified,
    mpmt_id_branch_modified,
    pmt_id_branch_modified,
    threshold_times,
    deleted_index_dict,
) = functions_spills.repeat_spills_nHits_with_channels(
    event_number_branch=event_number_branch,
    #times_branch_sorted_TOF=times_branch_sorted,
    times_branch_sorted_TOF=times_branch_sorted_TOF,
    times_branch_sorted=times_branch_sorted,
    charge_branch_sorted=charge_branch_sorted,
    mpmt_id_branch_sorted=mpmt_id_branch_sorted,
    pmt_id_branch_sorted=pmt_id_branch_sorted,
    threshold=nhits_threshold,
    window=nhits_window,
    death_window=death_time,
)

print("nHits filter applied.")

# Optionally still save deleted indices separately (kept for compatibility) ############################################
with open(f'{output_path}deleted_indices_nHits_{run_number}.pkl', 'wb') as f:
    pickle.dump(deleted_index_dict, f)

total_elements = sum(len(v) for v in threshold_times.values())
print(
    f"Total number of elements in threshold_times for run {run_number}:",
    total_elements,
    max(event_number_branch),
    total_elements / max(event_number_branch),
)

def flatten_values_and_offsets(list_of_arrays):
    if len(list_of_arrays) == 0:
        return np.array([], dtype=float), np.array([0], dtype=int)
    values = np.concatenate(list_of_arrays) if len(list_of_arrays) > 0 else np.array([], dtype=float)
    offsets = np.cumsum([0] + [len(arr) for arr in list_of_arrays])
    return values, offsets

# Build a single consolidated payload and save into ONE file ###########################################################
times_values_TOF, times_offsets_TOF = flatten_values_and_offsets(times_branch_modified_TOF)
times_values, times_offsets = flatten_values_and_offsets(times_branch_modified)
charge_values, charge_offsets = flatten_values_and_offsets(charge_branch_modified)
mpmt_values, mpmt_offsets = flatten_values_and_offsets(mpmt_id_branch_modified)
pmt_values, pmt_offsets = flatten_values_and_offsets(pmt_id_branch_modified)

payload = {
    "run_number": run_number,
    "nhits_threshold": nhits_threshold,
    "nhits_window": nhits_window,
    "death_time": death_time,
    "threshold_times": threshold_times,
    "deleted_indices_nHits": deleted_index_dict,
    "times_TOF": {
        "values": times_values_TOF,
        "offsets": times_offsets_TOF,
    },
    "times": {
        "values": times_values,
        "offsets": times_offsets,
    },
    "charge": {
        "values": charge_values,
        "offsets": charge_offsets,
    },
    "mpmt_id": {
        "values": mpmt_values,
        "offsets": mpmt_offsets,
    },
    "pmt_id": {
        "values": pmt_values,
        "offsets": pmt_offsets,
    },
}

if run_number == '2384':
    with open(f'{output_path}filtered_file_{run_number}_{version}_wTOF{drun}.pkl', 'wb') as f:
        pickle.dump(payload, f)
else:
    with open(f'{output_path}filtered_file_{run_number}_{version}_wTOF.pkl', 'wb') as f:
        pickle.dump(payload, f)

print(f"Saved single consolidated file:")
print(f" - {output_path}filtered_file_{run_number}_wTOF.pkl (contains times_TOF, charge, mpmt_id, metadata)")