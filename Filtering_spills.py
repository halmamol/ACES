print("Starting Spills Filter Algorithm")

import uproot
import numpy as np
import matplotlib.pyplot as plt
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
run_number = "1766"   # Run number
nhits_threshold = 300 # Threshold for nHits
nhits_window = 5000   # Window for nHits
death_time = 6000     # Death time for nHits

# Paths for files #############################################################################################
root_dir = f"/mnt/lustre/scratch/nlsas/home/usc/ie/dcr/hk/raw_data/{run_number}/"
root_file_path = f"{root_dir}WCTE_offline_R{run_number}S0P{partition}.root"
output_path = f"/mnt/lustre/scratch/nlsas/home/usc/ie/dcr/hk/processed_wcte_data"

# Showing results
if partition == "all":
    print("Analysing all partitions.")

    print(f"Loading Run {run_number}...")
    root_files = sorted(glob.glob(os.path.join(root_dir, "*.root")))
    root_files = sorted(root_files, key=lambda file_path: int(file_path.split("P")[-1].split(".")[0]))
    root_files = root_files[:-1]

    print(f"Found {len(root_files)} ROOT files.")

    times_branch_sorted, times_branch_sorted_TOF, charge_branch_sorted, mpmt_id_branch_sorted, pmt_id_branch_sorted, event_number_branch, _ = functions_spills.multiple_partition(root_files)

    print("Runs loaded.")
    N_events = max(event_number_branch) + 1

print(f"Total number of events in run {run_number}: {N_events}")

# Filter spills using nHits threshold ###################################################################################
print(f"Applying filter for Run {run_number}...")

(
    times_branch_modified_TOF,
    charge_branch_modified,
    mpmt_id_branch_modified,
    pmt_id_branch_modified,
    threshold_times,
    deleted_index_dict,
) = functions_spills.repeat_spills_nHits_with_channels(
    event_number_branch=event_number_branch,
    times_branch_sorted_TOF=times_branch_sorted,
    #times_branch_sorted_TOF=times_branch_sorted_TOF,
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
times_values, times_offsets = flatten_values_and_offsets(times_branch_modified_TOF)
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

with open(f'{output_path}filtered_file_{run_number}.pkl', 'wb') as f:
    pickle.dump(payload, f)

print(f"Saved single consolidated file:")
print(f" - {output_path}filtered_file_{run_number}.pkl (contains times_TOF, charge, mpmt_id, metadata)")