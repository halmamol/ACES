print("Starting Spills Filter Algorithm")

import uproot
import numpy as np
import pandas as pd
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
run_number = "2386"  # Run number
nhits_threshold = 300  # Threshold for nHits
nhits_window = 5000  # Window for nHits
death_time = 6000  # Death time for nHits

# Paths for files #############################################################################################
root_dir = f"/data/halmazan/WCTE/MC/"
root_file_path = f"{root_dir}wcte_ambe_mc_digidata_{run_number}.root" #wcte_ambe_mc_digidata_0.root
#root_file_path = f"{root_dir}WCTE_offline_R{run_number}S0P{partition}.root"
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
    print('Run not on list')

#root_files = glob.glob(os.path.join(root_dir, "*.root"))
# Showing results
#if len(root_files) == 0:
#    raise FileNotFoundError(f"No .root files found in {root_dir}")
#if len(root_files) > 1:
#    print(f"Warning: found {len(root_files)} files, using the first one")

file_path = root_file_path

times_branch_sorted = []
times_branch_sorted_TOF = []
charge_branch_sorted = []
mpmt_id_branch_sorted = []
pmt_id_branch_sorted = []
event_number_branch = []

# Contador global de eventos
event_offset = 0

print(f"Procesando archivo: {file_path}")
file = uproot.open(file_path)
tree = file["WCTEReadoutWindows"]

#times_branch_sorted_i, times_branch_sorted_TOF_i, charge_branch_sorted_i, mpmt_id_branch_sorted_i, pmt_id_branch_sorted_i, event_number_branch_i = functions_spills.initial_treatment(tree)
times_branch_sorted_i, times_branch_sorted_TOF_i, charge_branch_sorted_i, mpmt_id_branch_sorted_i, pmt_id_branch_sorted_i, event_number_branch_i = functions_spills.initial_treatment(tree, source_pos)

times_branch_sorted.extend(times_branch_sorted_i)
times_branch_sorted_TOF.extend(times_branch_sorted_TOF_i)
charge_branch_sorted.extend(charge_branch_sorted_i)
mpmt_id_branch_sorted.extend(mpmt_id_branch_sorted_i)
pmt_id_branch_sorted.extend(pmt_id_branch_sorted_i)
event_number_branch.extend(event_number_branch_i)

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

with open(f'{output_path}filtered_file_MC_{run_number}_digidata_noTOF.pkl', 'wb') as f:
    pickle.dump(payload, f)

print(f"Saved single consolidated file:")
print(f" - {output_path}filtered_file_MC_noTOF.pkl (contains times_TOF, charge, mpmt_id, metadata)")