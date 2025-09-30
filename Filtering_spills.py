print("Starting Spills Filter Algorithm")

import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import functions_spills
import glob
import os
import argparse

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
root_dir = f"/scratch/halmazan/WCTE/files/data/{run_number}/"
root_file_path = f"{root_dir}WCTE_offline_R{run_number}S0P{partition}.root"
output_path = f"//scratch/halmazan/WCTE/files/data/filtered_files/"

# Showing results
if partition == "all":
    print("Analysing all partitions.")

    print(f"Loading Run {run_number}...")
    root_files = sorted(glob.glob(os.path.join(root_dir, "*.root")))
    root_files = sorted(root_files, key=lambda file_path: int(file_path.split("P")[-1].split(".")[0]))
    root_files = root_files[:-1]

    print(f"Found {len(root_files)} ROOT files.")

    times_branch_sorted, times_branch_sorted_TOF, charge_branch_sorted, mpmt_id_branch_sorted, event_number_branch, _ = functions_spills.multiple_partition(root_files)

    print("Runs loaded.")
    N_events = max(event_number_branch) + 1

print(f"Total number of events in run {run_number}: {N_events}")


# Filter spills using nHits threshold ###################################################################################

print(f"Applying filter for Run {run_number}...")
times_branch_modified, threshold_times, deleted_index_dict = functions_spills.repeat_spills_nHits(event_number_branch, times_branch_sorted_TOF, nhits_threshold, nhits_window, death_time)

print("nHits filter applied.")

#Filter spills using charge threshold #############################################################################################

# Save filtered data to pickle files #############################################################################################

import pickle

# Guardar

with open(f'{output_path}deleted_indices_nHits_{run_number}.pkl', 'wb') as f:
    pickle.dump(deleted_index_dict, f)


total_elements = sum(len(v) for v in threshold_times.values())
print(f"Total number of elements in threshold_times for run {run_number}:", total_elements, max(event_number_branch), total_elements/max(event_number_branch))


def a_array_plano_y_indices(lista_arrays):
    plano = np.concatenate(lista_arrays)
    indices = np.cumsum([0] + [len(arr) for arr in lista_arrays])
    return plano, indices

valores, indices = a_array_plano_y_indices(times_branch_modified)

with open(f'{output_path}filtered_file_{run_number}.pkl', 'wb') as f:
    pickle.dump((valores, indices), f)

