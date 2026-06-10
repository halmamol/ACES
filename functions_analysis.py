import numpy as np
import glob
import os
#import uproot
import pandas as pd
import json
from collections import defaultdict

def load_pmt_positions(json_path):
    """
    Load PMT positions from the WCTE geometry JSON.
    Returns dict: pmt_id → np.array([x, y, z]) in mm.
    """
    with open(json_path, "r") as f:
        geo_data = json.load(f)

    mpmt_data = geo_data.get("mpmts", {})
    pmt_positions = {}

    for mpmt_idx in sorted(mpmt_data.keys(), key=int):
        mpmt = mpmt_data[mpmt_idx]
        pmts = mpmt.get("pmts", {})
        for pmt_idx in sorted(pmts.keys(), key=int):
            pmt = pmts[pmt_idx]
            location = pmt["placement"]["location"]
            pmt_id = int(mpmt_idx) * 19 + int(pmt_idx)
            pmt_positions[pmt_id] = np.array(location, dtype=np.float64)

    return pmt_positions


def build_tof_map(pmt_positions, source_position, n_water=1.33):
    """
    Precompute ToF [ns] from a source position to every PMT.

    Parameters
    ----------
    pmt_positions    : dict {pmt_id: np.array([x,y,z])} in mm
    source_position  : np.array([x,y,z]) in mm
    n_water          : refractive index of water

    Returns
    -------
    tof_map : dict {pmt_id: tof_ns}
    """
    c_water_mm_ns = (3e8 / n_water) * 1e-6  # mm/ns
    source = np.array(source_position, dtype=np.float64)

    tof_map = {}
    for pmt_id, pos in pmt_positions.items():
        dist = np.linalg.norm(pos - source)
        tof_map[pmt_id] = dist / c_water_mm_ns

    return tof_map


def read_mpmt_offsets(json_path):
    """
    Read the mPMT ToF offset map (legacy format, keyed by 'slot*100+pos').
    Returns dict: (slot, pos) → tof_ns.
    """
    with open(json_path, "r") as f:
        raw = json.load(f)

    mpmt_map = {}
    for key, val in raw.items():
        k = int(key)
        slot = k // 100
        pos  = k % 100
        mpmt_map[(slot, pos)] = val

    return mpmt_map


def a_lista_de_arrays(plano, indices):
    return [plano[indices[i]:indices[i+1]] for i in range(len(indices)-1)]

def count_nHits(times_branch_event_arg, window, threshold_inf, threshold_sup):
    times_branch_event = np.sort(times_branch_event_arg.copy()) #just to make sure, but it is supposed to be sorted

    if times_branch_event[0]!=times_branch_event_arg[0]:
        print("hey Carla, the intial list was not sorted : WARNING")

    threshold_times = []
    nHits_Distribution = []

    i = 0
    n = len(times_branch_event)

    while i < n:
        time_hit = times_branch_event[i]
        mask = (times_branch_event >= time_hit) & (times_branch_event < time_hit + window)

        count = mask.sum()

        if count > threshold_inf and count<threshold_sup:
            threshold_times.append(time_hit)
            nHits_Distribution.append(count)

            i = np.searchsorted(times_branch_event, time_hit + window, side='left')
        else:
            i += 1
            
    return threshold_times, nHits_Distribution


def get_partition_and_local_event(event_number, bordes):
    for i in range(len(bordes)-1):
        if bordes[i] <= event_number < bordes[i+1]:
            return i, event_number - bordes[i]

    print(f"Warning: Event {event_number} not found in any partition")
    return -1, -1  # Si no encuentra


# Function to look up partition info
def get_partition_info(event_number, df):
    row = df[df['event_number'] == event_number]
    if not row.empty:
        partition = row['partition'].values[0]
        event_number_partition = row['event_number_partition'].values[0]
        return partition, event_number_partition
    else:
        return None, None  # or raise an error
    


def deltaT_calculation(csv_file):

    df = pd.read_csv(csv_file)

    # Asegurar tipos consistentes (por si event_number o start_time eran strings o enteros)
    df['event_number'] = df['event_number'].astype(int)
    df['prompt_time'] = df['prompt_time'].astype(float)

    # Reconstruir el diccionario anidado
    neutron_dict = defaultdict(lambda: defaultdict(list))

    for _, row in df.iterrows():
        event_number = row['event_number']
        start_time = row['prompt_time']
        delayed_time = row['delayed_time']
        neutron_dict[event_number][start_time].append(delayed_time)
    # Para neutron_dict
    neutron_dict = {k: dict(v) for k, v in neutron_dict.items()}
    
    deltaT = []
    for event_number in neutron_dict:
        for start_time in neutron_dict[event_number]:
            delayed_time = neutron_dict[event_number][start_time]
            deltaT.append(min(delayed_time) - start_time)

    return deltaT

def time_RMS_fun_time(times_event, t_in, window):

    mask = (times_event >= t_in) & (times_event < t_in + window)
    times_bin = times_event[mask]
    mean_t = times_bin.mean()
    RMS = np.sqrt(np.mean((times_bin - mean_t) ** 2))
    mean_t_rel = (times_bin - t_in).mean()
    
    return RMS, mean_t_rel


if __name__ == "__main__":

    """num_entries_list = []
    root_dir_bkg = "/data/cgarcia_2002/WCTE/data/2385_calib_time/"
    root_files_bkg = sorted(glob.glob(os.path.join(root_dir_bkg, "*.root")))
    root_files_bkg = sorted(root_files_bkg, key=lambda file_path: int(file_path.split("P")[-1].split(".")[0]))

    for file_path in root_files_bkg:
        print(f"Procesando archivo: {file_path}")
        file = uproot.open(file_path)
        tree = file["WCTEReadoutWindows"]
        num_entries_list.append(tree.num_entries)

    np.savetxt("Filtered_data/num_entries_list_sig.csv", num_entries_list, delimiter=",", fmt="%d")
"""


    # Load the CSV
    df = pd.read_csv('TestPartition.csv')

    event = 1924
    partition, local_event = get_partition_info(event)
    print(f"Global event {event} → Partition {partition}, Local event {local_event}")

