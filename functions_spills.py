import uproot
import awkward as ak
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import json
import sys
import os
sys.path.append(os.path.abspath("/mnt/netapp2/Store_uni/home/usc/ie/dcr/software/hk"))

# import ACES.functions_bonsai


def spill_nHitsTT(times_branch_event_arg, threshold_inf, window, death_window, charge_branch_event = [], threshold_sup = np.inf):
    times_branch_event = np.sort(times_branch_event_arg.copy()) #just to make sure, but it is supposed to be sorted

    threshold_times = []
    indices_to_delete = []
    nhits_range = []

    i = 0
    n = len(times_branch_event)

    while i < n:
        time_hit = times_branch_event[i]
        mask = (times_branch_event >= time_hit) & (times_branch_event < time_hit + window)

        if len(charge_branch_event)!=0:
            count = max(charge_branch_event[mask])
        else:
            count = mask.sum()

        if count > threshold_inf and count<threshold_sup:
            threshold_times.append(time_hit)
            nhits_range.append(count)
            # Zero out the next death window ns after the hit window
            mask_2 = (times_branch_event >= time_hit) & (times_branch_event < time_hit + window + death_window)

            indices = np.where(mask_2)[0]
            indices_to_delete.extend(indices.tolist())

            # Jump ahead by death window
            i += 1
            while i < n and times_branch_event[i] < time_hit + window + death_window:
                i += 1
        else:
            i += 1

    return threshold_times, indices_to_delete, nhits_range


def repeat_spills_nHits(event_number_branch, times_branch_sorted, threshold, window, death_window):
    threshold_times = {}
    times_branch_modified = []
    deleted_indices_by_event = {}

    for event in event_number_branch:
        if event%1000 == 0:
            print(f"Filtering nHits event {event}...")
        remaining_times = times_branch_sorted[event].copy()
        all_thresholds = []
        all_deleted_indices = []
        pass_counter = 0

        while True:
            new_thresholds, indices_to_delete, _ = spill_nHitsTT(remaining_times, threshold, window, death_window)

            if not new_thresholds:
                break

            # Save deleted indices relative to the original event array
            deleted_indices = np.where(np.isin(times_branch_sorted[event], remaining_times[indices_to_delete]))[0]
            all_deleted_indices.extend(deleted_indices.tolist())

            remaining_times = np.delete(remaining_times, indices_to_delete)
            all_thresholds.extend(new_thresholds)
            pass_counter += 1

        if pass_counter > 1:
            print(f"Event {event}: spill_nHitsTT applied {pass_counter} times")

        times_branch_modified.append(remaining_times)

        if all_thresholds:
            threshold_times[event] = all_thresholds

        if all_deleted_indices:
            deleted_indices_by_event[event] = all_deleted_indices

    return times_branch_modified, threshold_times, deleted_indices_by_event

def repeat_spills_nHits_with_channels(
    event_number_branch,
    times_branch_sorted_TOF,
    charge_branch_sorted,
    mpmt_id_branch_sorted,
    pmt_id_branch_sorted,
    threshold,
    window,
    death_window,
):
    """
    Apply spill_nHitsTT repeatedly per event and delete hits across times (TOF), charge,
    and mpmt_id arrays in lockstep using the same indices.

    Returns:
        times_branch_modified_TOF: list[np.ndarray]
        charge_branch_modified: list[np.ndarray]
        mpmt_id_branch_modified: list[np.ndarray]
        threshold_times: dict[event_id -> list[thresholds]]
        deleted_indices_by_event: dict[event_id -> list[int]] (indices in original event arrays)
    """
    import numpy as np

    times_branch_modified_TOF = []
    charge_branch_modified = []
    mpmt_id_branch_modified = []
    pmt_id_branch_modified = []  
    threshold_times = {}
    deleted_indices_by_event = {}

    for event in event_number_branch:
        if event % 1000 == 0:
            print(f"Filtering nHits event {event}...")

        t_orig = times_branch_sorted_TOF[event]
        c_orig = charge_branch_sorted[event]
        m_orig = mpmt_id_branch_sorted[event]
        p_orig = pmt_id_branch_sorted[event]

        # Ensure arrays are aligned
        if not (len(t_orig) == len(c_orig) == len(m_orig) == len(p_orig)):
            raise ValueError(
                f"Event {event}: arrays have mismatched lengths "
                f"(times={len(t_orig)}, charge={len(c_orig)}, mpmt_id={len(m_orig)}, pmt_id={len(p_orig)})"
            )

        remaining_times = t_orig.copy()
        # idx_map maps indices in 'remaining_times' to indices in the original arrays
        idx_map = np.arange(len(t_orig))

        all_thresholds = []
        all_deleted_orig_indices = []
        pass_counter = 0

        while True:
            new_thresholds, indices_to_delete, _ = spill_nHitsTT(
                remaining_times, threshold, window, death_window
            )

            if not new_thresholds:
                break

            # Map deletions to original indices BEFORE deleting
            orig_indices = idx_map[indices_to_delete]
            all_deleted_orig_indices.extend(orig_indices.tolist())

            # Delete from the working arrays
            remaining_times = np.delete(remaining_times, indices_to_delete)
            idx_map = np.delete(idx_map, indices_to_delete)

            all_thresholds.extend(new_thresholds)
            pass_counter += 1

        if pass_counter > 1:
            print(f"Event {event}: spill_nHitsTT applied {pass_counter} times")

        # Build keep mask for original arrays
        if all_deleted_orig_indices:
            keep_mask = np.ones(len(t_orig), dtype=bool)
            keep_mask[all_deleted_orig_indices] = False
        else:
            keep_mask = np.ones(len(t_orig), dtype=bool)

        # Append filtered arrays
        times_branch_modified_TOF.append(t_orig[keep_mask])
        charge_branch_modified.append(c_orig[keep_mask])
        mpmt_id_branch_modified.append(m_orig[keep_mask])
        pmt_id_branch_modified.append(p_orig[keep_mask])  # Not returned but could be if needed

        # Metadata
        if all_thresholds:
            threshold_times[event] = all_thresholds
        if all_deleted_orig_indices:
            deleted_indices_by_event[event] = all_deleted_orig_indices

    return (
        times_branch_modified_TOF,
        charge_branch_modified,
        mpmt_id_branch_modified,
        pmt_id_branch_modified,
        threshold_times,
        deleted_indices_by_event,
    )

def repeat_spills_Charge(event_number_branch, times_branch_sorted, charge_branch_sorted, window, death_window, threshold = 5000):
    threshold_charges = {}
    times_branch_modified_chargesTT = []
    charge_branch_modified_chargesTT = []
    deleted_indices_by_event = {}

    for event in event_number_branch:
        if event % 1000 == 0:
            print(f"Filtering charge on event {event}...")
        remaining_times = times_branch_sorted[event].copy()
        remaining_charges = charge_branch_sorted[event].copy()
        all_thresholds = []
        all_deleted_indices = []
        pass_counter = 0

        while True:
            new_thresholds, indices_to_delete, _ = spill_nHitsTT(remaining_times, threshold, window, death_window, remaining_charges)

            if not new_thresholds:
                break

            # Save deleted indices relative to the original event array
            deleted_indices = np.where(np.isin(times_branch_sorted[event], remaining_times[indices_to_delete]))[0]
            all_deleted_indices.extend(deleted_indices.tolist())

            remaining_times = np.delete(remaining_times, indices_to_delete)
            remaining_charges = np.delete(remaining_charges, indices_to_delete)

            all_thresholds.extend(new_thresholds)
            pass_counter += 1

        if pass_counter > 1:
            print(f"Event {event}: spill_nHitsTT applied {pass_counter} times")

        times_branch_modified_chargesTT.append(remaining_times)
        charge_branch_modified_chargesTT.append(remaining_charges)


        if all_thresholds:
            threshold_charges[event] = all_thresholds

        if all_deleted_indices:
            deleted_indices_by_event[event] = all_deleted_indices


    return times_branch_modified_chargesTT, charge_branch_modified_chargesTT, threshold_charges, deleted_indices_by_event

def plot_TotalCharge_mPMT(mpmt_id, charge):

    charge_per_mpmt = {}

    for m, c in zip(mpmt_id, charge):
        if m in charge_per_mpmt:
            charge_per_mpmt[m] += c
        else:
            charge_per_mpmt[m] = c

    # Step 2: Prepare data for plotting
    mpmt_ids = sorted(charge_per_mpmt.keys())
    charges = [charge_per_mpmt[m] for m in mpmt_ids]

    # Step 3: Plot
    plt.figure(figsize=(8, 4))
    plt.bar(mpmt_ids, charges, width = 1, color='blue',  align='edge', edgecolor='navy')
    plt.xlim(0, 150)
    plt.xlabel("mPMT ID")
    plt.ylabel("Total Charge [u.a]")
    plt.title("Charge per mPMT in 50 ns window")
    plt.tight_layout()
    plt.show()

def plot_TotalCharge_Time(time, charge, bin_time, total_time):

    sum_charges = np.zeros(int(total_time/bin_time))
    div = (time-min(time))//bin_time

    for i, n in enumerate(div):
        sum_charges[int(n)] += charge[i]

    return sum_charges


def time_RMS_fun(times_event, window):

    RMS_list = []
    n = len(times_event)
    i = 0
    
    while i <n:
        t_in = times_event[i]
        mask = (times_event >= t_in) & (times_event < t_in + window)
        times_bin = times_event[mask]

        N = len(times_bin)
        if N == 0 or N==1:  #salen t_RMS bajos pero es mentira, no tiene ningun significado la rms de un unico valor, los evito
            i+=1  
            continue

        i = np.searchsorted(times_event, t_in + window, side='left')

        mean_t = times_bin.mean()
        RMS = np.sqrt(np.mean((times_bin - mean_t) ** 2))
        RMS_list.append(RMS)

    return RMS_list


def read_mpmt_offsets(path):
    with open(path) as f:
        mcc_map = json.load(f)

        d = {}
        for k,v in zip(mcc_map.keys(), mcc_map.values()):
            card, channel = [int(i) for i in str(int(k)/100).split(".")]
            d[(card, channel)] = v

    return d

def correction_TOF(mpmt_map, mpmt_slot_branch, pmt_position, max_slot = 106, max_pos = 19):
    lookup = np.zeros((max_slot, max_pos))
    for (card, chan), shift in mpmt_map.items():
        lookup[card, chan] = shift

        # Hit Times Correction
    flat_slot_ids     = ak.ravel(mpmt_slot_branch)
    flat_pos_ids      = ak.ravel(pmt_position)
    flat_corrections  = lookup[flat_slot_ids, flat_pos_ids]
    corrections       = ak.unflatten(flat_corrections, ak.num(mpmt_slot_branch))

    return corrections

def initial_treatment(tree):

    times_branch = tree["hit_pmt_calibrated_times"].array()
    charge_branch = tree["hit_pmt_charges"].array()
    mpmt_slot_branch = tree["hit_mpmt_slot_ids"].array()
    pmt_position = tree["hit_pmt_position_ids"].array()
    event_number_branch = tree["event_number"].array(library="np")
    hit_cards_id = tree["hit_mpmt_card_ids"].array()
    has_time_cte = tree["hit_pmt_has_time_constant"].array()

    # Crear una máscara booleana por evento donde ambos valores sean >= 0
    valid_mask = (mpmt_slot_branch >= 0) & (pmt_position >= 0) & (charge_branch < 1e4) & (hit_cards_id < 120) &  (has_time_cte != 0)

    # Aplicar la máscara a cada rama para eliminar los hits inválidos
    times_branch_clean = times_branch[valid_mask]
    charge_branch_clean = charge_branch[valid_mask]
    mpmt_slot_branch_clean = mpmt_slot_branch[valid_mask]
    pmt_position_clean = pmt_position[valid_mask]

    mpmt_map = read_mpmt_offsets("/mnt/lustre/scratch/nlsas/home/usc/ie/dcr/hk/mpmt_tof_map/mpmt_tof_pos1525.json")
    corrections = correction_TOF(mpmt_map, mpmt_slot_branch_clean, pmt_position_clean)
    corrected_times = times_branch_clean - corrections
    
    # Ordenar los tiempos por evento y obtener los índices
    sorted_idx = ak.argsort(corrected_times, axis=1)

    # Usar los índices para reordenar todos los branches
    times_sorted_TOF = corrected_times[sorted_idx]
    times_sorted = times_branch_clean[sorted_idx]
    charges_sorted = charge_branch_clean[sorted_idx]
    mpmt_sorted = mpmt_slot_branch_clean[sorted_idx]
    pmt_position_sorted = pmt_position_clean[sorted_idx]

    # Convertir a listas de NumPy arrays
    times_sorted_np = [np.array(evt) for evt in times_sorted]
    times_sorted_TOF_np = [np.array(evt) for evt in times_sorted_TOF]
    charges_sorted_np = [np.array(evt) for evt in charges_sorted]
    mpmt_sorted_np = [np.array(evt) for evt in mpmt_sorted]
    pmt_position_sorted_np = [np.array(evt) for evt in pmt_position_sorted]

    return times_sorted_np, times_sorted_TOF_np, charges_sorted_np, mpmt_sorted_np, pmt_position_sorted_np, event_number_branch

def delete_indices_list(list_to_delete, indices):
    """List to delete es la lista de listas incial
    Indices es un diccionario donde la clave es el evento y el valor los indices a eliminar"""

    new_list = []

    for event, values in enumerate(list_to_delete):
        if event in indices:
            indices_to_delete = set(indices[event])  # convert to set for faster lookup

            # Keep elements whose index is NOT in indices_to_delete
            filtered_list = [t for i, t in enumerate(values) if i not in indices_to_delete]
        else:
            filtered_list = values.copy()

        new_list.append(np.array(filtered_list))

    return new_list

def counting_nHits_window(event_number_branch, times_branch, bin_window):
    nHits = []

    for event_number in event_number_branch:
        times_branch_event = times_branch[event_number]
        n = len(times_branch_event)
        i=0

        while i<n:

            t_in = times_branch_event[i]
            mask = (times_branch_event >= t_in) & (times_branch_event < t_in + bin_window)
            count = mask.sum()
            nHits.append(count)
            i = np.searchsorted(times_branch_event, t_in + bin_window, side='left')

    return nHits

def multiple_partition(root_files):

    times_branch_sorted = []
    times_branch_sorted_TOF = []
    charge_branch_sorted = []
    mpmt_id_branch_sorted = []
    pmt_id_branch_sorted = []
    event_number_branch = []

    dir_event_partition = {}
    # Contador global de eventos
    event_offset = 0

    for file_path in root_files:
        print(f"Procesando archivo: {file_path}")
        file = uproot.open(file_path)
        tree = file["WCTEReadoutWindows"]

        times_branch_sorted_i, times_branch_sorted_TOF_i, charge_branch_sorted_i, mpmt_id_branch_sorted_i, pmt_id_branch_sorted_i, event_number_branch_i = initial_treatment(tree)
        
        dir_event_partition[int(file_path.split("P")[-1].split(".")[0])] = len(event_number_branch_i)
        # Ajustar los event_numbers con offset para que no se repitan
        new_event_numbers = [i + event_offset for i in event_number_branch_i]

        times_branch_sorted.extend(times_branch_sorted_i)
        times_branch_sorted_TOF.extend(times_branch_sorted_TOF_i)
        charge_branch_sorted.extend(charge_branch_sorted_i)
        mpmt_id_branch_sorted.extend(mpmt_id_branch_sorted_i)
        pmt_id_branch_sorted.extend(pmt_id_branch_sorted_i)
        event_number_branch.extend(new_event_numbers)

        # Actualizar offset para el siguiente archivo
        event_offset += tree.num_entries

    return times_branch_sorted, times_branch_sorted_TOF, charge_branch_sorted, mpmt_id_branch_sorted, pmt_id_branch_sorted, event_number_branch, dir_event_partition
