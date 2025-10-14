import numpy as np
import pandas as pd  
import json

import functions_bonsai
"""""""""""""""""
def nHitsTimeWindow(times_branch_event_arg, threshold_inf, window, death_window, charge_branch_event = [], threshold_sup = np.inf):
    times_branch_event = np.sort(np.array(times_branch_event_arg).copy()) #just to make sure, but it is supposed to be sorted

    threshold_times = []
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

            # Jump ahead by death window
            i += 1
            while i < n and times_branch_event[i] < time_hit + window + death_window:
                i += 1
        else:
            i += 1

    return threshold_times, nhits_range
"""""""""""""""""

def nHitsTimeWindow(times_branch_event_arg, threshold_inf, window, death_window, charge_branch_event=[], threshold_sup=np.inf):
    times = np.sort(np.array(times_branch_event_arg).copy())
    threshold_times = []
    nhits_range = []

    n = len(times)
    start = 0
    end = 0
    i = 0

    max_t = max(times)

    while i < n:
        t0 = times[i]
        
        if t0 % 1000 == 0:
                print(f"Current time stamp {t0}")
                print(f"Percentaje of total time evaluated: {t0*100/max_t}%")
                print(f"Number of nHits found {len(threshold_times)}")
        # Move 'end' to include all times within [t0, t0+window)
        while end < n and times[end] < t0 + window:
            end += 1

        # Compute count/charge in window
        if len(charge_branch_event) != 0:
            print("Charge")
            max_charge = np.max(charge_branch_event[i:end])
            count = max_charge
        else:
            count = end - i
        
        if count > threshold_inf and count < threshold_sup:
            #print(f"Time stamp {t0}")
            #print(f"Hits Count {count}")
            threshold_times.append(t0)
            nhits_range.append(count)
            # Skip forward by death window
            # Find first index after t0+window+death_window
            skip_to = t0 + window + death_window
            while i < n and times[i] < skip_to:
                i += 1
            # Move end as well (since i advanced)
            end = max(end, i)
        else:
            i += 1

    return threshold_times, nhits_range


def prompt_candidates(event_branch, times_branch_arg, window_sliding, window_clean, threshold_inf, threshold_sup):

    def prompt_candidates_event(times_branch_event_arg, window_sliding, window_clean, threshold_inf, threshold_sup):
        valid_thresholds= []
        threshold_list, n_hits = nHitsTimeWindow(times_branch_event_arg, threshold_inf, window_sliding, 0, threshold_sup = threshold_sup)

        for time_prompt in threshold_list:

            mask_clean_1 = (times_branch_event_arg >= time_prompt - window_clean) & (times_branch_event_arg < time_prompt)
            mask_clean_2 = (times_branch_event_arg > time_prompt + window_sliding) & (times_branch_event_arg < time_prompt + window_sliding + window_clean)
            
            if mask_clean_1.sum() != 0 :
                clean_list, _ = nHitsTimeWindow(times_branch_event_arg[mask_clean_1], threshold_inf, window_sliding, 0, threshold_sup=threshold_sup)
            else:
                clean_list = []
            
            if mask_clean_2.sum() != 0 :
                clean_list_2, _ = nHitsTimeWindow(times_branch_event_arg[mask_clean_2], threshold_inf, window_sliding, 0, threshold_sup=threshold_sup)
            else:
                clean_list_2 = []

            if len(clean_list) + len(clean_list_2) == 0:
                #valid_thresholds.append(time_prompt)
                mask_hits = (times_branch_event_arg >= time_prompt) & (times_branch_event_arg < time_prompt + window_sliding)
                #n_hits = mask_hits.sum()
                valid_thresholds.append((time_prompt, n_hits))
            #else:
            #    print(f"Not using trigger {time_prompt} because it has other signal_candidate too close in event {event}") 

        return valid_thresholds
    
    
    theshold_times = {}
    
    for event in event_branch:
        if event%1000 == 0:
            print(f"Searching prompt candidates on event {event}...")
        threshold_times_event = prompt_candidates_event(times_branch_arg[event], window_sliding, window_clean, threshold_inf, threshold_sup)
        
        if len(threshold_times_event)!=0:
            theshold_times[event] = threshold_times_event

    return theshold_times

def prompt_candidates_wBonsai(
    times_event: np.ndarray,
    charge_event: np.ndarray,
    mpmt_event: np.ndarray,
    pmt_event: np.ndarray,
    window_sliding: float,
    window_clean: float,
    threshold_inf: float,
    threshold_sup: float,
):

   # Sanity: arrays must be aligned
    if not (len(times_event) == len(charge_event) == len(mpmt_event) == len(pmt_event)):
        raise ValueError(
            f"Arrays have mismatched lengths (times={len(times_event)}, "
            f"charge={len(charge_event)}, mpmt={len(mpmt_event)}, pmt={len(pmt_event)})"
        )

    valid_thresholds = []
    print("Searching prompt nHits threshold")
    # Find candidate trigger times
    threshold_list,  _ = nHitsTimeWindow(
        times_event,
        threshold_inf,
        window_sliding,
        0,
        threshold_sup=threshold_sup,
    )
    print(len(threshold_list))
    print("Cleaning candidates on windows before and after")
    for time_prompt in threshold_list:
        # Clean windows before and after the candidate
        mask_clean_1 = (times_event >= time_prompt - window_clean) & (times_event < time_prompt)
        mask_clean_2 = (times_event > time_prompt + window_sliding) & (
            times_event < time_prompt + window_sliding + window_clean
        )

        if mask_clean_1.sum() != 0:
            clean_list_1, _ = nHitsTimeWindow(
                times_event[mask_clean_1], threshold_inf, window_sliding, 0, threshold_sup=threshold_sup
            )
        else:
            clean_list_1 = []

        if mask_clean_2.sum() != 0:
            clean_list_2, _ = nHitsTimeWindow(
                times_event[mask_clean_2], threshold_inf, window_sliding, 0, threshold_sup=threshold_sup
            )
        else:
            clean_list_2 = []

        #print("Evaluating candidate at time")
        # Accept only if no other signal candidate too close
        if len(clean_list_1) + len(clean_list_2) == 0:
            mask_hits = (times_event >= time_prompt) & (times_event < time_prompt + window_sliding)
            n_hits_here = int(mask_hits.sum())
            if n_hits_here == 0:
                continue
            
            times_in_prompt = times_event[mask_hits]
            charges_in_prompt = charge_event[mask_hits]
            mpmt_in_prompt = mpmt_event[mask_hits]
            pmt_in_prompt = pmt_event[mask_hits]

            vertex = functions_bonsai.run_BONSAI_candidate(times_in_prompt, charges_in_prompt, mpmt_in_prompt, pmt_in_prompt)

            valid_thresholds.append(
                {
                    "time": float(time_prompt),
                    "n_hits": n_hits_here,
                    "charge": charges_in_prompt,
                    "mpmt_id": mpmt_in_prompt,
                    "pmt_id": pmt_in_prompt,
                    "vertex_x": vertex["x"][0],
                    "vertex_y": vertex["y"][0],
                    "vertex_z": vertex["z"][0],
                    # If useful, you can also include:
                    # "indices": np.where(mask_hits)[0],
                    # "times": times_event[mask_hits],
                }
            )
            if len(valid_thresholds) % 10 == 0:
                print(f"Found Prompt Events {len(valid_thresholds)}...")

    return valid_thresholds


def neutron_detection(event_branch, times_branch_event_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup = np.inf, window_prompt = 100):

    def neutron_detection_event(times_branch_event_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt):
        
        dict_neutrons_event = {}
        last_prompt = None

        #for time_prompt in threshold_times:
        for time_prompt, _, _ in threshold_times:    
            if last_prompt is not None and (time_prompt-last_prompt) < (window_sliding + window_prompt):
                continue

            mask = (times_branch_event_arg >= time_prompt + window_prompt) & (times_branch_event_arg < time_prompt + window_prompt + window_sliding)
            #neutron_nhits = mask.sum()
            if mask.sum() == 0:
                continue
            
            neutron_candidates, neutron_nhits = nHitsTimeWindow(times_branch_event_arg[mask], threshold_inf, window_neutron, 0, threshold_sup=threshold_sup)
            if len(neutron_candidates) != 0:
                # Store as list of tuples (neutron_time, neutron_nhits)

                dict_neutrons_event[time_prompt] = list(zip(neutron_candidates, neutron_nhits))#[(nt, neutron_nhits) for nt in neutron_candidates]
                last_prompt = time_prompt

        return dict_neutrons_event
    
    dict_neutrons = {}
    for event in event_branch:
        if event % 1000 == 0:
            print(f"Searching neutron coincidences on event {event}...")
        if event in threshold_times:
            if len(threshold_times[event]) == 0:
                continue 
            dict_neutrons_event = neutron_detection_event(times_branch_event_arg[event], threshold_times[event], window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt)
            if dict_neutrons_event:
                dict_neutrons[event] = dict_neutrons_event
        #else:
            #print(f"Event {event} has no threshold times, skipping neutron detection.")

    return dict_neutrons

def neutron_detection_wBonsai(times_branch_event_arg, charge_branch_event_arg, mpmt_branch_event_arg, pmt_branch_event_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup = np.inf, window_prompt = 100):

    valid_thresholds = []
    last_prompt = None

    #for time_prompt in threshold_times:
    for time_prompt, nhits_prompt, trms_prompt,_ , _, _, x_prompt, y_prompt, z_prompt  in threshold_times:    
        if last_prompt is not None and (time_prompt-last_prompt) < (window_sliding + window_prompt):
            continue

        mask = (times_branch_event_arg >= time_prompt + window_prompt) & (times_branch_event_arg < time_prompt + window_prompt + window_sliding)
        
        if mask.sum() == 0:
            continue

        neutron_candidates,  neutron_nhits = nHitsTimeWindow(times_branch_event_arg[mask], threshold_inf, window_neutron, 0, threshold_sup=threshold_sup)
        
        #for time_delayed in neutron_candidates:
        for i in range(len(neutron_candidates)):
            time_delayed = neutron_candidates[i]
            delayed_nhits = neutron_nhits[i]
            mask_delayed = (times_branch_event_arg >= time_delayed) & (times_branch_event_arg < time_delayed + window_sliding)
            #neutron_nhits = mask_delayed.sum()

            times_in_delayed = times_branch_event_arg[mask_delayed]
            charges_in_delayed = charge_branch_event_arg[mask_delayed]
            mpmt_in_delayed = mpmt_branch_event_arg[mask_delayed]
            pmt_in_delayed = pmt_branch_event_arg[mask_delayed]
            
            vertex = functions_bonsai.run_BONSAI_candidate(times_in_delayed, charges_in_delayed, mpmt_in_delayed, pmt_in_delayed)

            valid_thresholds.append(
                {
                    'prompt_nhits': nhits_prompt,
                    'prompt_time': float(time_prompt),
                    'prompt_trms': trms_prompt,
                    'prompt_x': x_prompt,
                    'prompt_y': y_prompt,
                    'prompt_z': z_prompt,
                    "delayed_time": float(time_delayed),
                    "delayed_nhits": delayed_nhits,
                    "delayed_x": vertex["x"][0],
                    "delayed_y": vertex["y"][0],
                    "delayed_z": vertex["z"][0],
                }
            )

    return valid_thresholds
    
        #else:
            #print(f"Event {event} has no threshold times, skipping neutron detection.")
