import numpy as np
import pandas as pd  
import json

import functions_bonsai
import functions_multilateration
import functions_analysis


def nHitsTimeWindow(times_branch_event_arg, threshold_inf, window, death_window, threshold_sup = np.inf):
    times_branch_event = np.sort(times_branch_event_arg.copy()) #just to make sure, but it is supposed to be sorted

    threshold_times = []
    nhits_range = []

    i = 0
    n = len(times_branch_event)

    while i < n:
        time_hit = times_branch_event[i]
        mask = (times_branch_event >= time_hit) & (times_branch_event < time_hit + window)

        #if len(charge_branch_event)!=0:
        #    count = max(charge_branch_event[mask])
        #else:
        #    count = mask.sum()
        count = int(mask.sum())

        if threshold_inf <= count < threshold_sup:
        #if count > threshold_inf and count<threshold_sup:
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

def nHitstRMSTimeWindow(
    times_branch_event_arg,
    threshold_inf,
    window,
    death_window,
    time_rms_fun,             # e.g. functions_analysis.time_RMS_fun_time
    rms_cut_ns=10.0,
    threshold_sup=np.inf,     # reject if count >= threshold_sup
):
    """
    Sliding window finder:
      - find windows with count >= threshold_inf  (NO upper limit in search)
      - compute t_rms for that window
      - reject if t_rms > rms_cut_ns
      - reject if count >= threshold_sup
      - return surviving candidate start times + nhits + trms
    """
    times = np.sort(np.asarray(times_branch_event_arg, float))
    n = len(times)

    cand_times = []
    cand_nhits = []
    cand_trms  = []

    i = 0
    while i < n:
        t_start = times[i]
        in_win = (times >= t_start) & (times < t_start + window)
        count = int(in_win.sum())

        if count >= threshold_inf:
            # compute rms as soon as it passes the lower threshold
            t_rms, _ = time_rms_fun(times[in_win], t_start, window)
            t_rms = float(t_rms)

            # apply your rejection cuts
            if (t_rms <= rms_cut_ns) and (count < threshold_sup):
                cand_times.append(float(t_start))
                cand_nhits.append(count)
                cand_trms.append(t_rms)

            # lockout region
            t_skip_until = t_start + window + death_window
            i += 1
            while i < n and times[i] < t_skip_until:
                i += 1
        else:
            i += 1

    return cand_times, cand_nhits, cand_trms


def prompt_candidates(#event_branch, times_branch_arg, window_sliding, window_clean, threshold_inf, threshold_sup):
    event_branch,
    times_branch_arg,
    charge_branch_arg,
    mpmt_branch_arg,
    pmt_branch_arg,
    window_sliding,
    window_clean,
    threshold_inf,
    threshold_sup,
):
    """
    For each event, find "prompt" candidates using nHitsTimeWindow, enforce clean
    windows before/after, and return (per event) a list of prompt entries that
    include the trigger time, number of hits in the prompt window, and the charge
    and mPMT-id arrays for the hits inside that window.

    Returns:
        dict:
            {
              event_id: [
                {
                  "time": float,            # trigger time (prompt)
                  "n_hits": int,            # number of hits in [time, time+window_sliding)
                  "charge": np.ndarray,     # charges of hits in the prompt window
                  "mpmt_id": np.ndarray,    # mPMT ids of hits in the prompt window
                  # Optional: you may also add "indices" or "times" if needed
                },
                ...
              ],
              ...
            }
    """
    import numpy as np

    def prompt_candidates_event(
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

        # Find candidate trigger times
        threshold_list,  _ = nHitsTimeWindow(
            times_event,
            threshold_inf,
            window_sliding,
            0,
            threshold_sup=threshold_sup,
        )

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

                vertex = (0, 0, 0)

                valid_thresholds.append(
                    {
                        "time": float(time_prompt),
                        "n_hits": n_hits_here,
                        "charge": charges_in_prompt,
                        "mpmt_id": mpmt_in_prompt,
                        "pmt_id": pmt_in_prompt,
                        "vertex_x": vertex[0],
                        "vertex_y": vertex[1],
                        "vertex_z": vertex[2],
                        # If useful, you can also include:
                        # "indices": np.where(mask_hits)[0],
                        # "times": times_event[mask_hits],
                    }
                )

        return valid_thresholds

    threshold_times = {}

    for event in event_branch:
        if event % 1000 == 0:
            print(f"Searching prompt candidates on event {event}...")

        times_event = times_branch_arg[event]
        charge_event = charge_branch_arg[event]
        mpmt_event = mpmt_branch_arg[event]
        pmt_event = pmt_branch_arg[event]

        results = prompt_candidates_event(
            times_event,
            charge_event,
            mpmt_event,
            pmt_event,
            window_sliding,
            window_clean,
            threshold_inf,
            threshold_sup,
        )

        if len(results) != 0:
            threshold_times[event] = results

    return threshold_times

def prompt_candidates_wBonsai(
    event_branch,
    times_branch_arg,
    charge_branch_arg,
    mpmt_branch_arg,
    pmt_branch_arg,
    window_sliding,
    window_clean,
    threshold_inf,
    threshold_sup,
):
    """
    For each event, find "prompt" candidates using nHitsTimeWindow, enforce clean
    windows before/after, and return (per event) a list of prompt entries that
    include the trigger time, number of hits in the prompt window, and the charge
    and mPMT-id arrays for the hits inside that window.

    Returns:
        dict:
            {
              event_id: [
                {
                  "time": float,            # trigger time (prompt)
                  "n_hits": int,            # number of hits in [time, time+window_sliding)
                  "charge": np.ndarray,     # charges of hits in the prompt window
                  "mpmt_id": np.ndarray,    # mPMT ids of hits in the prompt window
                  # Optional: you may also add "indices" or "times" if needed
                },
                ...
              ],
              ...
            }
    """
    import numpy as np

    def prompt_candidates_event(
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

        # Find candidate trigger times
        threshold_list,  _ = nHitsTimeWindow(
            times_event,
            threshold_inf,
            window_sliding,
            0,
            threshold_sup=threshold_sup,
        )

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

        return valid_thresholds

    threshold_times = {}

    for event in event_branch:
        if event % 1000 == 0:
            print(f"Searching prompt candidates on event {event}...")

        times_event = times_branch_arg[event]
        charge_event = charge_branch_arg[event]
        mpmt_event = mpmt_branch_arg[event]
        pmt_event = pmt_branch_arg[event]

        results = prompt_candidates_event(
            times_event,
            charge_event,
            mpmt_event,
            pmt_event,
            window_sliding,
            window_clean,
            threshold_inf,
            threshold_sup,
        )

        if len(results) != 0:
            threshold_times[event] = results

    return threshold_times

def prompt_candidates_wScintillationLikelihood(
    event_branch,
    times_branch_arg,
    charge_branch_arg,
    mpmt_branch_arg,
    pmt_branch_arg,
    window_sliding,
    window_clean,
    threshold_inf,
    threshold_sup,
):
    """
    For each event, find "prompt" candidates using nHitsTimeWindow, enforce clean
    windows before/after, and return (per event) a list of prompt entries that
    include the trigger time, number of hits in the prompt window, and the charge
    and mPMT-id arrays for the hits inside that window.

    Returns:
        dict:
            {
              event_id: [
                {
                  "time": float,            # trigger time (prompt)
                  "n_hits": int,            # number of hits in [time, time+window_sliding)
                  "charge": np.ndarray,     # charges of hits in the prompt window
                  "mpmt_id": np.ndarray,    # mPMT ids of hits in the prompt window
                  # Optional: you may also add "indices" or "times" if needed
                },
                ...
              ],
              ...
            }
    """
    import numpy as np

    def prompt_candidates_event(
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

        # Find candidate trigger times
        threshold_list,  _ = nHitsTimeWindow(
            times_event,
            threshold_inf,
            window_sliding,
            0,
            threshold_sup=threshold_sup,
        )

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

                vertex = functions_bonsai.scintillation_candidates(times_in_prompt, charges_in_prompt, mpmt_in_prompt, pmt_in_prompt)

                valid_thresholds.append(
                    {
                        "time": float(time_prompt),
                        "n_hits": n_hits_here,
                        "charge": charges_in_prompt,
                        "mpmt_id": mpmt_in_prompt,
                        "pmt_id": pmt_in_prompt,
                        "vertex_x": vertex[0],
                        "vertex_y": vertex[1],
                        "vertex_z": vertex[2],
                        # If useful, you can also include:
                        # "indices": np.where(mask_hits)[0],
                        # "times": times_event[mask_hits],
                    }
                )

        return valid_thresholds

    threshold_times = {}

    for event in event_branch:
        if event % 1000 == 0:
            print(f"Searching prompt candidates on event {event}...")

        times_event = times_branch_arg[event]
        charge_event = charge_branch_arg[event]
        mpmt_event = mpmt_branch_arg[event]
        pmt_event = pmt_branch_arg[event]

        results = prompt_candidates_event(
            times_event,
            charge_event,
            mpmt_event,
            pmt_event,
            window_sliding,
            window_clean,
            threshold_inf,
            threshold_sup,
        )

        if len(results) != 0:
            threshold_times[event] = results

    return threshold_times

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

def neutron_detection_wBonsai(event_branch, times_branch_arg, charge_branch_arg, mpmt_branch_arg, pmt_branch_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup = np.inf, window_prompt = 100):
    #print('inside loop')
    def neutron_detection_event(event_number, times_branch_event_arg, charge_branch_event_arg, mpmt_branch_event_arg, pmt_branch_event_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt):
        #print('checking event', event_number)
        dict_neutrons_event = {}
        valid_thresholds = []
        last_prompt = None
        prompt_lockout_until = None
        print(f'Event {event_number} has {len(threshold_times)} prompt candidates')
        #for time_prompt in threshold_times:
        for time_prompt, nhits_prompt, trms_prompt, times_prompt, _ , mpmt_prompt, pmt_prompt, x_prompt, y_prompt, z_prompt  in threshold_times:
            print(f'Prompt time {time_prompt}')    
            if last_prompt is not None and (time_prompt-last_prompt) < (window_sliding + window_prompt):
                continue

            if prompt_lockout_until is not None and time_prompt < prompt_lockout_until:
                continue

            mask = (times_branch_event_arg >= time_prompt + window_prompt) & (times_branch_event_arg < time_prompt + window_prompt + window_sliding)
            #print('After isolating prompt')
            if mask.sum() == 0:
                continue

            neutron_candidates,  neutron_nhits = nHitsTimeWindow(times_branch_event_arg[mask], threshold_inf, window_neutron, 0, threshold_sup=threshold_sup)
            #print('Delayed found, checking now masks')
            #for time_delayed in neutron_candidates:

            ##########################################

            if neutron_candidates is None or neutron_nhits is None:
                continue
            if len(neutron_candidates) == 0 or len(neutron_nhits) == 0:
                continue

            # If they can ever be mismatched, either skip or truncate safely
            if len(neutron_candidates) != len(neutron_nhits):
                # safest behavior: skip this prompt (or you can handle it differently)
                print("Warning: mismatch candidates/nhits", len(neutron_candidates), len(neutron_nhits))
                continue

            # Highest nHits candidate
            i_best = int(np.argmax(neutron_nhits))
            time_delayed = neutron_candidates[i_best]
            delayed_nhits = neutron_nhits[i_best]

            # Earliest candidate
            #i_first = int(np.argmin(neutron_candidates))  # earliest time
            #time_delayed = neutron_candidates[i_first]
            #delayed_nhits = neutron_nhits[i_first]

            mask_delayed = (times_branch_event_arg >= time_delayed) & (times_branch_event_arg < time_delayed + window_sliding)

            times_in_delayed = times_branch_event_arg[mask_delayed]
            charges_in_delayed = charge_branch_event_arg[mask_delayed]
            mpmt_in_delayed = mpmt_branch_event_arg[mask_delayed]
            pmt_in_delayed = pmt_branch_event_arg[mask_delayed]

            #print('Running bonsai')
            if (event_number == 6621):
                #print(times_in_delayed, charges_in_delayed, mpmt_in_delayed, pmt_in_delayed)
                continue
            vertex = functions_bonsai.run_BONSAI_candidate(times_in_delayed, charges_in_delayed, mpmt_in_delayed, pmt_in_delayed)
            #print('vertex okay, saving info')
            valid_thresholds.append(
                {
                    'event_number': event_number,
                    'prompt_nhits': nhits_prompt,
                    'prompt_time': float(time_prompt),
                    'prompt_trms': trms_prompt,
                    'prompt_times': times_prompt,
                    'prompt_x': x_prompt,
                    'prompt_y': y_prompt,
                    'prompt_z': z_prompt,
                    'prompt_mpmt': mpmt_prompt,
                    'prompt_pmt': pmt_prompt,
                    "delayed_time": float(time_delayed),
                    "delayed_nhits": delayed_nhits,
                    "delayed_x": vertex["x"][0],
                    "delayed_y": vertex["y"][0],
                    "delayed_z": vertex["z"][0],
                }
            )
            print("  picked delayed", time_delayed, "-> lockout", time_delayed + window_neutron)

            prompt_lockout_until = time_delayed + window_neutron
            

        return valid_thresholds
    
    dict_neutrons = {}
    for event in event_branch:
        if event == 6621:
            continue
        if event % 1000 == 0:
            print(f"Searching neutron coincidences on event {event}...")
        if event in threshold_times:
            if len(threshold_times[event]) == 0:
                continue 
            results = neutron_detection_event(event, times_branch_arg[event], charge_branch_arg[event], mpmt_branch_arg[event], pmt_branch_arg[event], threshold_times[event], window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt)
            if len(results) != 0:
                dict_neutrons[event] = results

    return dict_neutrons
        #else:
            #print(f"Event {event} has no threshold times, skipping neutron detection.")

def neutron_detection_wMulti(event_branch, times_branch_arg, charge_branch_arg, mpmt_branch_arg, pmt_branch_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup = np.inf, window_prompt = 100):
    #print('inside loop')
    def neutron_detection_event(event_number, times_branch_event_arg, charge_branch_event_arg, mpmt_branch_event_arg, pmt_branch_event_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt):
        #print('checking event', event_number)
        dict_neutrons_event = {}
        valid_thresholds = []
        last_prompt = None
        prompt_lockout_until = None
        print(f'Event {event_number} has {len(threshold_times)} prompt candidates')
        #for time_prompt in threshold_times:
        for time_prompt, nhits_prompt, trms_prompt, times_prompt, _ , mpmt_prompt, pmt_prompt, x_prompt, y_prompt, z_prompt  in threshold_times:
            print(f'Prompt time {time_prompt}')    
            if last_prompt is not None and (time_prompt-last_prompt) < (window_sliding + window_prompt):
                continue

            if prompt_lockout_until is not None and time_prompt < prompt_lockout_until:
                continue

            mask = (times_branch_event_arg >= time_prompt + window_prompt) & (times_branch_event_arg < time_prompt + window_prompt + window_sliding)
            if mask.sum() == 0:
                continue

            neutron_candidates,  neutron_nhits = nHitsTimeWindow(times_branch_event_arg[mask], threshold_inf, window_neutron, window_neutron, threshold_sup=threshold_sup)

            ##########################################

            if neutron_candidates is None or neutron_nhits is None:
                continue
            if len(neutron_candidates) == 0 or len(neutron_nhits) == 0:
                continue

            # If they can ever be mismatched, either skip or truncate safely
            if len(neutron_candidates) != len(neutron_nhits):
                # safest behavior: skip this prompt (or you can handle it differently)
                print("Warning: mismatch candidates/nhits", len(neutron_candidates), len(neutron_nhits))
                continue

            # Highest nHits candidate
            i_best = int(np.argmax(neutron_nhits))
            time_delayed = neutron_candidates[i_best]
            delayed_nhits = neutron_nhits[i_best]


            # Earliest candidate
            #i_first = int(np.argmin(neutron_candidates))  # earliest time
            #time_delayed = neutron_candidates[i_first]
            #delayed_nhits = neutron_nhits[i_first]

            ###
            ##mask_delayed = (times_branch_event_arg >= time_delayed) & (times_branch_event_arg < time_delayed + window_sliding)
            mask_delayed = mask & (times_branch_event_arg >= time_delayed) & (times_branch_event_arg < time_delayed + window_neutron)
            
            times_in_delayed = times_branch_event_arg[mask_delayed]
            charges_in_delayed = charge_branch_event_arg[mask_delayed]
            mpmt_in_delayed = mpmt_branch_event_arg[mask_delayed]
            pmt_in_delayed = pmt_branch_event_arg[mask_delayed]
            #print('Running bonsai')
            if (event_number == 6621):
                #print(times_in_delayed, charges_in_delayed, mpmt_in_delayed, pmt_in_delayed)
                continue
            t_in = np.min(times_in_delayed)

            t_rms = functions_analysis.time_RMS_fun_time(times_in_delayed, t_in, window_neutron)
            #vertex = functions_multilateration.run_multilateration_candidate(times_in_delayed, mpmt_in_delayed, pmt_in_delayed, sigma_t=1.0)
            vertex = functions_multilateration.run_multilateration_timecal(times_in_delayed, mpmt_in_delayed, pmt_in_delayed, sigma_t=1.0)
            #print('vertex okay, saving info')
            valid_thresholds.append(
                {
                    'event_number': event_number,
                    'prompt_nhits': nhits_prompt,
                    'prompt_time': float(time_prompt),
                    'prompt_trms': trms_prompt,
                    'prompt_times': times_prompt,
                    'prompt_x': x_prompt,
                    'prompt_y': y_prompt,
                    'prompt_z': z_prompt,
                    'prompt_mpmt': mpmt_prompt,
                    'prompt_pmt': pmt_prompt,
                    "delayed_time": float(time_delayed),
                    "delayed_trms": float(t_rms),
                    "delayed_nhits": delayed_nhits,
                    "delayed_x": vertex["x"],
                    "delayed_y": vertex["y"],
                    "delayed_z": vertex["z"],
                    "vertex_chi2": vertex.get("chi2_ndof", None),
                }
            )
            print("  picked delayed", time_delayed, "-> lockout", time_delayed + window_neutron)

            prompt_lockout_until = time_delayed + window_neutron


        return valid_thresholds
    
    dict_neutrons = {}
    for event in event_branch:
        if event == 6621:
            continue
        if event % 1000 == 0:
            print(f"Searching neutron coincidences on event {event}...")
        if event in threshold_times:
            if len(threshold_times[event]) == 0:
                continue 
            results = neutron_detection_event(event, times_branch_arg[event], charge_branch_arg[event], mpmt_branch_arg[event], pmt_branch_arg[event], threshold_times[event], window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt)
            if len(results) != 0:
                dict_neutrons[event] = results

    return dict_neutrons
        #else:
            #print(f"Event {event} has no threshold times, skipping neutron detection.")

def neutron_detection_wMulti_wtRMS(event_branch, times_branch_arg, charge_branch_arg, mpmt_branch_arg, pmt_branch_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup = np.inf, window_prompt = 100):

    def neutron_detection_event(event_number, times_branch_event_arg, charge_branch_event_arg, mpmt_branch_event_arg, pmt_branch_event_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt):

        dict_neutrons_event = {}
        valid_thresholds = []
        last_prompt = None
        prompt_lockout_until = None
        print(f'Event {event_number} has {len(threshold_times)} prompt candidates')
        #for time_prompt in threshold_times:
        for time_prompt, nhits_prompt, trms_prompt, times_prompt, _ , mpmt_prompt, pmt_prompt, x_prompt, y_prompt, z_prompt  in threshold_times:
            print(f'Prompt time {time_prompt}')    
            if last_prompt is not None and (time_prompt-last_prompt) < (window_sliding + window_prompt):
                continue

            if prompt_lockout_until is not None and time_prompt < prompt_lockout_until:
                continue

            mask = (times_branch_event_arg >= time_prompt + window_prompt) & (times_branch_event_arg < time_prompt + window_prompt + window_sliding)
            if mask.sum() == 0:
                continue

            neutron_candidates,  neutron_nhits = nHitsTimeWindow(times_branch_event_arg[mask], threshold_inf, window_neutron, window_neutron, threshold_sup=threshold_sup)

            best = None  # will store dict with everything we need
            for time_delayed, delayed_nhits in zip(neutron_candidates, neutron_nhits):

                mask_delayed = mask & (times_branch_event_arg >= time_delayed) & (times_branch_event_arg < time_delayed + window_neutron)
                if mask_delayed.sum() < 6:
                    continue

                times_in_delayed   = times_branch_event_arg[mask_delayed]
                charges_in_delayed = charge_branch_event_arg[mask_delayed]
                mpmt_in_delayed    = mpmt_branch_event_arg[mask_delayed]
                pmt_in_delayed     = pmt_branch_event_arg[mask_delayed]

                # pre RMS (optional)
                t_in = float(np.min(times_in_delayed))
                t_rms_pre = functions_analysis.time_RMS_fun_time(times_in_delayed, t_in, window_neutron)

                # vertex fit
                vertex = functions_multilateration.run_multilateration_candidate(
                    times_in_delayed, mpmt_in_delayed, pmt_in_delayed, sigma_t=1.0
                )
                if isinstance(vertex, dict) and ("success" in vertex) and (not vertex["success"]):
                    continue

                # extract (x,y,z,t0)
                if isinstance(vertex, dict):
                    v_xyz = np.array([vertex["x"], vertex["y"], vertex["z"]], dtype=float)
                    t0_guess = float(vertex["t0"]) if "t0" in vertex else float(np.min(times_in_delayed))
                    chi2 = vertex.get("chi2_ndof", None)
                else:
                    v_arr = np.asarray(vertex, dtype=float).ravel()
                    v_xyz = v_arr[:3]
                    t0_guess = float(v_arr[3]) if v_arr.size >= 4 else float(np.min(times_in_delayed))
                    chi2 = None

                # geometry -> PMT positions (cm)
                x_p, y_p, z_p, _ = functions_bonsai.getxyz(functions_bonsai.geo, mpmt_in_delayed, pmt_in_delayed)
                pmt_pos_cm = np.column_stack([x_p, y_p, z_p])

                # dt residuals (cm/ns)
                c_water = 29.9792458 / 1.33  # cm/ns
                L_cm = np.linalg.norm(pmt_pos_cm - v_xyz[None, :], axis=1)
                tof_ns = L_cm / c_water
                dt_ns = times_in_delayed - t0_guess - tof_ns

                good = (np.abs(dt_ns) < 3.0)
                clean_nhits = int(good.sum())
                if clean_nhits < 6:
                    continue

                # cleaned hit lists
                times_cl   = times_in_delayed[good]
                charges_cl = charges_in_delayed[good]
                mpmt_cl    = mpmt_in_delayed[good]
                pmt_cl     = pmt_in_delayed[good]

                t_in_cl = float(np.min(times_cl))
                t_rms_aft = functions_analysis.time_RMS_fun_time(times_cl, t_in_cl, window_neutron)

                # (optional) refit vertex after cleaning
                vertex_cl = functions_multilateration.run_multilateration_candidate(
                    times_cl, mpmt_cl, pmt_cl, sigma_t=1.0
                )
                if isinstance(vertex_cl, dict) and ("success" in vertex_cl) and (not vertex_cl["success"]):
                    continue

                # get final coords + chi2 after cleaning
                if isinstance(vertex_cl, dict):
                    vx, vy, vz = float(vertex_cl["x"]), float(vertex_cl["y"]), float(vertex_cl["z"])
                    chi2_cl = vertex_cl.get("chi2_ndof", None)
                else:
                    v_arr2 = np.asarray(vertex_cl, dtype=float).ravel()
                    vx, vy, vz = float(v_arr2[0]), float(v_arr2[1]), float(v_arr2[2])
                    chi2_cl = None

                # -------- scoring: prefer more cleaned hits, then smaller RMS, then smaller chi2 --------
                # You can tune this easily.
                score_key = (
                    clean_nhits,          # maximize
                    -float(t_rms_aft),    # minimize RMS -> maximize negative
                    -(chi2_cl if (chi2_cl is not None) else 0.0),  # minimize chi2 if available
                    -float(delayed_nhits) # tie-breaker: original nhits
                )

                if (best is None) or (score_key > best["score_key"]):
                    best = {
                        "score_key": score_key,
                        "time_delayed": float(time_delayed),
                        "delayed_nhits_raw": int(delayed_nhits),
                        "clean_nhits": clean_nhits,
                        "t_rms_pre": float(t_rms_pre),
                        "t_rms_aft": float(t_rms_aft),
                        "vertex": vertex_cl,   # final vertex after cleaning
                        "chi2": chi2_cl,
                    }

            # After evaluating all candidates, either keep best or skip this prompt
            if best is None:
                continue

            time_delayed = best["time_delayed"]
            delayed_nhits = best["delayed_nhits_raw"]
            t_rms_aft = best["t_rms_aft"]
            vertex = best["vertex"]

            print("  picked BEST delayed", time_delayed, "clean_nhits", best["clean_nhits"], "-> lockout", time_delayed + window_neutron)

            # lockout and append to valid_thresholds (same as you already do)
            prompt_lockout_until = time_delayed + window_neutron

            valid_thresholds.append(
                {
                    'event_number': event_number,
                    'prompt_nhits': nhits_prompt,
                    'prompt_time': float(time_prompt),
                    'prompt_trms': trms_prompt,
                    'prompt_times': times_prompt,
                    'prompt_x': x_prompt,
                    'prompt_y': y_prompt,
                    'prompt_z': z_prompt,
                    'prompt_mpmt': mpmt_prompt,
                    'prompt_pmt': pmt_prompt,
                    "delayed_time": float(time_delayed),
                    "delayed_trms": float(t_rms_aft),
                    "delayed_nhits": int(delayed_nhits),
                    "delayed_x": float(vertex["x"]) if isinstance(vertex, dict) else float(np.asarray(vertex).ravel()[0]),
                    "delayed_y": float(vertex["y"]) if isinstance(vertex, dict) else float(np.asarray(vertex).ravel()[1]),
                    "delayed_z": float(vertex["z"]) if isinstance(vertex, dict) else float(np.asarray(vertex).ravel()[2]),
                    "vertex_chi2": (vertex.get("chi2_ndof", None) if isinstance(vertex, dict) else None),
                }
            )
            print("  picked delayed", time_delayed, "-> lockout", time_delayed + window_neutron)

            prompt_lockout_until = time_delayed + window_neutron

        return valid_thresholds
    
    dict_neutrons = {}
    for event in event_branch:
        if event == 6621:
            continue
        if event % 1000 == 0:
            print(f"Searching neutron coincidences on event {event}...")
        if event in threshold_times:
            if len(threshold_times[event]) == 0:
                continue 
            results = neutron_detection_event(event, times_branch_arg[event], charge_branch_arg[event], mpmt_branch_arg[event], pmt_branch_arg[event], threshold_times[event], window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt)
            if len(results) != 0:
                dict_neutrons[event] = results

    return dict_neutrons



def neutron_detection_wMulti_wtcut(
    event_branch,
    times_branch_tof_arg,
    times_branch_raw_arg,
    charge_branch_arg,
    mpmt_branch_arg,
    pmt_branch_arg,
    threshold_times,
    window_sliding,
    window_neutron,
    threshold_inf,
    threshold_sup=np.inf,
    window_prompt=100,
):


    def neutron_detection_event(
        event_number,
        times_branch_tof_event_arg,
        times_branch_raw_event_arg,
        charge_branch_event_arg,
        mpmt_branch_event_arg,
        pmt_branch_event_arg,
        threshold_times_event,
        window_sliding,
        window_neutron,
        threshold_inf,
        threshold_sup,
        window_prompt,
    ):
        dict_neutrons_event = {}
        valid_thresholds = []
        last_prompt = None
        prompt_lockout_until = None
        print(f'Event {event_number} has {len(threshold_times_event)} prompt candidates')
        
        #for time_prompt, nhits_prompt, trms_prompt, times_prompt, _ , mpmt_prompt, pmt_prompt, x_prompt, y_prompt, z_prompt  in threshold_times:
        for cand in threshold_times_event:
            if isinstance(cand, dict):
                time_prompt = cand.get("prompt_ti", cand.get("time"))
                nhits_prompt = cand.get("prompt_nhits", cand.get("n_hits"))
                trms_prompt = cand.get("prompt_trms")
                tmean_prompt = cand.get("prompt_tmean")
                times_prompt = cand.get("prompt_times_tof", cand.get("prompt_times"))
                mpmt_prompt = cand.get("prompt_mpmt")
                pmt_prompt = cand.get("prompt_pmt")
                x_prompt = cand.get("prompt_vertex_x", cand.get("vertex_x"))
                y_prompt = cand.get("prompt_vertex_y", cand.get("vertex_y"))
                z_prompt = cand.get("prompt_vertex_z", cand.get("vertex_z"))

            elif isinstance(cand, (tuple, list)):
                if len(cand) == 11:
                    (
                        time_prompt,
                        nhits_prompt,
                        trms_prompt,
                        tmean_prompt,
                        times_prompt,
                        _,
                        mpmt_prompt,
                        pmt_prompt,
                        x_prompt,
                        y_prompt,
                        z_prompt,
                    ) = cand

                elif len(cand) == 12:
                    (
                        time_prompt,
                        nhits_prompt,
                        trms_prompt,
                        tmean_prompt,
                        times_prompt_tof,
                        times_prompt_raw,
                        charge_prompt,
                        mpmt_prompt,
                        pmt_prompt,
                        x_prompt,
                        y_prompt,
                        z_prompt,
                    ) = cand
                    times_prompt = times_prompt_tof

                else:
                    print(f"Skipping event {event_number}: unsupported candidate format {cand}")
                    continue

            else:
                print(f"Skipping event {event_number}: candidate is scalar/non-iterable: {cand} ({type(cand)})")
                continue
            
            print(f'Prompt time {time_prompt}')    
            if last_prompt is not None and (time_prompt-last_prompt) < (window_sliding + window_prompt):
                continue

            if prompt_lockout_until is not None and time_prompt < prompt_lockout_until:
                continue

            mask = (times_branch_tof_event_arg >= time_prompt + window_prompt) & (times_branch_tof_event_arg < time_prompt + window_prompt + window_sliding)
            #mask = (times_branch_event_arg >= time_prompt + window_prompt) & (times_branch_event_arg < time_prompt + window_prompt + window_sliding)
            
            if mask.sum() == 0:
                continue

            neutron_candidates, neutron_nhits, neutron_trms = nHitstRMSTimeWindow(
                times_branch_tof_event_arg[mask],
                #times_branch_event_arg[mask],
                threshold_inf=threshold_inf,
                window=window_neutron,
                death_window=window_neutron,
                time_rms_fun=functions_analysis.time_RMS_fun_time,
                rms_cut_ns=10.0,
                threshold_sup=threshold_sup,
            )

            if neutron_candidates is None or neutron_nhits is None or len(neutron_candidates) == 0:
                continue

            if neutron_candidates is None or neutron_nhits is None:
                continue
            if len(neutron_candidates) == 0 or len(neutron_nhits) == 0:
                continue

            # If they can ever be mismatched, truncate safely (or skip)
            n_cand = min(len(neutron_candidates), len(neutron_nhits))
            neutron_candidates = neutron_candidates[:n_cand]
            neutron_nhits = neutron_nhits[:n_cand]

            # Loop over ALL delayed candidates (don’t just take argmax)
            for time_delayed, delayed_nhits in zip(neutron_candidates, neutron_nhits):

                #mask_delayed = (times_branch_event_arg >= time_delayed) & (times_branch_event_arg < time_delayed + window_sliding)
                mask_delayed = (times_branch_tof_event_arg >= time_delayed) & (times_branch_tof_event_arg < time_delayed + window_sliding)
                
                if mask_delayed.sum() == 0:
                    continue

                #times_in_delayed = times_branch_event_arg[mask_delayed]
                times_in_delayed_tof = times_branch_tof_event_arg[mask_delayed]
                times_in_delayed_raw = times_branch_raw_event_arg[mask_delayed]
                charges_in_delayed = charge_branch_event_arg[mask_delayed]
                mpmt_in_delayed = mpmt_branch_event_arg[mask_delayed]
                pmt_in_delayed = pmt_branch_event_arg[mask_delayed]

                #t_in = np.min(times_in_delayed)
                #t_rms = functions_analysis.time_RMS_fun_time(times_in_delayed, t_in, window_neutron)

                t_in = np.min(times_in_delayed_tof)
                t_rms, _ = functions_analysis.time_RMS_fun_time(times_in_delayed_tof, t_in, window_neutron)
                # avoid pathological event if you still need this special-case
                if event_number == 6621:
                    continue

                # Reconstruct delayed vertex
                #vertex = functions_multilateration.run_multilateration_candidate(times_in_delayed_raw, mpmt_in_delayed, pmt_in_delayed, sigma_t=1.0)
                vertex = functions_multilateration.run_multilateration_timecal(times_in_delayed_raw, mpmt_in_delayed, pmt_in_delayed, sigma_t=1.0)
                # If recon fails, skip
                if not vertex.get("success", True):
                    continue

                valid_thresholds.append(
                    {
                        "event_number": event_number,

                        "prompt_nhits": nhits_prompt,
                        "prompt_time": float(time_prompt),
                        "prompt_trms": None if trms_prompt is None else float(trms_prompt),
                        "prompt_tmean": None if tmean_prompt is None else float(tmean_prompt),
                        #"prompt_times_tof": times_prompt_tof,
                        #"prompt_times_raw": times_prompt_raw,
                        #"prompt_x": x_prompt,
                        #"prompt_y": y_prompt,
                        #"prompt_z": z_prompt,
                        #"prompt_mpmt": mpmt_prompt,
                        #"prompt_pmt": pmt_prompt,

                        "delayed_time_tof": float(time_delayed),
                        "delayed_nhits": int(delayed_nhits),
                        "delayed_trms_tof": float(t_rms),

                        #"delayed_times_tof": times_in_delayed_tof.tolist(),
                        #"delayed_times_raw": times_in_delayed_raw.tolist(),
                        #"delayed_charge": charges_in_delayed.tolist(),
                        #"delayed_mpmt": mpmt_in_delayed.tolist(),
                        #"delayed_pmt": pmt_in_delayed.tolist(),

                        "delayed_x": float(vertex["x"]),
                        "delayed_y": float(vertex["y"]),
                        "delayed_z": float(vertex["z"]),
                        "vertex_chi2": float(vertex.get("chi2_ndof", np.nan)),
                    }
                )

            print("  picked delayed", time_delayed, "-> lockout", time_delayed + window_neutron)

        return valid_thresholds
    
    dict_neutrons = {}
    for event in event_branch:
        if event == 6621:
            continue
        if event % 1000 == 0:
            print(f"Searching neutron coincidences on event {event}...")
        if event in threshold_times:
            if len(threshold_times[event]) == 0:
                continue 
            #results = neutron_detection_event(event, times_branch_arg[event], charge_branch_arg[event], mpmt_branch_arg[event], pmt_branch_arg[event], threshold_times[event], window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt)
            results = neutron_detection_event(
            event,
            times_branch_tof_arg[event],
            times_branch_raw_arg[event],
            charge_branch_arg[event],
            mpmt_branch_arg[event],
            pmt_branch_arg[event],
            threshold_times[event],
            window_sliding,
            window_neutron,
            threshold_inf,
            threshold_sup,
            window_prompt,
            )
            if len(results) != 0:
                dict_neutrons[event] = results

    return dict_neutrons

def neutron_detection_wMulti_wtcut_clean(
    event_branch,
    times_branch_arg,
    charge_branch_arg,
    mpmt_branch_arg,
    pmt_branch_arg,
    threshold_times,
    window_sliding,
    window_neutron,
    threshold_inf,
    threshold_sup=np.inf,
    window_prompt=100,
    # cleaning cut params
    clean_half_window=100,     # +/- 100 ns around delayed (total extra 200 ns)
    clean_max_hits=200,         # reject if >200 hits in that widened window
):

    def neutron_detection_event_clean(
        event_number,
        times_branch_event_arg,
        charge_branch_event_arg,
        mpmt_branch_event_arg,
        pmt_branch_event_arg,
        threshold_times,
        window_sliding,
        window_neutron,
        threshold_inf,
        threshold_sup,
        window_prompt,
        clean_half_window,
        clean_max_hits,
    ):
        valid_thresholds = []
        last_prompt = None
        prompt_lockout_until = None

        print(f"Event {event_number} has {len(threshold_times)} prompt candidates")

        for (
            time_prompt, nhits_prompt, trms_prompt, times_prompt, _,
            mpmt_prompt, pmt_prompt, x_prompt, y_prompt, z_prompt
        ) in threshold_times:

            print(f"Prompt time {time_prompt}")

            if last_prompt is not None and (time_prompt - last_prompt) < (window_sliding + window_prompt):
                continue
            if prompt_lockout_until is not None and time_prompt < prompt_lockout_until:
                continue

            mask = (
                (times_branch_event_arg >= time_prompt + window_prompt) &
                (times_branch_event_arg <  time_prompt + window_prompt + window_sliding)
            )
            if mask.sum() == 0:
                continue

            neutron_candidates, neutron_nhits, neutron_trms = nHitstRMSTimeWindow(
                times_branch_event_arg[mask],
                threshold_inf=threshold_inf,
                window=window_neutron,
                death_window=window_neutron,
                time_rms_fun=functions_analysis.time_RMS_fun_time,
                rms_cut_ns=10.0,
                threshold_sup=threshold_sup,
            )

            if neutron_candidates is None or neutron_nhits is None or len(neutron_candidates) == 0:
                continue

            n_cand = min(len(neutron_candidates), len(neutron_nhits))
            neutron_candidates = neutron_candidates[:n_cand]
            neutron_nhits = neutron_nhits[:n_cand]

            # --- cleaning cut state (per prompt) ---
            # candidates earlier than this are auto-skipped (like your t0_delayed logic)
            t0_delayed = -float("inf")

            for time_delayed, delayed_nhits in zip(neutron_candidates, neutron_nhits):

                # 1) rolling veto / lockout for delayed candidates
                if time_delayed < t0_delayed:
                    continue

                # 2) 200 ns cleaning window around delayed (±half_window, plus your delayed window)
                #    You used: [time_delayed-half_window, time_delayed+window_neutron+half_window]
                mask_200ns = (
                    (times_branch_event_arg >= time_delayed - clean_half_window) &
                    (times_branch_event_arg <= time_delayed + window_neutron + clean_half_window)
                )
                if int(mask_200ns.sum()) > clean_max_hits:
                    t0_delayed = time_delayed + window_neutron + clean_half_window
                    continue

                # Your original delayed-window content
                mask_delayed = (
                    (times_branch_event_arg >= time_delayed) &
                    (times_branch_event_arg <  time_delayed + window_sliding)
                )
                if mask_delayed.sum() == 0:
                    continue

                times_in_delayed = times_branch_event_arg[mask_delayed]
                charges_in_delayed = charge_branch_event_arg[mask_delayed]
                mpmt_in_delayed = mpmt_branch_event_arg[mask_delayed]
                pmt_in_delayed = pmt_branch_event_arg[mask_delayed]

                t_in = float(times_in_delayed.min())
                t_rms = float(functions_analysis.time_RMS_fun_time(times_in_delayed, t_in, window_neutron))

                # preserve your special-case skip
                if event_number == 6621:
                    continue

                vertex = functions_multilateration.run_multilateration_candidate(
                    times_in_delayed, mpmt_in_delayed, pmt_in_delayed, sigma_t=1.0
                )
                if not vertex.get("success", True):
                    continue

                valid_thresholds.append(
                    {
                        "event_number": event_number,

                        "prompt_nhits": nhits_prompt,
                        "prompt_time": float(time_prompt),
                        "prompt_trms": trms_prompt,
                        "prompt_times": times_prompt,
                        "prompt_x": x_prompt,
                        "prompt_y": y_prompt,
                        "prompt_z": z_prompt,
                        "prompt_mpmt": mpmt_prompt,
                        "prompt_pmt": pmt_prompt,

                        "delayed_time": float(time_delayed),
                        "delayed_nhits": int(delayed_nhits),
                        "delayed_trms": float(t_rms),
                        "delayed_x": float(vertex["x"]),
                        "delayed_y": float(vertex["y"]),
                        "delayed_z": float(vertex["z"]),
                        "vertex_chi2": float(vertex.get("chi2_ndof", np.nan)),
                    }
                )

            # if you want prompt-level lockout based on the *last accepted delayed*,
            # you can set it here (optional). Example: lock out until end of last veto.
            # prompt_lockout_until = t0_delayed

        return valid_thresholds

    dict_neutrons = {}
    for event in event_branch:
        if event == 6621:
            continue
        if event % 1000 == 0:
            print(f"Searching neutron coincidences on event {event}...")

        if event in threshold_times:
            if len(threshold_times[event]) == 0:
                continue

            results = neutron_detection_event_clean(
                event,
                times_branch_arg[event],
                charge_branch_arg[event],
                mpmt_branch_arg[event],
                pmt_branch_arg[event],
                threshold_times[event],
                window_sliding,
                window_neutron,
                threshold_inf,
                threshold_sup,
                window_prompt,
                clean_half_window,
                clean_max_hits,
            )
            if len(results) != 0:
                dict_neutrons[event] = results

    return dict_neutrons

def neutron_detection_wMulti_clean(event_branch, times_branch_arg, charge_branch_arg, mpmt_branch_arg, pmt_branch_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup = np.inf, window_prompt = 100):
    #print('inside loop')
    def neutron_detection_event(event_number, times_branch_event_arg, charge_branch_event_arg, mpmt_branch_event_arg, pmt_branch_event_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt):
        #print('checking event', event_number)

        half_window = 200  # ns

        dict_neutrons_event = {}
        valid_thresholds = []
        last_prompt = None
        prompt_lockout_until = None
        print(f'Event {event_number} has {len(threshold_times)} prompt candidates')
        #for time_prompt in threshold_times:
        for time_prompt, nhits_prompt, trms_prompt, times_prompt, _ , mpmt_prompt, pmt_prompt, x_prompt, y_prompt, z_prompt  in threshold_times:
            print(f'Prompt time {time_prompt}')    
            if last_prompt is not None and (time_prompt-last_prompt) < (window_sliding + window_prompt):
                continue

            if prompt_lockout_until is not None and time_prompt < prompt_lockout_until:
                continue

            mask = (times_branch_event_arg >= time_prompt + window_prompt) & (times_branch_event_arg < time_prompt + window_prompt + window_sliding)
            if mask.sum() == 0:
                continue

            times_event = times_branch_event_arg[mask]
            charge_event = charge_branch_event_arg[mask]
            mpmt_event   = mpmt_branch_event_arg[mask]
            pmt_event    = pmt_branch_event_arg[mask]

            neutron_candidates,  _ = nHitsTimeWindow(times_event, threshold_inf,  window_neutron, 0,threshold_sup=threshold_sup,)

            t0_delayed = time_prompt

            for time_delayed in neutron_candidates:
                # Reject candidate if there are >200 hits in a 200 ns window around it (±100 ns)
                if time_delayed < t0_delayed:
                    continue

                mask_200ns = (times_event >= time_delayed - half_window) & (times_event <= time_delayed + window_neutron + half_window)
                if int(mask_200ns.sum()) > 50:
                    t0_delayed = time_delayed + window_neutron + half_window
                    continue

                # ...keep the rest of your candidate logic here (no clean_1 / clean_2)
                mask_delayed = (times_event >= time_delayed) & (times_event < time_delayed + window_neutron)
                n_hits_here = int(mask_delayed.sum())
                if n_hits_here == 0:
                    continue

                times_in_delayed   = times_event[mask_delayed]
                charges_in_delayed = charge_event[mask_delayed]
                mpmt_in_delayed    = mpmt_event[mask_delayed]
                pmt_in_delayed     = pmt_event[mask_delayed]

                t_in = np.min(times_in_delayed)
                t_rms = functions_analysis.time_RMS_fun_time(times_in_delayed, t_in, window_neutron)
                vertex = functions_multilateration.run_multilateration_candidate(times_in_delayed, mpmt_in_delayed, pmt_in_delayed, sigma_t=1.0)
                #print('vertex okay, saving info')
                valid_thresholds.append(
                    {
                        'event_number': event_number,
                        'prompt_nhits': nhits_prompt,
                        'prompt_time': float(time_prompt),
                        'prompt_trms': trms_prompt,
                        'prompt_times': times_prompt,
                        'prompt_x': x_prompt,
                        'prompt_y': y_prompt,
                        'prompt_z': z_prompt,
                        'prompt_mpmt': mpmt_prompt,
                        'prompt_pmt': pmt_prompt,
                        "delayed_time": float(time_delayed),
                        "delayed_trms": float(t_rms),
                        "delayed_nhits": n_hits_here,
                        "delayed_x": vertex["x"],
                        "delayed_y": vertex["y"],
                        "delayed_z": vertex["z"],
                        "vertex_chi2": vertex.get("chi2_ndof", None),
                    }
                )
                print("  picked delayed", time_delayed, "-> lockout", time_delayed + window_neutron)

                prompt_lockout_until = time_delayed + window_neutron

        return valid_thresholds
    
    dict_neutrons = {}
    for event in event_branch:
        if event == 6621:
            continue
        if event % 1000 == 0:
            print(f"Searching neutron coincidences on event {event}...")
        if event in threshold_times:
            if len(threshold_times[event]) == 0:
                continue 
            results = neutron_detection_event(event, times_branch_arg[event], charge_branch_arg[event], mpmt_branch_arg[event], pmt_branch_arg[event], threshold_times[event], window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt)
            if len(results) != 0:
                dict_neutrons[event] = results

    return dict_neutrons
        #else:
            #print(f"Event {event} has no threshold times, skipping neutron detection.")

            

def neutron_detection_wMulti_clean_prompts(event_branch, times_branch_arg, charge_branch_arg, mpmt_branch_arg, pmt_branch_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup = np.inf, window_prompt = 100):
    def neutron_detection_event(
        event_number,
        times_branch_event_arg,
        charge_branch_event_arg,
        mpmt_branch_event_arg,
        pmt_branch_event_arg,
        threshold_times,
        window_sliding,
        window_neutron,
        threshold_inf,
        threshold_sup,
        window_prompt,
    ):
        # --- Tuning parameters ---
        half_window = 100          # ns  -> 200 ns total window around delayed candidate
        max_hits_200ns = 200       # reject if hits in ±100 ns exceeds this

        single_prompt_search = window_sliding  # ns (150 us)
        event_end_window = 270000      # ns (270 us after prompt integration window)

        valid_thresholds = []
        last_prompt = None
        prompt_lockout_until = None

        # Ensure prompts are time-ordered
        threshold_times = sorted(threshold_times, key=lambda x: x[0])

        n_prompts = len(threshold_times)
        print(f"Event {event_number} has {n_prompts} prompt candidates")

        for i, (time_prompt, nhits_prompt, trms_prompt, times_prompt, _,
                mpmt_prompt, pmt_prompt, x_prompt, y_prompt, z_prompt) in enumerate(threshold_times):

            print(f"Prompt time {time_prompt}")

            # Optional: keep your existing skip/lockout logic
            if last_prompt is not None and (time_prompt - last_prompt) < (window_sliding + window_prompt):
                continue
            if prompt_lockout_until is not None and time_prompt < prompt_lockout_until:
                continue

            # Start searching after the prompt integration window
            t0 = time_prompt + window_prompt

            if n_prompts == 1:
                # Only one prompt: search a fixed 150 us region after it
                t1 = t0 + single_prompt_search
            else:
                # More than one prompt:
                if i + 1 < n_prompts:
                    # search between this prompt and the next prompt time
                    next_time_prompt = threshold_times[i + 1][0]
                    t1 = next_time_prompt
                else:
                    # last prompt: search until end of the (270 us) window
                    t1 = t0 + event_end_window

            # Guard against empty/inverted windows
            if t1 <= t0:
                last_prompt = time_prompt
                continue

            mask = (times_branch_event_arg >= t0) & (times_branch_event_arg < t1)
            if mask.sum() == 0:
                last_prompt = time_prompt
                continue

            times_event  = times_branch_event_arg[mask]
            charge_event = charge_branch_event_arg[mask]
            mpmt_event   = mpmt_branch_event_arg[mask]
            pmt_event    = pmt_branch_event_arg[mask]

            neutron_candidates, _ = nHitsTimeWindow(
                times_event,
                threshold_inf,
                window_neutron,
                0,
                threshold_sup=threshold_sup,
            )

            for time_delayed in neutron_candidates:
                # Reject if too many hits in a 200 ns window around delayed candidate
                mask_200ns = (times_event >= time_delayed - half_window) & (times_event <= time_delayed + half_window)
                if int(mask_200ns.sum()) > max_hits_200ns:
                    continue

                # Count hits in the delayed candidate window
                mask_delayed = (times_event >= time_delayed) & (times_event < time_delayed + window_neutron)
                n_hits_here = int(mask_delayed.sum())
                if n_hits_here == 0:
                    continue

                times_in_delayed   = times_event[mask_delayed]
                charges_in_delayed = charge_event[mask_delayed]
                mpmt_in_delayed    = mpmt_event[mask_delayed]
                pmt_in_delayed     = pmt_event[mask_delayed]

                t_in = float(np.min(times_in_delayed))
                t_rms = functions_analysis.time_RMS_fun_time(times_in_delayed, t_in, window_neutron)

                vertex = functions_multilateration.run_multilateration_candidate(
                    times_in_delayed,
                    mpmt_in_delayed,
                    pmt_in_delayed,
                    sigma_t=1.0,
                )

                valid_thresholds.append(
                    {
                        "event_number": event_number,
                        "prompt_nhits": nhits_prompt,
                        "prompt_time": float(time_prompt),
                        "prompt_trms": trms_prompt,
                        "prompt_times": times_prompt,
                        "prompt_x": x_prompt,
                        "prompt_y": y_prompt,
                        "prompt_z": z_prompt,
                        "prompt_mpmt": mpmt_prompt,
                        "prompt_pmt": pmt_prompt,
                        "delayed_time": float(time_delayed),
                        "delayed_trms": float(t_rms),
                        "delayed_nhits": n_hits_here,
                        "delayed_x": vertex["x"],
                        "delayed_y": vertex["y"],
                        "delayed_z": vertex["z"],
                        "vertex_chi2": vertex.get("chi2_ndof", None),
                    }
                )

                print("  picked delayed", time_delayed, "-> lockout", time_delayed + window_neutron)
                prompt_lockout_until = time_delayed + window_neutron

            last_prompt = time_prompt

        return valid_thresholds

    dict_neutrons = {}
    for event in event_branch:
        if event == 6621:
            continue
        if event % 1000 == 0:
            print(f"Searching neutron coincidences on event {event}...")

        if event in threshold_times and len(threshold_times[event]) > 0:
            results = neutron_detection_event(
                event,
                times_branch_arg[event],
                charge_branch_arg[event],
                mpmt_branch_arg[event],
                pmt_branch_arg[event],
                threshold_times[event],
                window_sliding,
                window_neutron,
                threshold_inf,
                threshold_sup,
                window_prompt,
            )
            if len(results) != 0:
                dict_neutrons[event] = results

    return dict_neutrons


def neutron_detection_wMulti_alln(event_branch, times_branch_arg, charge_branch_arg, mpmt_branch_arg, pmt_branch_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup = np.inf, window_prompt = 100):
    #print('inside loop')
    def neutron_detection_event(event_number, times_branch_event_arg, charge_branch_event_arg, mpmt_branch_event_arg, pmt_branch_event_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt):
        #print('checking event', event_number)
        dict_neutrons_event = {}
        valid_thresholds = []
        last_prompt = None
        prompt_lockout_until = None
        print(f'Event {event_number} has {len(threshold_times)} prompt candidates')
        #for time_prompt in threshold_times:
        for time_prompt, nhits_prompt, trms_prompt, times_prompt, _ , mpmt_prompt, pmt_prompt, x_prompt, y_prompt, z_prompt  in threshold_times:
            print(f'Prompt time {time_prompt}')    
            #if last_prompt is not None and (time_prompt-last_prompt) < (window_sliding + window_prompt):
            #    continue

            #if prompt_lockout_until is not None and time_prompt < prompt_lockout_until:
            #    continue

            mask = (times_branch_event_arg >= time_prompt + window_prompt) & (times_branch_event_arg < time_prompt + window_prompt + window_sliding)
            #print('After isolating prompt')
            if mask.sum() == 0:
                continue

            neutron_candidates,  neutron_nhits = nHitsTimeWindow(times_branch_event_arg[mask], threshold_inf, window_neutron, 0, threshold_sup=threshold_sup)
            #print('Delayed found, checking now masks')
            #for time_delayed in neutron_candidates:

            ##########################################
            # Find delayed candidates in the post-prompt region
            if neutron_candidates is None or neutron_nhits is None:
                continue
            if len(neutron_candidates) == 0 or len(neutron_nhits) == 0:
                continue

            # If they can ever be mismatched, truncate safely (or skip)
            n_cand = min(len(neutron_candidates), len(neutron_nhits))
            neutron_candidates = neutron_candidates[:n_cand]
            neutron_nhits = neutron_nhits[:n_cand]

            # Loop over ALL delayed candidates (don’t just take argmax)
            for time_delayed, delayed_nhits in zip(neutron_candidates, neutron_nhits):

                # optional: apply your delayed nhits cuts here
                # if delayed_nhits < delayed_nhits_min or delayed_nhits > delayed_nhits_max:
                #     continue

                mask_delayed = (times_branch_event_arg >= time_delayed) & (times_branch_event_arg < time_delayed + window_sliding)
                if mask_delayed.sum() == 0:
                    continue

                times_in_delayed = times_branch_event_arg[mask_delayed]
                charges_in_delayed = charge_branch_event_arg[mask_delayed]
                mpmt_in_delayed = mpmt_branch_event_arg[mask_delayed]
                pmt_in_delayed = pmt_branch_event_arg[mask_delayed]

                t_in = np.min(times_in_delayed)
                t_rms = functions_analysis.time_RMS_fun_time(times_in_delayed, t_in, window_neutron)
                # avoid pathological event if you still need this special-case
                if event_number == 6621:
                    continue

                # Reconstruct delayed vertex
                vertex = functions_multilateration.run_multilateration_candidate(times_in_delayed, mpmt_in_delayed, pmt_in_delayed, sigma_t=1.0)

                # If recon fails, skip
                if not vertex.get("success", True):
                    continue

                # Save *this* delayed candidate (note: add candidate index/time)
                valid_thresholds.append(
                    {
                        "event_number": event_number,

                        "prompt_nhits": nhits_prompt,
                        "prompt_time": float(time_prompt),
                        "prompt_trms": trms_prompt,
                        "prompt_times": times_prompt,
                        "prompt_x": x_prompt,
                        "prompt_y": y_prompt,
                        "prompt_z": z_prompt,
                        "prompt_mpmt": mpmt_prompt,
                        "prompt_pmt": pmt_prompt,
                        "delayed_time": float(time_delayed),
                        "delayed_nhits": int(delayed_nhits),
                        "delayed_trms": float(t_rms),
                        "delayed_x": float(vertex["x"]),
                        "delayed_y": float(vertex["y"]),
                        "delayed_z": float(vertex["z"]),
                        "vertex_chi2": float(vertex.get("chi2_ndof", np.nan)),
                    }
                )

            print("  picked delayed", time_delayed, "-> lockout", time_delayed + window_neutron)

            #prompt_lockout_until = time_delayed + window_neutron

        return valid_thresholds
    
    dict_neutrons = {}
    for event in event_branch:
        if event == 6621:
            continue
        if event % 1000 == 0:
            print(f"Searching neutron coincidences on event {event}...")
        if event in threshold_times:
            if len(threshold_times[event]) == 0:
                continue 
            results = neutron_detection_event(event, times_branch_arg[event], charge_branch_arg[event], mpmt_branch_arg[event], pmt_branch_arg[event], threshold_times[event], window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt)
            if len(results) != 0:
                dict_neutrons[event] = results

    return dict_neutrons
        #else:
            #print(f"Event {event} has no threshold times, skipping neutron detection.")

def accidentals_wBonsai(event_branch, times_branch_event_arg, charge_branch_event_arg, mpmt_branch_event_arg, pmt_branch_event_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup=None, window_prompt = 100):

    def accidentals_wBonsai_evt(event_number, times_branch_event_arg, charge_branch_event_arg, mpmt_branch_event_arg, pmt_branch_event_arg, threshold_times, window_sliding, window_neutron, threshold_inf, threshold_sup = None, window_prompt = 10):
        valid_thresholds = []
        #last_prompt = None
        

        # Windows Method (1)###################################################################
        for time_prompt in threshold_times:
            mask = (times_branch_event_arg >= time_prompt + window_prompt) & (times_branch_event_arg < time_prompt + window_prompt + window_sliding)
            if mask.sum() == 0:
                continue
            neutron_candidates,  neutron_nhits = nHitsTimeWindow(times_branch_event_arg[mask], threshold_inf, window_neutron, window_neutron, threshold_sup=threshold_sup)
        #########################################################################################            

        # No Windows Method (2)########################################################################################
        #time_prompt = threshold_times
        #
        #mask = (times_branch_event_arg >= time_prompt + window_prompt) & (times_branch_event_arg < time_prompt + window_prompt + window_sliding)
        #
        #if mask.sum() != 0:
        #    neutron_candidates,  neutron_nhits = nHitsTimeWindow(times_branch_event_arg[mask], threshold_inf, window_neutron, window_neutron, threshold_sup=threshold_sup)
    
        #########################################################################################
            

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
                        'event_number': event_number,
                        'vp_prompt_time': float(time_prompt),
                        "vp_delayed_time": float(time_delayed),
                        "vp_delayed_nhits": delayed_nhits,
                        "vp_delayed_x": vertex["x"][0],
                        "vp_delayed_y": vertex["y"][0],
                        "vp_delayed_z": vertex["z"][0],
                    }
                )
                if len(valid_thresholds) % 10 == 0:
                    print(f"Found Cand Events {len(valid_thresholds)}...")
        #for time_prompt in threshold_times:

        #for time_prompt in threshold_times:    
        #    #if last_prompt is not None and (time_prompt-last_prompt) < (window_sliding + window_prompt):
        #    #    continue
        #
        #    mask = (times_branch_event_arg >= time_prompt + window_prompt) & (times_branch_event_arg < time_prompt + window_prompt + window_sliding)
        #    
        #    if mask.sum() == 0:
        #        continue
        #
        #    neutron_candidates,  neutron_nhits = nHitsTimeWindow(times_branch_event_arg[mask], threshold_inf, window_neutron, window_neutron, threshold_sup=threshold_sup)
        #    
        #    #for time_delayed in neutron_candidates:
        #    for i in range(len(neutron_candidates)):
        #        time_delayed = neutron_candidates[i]
        #        delayed_nhits = neutron_nhits[i]
        #        mask_delayed = (times_branch_event_arg >= time_delayed) & (times_branch_event_arg < time_delayed + window_sliding)
        #        #neutron_nhits = mask_delayed.sum()
        #
        #        times_in_delayed = times_branch_event_arg[mask_delayed]
        #        charges_in_delayed = charge_branch_event_arg[mask_delayed]
        #        mpmt_in_delayed = mpmt_branch_event_arg[mask_delayed]
        #        pmt_in_delayed = pmt_branch_event_arg[mask_delayed]
        #        
        #        vertex = functions_bonsai.run_BONSAI_candidate(times_in_delayed, charges_in_delayed, mpmt_in_delayed, pmt_in_delayed)
        #
        #        valid_thresholds.append(
        #            {
        #                'event_number': event_number,
        #                'vp_prompt_time': float(time_prompt),
        #                "vp_delayed_time": float(time_delayed),
        #                "vp_delayed_nhits": delayed_nhits,
        #                "vp_delayed_x": vertex["x"][0],
        #                "vp_delayed_y": vertex["y"][0],
        #                "vp_delayed_z": vertex["z"][0],
        #            }
        #        )
        #        if len(valid_thresholds) % 10 == 0:
        #            print(f"Found Cand Events {len(valid_thresholds)}...")
        return valid_thresholds

    #dict_neutrons = {}
    #for event in event_branch:
    #    if event % 1000 == 0:
    #        print(f"Searching accidental coincidences on event {event}...")
    #    if event in threshold_times:
    #        if len(threshold_times[event]) == 0:
    #            continue 
    #        results = accidentals_wBonsai_evt(event, times_branch_event_arg[event], charge_branch_event_arg[event], mpmt_branch_event_arg[event], pmt_branch_event_arg[event], threshold_times[event], window_sliding, window_neutron, threshold_inf, threshold_sup, window_prompt)
    #        if len(results) != 0:
    #            dict_neutrons[event] = results

    dict_neutrons = {}

    for event, times in threshold_times.items():
        if event % 1000 == 0:
            print(f"Searching accidental coincidences on event {event}...")
        results = accidentals_wBonsai_evt(
            event,
            times_branch_event_arg[event],
            charge_branch_event_arg[event],
            mpmt_branch_event_arg[event],
            pmt_branch_event_arg[event],
            times,
            window_sliding,
            window_neutron,
            threshold_inf,
            threshold_sup,
            window_prompt
        )
        if results:
            dict_neutrons[event] = results
            
    return dict_neutrons

