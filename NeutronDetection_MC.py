print("Iniciando el script de detección de neutrones...")
print("Importando librerias necesarias...")

#import uproot
import numpy as np
import pandas as pd
import functions_coincidence
import functions_analysis
import glob
import os
import argparse
import pickle

class Numpy2to1Unpickler(pickle.Unpickler):
    MAP = {
        "numpy._core": "numpy.core",
        # If needed, add more granular mappings here:
        # "numpy._core.multiarray": "numpy.core.multiarray",
        # "numpy._core.overrides": "numpy.core.overrides",
        # "numpy._core._multiarray_umath": "numpy.core._multiarray_umath",
    }
    def find_class(self, module, name):
        for old, new in self.MAP.items():
            if module == old or module.startswith(old + "."):
                module = module.replace(old, new, 1)
                break
        return super().find_class(module, name)

def unfold(values: np.ndarray, offsets: np.ndarray):
    # Returns a list of per-event arrays using the offsets
    return [values[int(offsets[i]):int(offsets[i+1])] for i in range(len(offsets) - 1)]

# Crear el parser
parser = argparse.ArgumentParser(description="Window size")

# Agregar argumento opcional con valor por defecto
parser.add_argument(
    "--window_size",
    type = int, 
    help="Window size neutron",
    default=100
)

# Arguments for Analysis 
run_number = "2386"  # Run number
output_path = "/scratch/halmazan/WCTE/files/"

input_file = f'{output_path}filtered_files/filtered_file_MC_{run_number}_digidata_noTOF.pkl'
title = f'_{run_number}_prompt50threshold_prompt1000window_neutronsel_mult_wtfilter_w10ns_mtrms10_window100mus'
output_file = f'/scratch/halmazan/WCTE/files/AmBeCandidates/neutron_candidates_MC{title}.csv'

prompt_window = 1000  # Window for prompt candidates
prompt_dead_time = 200  # Death time for prompt candidates
prompt_t_rms_min = 200 # Minimum RMS time for prompt candidates
prompt_t_rms_max = 500 # Maximum RMS time for prompt candidates #400
prompt_nhits_min = 50 # Minimum number of hits for prompt candidates #150
prompt_nhits_max = 300 # Maximum number of hits for prompt candidates
afterpulse_time = 5000
coincidence_window = 150000  # Window for coincidence search
delayed_window = 10  # Window for delayed candidates
delayed_nhits_min = 10  # Minimum number of hits for delayed candidates
delayed_nhits_max = 50  # Maximum number of hits for delayed candidates

source_pos = []
sources_pos = {
    '2386': [0.0, +30.5, 0.0],
    '2387': [0.0, -30.5, 0.0],
    '2388': [0.0, -58.5, 0.0],
    '2389': [-47.64, +30.5, 29.17],
    '2390': [-47.64, -30.5, 29.17], 
}
if run_number in sources_pos.keys():
    source_pos = sources_pos[run_number]
else:
    print('Run not on list')

print("Opening: ", input_file)
print("Output saved in: ", output_file)

with open(input_file, 'rb') as f:
#    data = pickle.load(f)
    data = Numpy2to1Unpickler(f).load()

times_vals_TOF = data["times_TOF"]["values"]
times_offs_TOF = data["times_TOF"]["offsets"]
times_vals = data["times"]["values"]
times_offs = data["times"]["offsets"]
charge_vals = data["charge"]["values"]
charge_offs = data["charge"]["offsets"]
mpmt_vals = data["mpmt_id"]["values"]
mpmt_offs = data["mpmt_id"]["offsets"]
pmt_vals = data["pmt_id"]["values"]
pmt_offs = data["pmt_id"]["offsets"]

times_per_event_TOF  = unfold(times_vals_TOF,  times_offs_TOF)
times_per_event  = unfold(times_vals,  times_offs)
charge_per_event = unfold(charge_vals, charge_offs)
mpmt_per_event   = unfold(mpmt_vals,   mpmt_offs)
pmt_per_event   = unfold(pmt_vals,   pmt_offs)

print("Filtered data loaded.")
N_events = len(times_per_event)


print(f"Number of events in run {run_number}", N_events)


event_number_branch = np.arange(0, N_events, 1)

# Prompt candidates detection ###########################################################################################################

print("Searching prompt candidate events...")
threshold_times_prompt = functions_coincidence.prompt_candidates_wBonsai(event_number_branch, times_per_event_TOF, charge_per_event, mpmt_per_event, pmt_per_event, prompt_window, prompt_dead_time, prompt_nhits_min, prompt_nhits_max)
#threshold_times_prompt = functions_coincidence.prompt_candidates(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, prompt_window, prompt_dead_time, prompt_nhits_min, prompt_nhits_max)
#threshold_times_prompt = functions_coincidence.prompt_candidates(event_number_branch, times_per_event, prompt_window, prompt_dead_time, prompt_nhits_min, prompt_nhits_max)
gamma_info_prompt = {}  # separate dict

for event, candidates in threshold_times_prompt.items():
    times_event_tof = times_per_event_TOF[event]
    times_event_raw = times_per_event[event]
    charge_event = charge_per_event[event]
    mpmt_event = mpmt_per_event[event]
    pmt_event = pmt_per_event[event]

    filtered = []
    gamma_info = []

    for cand in candidates:
        if isinstance(cand, dict):
            t_in = cand["time"]
            n_hits = cand["n_hits"]
            x_arr = cand.get("vertex_x")
            y_arr = cand.get("vertex_y")
            z_arr = cand.get("vertex_z")
        else:
            t_in, n_hits = cand
            x_arr = y_arr = z_arr = None

        mask_prompt = (
            (times_event_tof >= t_in) &
            (times_event_tof <= t_in + prompt_window)
        )

        times_prompt_tof = times_event_tof[mask_prompt]
        times_prompt_raw = times_event_raw[mask_prompt]
        charge_prompt = charge_event[mask_prompt]
        mpmt_prompt = mpmt_event[mask_prompt]
        pmt_prompt = pmt_event[mask_prompt]

        t_rms, t_mean = functions_analysis.time_RMS_fun_time(
            times_event_tof, t_in, prompt_window
        )

        if prompt_t_rms_min <= t_rms <= prompt_t_rms_max:
            filtered.append((
                t_in,
                n_hits,
                t_rms,
                t_mean,
                times_prompt_tof.tolist(),
                times_prompt_raw.tolist(),
                charge_prompt.tolist(),
                mpmt_prompt.tolist(),
                pmt_prompt.tolist(),
                x_arr,
                y_arr,
                z_arr,
            ))

            gamma_info.append({
                "prompt_ti": t_in,
                "prompt_nhits": n_hits,
                "prompt_trms": t_rms,
                "prompt_tmean": t_mean,
                "prompt_times_tof": times_prompt_tof.tolist(),
                "prompt_times_raw": times_prompt_raw.tolist(),
                "prompt_charge": charge_prompt.tolist(),
                "prompt_mpmt": mpmt_prompt.tolist(),
                "prompt_pmt": pmt_prompt.tolist(),
                "prompt_vertex_x": x_arr,
                "prompt_vertex_y": y_arr,
                "prompt_vertex_z": z_arr,
            })

    threshold_times_prompt[event] = filtered
    gamma_info_prompt[event] = gamma_info
    
print("Prompt candidates found in run.")

df_gamma_candidates = pd.concat(
    [pd.DataFrame(recs).assign(event_number=int(event)) for event, recs in gamma_info_prompt.items()],
    ignore_index=True
)

df_gamma_candidates.to_csv(f'/scratch/halmazan/WCTE/files/AmBeCandidates/gamma_candidates_MC{title}.csv', index=False)

# Neutron detection ###########################################################################################################
prompt_candidates = sum(len(v) for v in threshold_times_prompt.values())
print("Prompt candidates before neutron search", prompt_candidates)
print("Searching for neutron events...")
#neutron_dict = functions_coincidence.neutron_detection_wBonsai(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, threshold_times_prompt, coincidence_window, delayed_window, delayed_nhits_min, delayed_nhits_max, prompt_window)
#neutron_dict = functions_coincidence.neutron_detection_wMulti_alln(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, threshold_times_prompt, coincidence_window, delayed_window, delayed_nhits_min, delayed_nhits_max, prompt_window)

#neutron_dict = functions_coincidence.neutron_detection_wMulti_wtcut(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, threshold_times_prompt, coincidence_window, delayed_window, delayed_nhits_min, delayed_nhits_max, afterpulse_time)
#neutron_dict = functions_coincidence.neutron_detection_wMulti_clean(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, threshold_times_prompt, coincidence_window, delayed_window, delayed_nhits_min, delayed_nhits_max, afterpulse_time)

neutron_dict = functions_coincidence.neutron_detection_wMulti_wtcut(
    event_branch=event_number_branch,
    times_branch_tof_arg=times_per_event_TOF,
    times_branch_raw_arg=times_per_event,
    charge_branch_arg=charge_per_event,
    mpmt_branch_arg=mpmt_per_event,
    pmt_branch_arg=pmt_per_event,
    threshold_times=threshold_times_prompt,
    window_sliding=coincidence_window,
    window_neutron=delayed_window,
    threshold_inf=delayed_nhits_min,
    threshold_sup=delayed_nhits_max,
    window_prompt=afterpulse_time,
)

prompt_candidates_wpair = sum(len(v) for v in threshold_times_prompt.values())
delayed_candidates =  sum(len(v) for v in neutron_dict.values())
print("Prompt candidates", prompt_candidates_wpair)
print("Neutron candidates", delayed_candidates)


# Save neutron candidates to CSV files ###########################################################################################################

print("Saving candidate neutron events on CSV...")

#df_neutron_candidates = pd.DataFrame(neutron_dict)
#df_neutron_candidates = pd.DataFrame(neutron_candidates)
df_neutron_candidates = pd.concat(
    [pd.DataFrame(recs).assign(event_number=int(event)) for event, recs in neutron_dict.items()],
    ignore_index=True
)

df_neutron_candidates.to_csv(f'{output_file}', index=False)

# Summary dataframe  ###########################################################################################################
df_summary = pd.DataFrame([{
    "run_number": run_number,
    "prompt_candidates": prompt_candidates,
    "prompt_candidates_wpair": prompt_candidates_wpair,
    "delayed_candidates": delayed_candidates,
}])

# Save summary to its own CSV
df_summary.to_csv(f'/scratch/halmazan/WCTE/files/AmBeCandidates/summary_MC{title}.csv', index=False)


print("CSV files saved.")
print("End of code.")
