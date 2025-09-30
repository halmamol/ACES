print("Iniciando el script de detección de neutrones...")
print("Importando librerias necesarias...")

#import uproot
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import functions_spills
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
output_path = "/scratch/halmazan/WCTE/files/data/"

prompt_window = 1500  # Window for prompt candidates
prompt_dead_time = 200  # Death time for prompt candidates
prompt_t_rms_min = 200
prompt_t_rms_max = 400
prompt_nhits_min = 150
prompt_nhits_max = 300
coincidence_window = 150000  # Window for coincidence search
delayed_window = 100  # Window for delayed candidates
delayed_nhits_min = 10  # Minimum number of hits for delayed candidates
delayed_nhits_max = 30  # Maximum number of hits for delayed candidates

print("Opening: ", f'{output_path}filtered_files/filtered_file_{run_number}.pkl')
print("Output saved in: ", f'{output_path}AmBeCandidates/neutron_candidates_{run_number}.csv')

with open(f'{output_path}filtered_files/filtered_file_{run_number}.pkl', 'rb') as f:
#    data = pickle.load(f)
    data = Numpy2to1Unpickler(f).load()

times_vals = data["times_TOF"]["values"]
times_offs = data["times_TOF"]["offsets"]
charge_vals = data["charge"]["values"]
charge_offs = data["charge"]["offsets"]
mpmt_vals = data["mpmt_id"]["values"]
mpmt_offs = data["mpmt_id"]["offsets"]
pmt_vals = data["pmt_id"]["values"]
pmt_offs = data["pmt_id"]["offsets"]

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
threshold_times_prompt = functions_spills.prompt_candidates_wBonsai(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, prompt_window, prompt_dead_time, prompt_nhits_min, prompt_nhits_max)

#for event in threshold_times_prompt:
#    t_in_list = threshold_times_prompt[event]  # List of (time_prompt, n_hits)
#    threshold_times_prompt[event] = [
#        (t_in, n_hits, functions_analysis.time_RMS_fun_time(times_per_event[event], t_in, prompt_window))
#        for (t_in, n_hits) in t_in_list
#        if prompt_t_rms_min <= functions_analysis.time_RMS_fun_time(times_per_event[event], t_in, prompt_window) <= prompt_t_rms_max
#    ]

for event, candidates in threshold_times_prompt.items():
    times_event = times_per_event[event]
    filtered = []
    for cand in candidates:
        if isinstance(cand, dict):
            t_in = cand["time"]
            n_hits = cand["n_hits"]
            charge_arr = cand["charge"]
            mpmt_arr = cand["mpmt_id"]
            pmt_arr = cand["pmt_id"]
            x_arr = cand["vertex_x"]  
            y_arr = cand["vertex_y"]
            z_arr = cand["vertex_z"]
        else:
            # If you only ever have tuples, you can't include charge/mpmt here
            t_in, n_hits = cand
            charge_arr = None
            mpmt_arr = None
            pmt_arr = None

        t_rms = functions_analysis.time_RMS_fun_time(times_event, t_in, prompt_window)

        if prompt_t_rms_min <= t_rms <= prompt_t_rms_max:
            filtered.append((t_in, n_hits, t_rms, charge_arr, mpmt_arr, pmt_arr, x_arr, y_arr, z_arr))

    threshold_times_prompt[event] = filtered
print("Prompt candidates found in run.")


# Neutron detection ###########################################################################################################

print("Searching for neutron events...")
neutron_dict = functions_spills.neutron_detection_wBonsai(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, threshold_times_prompt, coincidence_window, delayed_window, delayed_nhits_min, delayed_nhits_max, prompt_window)

print("Prompt candidates", sum(len(v) for v in threshold_times_prompt.values()))
print("Neutron candidates", sum(len(v) for v in neutron_dict.values()))


# Save neutron candidates to CSV files ###########################################################################################################

print("Saving candidate neutron events on CSV...")

"""
neutron_candidates = []

for event_number, times in neutron_dict.items():
    
    prompt_candidates = threshold_times_prompt.get(event_number, [])  # List of (time_prompt, n_hits, prompt_trms)

    prompt_nhits_map = {
    float(t): (
        float(n[0]) if isinstance(n, list) else float(n),
        float(trms[0]) if isinstance(trms, list) else float(trms),
        float(x[0]) if isinstance(x, list) else float(x),
        float(y[0]) if isinstance(y, list) else float(y),
        float(z[0]) if isinstance(z, list) else float(z)
    )
    for t, n, trms, _, _, x, y, z in prompt_candidates
    }
    for start_time, neutron_pairs in times.items():
        prompt_data = prompt_nhits_map.get(float(start_time), (None, None))
        prompt_nhits_val, prompt_trms_val, prompt_x_val, prompt_y_val, prompt_z_val = prompt_data
        for delayed_time, delayed_nhits, delayed_x, delayed_y, delayed_z in neutron_pairs:
            neutron_candidates.append({
                'event_number': event_number,
                'prompt_time': float(start_time),
                'prompt_nhits': prompt_nhits_val,
                'prompt_trms': prompt_trms_val,
                'prompt_x': prompt_x_val,
                'prompt_y': prompt_y_val,
                'prompt_z': prompt_z_val,
                'delayed_time': float(delayed_time),
                'delayed_nhits': float(delayed_nhits) if not isinstance(delayed_nhits, list) else float(delayed_nhits[0]),
                'delayed_x': float(delayed_x) if not isinstance(delayed_x, list) else float(delayed_x[0]),
                'delayed_y': float(delayed_y) if not isinstance(delayed_y, list) else float(delayed_y[0]),
                'delayed_z': float(delayed_z) if not isinstance(delayed_z, list) else float(delayed_z[0])
            })
"""
#df_neutron_candidates = pd.DataFrame(neutron_dict)
#df_neutron_candidates = pd.DataFrame(neutron_candidates)
df_neutron_candidates = pd.concat(
    [pd.DataFrame(recs).assign(event_number=int(event)) for event, recs in neutron_dict.items()],
    ignore_index=True
)

df_neutron_candidates.to_csv(f'{output_path}/AmBeCandidates/neutron_candidates_{run_number}.csv', index=False)


print("CSV files saved.")
print("End of code.")
