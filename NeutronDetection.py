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
run_number = "2384"  # Run number
output_path = "/data/halmazan/WCTE/data/"

input_file = f'{output_path}filtered_files/filtered_file_{run_number}.pkl'
output_file = f'/scratch/halmazan/WCTE/files/AmBeCandidates/neutron_candidates_{run_number}_promptvar_neutronsel_test2_mult.csv'

prompt_window = 1500  # Window for prompt candidates
prompt_dead_time = 200  # Death time for prompt candidates
prompt_t_rms_min = 200 # Minimum RMS time for prompt candidates
prompt_t_rms_max = 500 # Maximum RMS time for prompt candidates #400
prompt_nhits_min = 50 # Minimum number of hits for prompt candidates #150
prompt_nhits_max = 300 # Maximum number of hits for prompt candidates
coincidence_window = 150000  # Window for coincidence search
delayed_window = 40  # Window for delayed candidates
delayed_nhits_min = 10  # Minimum number of hits for delayed candidates
delayed_nhits_max = 30  # Maximum number of hits for delayed candidates


print("Opening: ", input_file)
print("Output saved in: ", output_file)

with open(input_file, 'rb') as f:
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
#threshold_times_prompt = functions_coincidence.prompt_candidates_wScintillationLikelihood(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, prompt_window, prompt_dead_time, prompt_nhits_min, prompt_nhits_max)
threshold_times_prompt = functions_coincidence.prompt_candidates_wBonsai(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, prompt_window, prompt_dead_time, prompt_nhits_min, prompt_nhits_max)

gamma_info_prompt = {}  # separate dict

for event, candidates in threshold_times_prompt.items():
    times_event = times_per_event[event]
    filtered = []
    gamma_info = []
    for cand in candidates:
        if isinstance(cand, dict):
            t_in = cand["time"]
            n_hits = cand["n_hits"]
            charge_arr = cand.get("charge")
            mpmt_arr = cand.get("mpmt_id")
            pmt_arr = cand.get("pmt_id")
            x_arr = cand.get("vertex_x")
            y_arr = cand.get("vertex_y")
            z_arr = cand.get("vertex_z")
        else:
            t_in, n_hits = cand
            charge_arr = mpmt_arr = pmt_arr = None
            x_arr = y_arr = z_arr = None  # important

        t_rms = functions_analysis.time_RMS_fun_time(times_event, t_in, prompt_window)
        times_prompt = times_event[(times_event >= t_in) & (times_event <= t_in + prompt_window)]

        if prompt_t_rms_min <= t_rms <= prompt_t_rms_max:
            filtered.append((t_in, n_hits, t_rms, times_prompt.tolist(),
                             charge_arr, mpmt_arr, pmt_arr, x_arr, y_arr, z_arr))

            gamma_info.append({
                "prompt_ti": t_in,
                "prompt_nhits": n_hits,
                "prompt_trms": t_rms,
                "prompt_times": times_prompt.tolist(),
                "prompt_charge": charge_arr,
            })

    threshold_times_prompt[event] = filtered
    gamma_info_prompt[event] = gamma_info
    
print("Prompt candidates found in run.")

df_gamma_candidates = pd.concat(
    [pd.DataFrame(recs).assign(event_number=int(event)) for event, recs in gamma_info_prompt.items()],
    ignore_index=True
)

df_gamma_candidates.to_csv(f'/scratch/halmazan/WCTE/files/AmBeCandidates/gamma_candidates_{run_number}_promptvar_neutronsel_test2_mult.csv', index=False)

# Neutron detection ###########################################################################################################
prompt_candidates = sum(len(v) for v in threshold_times_prompt.values())
print("Prompt candidates before neutron search", prompt_candidates)
print("Searching for neutron events...")
#neutron_dict = functions_coincidence.neutron_detection_wBonsai(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, threshold_times_prompt, coincidence_window, delayed_window, delayed_nhits_min, delayed_nhits_max, prompt_window)
neutron_dict = functions_coincidence.neutron_detection_wMulti(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, threshold_times_prompt, coincidence_window, delayed_window, delayed_nhits_min, delayed_nhits_max, prompt_window)

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
df_summary.to_csv(f'/scratch/halmazan/WCTE/files/AmBeCandidates/summary_{run_number}_promptvar_neutronsel_test2_mult.csv', index=False)


print("CSV files saved.")
print("End of code.")
