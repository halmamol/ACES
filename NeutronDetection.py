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
run_number = "2385"  # Run number
output_path = "/scratch/halmazan/WCTE/files/data/"

prompt_window = 1500  # Window for prompt candidates
prompt_dead_time = 200  # Death time for prompt candidates
prompt_t_rms_min = 200 # Minimum RMS time for prompt candidates
prompt_t_rms_max = 400 # Maximum RMS time for prompt candidates
prompt_nhits_min = 150 # Minimum number of hits for prompt candidates
prompt_nhits_max = 300 # Maximum number of hits for prompt candidates
coincidence_window = 150000  # Window for coincidence search
delayed_window = 100  # Window for delayed candidates
delayed_nhits_min = 10  # Minimum number of hits for delayed candidates
delayed_nhits_max = 30  # Maximum number of hits for delayed candidates

print("Opening: ", f'{output_path}filtered_files/filtered_file_{run_number}.pkl')
print("Output saved in: ", f'{output_path}AmBeCandidates/neutron_candidates_{run_number}.csv')

with open(f'{output_path}filtered_files/filtered_file_daqtime_{run_number}.pkl', 'rb') as f:
#    data = pickle.load(f)
    data = Numpy2to1Unpickler(f).load()

times_vals = data["times_TOF"]["values"]
times_offs = data["times_TOF"]["offsets"]
times_daq_vals = data["times_daq_TOF"]["values"]
times_daq_offs = data["times_daq_TOF"]["offsets"]
charge_vals = data["charge"]["values"]
charge_offs = data["charge"]["offsets"]
mpmt_vals = data["mpmt_id"]["values"]
mpmt_offs = data["mpmt_id"]["offsets"]
pmt_vals = data["pmt_id"]["values"]
pmt_offs = data["pmt_id"]["offsets"]

times_per_event  = unfold(times_vals,  times_offs)
times_daq_per_event  = unfold(times_daq_vals,  times_daq_offs)
charge_per_event = unfold(charge_vals, charge_offs)
mpmt_per_event   = unfold(mpmt_vals,   mpmt_offs)
pmt_per_event   = unfold(pmt_vals,   pmt_offs)

print("Filtered data loaded.")
N_events = len(times_per_event)


print(f"Number of events in run {run_number}", N_events)
print(f"Total suppousdely lenght of run ", {N_events*270000})


# Concatenate all per-event arrays
all_times = np.concatenate(times_daq_per_event)
all_charges = np.concatenate(charge_per_event)
all_mpmts = np.concatenate(mpmt_per_event)
all_pmt = np.concatenate(pmt_per_event)

# To add offsets to each time, build an array of per-hit offsets
# For each event, repeat the offset for the number of hits in that event
event_offsets = np.concatenate([
    np.full(len(ev), offset) for ev, offset in zip(times_per_event, times_offs)
])

# Now just add the offsets
times_run = all_times + event_offsets
charge_run = all_charges  # Already flat, nothing more is needed
mpmt_run = all_mpmts  # Already flat, nothing more is needed
pmt_run = all_pmt  # Already flat, nothing more is needed

print(f"Total suppousdely lenght of run ", max(times_run))
print("Event times concatenated.")

#times_run = times_run[:10000000]
#charge_run = charge_run[:10000000]
#mpmt_run = mpmt_run[:10000000]
#pmt_run = pmt_run[:10000000]



event_number_branch = np.arange(0, N_events, 1)

# Prompt candidates detection ###########################################################################################################

print("Searching prompt candidate events...")
threshold_times_prompt = functions_coincidence.prompt_candidates_wBonsai(times_run, charge_run, mpmt_run, pmt_run, prompt_window, prompt_dead_time, prompt_nhits_min, prompt_nhits_max)
print(f"Total prompt candidate events found: {len(threshold_times_prompt)}")
#print(f"Total promp candidates found in run: {sum(len(v) for v in threshold_times_prompt.values())}")

filtered = []
for cand in threshold_times_prompt:
    
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

    t_rms = functions_analysis.time_RMS_fun_time(times_run, t_in, prompt_window)
    #print(f'Value of t_rms: {t_rms}')
    if prompt_t_rms_min <= t_rms <= prompt_t_rms_max:
        #print('Candidate passed the RMS cut')
        filtered.append((t_in, n_hits, t_rms, charge_arr, mpmt_arr, pmt_arr, x_arr, y_arr, z_arr))
    
threshold_times_prompt = filtered
print(f"Total promp candidates after RMS cut: {len(threshold_times_prompt)}")

"""""""""""""""
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
"""""""""""""""
print("Prompt candidates found in run.")


# Neutron detection ###########################################################################################################

print("Searching for neutron events...")
neutron_dict = functions_coincidence.neutron_detection_wBonsai(times_run, charge_run, mpmt_run, pmt_run, threshold_times_prompt, coincidence_window, delayed_window, delayed_nhits_min, delayed_nhits_max, prompt_window)

print("Prompt candidates", len(threshold_times_prompt))
print("Neutron candidates", len(neutron_dict))


# Save neutron candidates to CSV files ###########################################################################################################

print("Saving candidate neutron events on CSV...")

df_neutron_candidates = pd.DataFrame(neutron_dict)
#df_neutron_candidates = pd.DataFrame(neutron_candidates)
#df_neutron_candidates = pd.concat(
#    [pd.DataFrame(recs).assign(event_number=int(event)) for event, recs in neutron_dict.items()],
#    ignore_index=True
#)

df_neutron_candidates.to_csv(f'{output_path}/AmBeCandidates/neutron_candidates_{run_number}_timestest.csv', index=False)


print("CSV files saved.")
print("End of code.")
