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
import random

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

# Parameters for Candidate Detection
prompt_window = 1500  # Window for prompt candidates
prompt_dead_time = 200  # Death time for prompt candidates
prompt_t_rms_min = 200 # Minimum RMS time for prompt candidates
prompt_t_rms_max = 400 # Maximum RMS time for prompt candidates
prompt_nhits_min = 150 # Minimum number of hits for prompt candidates
prompt_nhits_max = 300 # Maximum number of hits for prompt candidates
coincidence_window = 150000  # Window for coincidence search
delayed_window = 40  # Window for delayed candidates
delayed_nhits_min = 10  # Minimum number of hits for delayed candidates
delayed_nhits_max = 30  # Maximum number of hits for delayed candidates
### Accidentals
nBG = 10 # Number of accidental windows to search
deltat_vp = 50000 #ms Time difference for accidental search

print("Opening Files...")

with open(f'{output_path}filtered_files/filtered_file_{run_number}.pkl', 'rb') as f:
#    data = pickle.load(f)
    data = Numpy2to1Unpickler(f).load()

neutron_df = pd.read_csv(f'{output_path}AmBeCandidates/neutron_candidates_{run_number}_test.csv')

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
print(f"Total suppousdely lenght of run ", {N_events*270000})


# Concatenate all per-event arrays
all_charges = np.concatenate(charge_per_event)
all_mpmts = np.concatenate(mpmt_per_event)
all_pmt = np.concatenate(pmt_per_event)

# To add offsets to each time, build an array of per-hit offsets
# For each event, repeat the offset for the number of hits in that event
event_offsets = np.concatenate([
    np.full(len(ev), offset) for ev, offset in zip(times_per_event, times_offs)
])


event_number_branch = np.arange(0, N_events, 1)


# Running Accidental Search ###########################################################################################################
neutron_df = neutron_df.sort_values('prompt_time').reset_index(drop=True)

# compute differences between consecutive prompt_time values
#neutron_df['dt'] = neutron_df['prompt_time'].diff()

# Windows Method (1)###################################################################
neutron_df['dt'] = 270000 - neutron_df['prompt_time']


mask = neutron_df['dt'] < (deltat_vp * nBG)  

vp_times_by_event = {}
for event, group in neutron_df.loc[mask].groupby('event_number'):
    t0 = group['prompt_time'].iloc[0]   # or use .min() if preferred
    vp_times_by_event[event] = [t0 + j * deltat_vp for j in range(nBG)]
########################################################################################

# No Windows Method (2)########################################################################################
#neutron_events = set(neutron_df['event_number'].values)
#all_events = set(range(N_events))

#missing_events = list(all_events - neutron_events)
#acc_windows = len(neutron_events)*nBG
#sampled_missing = random.sample(missing_events, acc_windows)

#vp_times_by_event = {}
#for event in sampled_missing:
#    vp_times_by_event[event] = np.random.uniform(0, 120000)
########################################################################################
    
print(f"Total vp windows found: {len(vp_times_by_event)}")
print("Running accidental search...")
vp_neutron_dict = functions_coincidence.accidentals_wBonsai(event_number_branch, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, vp_times_by_event, coincidence_window, delayed_window, delayed_nhits_min, delayed_nhits_max, prompt_window)
#vp_neutron_dict = functions_coincidence.accidentals_wBonsai(sampled_missing, times_per_event, charge_per_event, mpmt_per_event, pmt_per_event, vp_times_by_event, coincidence_window, delayed_window, delayed_nhits_min, delayed_nhits_max, prompt_window)

print("VP Neutron candidates", len(vp_neutron_dict))

# Save accidental candidates to CSV files ###########################################################################################################

print("Saving accidental events on CSV...")

#df_vp_neutron_candidates = pd.DataFrame({
#    'event': list(vp_neutron_dict.keys()),
#    'results': list(vp_neutron_dict.values())
#})
df_vp_neutron_candidates = pd.concat(
    [pd.DataFrame(recs).assign(event_number=int(event)) for event, recs in vp_neutron_dict.items()],
    ignore_index=True
)
#df_vp_neutron_candidates = pd.DataFrame(vp_neutron_dict)
#df_neutron_candidates = pd.DataFrame(neutron_candidates)
#df_neutron_candidates = pd.concat(
#    [pd.DataFrame(recs).assign(event_number=int(event)) for event, recs in neutron_dict.items()],
#    ignore_index=True
#)

df_vp_neutron_candidates.to_csv(f'{output_path}/AmBeCandidates/accidental_candidates_{run_number}_40ns.csv', index=False)
#df_vp_neutron_candidates.to_csv(f'{output_path}/AmBeCandidates/accidental_candidates_noWindow_{run_number}.csv', index=False)


print("CSV files saved.")
print("End of code.")