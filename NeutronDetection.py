print("Iniciando el script de detección de neutrones...")
print("Importando librerias necesarias...")

import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import functions_spills
import functions_analysis
import glob
import os
import argparse
import pickle

from matplotlib import rcParams
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'
rcParams['figure.figsize'] = [10, 8]
rcParams['font.size'] = 18

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
output_path = "/scratch/halmazan/WCTE/files/data/"

prompt_window = 1500  # Window for prompt candidates
prompt_dead_time = 200  # Death time for prompt candidates
prompt_t_rms_min = 200
prompt_t_rms_max = 400
prompt_nhits_min = 150
prompt_nhits_max = 300
coincidence_window = 150000  # Window for coincidence search
delayed_window = 100  # Window for delayed candidates
delayed_nhits_min = 5  # Minimum number of hits for delayed candidates
delayed_nhits_max = 50  # Maximum number of hits for delayed candidates

print("Opening: ", f'{output_path}filtered_files/filtered_file_{run_number}.pkl')
print("Output saved in: ", f'{output_path}AmBeCandidates/neutron_candidates_{run_number}.csv')

with open(f'{output_path}filtered_files/filtered_file_{run_number}.pkl', 'rb') as f:
    valores_read, indices_read = pickle.load(f)

#with open('Filtered_data/datos_filtrados_sig.pkl', 'rb') as f:
#    valores_read_sig, indices_read_sig = pickle.load(f)

times_branch_filtered = functions_analysis.a_lista_de_arrays(valores_read, indices_read)
#times_branch_filtered_sig = functions_analysis.a_lista_de_arrays(valores_read_sig, indices_read_sig)

print("Filtered data loaded.")
N_events = len(times_branch_filtered)
#N_events_sig = len(times_branch_filtered_sig)

print(f"Number of events in run {run_number}", N_events)
#print("Numero de eventos señal", N_events_sig)

event_number_branch = np.arange(0, N_events, 1)
#event_number_branch_sig = np.arange(0, N_events_sig, 1)
# Prompt candidates detection ###########################################################################################################

print("Searching prompt candidate events...")
threshold_times_prompt = functions_spills.prompt_candidates(event_number_branch, times_branch_filtered, prompt_window, prompt_dead_time, prompt_nhits_min, prompt_nhits_max)

#for event in threshold_times_prompt:
#    t_in_list = threshold_times_prompt[event]
#    threshold_times_prompt[event] = [
#        t_in for t_in in t_in_list
#        if prompt_t_rms_min <= functions_analysis.time_RMS_fun_time(times_branch_filtered[event], t_in, prompt_window) <= prompt_t_rms_max]

#for event in threshold_times_prompt:
#    t_in_list = threshold_times_prompt[event]  # List of (time_prompt, n_hits)
#    threshold_times_prompt[event] = [
#        (t_in, n_hits)
#        for (t_in, n_hits) in t_in_list
#        if prompt_t_rms_min <= functions_analysis.time_RMS_fun_time(times_branch_filtered[event], t_in, prompt_window) <= prompt_t_rms_max
#    ]

for event in threshold_times_prompt:
    t_in_list = threshold_times_prompt[event]  # List of (time_prompt, n_hits)
    threshold_times_prompt[event] = [
        (t_in, n_hits, functions_analysis.time_RMS_fun_time(times_branch_filtered[event], t_in, prompt_window))
        for (t_in, n_hits) in t_in_list
        if prompt_t_rms_min <= functions_analysis.time_RMS_fun_time(times_branch_filtered[event], t_in, prompt_window) <= prompt_t_rms_max
    ]
print("Prompt candidates found in run.")

#threshold_times_prompt_sig = functions_spills.prompt_candidates(event_number_branch_sig, times_branch_filtered_sig, 1500, 200, 150, 300)

#for event in threshold_times_prompt_sig:
#    t_in_list = threshold_times_prompt_sig[event]
#    threshold_times_prompt_sig[event] = [
#        t_in for t_in in t_in_list
#        if 200 <= functions_analysis.time_RMS_fun_time(times_branch_filtered_sig[event], t_in, 1500) <= 400]

#print("Candidatos prompt encontrados en data run.")

# Neutron detection ###########################################################################################################

print("Searching for neutron events...")
neutron_dict = functions_spills.neutron_detection(event_number_branch, times_branch_filtered,  threshold_times_prompt, coincidence_window, delayed_window, delayed_nhits_min, delayed_nhits_max, prompt_window)
#print("Detectando neutrones en señal...")
#neutron_dict_sig = functions_spills.neutron_detection(event_number_branch_sig, times_branch_filtered_sig,  threshold_times_prompt_sig, 150000, window_size_neutron, 15, threshold_sup=21, window_prompt=1500)

print("Prompt candidates", sum(len(v) for v in threshold_times_prompt.values()))
print("Neutron candidates", sum(len(v) for v in neutron_dict.values()))

#print("Prompt candidates signal", sum(len(v) for v in threshold_times_prompt_sig.values()))
#print("Neutron candidates signal", sum(len(v) for v in neutron_dict_sig.values()))


# Save neutron candidates to CSV files ###########################################################################################################

print("Saving candidate neutron events on CSV...")

#df = pd.read_csv(f'{output_path}/TestPartition.csv')
#df_sig = pd.read_csv('csv_saveData/TestPartition_sig.csv')

neutron_candidates = []
#for event_number, times in neutron_dict.items():
#    for start_time, neutron_times in times.items():
#        for neutron_time in neutron_times:
#            #partition, event_number_partition = functions_analysis.get_partition_info(event_number, df)
#            neutron_candidates.append({
#                #'partition': partition, 
#                #'event_number_partition': event_number_partition,
#                'event_number': event_number,
#                'start_time': start_time,
#                'neutron_time': neutron_time      
#            })

for event_number, times in neutron_dict.items():
    #prompt_candidates = threshold_times_prompt.get(event_number, [])  # List of (time_prompt, n_hits)
    ## Optional: Split prompt_candidates into separate lists if you want
    #prompt_times = [t for t, n in prompt_candidates]
    #prompt_nhits = [n for t, n in prompt_candidates]
    
    #for start_time, neutron_pairs in times.items():
    #    prompt_times_py = [float(t) for t in prompt_times]
    #    prompt_nhits_py = [int(n) for n in prompt_nhits]
    #    for neutron_time, neutron_nhits in neutron_pairs: 
    #        neutron_candidates.append({
    #            'event_number': event_number,
    #            'prompt_time': start_time,
    #            'prompt_nhits': float(prompt_nhits_py[0]),      # List of nhits
    #            'delayed_time': neutron_time,
    #            'delayed_nhits': neutron_nhits
    #        })
    #prompt_candidates = threshold_times_prompt.get(event_number, [])  # List of (time_prompt, n_hits)
    ##prompt_nhits_map = {float(t): float(n) for t, n in prompt_candidates}
    #prompt_nhits_map = {float(t): float(n[0]) if isinstance(n, list) else float(n) for t, n in prompt_candidates}
    #for start_time, neutron_pairs in times.items():
    #    prompt_nhits_val = prompt_nhits_map.get(float(start_time), None)
    #    for neutron_time, neutron_nhits in neutron_pairs: 
    #        neutron_candidates.append({
    #            'event_number': event_number,
    #            'prompt_time': float(start_time),
    #            'prompt_nhits': prompt_nhits_val,
    #            'delayed_time': float(neutron_time),
    #            'delayed_nhits': neutron_nhits
    #        })
    
    #prompt_candidates = threshold_times_prompt.get(event_number, [])  # List of (time_prompt, n_hits)
    #prompt_nhits_map = {float(t): float(n[0]) if isinstance(n, list) else float(n) for t, n in prompt_candidates}
    #for start_time, neutron_pairs in times.items():
    #    prompt_nhits_val = prompt_nhits_map.get(float(start_time), None)
    #    # neutron_pairs is a list of (delayed_time, delayed_nhits)
    #    for delayed_time, delayed_nhits in neutron_pairs:
    #        neutron_candidates.append({
    #            'event_number': event_number,
    #            'prompt_time': float(start_time),
    #            'prompt_nhits': prompt_nhits_val,
    #            'delayed_time': float(delayed_time),
    #            'delayed_nhits': float(delayed_nhits)
    #        })
    
    prompt_candidates = threshold_times_prompt.get(event_number, [])  # List of (time_prompt, n_hits, prompt_trms)
    #prompt_nhits_map = {float(t): (float(n), float(trms)) for t, n, trms in prompt_candidates}
    prompt_nhits_map = {
    float(t): (
        float(n[0]) if isinstance(n, list) else float(n),
        float(trms[0]) if isinstance(trms, list) else float(trms)
    )
    for t, n, trms in prompt_candidates
    }
    for start_time, neutron_pairs in times.items():
        prompt_data = prompt_nhits_map.get(float(start_time), (None, None))
        prompt_nhits_val, prompt_trms_val = prompt_data
        for delayed_time, delayed_nhits in neutron_pairs:
            neutron_candidates.append({
                'event_number': event_number,
                'prompt_time': float(start_time),
                'prompt_nhits': prompt_nhits_val,
                'prompt_trms': prompt_trms_val,
                'delayed_time': float(delayed_time),
                'delayed_nhits': float(delayed_nhits) if not isinstance(delayed_nhits, list) else float(delayed_nhits[0])
            })
#neutron_candidates_sig = []
#for event_number, times in neutron_dict_sig.items():
#    for start_time, neutron_times in times.items():
#        for neutron_time in neutron_times:
#            partition, event_number_partition = functions_analysis.get_partition_info(event_number, df_sig)
#            neutron_candidates_sig.append({
#                'partition': partition, 
#                'event_number_partition': event_number_partition,
#                'event_number': event_number,
#                'start_time': start_time,
#                'neutron_time': neutron_time
#            })
            
df_neutron_candidates = pd.DataFrame(neutron_candidates)
#df_neutron_candidates_sig = pd.DataFrame(neutron_candidates_sig)
df_neutron_candidates.to_csv(f'{output_path}/AmBeCandidates/neutron_candidates_{run_number}.csv', index=False)
#df_neutron_candidates_sig.to_csv(f'/scratch/cgarcia_2002/Complete_analysis/csv_saveData/Neutron_candidates/neutron_nHitsTest/neutron_candidates_sig_150-300_15-21_200-400.csv', index=False)

print("CSV files saved.")
print("End of code.")

#Plot filtered data #################################################################################################
"""
print("Generando histogramas antes y después del filtrado...")
nDetections_event_in = []
nDetections_event_fin = []

nDetections_event_in_sig = []
nDetections_event_fin_sig = []

for x in times_branch_sorted_TOF:
    nDetections_event_in.append(len(x))

for x in times_branch_filtered:
    nDetections_event_fin.append(len(x))

for x in times_branch_sorted_TOF_sig:
    nDetections_event_in_sig.append(len(x))

for x in times_branch_filtered_sig:
    nDetections_event_fin_sig.append(len(x))

n_bins = 600

hist_in, bin_edges = np.histogram(nDetections_event_in, bins=n_bins)
hist_in_sig, _ = np.histogram(nDetections_event_in_sig, bins=bin_edges)  # usa los mismos bordes

hist_filtered, _ = np.histogram(nDetections_event_fin, bins = bin_edges)
hist_filtered_sig, _ = np.histogram(nDetections_event_fin_sig, bins = bin_edges)

plt.figure(facecolor='white')
plt.step(bin_edges[:-1], hist_in, where='post', color='crimson', linestyle='--', linewidth=1.5, label='Background')
plt.step(bin_edges[:-1], hist_in_sig *N_events/ N_events_sig, where='post', color='navy', linestyle='-', linewidth=1.5, label='Data')
plt.legend()
plt.xlabel('Number of Hits')
plt.ylabel('Number of Events')
plt.title('Histogram hits/event')
plt.tight_layout()
plt.xlim(0, 4000)
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/Initial_Hist_hitsEvent.png")


plt.figure(facecolor='white')
plt.fill_between(bin_edges[:-1], hist_filtered,  hatch='//////',  step='post', color='white', edgecolor='red', alpha=0.55, linestyle= '--', label='Bkg (after filter)')
plt.fill_between(bin_edges[:-1], hist_filtered_sig *N_events / N_events_sig, hatch='\\\\\\\\', step='post', color='white', edgecolor='blue', alpha=0.55, label='Data (after filter)')
plt.step(bin_edges[:-1], hist_in, where='post', color='crimson', linestyle='--', linewidth=1.5, label='Bkg (before filter)')
plt.step(bin_edges[:-1], hist_in_sig *N_events/ N_events_sig, where='post', color='navy', linestyle='-', linewidth=1.5, label='Data (before filter)')
plt.legend()
plt.xlabel('Number of Hits')
plt.ylabel('Number of Events')
plt.title('Histograms Before and After Filtering')
plt.tight_layout()
plt.xlim(0, 4000)
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/Comparision_SigBkg_Filtering.png")

bin_window = 4000

nHits_tot = functions_spills.counting_nHits_window(event_number_branch, times_branch_filtered, bin_window)
nHits_in = functions_spills.counting_nHits_window(event_number_branch, times_branch_sorted_TOF, bin_window)

nHits_tot_sig = functions_spills.counting_nHits_window(event_number_branch_sig, times_branch_filtered_sig, bin_window)
nHits_in_sig = functions_spills.counting_nHits_window(event_number_branch_sig, times_branch_sorted_TOF_sig, bin_window)

hist_in, bin_edges = np.histogram(nHits_in, bins=100, range=(0, 500))
hist_in_sig, _ = np.histogram(nHits_in_sig, bins=bin_edges)  # usa los mismos bordes

hist_filtered, _ = np.histogram(nHits_tot, bins = bin_edges)
hist_filtered_sig, _ = np.histogram(nHits_tot_sig, bins = bin_edges)

n_windows_ev = 270000 / bin_window

fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})

# Top plot: Background and Signal
axs[0].step(bin_edges[:-1], hist_in / (N_events * n_windows_ev), linewidth = 1, where='post', label='Background', color='red')
axs[0].step(bin_edges[:-1], hist_in_sig / (N_events_sig * n_windows_ev), linewidth = 1, where='post', label='Data', color='blue')
axs[0].set_ylabel("Fraction of windows")
axs[0].set_xlabel(f"hits in {bin_window} ns")
axs[0].set_title("Initial data")
axs[0].set_yscale('log')
axs[0].legend()

# Bottom plot: Signal/Background Ratio
ratio = np.divide(
    hist_in_sig / (N_events_sig * n_windows_ev),
    hist_in / (N_events * n_windows_ev),
    out=np.full_like(hist_in, 0, dtype=float),
    where=hist_in > 0)
axs[1].step(bin_edges[:-1], ratio, linewidth = 1, where='post', color='green', label='Data / Background')
axs[1].set_xlabel(f"hits in {bin_window} ns")
axs[1].set_ylabel("S/B ratio")
axs[1].legend()
plt.tight_layout()
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/Comparision_SigBkg_Windows_BeforeFiltering.png")

fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})

# Top plot: Background and Signal
axs[0].step(bin_edges[:-1], hist_filtered / (N_events * n_windows_ev), linewidth = 1, where='post', label='Background', color='red')
axs[0].step(bin_edges[:-1], hist_filtered_sig / (N_events_sig * n_windows_ev), linewidth = 1, where='post', label='Data', color='blue')
axs[0].set_ylabel("Fraction of windows")
axs[0].set_xlabel(f"hits in {bin_window} ns")
axs[0].set_title("After Filtering total data")
axs[0].set_yscale('log')
axs[0].legend()
# Bottom plot: Signal/Background Ratio
ratio_2 = np.divide(
    hist_filtered_sig / (N_events_sig * n_windows_ev),
    hist_filtered / (N_events * n_windows_ev),
    out=np.full_like(hist_filtered, 0, dtype=float),
    where=hist_filtered > 0)
axs[1].step(bin_edges[:-1], ratio_2, linewidth = 1, where='post', color='green', label='Data / Background')
axs[1].set_xlabel(f"hits in {bin_window} ns")
axs[1].set_ylabel("S/B ratio")
axs[1].legend()
plt.tight_layout()
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/Comparision_SigBkg_Windows_AfterFiltering.png")

print("Histogramas guardados.")"""




"""window_ns = 500  # podrias hacer plots con esto pero en este script no sirve de nada

times_branch_prompt = []
for event in event_number_branch:

    if event in threshold_times_prompt.keys():
    
        all_hits = [t for ref_time in threshold_times_prompt[event]
            for t in times_branch_filtered[event]
            if ref_time <= t <= ref_time + window_ns]
    else:
        all_hits = []
    
    times_branch_prompt.append(np.array(all_hits))

times_branch_prompt_sig = []
for event in event_number_branch_sig:
    if event in  threshold_times_prompt_sig.keys():

        all_hits = [t for ref_time in threshold_times_prompt_sig[event]
            for t in times_branch_filtered_sig[event]
            if ref_time <= t <= ref_time + window_ns]
    else:
        all_hits= []
    times_branch_prompt_sig.append(np.array(all_hits))"""



# Delta T calculation ###########################################################################################################

"""print("Calculando Delta T entre prompt y neutrones...")

deltaT = []
for event_number in neutron_dict:
    for start_time in neutron_dict[event_number]:
        neutron_times = neutron_dict[event_number][start_time]
        deltaT.append(min(neutron_times) - start_time)

deltaT_sig = []
for event_number in neutron_dict_sig:
    for start_time in neutron_dict_sig[event_number]:
        neutron_times = neutron_dict_sig[event_number][start_time]
        deltaT_sig.append(min(neutron_times) - start_time)


hist, bins_edges = np.histogram(deltaT, bins=100, range=(0, 100000))
hist_sig, _ = np.histogram(deltaT_sig, bins=bins_edges)

plt.figure()

plt.step(bins_edges[:-1], hist, where='post', linewidth=1, label='bkg', color='red')
plt.step(bins_edges[:-1], hist_sig*N_events/N_events_sig, where='post', linewidth=1, label='signal', color='blue')
plt.xlabel('Delta T (ns)')
plt.ylabel('Número de signal candidates')
plt.title('Delta T between neutron and prompt candidate')
plt.legend()
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/DeltaT_Neutron_Prompt_100-300_22-30_TestPartitions.png")

print("Delta T calculado y gráfico guardado.")"""

