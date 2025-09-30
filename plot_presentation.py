print("Iniciando el script de filtrado de spills")
print("Importando librerias necesarias")

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import functions_spills
import glob
import os
import functions_spills
import functions_analysis
import pickle

from matplotlib import rcParams
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'
rcParams['figure.figsize'] = [10, 8]
rcParams['font.size'] = 18

print("Cargando datos de bkg...")

root_dir_bkg = "/data/cgarcia_2002/WCTE/data/2384_calib_time/"
root_files_bkg = sorted(glob.glob(os.path.join(root_dir_bkg, "*.root")))
root_files_bkg = sorted(root_files_bkg, key=lambda file_path: int(file_path.split("P")[-1].split(".")[0]))

print(f"Found {len(root_files_bkg)} background ROOT files.")

times_branch_sorted, times_branch_sorted_TOF, charge_branch_sorted, mpmt_id_branch_sorted, event_number_branch, _ = functions_spills.multiple_partition(root_files_bkg)

print("Datos de background cargados.")
N_events = max(event_number_branch) + 1

#Signal data download #############################################################################################

print("Cargando datos de signal...")
root_dir_sig = "/data/cgarcia_2002/WCTE/data/2385_calib_time/"
root_files_sig = sorted(glob.glob(os.path.join(root_dir_sig, "*.root")))
root_files_sig = sorted(root_files_sig, key=lambda file_path: int(file_path.split("P")[-1].split(".")[0]))

print(f"Found {len(root_files_sig)} signal ROOT files.")

times_branch_sorted_sig, times_branch_sorted_TOF_sig, charge_branch_sorted_sig, mpmt_id_branch_sorted_sig, event_number_branch_sig, _ = functions_spills.multiple_partition(root_files_sig)

print("Datos de signal cargados.")
N_events_sig = max(event_number_branch_sig) + 1


with open('Filtered_data/datos_filtrados.pkl', 'rb') as f:
    valores_read, indices_read = pickle.load(f)

with open('Filtered_data/datos_filtrados_sig.pkl', 'rb') as f:
    valores_read_sig, indices_read_sig = pickle.load(f)

times_branch_filtered = functions_analysis.a_lista_de_arrays(valores_read, indices_read)
times_branch_filtered_sig = functions_analysis.a_lista_de_arrays(valores_read_sig, indices_read_sig)

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


n_bins = 1500

hist_in, bin_edges = np.histogram(nDetections_event_in, bins=n_bins)
hist_in_sig, _ = np.histogram(nDetections_event_in_sig, bins=bin_edges)  # usa los mismos bordes

hist_filtered, _ = np.histogram(nDetections_event_fin, bins = bin_edges)
hist_filtered_sig, _ = np.histogram(nDetections_event_fin_sig, bins = bin_edges)

plt.figure(facecolor='white')
plt.step(bin_edges[:-1], hist_in_sig *N_events/ N_events_sig, where='post', color='navy', linestyle='-', linewidth=1.5, label='Data')
plt.step(bin_edges[:-1], hist_in, where='post', color='crimson', linewidth=1.5, linestyle= '--', label='Background')
plt.legend()
plt.xlabel('Number of Hits')
plt.ylabel('Number of Events')
plt.title('Histogram hits/event')
plt.tight_layout()
plt.xlim(0, 4000)
plt.yscale('log')
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/Initial_Hist_hitsEvent_log.png")


plt.figure(facecolor='white')
plt.step(bin_edges[:-1], hist_in_sig *N_events/ N_events_sig, where='post', color='navy', linestyle='-', linewidth=1.5, label='Data (before filter)')
plt.fill_between(bin_edges[:-1], hist_filtered_sig *N_events / N_events_sig, hatch='\\\\\\\\', step='post', color='white', edgecolor='blue', alpha=0.55, label='Data (after filter)')
plt.step(bin_edges[:-1], hist_in, where='post', color='crimson', linewidth=1.5, label='Bkg (before filter)')
plt.fill_between(bin_edges[:-1], hist_filtered,  hatch='//////',  step='post', color='white', edgecolor='red', alpha=0.55, label='Bkg (after filter)')
plt.legend()
plt.xlabel('Number of Hits')
plt.ylabel('Number of Events')
plt.title('Histograms Before and After Filtering')
plt.tight_layout()
plt.xlim(0, 4000)
plt.yscale('log')
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/Comparision_SigBkg_Filtering_log.png")