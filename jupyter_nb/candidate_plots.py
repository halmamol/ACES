import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from collections import defaultdict
import matplotlib.ticker as ticker
import glob
import os
import pickle
import sys


def deltaR(df):
    # Ensure numeric columns (coerce any stray strings/lists)
    cols = ["prompt_x","prompt_y","prompt_z","delayed_x","delayed_y","delayed_z"]
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Build position arrays
    p = df[["prompt_x","prompt_y","prompt_z"]].to_numpy()
    d = df[["delayed_x","delayed_y","delayed_z"]].to_numpy()

    # Radii from origin
    df["prompt_r"] = np.linalg.norm(p, axis=1)
    df["delayed_r"] = np.linalg.norm(d, axis=1)

    # DeltaR: 3D distance between prompt and delayed vertices
    df["DeltaR"] = np.abs(df["delayed_r"] - df["prompt_r"])

def deltaT(df):
    # DeltaT: Time difference between prompt and delayed vertices
    df["DeltaT"] = np.abs(df["delayed_time"]/1000 - df["prompt_time"]/1000)

def exp_func(x, A, tau):
    return A * np.exp(-x / tau)

def deltaR_to_fixed(df, pos):
    # Ensure numeric columns
    cols = ["delayed_x", "delayed_y", "delayed_z"]
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    
    # Build delayed position array
    d = df[["delayed_x", "delayed_y", "delayed_z"]].to_numpy()

    # Calculate DeltaR: Euclidean distance between pos and delayed positions
    df["DeltaR"] = np.linalg.norm(d - pos, axis=1)

def pl_nHits_pos(df, run):
    nhits_bins = np.linspace(df['delayed_nhits'].min(), df['delayed_nhits'].max(), 10)
    x_bins = np.linspace(df['delayed_x'].min(), df['delayed_x'].max(), 20)
    y_bins = np.linspace(df['delayed_y'].min(), df['delayed_y'].max(), 20)
    z_bins = np.linspace(df['delayed_z'].min(), df['delayed_z'].max(), 20)

    fig, axes = plt.subplots(1, 3, figsize=(22, 5), sharey=True)

    # X
    hist_x, xedges, nhits_edges = np.histogram2d(df['delayed_x'], df['delayed_nhits'], bins=[x_bins, nhits_bins])
    im0 = axes[0].imshow(hist_x.T, origin='lower', aspect='auto',
                        extent=[xedges[0], xedges[-1], nhits_edges[0], nhits_edges[-1]],
                        cmap='plasma')
    axes[0].set_xlabel('Reconstructed X (cm)')
    axes[0].set_ylabel('nHits')

    # Y
    hist_y, yedges, nhits_edges = np.histogram2d(df['delayed_y'], df['delayed_nhits'], bins=[y_bins, nhits_bins])
    im1 = axes[1].imshow(hist_y.T, origin='lower', aspect='auto',
                        extent=[yedges[0], yedges[-1], nhits_edges[0], nhits_edges[-1]],
                        cmap='plasma')
    axes[1].set_xlabel('Reconstructed Y (cm)')

    # Z
    hist_z, zedges, nhits_edges = np.histogram2d(df['delayed_z'], df['delayed_nhits'], bins=[z_bins, nhits_bins])
    im2 = axes[2].imshow(hist_z.T, origin='lower', aspect='auto',
                        extent=[zedges[0], zedges[-1], nhits_edges[0], nhits_edges[-1]],
                        cmap='plasma')
    axes[2].set_xlabel('Reconstructed Z (cm)')

    # Single colorbar for all
    fig.colorbar(im2, ax=axes, orientation='vertical', label='Events (a.u.)', pad=0.01, fraction=0.05)

    fig.suptitle(f'Run {run}')
    plt.show()

def pl_deltaT_pos(df, run):
    # Define bins for all axes
    deltaT_bins = np.linspace(df['DeltaT'].min(), df['DeltaT'].max(), 20)
    x_bins = np.linspace(df['delayed_x'].min(), df['delayed_x'].max(), 20)
    y_bins = np.linspace(df['delayed_y'].min(), df['delayed_y'].max(), 20)
    z_bins = np.linspace(df['delayed_z'].min(), df['delayed_z'].max(), 20)

    fig, axes = plt.subplots(1, 3, figsize=(22, 5), sharey=True)

    # X
    hist_x, xedges, nhits_edges = np.histogram2d(df['delayed_x'], df['DeltaT'], bins=[x_bins, deltaT_bins])
    im0 = axes[0].imshow(hist_x.T, origin='lower', aspect='auto',
                        extent=[xedges[0], xedges[-1], nhits_edges[0], nhits_edges[-1]],
                        cmap='plasma')
    axes[0].set_xlabel('Reconstructed X (cm)')
    axes[0].set_ylabel(r'$\Delta t$ ($\mu$s)')

    # Y
    hist_y, yedges, nhits_edges = np.histogram2d(df['delayed_y'], df['DeltaT'], bins=[y_bins, deltaT_bins])
    im1 = axes[1].imshow(hist_y.T, origin='lower', aspect='auto',
                        extent=[yedges[0], yedges[-1], nhits_edges[0], nhits_edges[-1]],
                        cmap='plasma')
    axes[1].set_xlabel('Reconstructed Y (cm)')

    # Z
    hist_z, zedges, nhits_edges = np.histogram2d(df['delayed_z'], df['DeltaT'], bins=[z_bins, deltaT_bins])
    im2 = axes[2].imshow(hist_z.T, origin='lower', aspect='auto',
                        extent=[zedges[0], zedges[-1], nhits_edges[0], nhits_edges[-1]],
                        cmap='plasma')
    axes[2].set_xlabel('Reconstructed Z (cm)')

    # Single colorbar for all
    fig.colorbar(im2, ax=axes, orientation='vertical', label='Events (a.u.)', pad=0.01, fraction=0.05)

    fig.suptitle(f'Run {run}')
    plt.show()

def pl_prompt(N_events_bkg, N_events_sig, df, df_bkg, sig_run, bkg_run, prompt_nhits_min, prompt_nhits_max, prompt_t_rms_min, prompt_t_rms_max):
    prompt_nhits = df.prompt_nhits.values
    prompt_nhits_bkg = df_bkg.prompt_nhits.values

    hist, bins_edges = np.histogram(prompt_nhits_bkg, bins=30, range=(prompt_nhits_min, prompt_nhits_max))
    hist_sig, _ = np.histogram(prompt_nhits, bins = bins_edges)

    bin_centers = (bins_edges[:-1] + bins_edges[1:]) / 2  
    data = hist_sig * N_events_bkg / N_events_sig - hist


    plt.figure()
    plt.step(bins_edges[:-1], hist *N_events_sig/N_events_bkg, where='post', label=f'background - run {bkg_run}', color='red', lw=2)
    plt.step(bins_edges[:-1], hist_sig, where='post', label=f'signal - run {sig_run}', color='blue', lw=2, linestyle='--')
    #plt.step(bin_centers, data, where='mid', linewidth=2, label='signal - bkg', color='green', alpha=0.6,) 
    plt.xlabel('nHits')
    plt.ylabel('events (a.u.)')
    plt.title('Prompt signal candidates')
    plt.legend()
    plt.grid(alpha=0.3)

    prompt_nhits = df.prompt_trms.values
    prompt_nhits_bkg = df_bkg.prompt_trms.values

    hist, bins_edges = np.histogram(prompt_nhits_bkg, bins=30, range=(prompt_t_rms_min, prompt_t_rms_max))
    hist_sig, _ = np.histogram(prompt_nhits, bins = bins_edges)

    bin_centers = (bins_edges[:-1] + bins_edges[1:]) / 2  
    data = hist_sig * N_events_bkg / N_events_sig - hist

    plt.figure()
    plt.step(bins_edges[:-1], hist*N_events_sig/N_events_bkg, where='post', label=f'background - run {bkg_run}', color='red', lw=2)
    plt.step(bins_edges[:-1], hist_sig, where='post', label=f'signal - run {sig_run}', color='blue', lw=2, linestyle='--') 
    #plt.step(bin_centers, data, where='mid', linewidth=2, label='signal - bkg', color='green', alpha=0.6,) 
    plt.xlabel('tRMS [ns]')
    plt.ylabel('events (a.u.)')
    plt.title('Prompt signal candidates')
    plt.legend()
    plt.grid(alpha=0.3)

def pl_delayed(N_events_bkg, N_events_sig, df, df_bkg, sig_run, bkg_run, delayed_nhits_min, delayed_nhits_max):
    prompt_nhits = df.delayed_nhits.values
    prompt_nhits_bkg = df_bkg.delayed_nhits.values

    hist, bins_edges = np.histogram(prompt_nhits_bkg, bins=20, range=(delayed_nhits_min, delayed_nhits_max))
    hist_sig, _ = np.histogram(prompt_nhits, bins = bins_edges)

    plt.figure()
    plt.step(bins_edges[:-1], hist*N_events_sig/N_events_bkg, where='post', label=f'background - run {bkg_run}', color='red', lw=2)
    plt.step(bins_edges[:-1], hist_sig, where='post', label=f'signal - run {sig_run}', color='blue', lw=2, linestyle='--') 
    plt.xlabel('nHits')
    plt.ylabel('events (a.u.)')
    plt.title('Delayed signal candidates')
    plt.legend()
    plt.grid(alpha=0.3)

    nhits_bins = np.linspace(df['delayed_nhits'].min(), df['delayed_nhits'].max(), 10)
    deltaT_bins = np.linspace(df['DeltaT'].min(), df['DeltaT'].max(), 10)
    hist, xedges, yedges = np.histogram2d(df['DeltaT'], df['delayed_nhits'], bins=[deltaT_bins, nhits_bins])

    plt.figure(figsize=(7,6))
    plt.imshow(hist.T, origin='lower', aspect='auto',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            cmap='plasma')

    plt.xlabel(r'$\Delta t$ ($\mu$s)')
    plt.ylabel('nHits')
    plt.title(f'Run {sig_run}')

    cb = plt.colorbar()
    cb.set_label('Events (a.u.)')
    plt.show()

def pl_delayed_vp(N_events_bkg, N_events_sig, df, df_bkg, sig_run, bkg_run, delayed_nhits_min, delayed_nhits_max):
    prompt_nhits = df.delayed_nhits.values
    prompt_nhits_bkg = df_bkg.vp_delayed_nhits.values

    hist, bins_edges = np.histogram(prompt_nhits_bkg, bins=20, range=(delayed_nhits_min, delayed_nhits_max))
    hist_sig, _ = np.histogram(prompt_nhits, bins = bins_edges)

    plt.figure()
    plt.step(bins_edges[:-1], hist*N_events_sig/N_events_bkg, where='post', label=f'background - run {bkg_run}', color='red', lw=2)
    plt.step(bins_edges[:-1], hist_sig, where='post', label=f'signal - run {sig_run}', color='blue', lw=2, linestyle='--') 
    plt.xlabel('nHits')
    plt.ylabel('events (a.u.)')
    plt.title('Delayed signal candidates')
    plt.legend()
    plt.grid(alpha=0.3)

    nhits_bins = np.linspace(df['delayed_nhits'].min(), df['delayed_nhits'].max(), 10)
    deltaT_bins = np.linspace(df['DeltaT'].min(), df['DeltaT'].max(), 10)
    hist, xedges, yedges = np.histogram2d(df['DeltaT'], df['delayed_nhits'], bins=[deltaT_bins, nhits_bins])

    plt.figure(figsize=(7,6))
    plt.imshow(hist.T, origin='lower', aspect='auto',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            cmap='plasma')

    plt.xlabel(r'$\Delta t$ ($\mu$s)')
    plt.ylabel('nHits')
    plt.title(f'Run {sig_run}')

    cb = plt.colorbar()
    cb.set_label('Events (a.u.)')
    plt.show()

def pl_bonsai(N_events_bkg, N_events_sig, df, df_bkg, pos, sig_run):
    cols = ["delayed_x", "delayed_y", "delayed_z"]
    labels = ["X [cm]", "Y [cm]", "Z [cm]"]
    source_pos = {"prompt_x": pos[0], "prompt_y": pos[1], "prompt_z": pos[2]}

    alpha = N_events_bkg / N_events_sig   # scale signal to background exposure

    # Clean arrays
    def clean(arr):
        arr = np.asarray(arr)
        return arr[np.isfinite(arr)]

    sx, sy, sz = [clean(df[c].to_numpy()) for c in cols]
    bx, by, bz = [clean(df_bkg[c].to_numpy()) for c in cols]

    n_bins = 25

    global_min = 0 
    global_max = 25

    # Shared bin ranges for consistency
    xmin, xmax = -350, 350
    ymin, ymax = -350, 350
    zmin, zmax = -350, 350
    xbins = np.linspace(xmin, xmax, n_bins+1)
    ybins = np.linspace(ymin, ymax, n_bins+1)
    zbins = np.linspace(zmin, zmax, n_bins+1)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)

    # SIGNAL: X vs Y
    h_sig_xy, xedges, yedges, img_sig_xy = axes[0,0].hist2d(
        sx, sy, bins=[xbins, ybins], weights=np.full(sx.shape, alpha),
        cmap='plasma', vmin=global_min, vmax=global_max, cmin=0.0001
    )
    axes[0,0].plot(source_pos["prompt_x"], source_pos["prompt_y"], 'rx', ms=8, mew=4, label='Source position')
    axes[0,0].set_xlabel(labels[0]); axes[0,0].set_ylabel(labels[1])
    axes[0,0].set_title("Signal")
    axes[0,0].legend()
    fig.colorbar(img_sig_xy, ax=axes[0,0], label="Events (a.u.)")

    # SIGNAL: X vs Z
    h_sig_xz, xedges, zedges, img_sig_xz = axes[0,1].hist2d(
        sx, sz, bins=[xbins, zbins], weights=np.full(sx.shape, alpha),
        cmap='plasma', vmin=global_min, vmax=global_max, cmin=0.0001
    )
    axes[0,1].plot(source_pos["prompt_x"], source_pos["prompt_z"], 'rx', ms=8, mew=4, label='Source position')
    axes[0,1].set_xlabel(labels[0]); axes[0,1].set_ylabel(labels[2])
    axes[0,1].set_title("Signal")
    #axes[0,1].legend()
    fig.colorbar(img_sig_xz, ax=axes[0,1], label="Events (a.u.)")

    # BACKGROUND: X vs Y
    h_bkg_xy, xedges, yedges, img_bkg_xy = axes[1,0].hist2d(
        bx, by, bins=[xbins, ybins],
        cmap='plasma', vmin=global_min, vmax=global_max, cmin=0.0001
    )
    axes[1,0].plot(source_pos["prompt_x"], source_pos["prompt_y"], 'rx', ms=8, mew=4,  label='Source position')
    axes[1,0].set_xlabel(labels[0]); axes[1,0].set_ylabel(labels[1])
    axes[1,0].set_title("Background")
    #axes[1,0].legend()
    fig.colorbar(img_bkg_xy, ax=axes[1,0], label="Events (a.u.)")

    # BACKGROUND: X vs Z
    h_bkg_xz, xedges, zedges, img_bkg_xz = axes[1,1].hist2d(
        bx, bz, bins=[xbins, zbins],
        cmap='plasma', vmin=global_min, vmax=global_max, cmin=0.0001
    )
    axes[1,1].plot(source_pos["prompt_x"], source_pos["prompt_z"], 'rx', ms=8, mew=4, label='Source position')
    axes[1,1].set_xlabel(labels[0]); axes[1,1].set_ylabel(labels[2])
    axes[1,1].set_title("Background")
    #axes[1,1].legend()
    fig.colorbar(img_bkg_xz, ax=axes[1,1], label="Events (a.u.)")

    fig.suptitle(f"AmBe Run {sig_run} - Delayed Candidates")
    plt.show()

    cols = ["delayed_x", "delayed_y", "delayed_z"]
    labels = ["X [cm]", "Y [cm]", "Z [cm]"]

    source_pos = {"delayed_x": pos[0], "delayed_y": pos[1], "delayed_z": pos[2]}

    alpha = N_events_bkg / N_events_sig  # scale signal to background exposure

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True, constrained_layout=True)

    for ax, col, label in zip(axes, cols, labels):
        # Drop non-finite
        s = df[col].to_numpy()
        b = df_bkg[col].to_numpy()
        s = s[np.isfinite(s)]
        b = b[np.isfinite(b)]

        # Common bins from combined range
        xmin = np.min([s.min(), b.min()]) if len(s) and len(b) else (s.min() if len(s) else b.min())
        xmax = np.max([s.max(), b.max()]) if len(s) and len(b) else (s.max() if len(s) else b.max())
        bins = np.linspace(xmin, xmax, 26)  # 50 bins

        # Histograms: scale signal counts by alpha via weights
        ax.hist(s, bins=bins, weights=np.full(s.shape, alpha), color="blue", alpha=0.6, label=f'signal')
        ax.hist(b, bins=bins, color="red", alpha=0.6, label=f'background')

        # Source marker line
        vline = source_pos[col]
        ax.axvline(vline, color="k", linestyle="--", linewidth=1.5, label="Source position")

        ax.set_xlabel(label)
        ax.grid(alpha=0.3)
        if(col =='delayed_x'):
            ax.legend()

    axes[0].set_ylabel("Events (a.u.)")
    fig.suptitle(f"AmBe Run {sig_run} - Delayed Candidates")
    plt.show()

def exp_func(x, A, tau):
    return A * np.exp(-x / tau)

def pl_deltaT(N_events_sig, N_events_bkg, df, df_bkg, sig_run, bkg_run):
    prompt_nhits = df.DeltaT.values
    prompt_nhits_bkg = df_bkg.DeltaT.values

    hist, bins_edges = np.histogram(prompt_nhits_bkg, bins=30, range=(0, 150))
    hist_sig, _ = np.histogram(prompt_nhits, bins = bins_edges)

    plt.figure()
    plt.step(bins_edges[:-1], hist*N_events_sig/N_events_bkg, where='post', label=f'background - run {bkg_run}', color='red', lw=2)
    plt.step(bins_edges[:-1], hist_sig, where='post', label=f'signal - run {sig_run}', color='blue', lw=2)
    plt.xlabel(r'$\Delta t$ ($\mu$s)')
    plt.ylabel('events (a.u.)')
    plt.title(f'Delayed signal candidates - Run {sig_run}')
    plt.legend()
    plt.grid(alpha=0.3)

    plt.figure()
    plt.step(bins_edges[:-1], hist_sig - hist*N_events_sig/N_events_bkg, where='post', linewidth=2, color='green', label=f'signal - bkg')
    plt.xlabel(r'$\Delta t$ ($\mu$s)')
    plt.ylabel('Number of signal candidates')
    plt.grid(alpha=0.3)
    plt.legend()
    plt.show()

    bin_centers = 0.5 * (bins_edges[:-1] + bins_edges[1:])
    values = hist_sig - hist*N_events_sig/N_events_bkg

    np.savez(f'/scratch/halmazan/WCTE/files/data/AmBeCandidates/histograms_run{sig_run}.npz', deltat=values, bins_deltat=bin_centers)

    # Prepare bin centers and data for full plot
    bin_centers = (bins_edges[:-1] + bins_edges[1:]) / 2  
    data = hist_sig  - hist*N_events_sig/N_events_bkg
    errors = np.sqrt(hist_sig  + hist*N_events_sig/N_events_bkg)

    fit_min = 10
    fit_max = 150
    # Mask for fit region
    region_mask = (bin_centers >= fit_min) & (bin_centers <= fit_max)
    fit_x = bin_centers[region_mask]
    fit_y = data[region_mask]

    # Fit only in region
    popt, pcov = curve_fit(exp_func, fit_x, fit_y, p0=(max(fit_y), 500))
    A_fit, tau_fit = popt
    A_err, tau_err = np.sqrt(np.diag(pcov))

    # Plot: full data, fit only in region
    plt.figure()
    plt.step(bin_centers, data, where='mid', linewidth=2, label='Data', color='green')
    plt.errorbar(bin_centers, data, yerr=errors, fmt='none', ecolor='black', elinewidth=1, capsize=2)

    # Overlay fit on selected region only
    x_fit = np.linspace(fit_min, fit_max, 500)
    plt.plot(x_fit, exp_func(x_fit, *popt), 'r-', label=f'Fit region: τ ≈ {tau_fit:.0f} ± {tau_err:.0f} μs', color="black")

    plt.xlabel(r'Time ($\mu$s)')
    plt.ylabel('Number of signal candidates')
    plt.grid(alpha=0.3)
    plt.title(f'Delayed signal candidates - Run {sig_run}')
    plt.legend()
    plt.tight_layout()
    plt.ylim(bottom=0)
    plt.show()