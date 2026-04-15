from collections import defaultdict
import numpy as np
from scipy.optimize import least_squares
import functions_bonsai

def run_multilateration_candidate(
    times, mpmt_ids, pmt_ids,
    *,
    sigma_t=1.0,
    n=1.33,
    c_cm_per_ns=29.9792458,   # cm/ns
    guess=(0., 0., 0., 0.),
    mins=(-300., -300., -300., -300.),
    maxs=(300., 300., 300., 300.),
    drop_invalid_geo=True,
    earliest_per_channel=True,
    early_window_ns=100.0,    # set None to disable extra windowing
    robust_loss="soft_l1",
    f_scale=2.0,
    **kwargs
):
    times = np.asarray(times, dtype=float)
    mpmt_ids = np.asarray(mpmt_ids, dtype=int)
    pmt_ids = np.asarray(pmt_ids, dtype=int)

    x_pmt, y_pmt, z_pmt, _ = functions_bonsai.getxyz(functions_bonsai.geo, mpmt_ids, pmt_ids)

    if drop_invalid_geo:
        good = (x_pmt > -900) & (y_pmt > -900) & (z_pmt > -900) & np.isfinite(times)
        times = times[good]
        mpmt_ids = mpmt_ids[good]
        pmt_ids = pmt_ids[good]
        x_pmt, y_pmt, z_pmt = x_pmt[good], y_pmt[good], z_pmt[good]

    # --- NEW: keep earliest hit per (mpmt_id,pmt_id) ---
    if earliest_per_channel and len(times) > 0:
        # Optional: reject very late hits globally (helps a lot for delayed candidates)
        if early_window_ns is not None:
            tmin = float(np.min(times))
            wmask = times <= (tmin + float(early_window_ns))
            times = times[wmask]
            mpmt_ids = mpmt_ids[wmask]
            pmt_ids = pmt_ids[wmask]
            x_pmt, y_pmt, z_pmt = x_pmt[wmask], y_pmt[wmask], z_pmt[wmask]

        # Find earliest time index for each channel
        # Approach: sort by time, then take first occurrence of each (mpmt,pmt).
        order = np.argsort(times)
        times_s = times[order]
        mpmt_s = mpmt_ids[order]
        pmt_s = pmt_ids[order]
        x_s, y_s, z_s = x_pmt[order], y_pmt[order], z_pmt[order]

        seen = set()
        keep_idx = []
        for i, key in enumerate(zip(mpmt_s, pmt_s)):
            if key in seen:
                continue
            seen.add(key)
            keep_idx.append(i)

        keep_idx = np.asarray(keep_idx, dtype=int)

        times = times_s[keep_idx]
        mpmt_ids = mpmt_s[keep_idx]
        pmt_ids = pmt_s[keep_idx]
        x_pmt, y_pmt, z_pmt = x_s[keep_idx], y_s[keep_idx], z_s[keep_idx]
    # --- END NEW ---

    if len(times) < 6:
        return {"x": np.nan, "y": np.nan, "z": np.nan, "eps": np.nan,
                "success": False, "n_hits_used": int(len(times)), "result": None}

    pmt_locs = np.column_stack([x_pmt, y_pmt, z_pmt]).astype(float)
    vc = float(c_cm_per_ns) / float(n)  # cm/ns
    sigma_ts = np.full(times.shape, float(sigma_t), dtype=float)

    # shift time origin
    t0 = float(np.min(times))
    times0 = times - t0

    # better eps guess (important)
    loc0 = np.array([0., 0., 0.])
    tof0 = np.linalg.norm(pmt_locs - loc0, axis=1) / vc
    eps_guess = float(np.median(times0 - tof0))

    x0 = np.array([0., 0., 0., eps_guess], dtype=float)

    def rho(pars):
        loc = pars[0:3]
        eps = pars[3]
        dists = np.linalg.norm(pmt_locs - loc, axis=1)
        tofs = dists / vc
        return (times0 - eps - tofs) / sigma_ts

    def jac(pars):
        loc = pars[0:3]
        light_vecs = pmt_locs - loc
        dists = np.linalg.norm(light_vecs, axis=1)
        dists = np.where(dists == 0, 1e-12, dists)
        jac_xyz = light_vecs / dists.reshape(-1, 1) / vc / sigma_ts.reshape(-1, 1)
        jac_eps = -1.0 / sigma_ts
        return np.column_stack([jac_xyz, jac_eps])

    loss = robust_loss if robust_loss is not None else "linear"

    result = least_squares(
        rho, x0, jac,
        bounds=(np.array(mins, dtype=float), np.array(maxs, dtype=float)),
        loss=loss,
        f_scale=float(f_scale),
        **kwargs
    )

    if not result.success:
        return {"x": np.nan, "y": np.nan, "z": np.nan, "eps": np.nan,
                "success": False, "n_hits_used": int(len(times)), "result": result}

    pulls = rho(result.x)                 # 1D array of pulls
    chi2 = float(np.sum(pulls**2))        # TimeCal style
    ndof = int(len(pulls) - 4)            # 4 fit params: x,y,z,eps
    chi2_ndof = float(chi2 / ndof) if ndof > 0 else np.inf

    x, y, z, eps0 = result.x
    eps_abs = float(eps0 + t0)

    return {
        "x": float(x),
        "y": float(y),
        "z": float(z),
        "eps": float(eps_abs),
        "success": True,
        "n_hits_used": int(len(pulls)),
        "chi2": chi2,
        "ndof": ndof,
        "chi2_ndof": chi2_ndof,
        "result": result,
    }

    import numpy as np
import functions_bonsai

def run_grid_vertex_candidate(
    times, mpmt_ids, pmt_ids,
    *,
    n=1.33,
    c_cm_per_ns=29.9792458,      # cm/ns
    xyz_bounds_cm=300.0,         # global search volume
    coarse_step_cm=10.0,
    fine_step_cm=1.0,
    refine_halfwidth_cm=20.0,    # +/- around best coarse
    dt_cut_ns=3.0,               # cut used in fine stage (set None to disable)
    earliest_per_channel=True,   # recommended
    drop_invalid_geo=True,
):
    """
    Grid-search vertex:
      1) coarse scan over cube [-bounds,+bounds] with coarse_step_cm
      2) compute corrected times t_corr = t - TOF
         choose t0 = median(t_corr)
         dt = t_corr - t0
         metric = TRMS = sqrt(mean(dt^2))
      3) refine in smaller box around coarse best with fine_step_cm
         optionally apply |dt|<dt_cut in the metric
    Returns dict with best vertex + quality metrics.
    """

    times = np.asarray(times, dtype=float)
    mpmt_ids = np.asarray(mpmt_ids, dtype=int)
    pmt_ids = np.asarray(pmt_ids, dtype=int)

    # geometry lookup
    x_pmt, y_pmt, z_pmt, _ = functions_bonsai.getxyz(functions_bonsai.geo, mpmt_ids, pmt_ids)

    # drop invalid geometry
    if drop_invalid_geo:
        good = (x_pmt > -900) & (y_pmt > -900) & (z_pmt > -900) & np.isfinite(times)
        times = times[good]
        mpmt_ids = mpmt_ids[good]
        pmt_ids = pmt_ids[good]
        x_pmt, y_pmt, z_pmt = x_pmt[good], y_pmt[good], z_pmt[good]

    if len(times) < 6:
        return {"x": np.nan, "y": np.nan, "z": np.nan, "t0": np.nan,
                "trms": np.nan, "chi2": np.nan, "ndof": -1,
                "success": False, "n_hits_used": int(len(times))}

    # earliest hit per channel (optional but usually helps a lot)
    if earliest_per_channel:
        order = np.argsort(times)
        times_s = times[order]
        mpmt_s = mpmt_ids[order]
        pmt_s = pmt_ids[order]
        x_s, y_s, z_s = x_pmt[order], y_pmt[order], z_pmt[order]

        seen = set()
        keep = []
        for i, key in enumerate(zip(mpmt_s, pmt_s)):
            if key in seen:
                continue
            seen.add(key)
            keep.append(i)
        keep = np.asarray(keep, dtype=int)

        times = times_s[keep]
        x_pmt, y_pmt, z_pmt = x_s[keep], y_s[keep], z_s[keep]

    # shift time origin (like bonsai) for numeric stability; t0 returned in shifted basis
    tmin = float(np.min(times))
    times0 = times - tmin

    pmt_locs = np.column_stack([x_pmt, y_pmt, z_pmt]).astype(float)
    vc = float(c_cm_per_ns) / float(n)  # cm/ns

    def score_vertex(vxyz, use_cut=False):
        """Return (TRMS, t0, chi2, ndof, n_used) for a trial vertex."""
        vxyz = np.asarray(vxyz, dtype=float)
        dists = np.linalg.norm(pmt_locs - vxyz, axis=1)
        tofs = dists / vc
        t_corr = times0 - tofs

        # best t0 for L2-like alignment is mean; for robustness use median:
        t0 = float(np.median(t_corr))
        dt = t_corr - t0

        if use_cut and (dt_cut_ns is not None):
            mask = np.abs(dt) < float(dt_cut_ns)
            if np.count_nonzero(mask) < 6:
                # too few hits survive: penalize heavily
                return (np.inf, t0, np.inf, -1, int(np.count_nonzero(mask)))
            dt_use = dt[mask]
        else:
            dt_use = dt

        trms = float(np.sqrt(np.mean(dt_use**2)))

        # chi2-like: sum(dt^2) with implicit sigma=1 ns
        chi2 = float(np.sum(dt_use**2))
        ndof = int(len(dt_use) - 4)  # keep same convention as multilateration (x,y,z,t0)
        return (trms, t0, chi2, ndof, int(len(dt_use)))

    def make_grid(xmin, xmax, step):
        xs = np.arange(xmin, xmax + 0.5*step, step)
        return xs

    # --- Stage 1: coarse global grid ---
    xs = make_grid(-xyz_bounds_cm, xyz_bounds_cm, coarse_step_cm)
    best = (np.inf, None, None, None, None)  # trms, vxyz, t0, chi2, info
    for x in xs:
        for y in xs:
            for z in xs:
                trms, t0, chi2, ndof, n_used = score_vertex((x, y, z), use_cut=False)
                if trms < best[0]:
                    best = (trms, (float(x), float(y), float(z)), t0, chi2, (ndof, n_used))

    coarse_trms, (bx, by, bz), coarse_t0, coarse_chi2, (coarse_ndof, coarse_nused) = best

    # --- Stage 2: fine grid around coarse best ---
    fxmin, fxmax = bx - refine_halfwidth_cm, bx + refine_halfwidth_cm
    fymin, fymax = by - refine_halfwidth_cm, by + refine_halfwidth_cm
    fzmin, fzmax = bz - refine_halfwidth_cm, bz + refine_halfwidth_cm

    fxs = make_grid(fxmin, fxmax, fine_step_cm)
    fys = make_grid(fymin, fymax, fine_step_cm)
    fzs = make_grid(fzmin, fzmax, fine_step_cm)
    
    best2 = (np.inf, None, None, None, None)
    for x in fxs:
        for y in fys:
            for z in fzs:
                trms, t0, chi2, ndof, n_used = score_vertex((x, y, z), use_cut=True)
                if trms < best2[0]:
                    best2 = (trms, (float(x), float(y), float(z)), t0, chi2, (ndof, n_used))

    # Option A fallback: if dt_cut rejects everything, redo without the cut
    if best2[1] is None:
        best2 = (np.inf, None, None, None, None)
        for x in fxs:
            for y in fys:
                for z in fzs:
                    trms, t0, chi2, ndof, n_used = score_vertex((x, y, z), use_cut=False)
                    if trms < best2[0]:
                        best2 = (trms, (float(x), float(y), float(z)), t0, chi2, (ndof, n_used))

    fine_trms, (fx, fy, fz), fine_t0, fine_chi2, (fine_ndof, fine_nused) = best2

    # Convert t0 back to absolute basis if desired:
    # t0_abs = fine_t0 + tmin
    # Here I return t0 in absolute time basis:
    t0_abs = float(fine_t0 + tmin)

    return {
        "x": float(fx),
        "y": float(fy),
        "z": float(fz),
        "t0": t0_abs,
        "trms": float(fine_trms),
        "chi2": float(fine_chi2),
        "ndof": int(fine_ndof),
        "chi2_ndof": float(fine_chi2 / fine_ndof) if fine_ndof > 0 else np.inf,
        "success": np.isfinite(fine_trms),
        "n_hits_used": int(fine_nused),
        "coarse": {
            "x": float(bx), "y": float(by), "z": float(bz),
            "trms": float(coarse_trms),
            "t0": float(coarse_t0 + tmin),
            "chi2": float(coarse_chi2),
            "ndof": int(coarse_ndof),
            "n_hits_used": int(coarse_nused),
        }
    }