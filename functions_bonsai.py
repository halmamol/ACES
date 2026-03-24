import ROOT 
import glob 
import numpy as np 
import pandas as pd 
import sys 
import array
import pandas as pd
import os
import cppyy

cppyy.add_include_path(os.environ["WCSIM_BUILD_DIR"] + "/include/")
cppyy.load_library(os.environ["WCSIM_BUILD_DIR"] + "/lib/libWCSimRoot.so")

cppyy.add_include_path(os.environ["BONSAIDIR"] + "/bonsai/")
cppyy.load_library(os.environ["BONSAIDIR"] + "/libWCSimBonsai.so")

# # Setup HKBONSAI with WCTE geo
simfile = ROOT.TFile(os.environ["BONSAIDIR"]+"/NiCf/wcsim_dummy.root")
simtree = simfile.Get("wcsimGeoT")

def get_geo_mapping():
    geo = pd.read_csv(os.environ["BONSAIDIR"]+"/NiCf/geofile_NuPRISMBeamTest_16cShort_mPMT.txt", 
                      index_col=False,      # Do not use any column as the index
                      sep='\s+',            # Use whitespace as separator
                      skiprows=5,           # Skip the first 5 header lines
                      names=["id","mpmtid","spmtid",
                             "x","y","z","dx","dy","dz", "cyloc"])  # Explicit column names

    return geo
geo = get_geo_mapping()

geotree = None
for geoevent in simtree:
    geotree = geoevent.wcsimrootgeom
    break
    
bonsai = cppyy.gbl.WCSimBonsai()
bonsai.Init(geotree)

def getxyz(geo, mpmtids, posids):
    # Build a single lookup dictionary: {(mpmtid, spmtid): (x, y, z, id)}
    lookup = {
        (row.mpmtid, row.spmtid): (row.x, row.y, row.z, row.id)
        for row in geo.itertuples(index=False)
    }

    # Adjust input IDs to match geometry convention
    keys = zip((mid for mid in mpmtids), (sid + 1 for sid in posids))

    # Use the lookup dictionary to retrieve values efficiently
    results = [lookup.get(k, (-999.9, -999.9, -999.9, -999)) for k in keys]

    # Unpack results into separate arrays
    if len(results) == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])
    
    x, y, z, c = map(np.array, zip(*results))

    return x, y, z, c

def run_BONSAI_candidate(times, charges, mpmt, pmt, g=geo):
    
    # Start Filter
    ns = 1

    vertex = {
    "nhits": [],
    "nhitso": [],
    "x": [],
    "y": [],
    "z": [],
    "result0": [],
    "result1": [],
    "result2": [],
    "result3": [],
    "result4": [],
    "result5": [],
    "good0":[],
    "good1":[],
    "good2":[] 
    }

    _, _ , _, cables = getxyz(g, mpmt, pmt)

    # Run Bonsai
    bsVertex = array.array('f',3*[0.0])
    bsResult = array.array('f',6*[0.0])
    bsGood = array.array('f',3*[0.0])
    bsNhit = array.array('i',[len(cables)])
    bsNsel = array.array('i',[0])

    # Generate hit collection for this triggger
    bsCAB_a = array.array('i', cables)
    bsT_a = array.array('f', times - np.min(times) + 200)
    bsQ_a = array.array('f', charges)

    # Run Bonsai
    try:
        nhits = bonsai.BonsaiFit(bsVertex, bsResult, bsGood, bsNsel, bsNhit, bsCAB_a, bsT_a, bsQ_a);
    except:
        print("BONSAIFAILED");
        pass
#         print(nhits, bsVertex, bsResult, bsGood, bsNsel, bsNhit)

    vertex["nhits"].append(nhits)
    vertex["nhitso"].append(len(times))
    
    vertex["x"].append(bsVertex[0])
    vertex["y"].append(bsVertex[1])
    vertex["z"].append(bsVertex[2])
    vertex["result0"].append(bsResult[0])
    vertex["result1"].append(bsResult[1])
    vertex["result2"].append(bsResult[2])
    vertex["result3"].append(bsResult[3])
    vertex["result4"].append(bsResult[4])
    vertex["result5"].append(bsResult[5])
    vertex["good0"].append(bsGood[0])
    vertex["good1"].append(bsGood[1])
    vertex["good2"].append(bsGood[2])
    
    
    return vertex

#------------------------------------------
from scipy.integrate import quad
from scipy.optimize import minimize

# Constants
speed_of_light = 299792458  # m/s

# Function to calculate time-of-flight (t_TOF)
def time_of_flight(r_s, r_PMT):
    """
    Calculate time-of-flight for light traveling from source vertex to PMT.
    :param r_s: Source vertex [x, y, z] (numpy array)
    :param r_PMT: PMT position [x, y, z] (numpy array)
    :return: Time-of-flight in seconds
    """
    distance = np.linalg.norm(r_s - r_PMT)
    return distance / speed_of_light

# Likelihood function for a single PMT hit
def likelihood(t_hit, r_s, T0, r_PMT, lambda_decay):
    """
    Compute likelihood for a single hit given PMT time, position, and source vertex.
    :param t_hit: Time of the hit (float, in seconds)
    :param r_s: Source vertex [x, y, z] (numpy array)
    :param T0: First scintillation photon time (float, in seconds)
    :param r_PMT: PMT position [x, y, z] (numpy array)
    :param lambda_decay: Decay constant (float)
    :return: Likelihood value (float)
    """
    # Compute time-of-flight
    t_TOF = time_of_flight(r_s, r_PMT)
    
    # Integrand for unknown emission time
    def integrand(t_emission):
        return lambda_decay * np.exp(-lambda_decay * (t_hit - T0 - t_TOF - t_emission))
    
    # Perform numerical integration over emission times
    likelihood_value, _ = quad(integrand, 0, np.inf)
    return likelihood_value

# Objective function for optimization
def total_likelihood(params, times_in_delayed, charges_in_delayed, x, y, z):
    """
    Total likelihood (negative log-likelihood for optimization).
    :param params: Array of [x, y, z, T0, lambda_decay]
    :param times_in_delayed: Hit times [t_hit^(i)] (1D numpy array)
    :param charges_in_delayed: Hit charges [q_hit^(i)] (1D numpy array)
    :param x: X coordinates of PMT positions (1D numpy array)
    :param y: Y coordinates of PMT positions (1D numpy array)
    :param z: Z coordinates of PMT positions (1D numpy array)
    :return: Negative log-likelihood (float)
    """
    # Extract parameters
    r_s = np.array([params[0], params[1], params[2]])  # Source vertex [x, y, z]
    T0 = params[3]  # First photon time
    lambda_decay = params[4]  # Decay constant
    
    # Compute total likelihood
    total_log_likelihood = 0
    for i, t_hit in enumerate(times_in_delayed):
        r_PMT = np.array([x[i], y[i], z[i]])  # PMT position
        likelihood_value = likelihood(t_hit, r_s, T0, r_PMT, lambda_decay)
        weighted_likelihood = likelihood_value * charges_in_delayed[i]  # Weight by charge
        total_log_likelihood -= np.log(weighted_likelihood)  # Negative log-likelihood
    
    return total_log_likelihood

# Optimization function
def optimize_parameters(times_in_delayed, charges_in_delayed, x, y, z, initial_guess):
    """
    Optimize parameters (r_s, T0, lambda_decay) to maximize likelihood.
    :param times_in_delayed: Hit times [t_hit^(i)] (1D numpy array)
    :param charges_in_delayed: Hit charges [q_hit^(i)] (1D numpy array)
    :param x: X coordinates of PMT positions (1D numpy array)
    :param y: Y coordinates of PMT positions (1D numpy array)
    :param z: Z coordinates of PMT positions (1D numpy array)
    :param initial_guess: Initial guess for parameters [x, y, z, T0, lambda_decay]
    :return: Optimized parameters (numpy array)
    """
    result = minimize(
        total_likelihood,
        initial_guess,
        args=(times_in_delayed, charges_in_delayed, x, y, z),
        method='L-BFGS-B',  # Bound optimization
        bounds=[(-10, 10), (-10, 10), (-10, 10), (0, 1), (1e-4, 10)]  # Example bounds for params
    )
    return result.x  # Return optimized parameters

def scintillation_candidates(times, charges, mpmt, pmt, g=geo):

    x_pmt, y_pmt , z_pmt, cables = getxyz(g, mpmt, pmt)

    initial_guess = [0, 0, 0, 0.001, 0.1]  # Initial guess for [x, y, z, T0, lambda_decay]

    # Perform optimization
    optimized_parameters = optimize_parameters(times, charges, x_pmt, y_pmt, z_pmt, initial_guess)

    print("Optimized Parameters:", optimized_parameters)
    print("Optimized Source Vertex (r_s):", optimized_parameters[:3])
    print("Optimized T0:", optimized_parameters[3])
    print("Optimized Lambda Decay:", optimized_parameters[4])
    
    return optimized_parameters
            

