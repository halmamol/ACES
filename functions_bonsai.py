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
            

