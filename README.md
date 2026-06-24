# ACES: AmBe Coincidence Event Selection for WCTE Experiment 

Requirements of:
- ROOT (works with ROOT 6.26.04)
- BONSAI (using WCTE version from https://github.com/WCTE/hk-BONSAI) -> Only used for running BONSAI functions
- WCTE/Geometry (using WCTE version from https://github.com/WCTE/Geometry) -> Needed to run multilaterator functions. Note there are some paths pointing to this repository hardcoded on the code. 

The code consists on two parts: 

- Filtering Spills Code: It uses <Filtering_spills.py> to filter spills on the selected `run_number`. Filtered data is saved on a pkl to analised later.
- Neutron Selection Code: Takes the pkl after filtering the spills, and select coincidence events for selected `run_number`. Output is on a DataFrame with important information. Selected events can be checked using <jupyter_nb/CandidateEvents.ipynb>. Same type of files are created using *_MC* for simulations. 

*`functions_spills.py`: should include all the functions for filtering, in reality its all the functions I created in a first period of time

*`functions_analysis.py`: should include all the functions realted to the search of prompt-neutron candidate. In reality its all the functions I created in a second period of time

*`functions_bonsai.py`: include all the functions related to BONSAI implementation of the code

*`functions_multilateration.py`: include all the functions related to multilateration implementation of the code from Dean (https://github.com/WCTE/TimeCal/blob/main/TimeCal/TC_Multilaterator.py), some functions use the BONSAI file (`functions_bonsai.py`)

Code based on: https://github.com/carla-garcia/DIPC_HK
