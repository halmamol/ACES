ACES: AmBe Coincidence Event Selection for WCTE Experiment 

1) Filtering_spills_update.py: save a pkl with the data filtered as datos_filtrados.pkl inside the directory Filtered_data

2) NeutronDetection_update2.py: search for promtp + neutron candidates and save its times in the directory csv_saveData/Neutron_candidates

3) jupyter_nb/CandidateEvents.ipynb: calculate deltaT and make the fit in the notebook reading the neutron candidates csv files


*functions_spills.py: should include all the functions for filtering, in reality its all the functions I created in a first period of time
*functions_analysis.py: should include all the functions realted to the search of prompt-neutron candidate. In reality its all the functions I created in a second period of time
*functions_bonsai.py: include all the functions related to BONSAI implementation of the code

Code based on: https://github.com/carla-garcia/DIPC_HK
