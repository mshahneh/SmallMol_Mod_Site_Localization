args = {
    'adduct': {"name": "adduct", "short":'add', "defult":"M+H"},
    'USI1': {"name": "USI1", "short":'USI1', "defult":"mzspec:GNPS:TASK-5700dee92610412ea452a4262add2b93-f.MSV000086107/ccms_peak/VVP3-2_EtAc_MeOh.mzML:scan:2303"},
    'USI2': {"name": "USI2", "short":'USI2', "defult":"mzspec:GNPS:TASK-5700dee92610412ea452a4262add2b93-f.MSV000086107/ccms_peak/VVP77_MeOh.mzML:scan:3024"},
    'SMILES1': {"name": "SMILES1", "short":'SMILES1', "defult":"OC1=CC=CC=C1C2=NC(C3N(C)C(C(O)=O)CS3)CS2"}, 
    'SMILES2': {"name": "SMILES2", "short":'SMILES2', "defult":None},
    'mz_tolerance': {"name": "mz_tolerance", "short":'mz', "defult":0.05},
    'ppm': {"name": "ppm", "short":'ppm', "defult":1.01},
    'filter_peaks_method': {"name": "filter_peaks_method", "short":'filter', "defult":"top_k"},
    'filter_peaks_variable': {"name": "filter_peaks_variable", "short":'filter_var', "defult":50},
    'presence_only': {"name": "presence_only", "short":'presence', "defult":True},
    'consider_intensity': {"name": "consider_intensity", "short":'intensity', "defult":True},
    'shifted_only': {"name": "shifted_only", "short":'shifted', "defult":True},
    'USIhelp': {"name": "USIhelp", "short":'USIhelp', "defult":None},
    'SMILEShelp': {"name": "SMILEShelp", "short":'SMILEShelp', "defult":None},
}