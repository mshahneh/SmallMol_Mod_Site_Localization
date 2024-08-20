from pyteomics import mgf
import pandas as pd
import numpy as np

def is_shifted(val1, val2, ppm = None, mz_tol = None):
    """
    Check if two values are shifted by a given ppm value or mz_tol value
    """
    diff = abs(val1 - val2)
    if ppm is None and mz_tol is None:
        raise ValueError("Either ppm or mz_tol must be provided")
    if mz_tol is not None and diff > mz_tol:
        return True
    if ppm is not None and diff > max(val1, val2) * ppm / 1e6:
        return True
    return False

def _read_mgf(mgf_path: str) -> pd.DataFrame:
    msms_df = []
    with mgf.MGF(mgf_path) as reader:
        for spectrum in reader:
            try:
                d = spectrum['params']
                d['spectrum'] = np.array([spectrum['m/z array'],
                                        spectrum['intensity array']])
                if 'precursor_mz' not in d:
                    d['precursor_mz'] = d['pepmass'][0]
                else:
                    d['precursor_mz'] = float(d['precursor_mz'])
                msms_df.append(d)
            except Exception as e:
                print(e)

    msms_df = pd.DataFrame(msms_df)
    if 'precursor_mz' in msms_df.columns and 'scans' in msms_df.columns:
        msms_df['precursor_mz'] = msms_df['precursor_mz'].astype(float)
        msms_df['scans'] = msms_df['scans'].astype(int)
    return msms_df

def _write_mgf(msms_df: pd.DataFrame, mgf_path: str):
    specs = []
    for i, row in msms_df.iterrows():
        spectrum = {
            'params': row.drop('spectrum').to_dict(),
            'm/z array': row['spectrum'][0],
            'intensity array': row['spectrum'][1]
        }
        specs.append(spectrum)
    with open(mgf_path, 'w') as out:
        mgf.write(specs, out)