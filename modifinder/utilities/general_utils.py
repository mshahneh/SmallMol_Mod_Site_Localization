"""
General utility functions
"""
from pyteomics import mgf
import pandas as pd
import numpy as np

def is_shifted(val1:float, val2:float, ppm:float=None, mz_tol:float=None) -> bool:
    """
    Determine if two values differ by more than a specified tolerance.
    
    The function checks if the absolute difference between two values exceeds either a given parts per million (ppm) value or a mass/charge (m/z) tolerance. 
    If only one of ppm or mz_tol is provided, the function uses that value for comparison. 
    If both are provided, the function checks both conditions and returns True if either condition is satisfied.
    
    Parameters:
        :val1 (float): The first value to compare.
        :val2 (float): The second value to compare.
        :ppm (float, optional): The parts per million tolerance. Default is None.
        :mz_tol (float, optional): The m/z tolerance. Default is None.
    
    Returns:
        :bool: True if the values differ by more than the specified tolerance, False otherwise.
    
    Raises:
        :ValueError: If neither ppm nor mz_tol is provided.
    """
    diff = abs(val1 - val2)
    if ppm is None:
        if mz_tol is None:
            raise ValueError("Either ppm or mz_tol must be provided")
        return diff > mz_tol
    else:
        if mz_tol is None:
            return diff > max(val1, val2) * ppm / 1e6
        else:
            return diff > mz_tol or diff > max(val1, val2) * ppm / 1e6


def read_mgf(mgf_path: str) -> pd.DataFrame:
    """
    Read an MGF file into a pandas DataFrame
    
    input:
        :mgf_path: path to the MGF file
    return: 
        :pd.DataFrame: pandas DataFrame with columns as metadata and 'spectrum' as the m/z and intensity values
    """

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


def write_mgf(msms_df: pd.DataFrame, mgf_path: str):
    """
    Writes a pandas DataFrame to an MGF file

    input:
        :msms_df: pandas DataFrame with column 'spectrum' and other columns as metadata
        :mgf_path: path to write the MGF file
    Returns:
      None
    """

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