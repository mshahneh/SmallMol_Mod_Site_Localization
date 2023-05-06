
from base_tester import multiprocessing_wrapper, create_modifiSiteLoc

args = {}

def compare_score_row(index, data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers, matches_array):
    """
    Compares the different scoring functions.
    """
    # site_locator = create_modifiSiteLoc(matches_array[index][0], matches_array[index][1], args, data_dict_filtered, cachedStructures_filtered)
    return len(data_dict_filtered[matches_array[index][0]]['peaks_json'])

if __name__ == '__main__':
    res = multiprocessing_wrapper(compare_score_row, library = "BERKELEY-LAB", accepted_adduct = "M+H", count = 200, processes = 8)
    print(res)
