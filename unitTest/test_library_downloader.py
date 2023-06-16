import unittest
import json
import os
import shutil
from library_downloader import json_to_dict, download, calculate_matches, get_gnps_library


class TestLibraryDownloader(unittest.TestCase):

    # test json_to_dict
    def test_json_to_dict(self):
        # create json variable:
        Jdata = json.loads('[{"key1": "value1", "key2": "value2"}, {"key1": "value3", "key4": "value4"}, {"key1": "value5", "key6": "value6"}]')
        expected = {"value1": {"key1": "value1", "key2": "value2"}, "value3": {"key1": "value3", "key4": "value4"}, "value5": {"key1": "value5", "key6": "value6"}}
        self.assertEqual(json_to_dict(Jdata, "key1"), expected)
    
    # test get_gnps_library
    def test_get_gnps_library(self):
        url = "https://gnps-external.ucsd.edu/gnpslibrary/GNPS-MSMLS.json"
        self.assertIsNotNone(get_gnps_library(url))
    
    # test download
    def test_download(self):
        url = "https://gnps-external.ucsd.edu/gnpslibrary/GNPS-MSMLS.json"
        location = "../data/temp/"
        download(url, location, 0.5, 0.1)
        self.assertTrue(os.path.exists(os.path.join(location,"data_dict_filtered.pkl")))
        self.assertTrue(os.path.exists(os.path.join(location,"matches.pkl")))
        self.assertTrue(os.path.exists(os.path.join(location,"cachedStructures.pkl")))

        #delete location folder
        shutil.rmtree(location)

    
    # # test calculate_matches
    # def test_calculate_matches(self):
    #     # create json variable:
    #     Jdata = json.loads('[{"key1": "value1", "key2": "value2"}, {"key1": "value3", "key4": "value4"}, {"key1": "value5", "key6": "value6"}]')
    #     data_dict_filtered = json_to_dict(Jdata, "key1")
    #     cachedStructures_filtered = {}
    #     matches = calculate_matches(data_dict_filtered, cachedStructures_filtered)
    #     self.assertEqual(matches, {})
    

if __name__ == '__main__':
    print('Running unit tests')
    unittest.main()


