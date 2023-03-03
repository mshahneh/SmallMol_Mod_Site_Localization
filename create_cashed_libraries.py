import library_downloader as ld

libraries = {
    "GNPS-MSMLS": "https://external.gnps2.org/gnpslibrary/GNPS-MSMLS.json",
    "GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE": "https://external.gnps2.org/gnpslibrary/GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE.json",
    "GNPS-LIBRARY": "https://gnps-external.ucsd.edu/gnpslibrary/GNPS-LIBRARY.json"
}

for library_name in libraries:
    print("Downloading library: " + library_name)
    data_dict_filtered, matches, cachedStructures_filtered = ld.download(libraries[library_name], "data/libraries/"+library_name+"/", 0.5, 0.1)
    print("")