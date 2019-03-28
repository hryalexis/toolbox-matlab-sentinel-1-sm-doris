# The Toolbox to compute the Sentinel-1 A/B in StripMap mode with Matlab and Doris 5.0

  The Sentinel-1 SM Toolbox is a set of Matlab functions to easly compute the Sentinel-1 A/B interferograms in StripMap mode. The objective is to have a simpliest file of parameters to run Doris 5.0. 
  
The Toolbox contains: 
  - computation_Sentinel_SM.m: main function
  - read_data_Sentinel_SM.m
  - readSent1Data.m 
  - getTIFFinfo.m
  - detection_orbit_files.m 
  - cmap_sar.csv
  - card_raw.input
  
The others requirements:
  - xml2struct.m
  - ToolBox Envi for Matlab to create the UTM products
  - ToolBox Mapping of Matlab to create the DEM geometry products
  
System requirements: 
  - Linux / Mac (using of the bash language)

# Installation
```sh
Installation of Doris 5.0: 
```


# Authors / developers

The contributors are:
  - Alexis Hrysiewicz (Laboratoire Magmas et Volcans)
  - Delphine Smittarello (IsTerre)
  
# License
The current license of the software is LGPL v3.0.
