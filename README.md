# The Toolbox to compute the Sentinel-1 A/B in StripMap mode with Matlab and Doris 5.0

  The Sentinel-1 SM Toolbox is a set of Matlab functions to easly compute the Sentinel-1 A/B interferograms in StripMap mode. The objective is to have a simpliest file of parameters to run Doris 5.0. 
  
The Toolbox contains: 
  - computation_Sentinel_SM.m: main function ot run the interferogram computation
  - read_data_Sentinel_SM.m: function to read the S1 data and the orbit files, to create the file for Doris 5.0
  - detection_orbit_files.m: function to detect the good orbit file in one directory 
  - readSent1Data.m: function to read the SLC images by Louis-Philippe Rousseau
  - getTIFFinfo.m: function to read the tiff informations by Louis-Philippe Rousseau
  - cmap_sar.csv: colormap for the interfergrams
  - sar.m colormap for the interfergrams (in the Matlab format)
  - card_raw.input: input card to run the computation_Sentinel_SM.m function. 
  
The others requirements:
  - xml2struct.m
  - ToolBox Envi for Matlab to create the UTM products
  - ToolBox Mapping of Matlab to create the DEM geometry products
  - The readSent1Data.m function to read the SLC images by Louis-Philippe Rousseau
  - The getTIFFinfo.m function to read the tiff informations by Louis-Philippe Rousseau

The links to find the Matlab files are: 
  - https://fr.mathworks.com/matlabcentral/fileexchange/28518-xml2struct
  - https://fr.mathworks.com/matlabcentral/fileexchange/27172-envi-file-reader-writer
  
The functions by Louis-Philippe Rousseau can be fund in his Github: 
  - https://github.com/lprouss/sent1-L1-utilities
  
System requirements: 
  - Linux / Mac (using of the bash language)

# Installation 

1) Installation of Doris 5.0: https://github.com/TUDelftGeodesy/Doris
2) Download of the Toolbox
3) Download of the files in Mathworks website

# How to compute a S1 StripMap interferogram

1) Download of the two S1 StripMap SLC
2) Download of the orbit files (optionnal)
3) Creating of the full input card or the incomplet input card
4) Running of Matlab 

```sh
computation_computation_Sentinel_SM('full_path_of_input_car')
```
or
```sh
computation_Sentinel_SM('path_input_card.input',path_EXEC,path_MASTER,path_SLAVE,pol)
```

# Authors / developers

The contributors are:
  - Alexis Hrysiewicz (Laboratoire Magmas et Volcans)
  - Delphine Smittarello (IsTerre)
  
# License
The current license of the software is LGPL v3.0.
