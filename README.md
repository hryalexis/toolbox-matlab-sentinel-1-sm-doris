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
  
The other requirements:
  - xml2struct.m
  - ToolBox Envi for Matlab to create the UTM products (optionnal)
  - ToolBox Mapping of Matlab to create the DEM geometry products (optionnal)
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
3) Download of the files from the Mathworks website

# How to compute a S1 StripMap interferogram

1) Download of two S1 StripMap images 
2) Download of the orbit files (optionnal)
3) Creating of the full input card or the incomplet input card
4) In Matlab:

```sh
computation_computation_Sentinel_SM('full_path_of_input_car')
```
or
```sh
computation_Sentinel_SM('path_input_card.input',path_EXEC,path_MASTER,path_SLAVE,pol)
```

# An input card example:

```sh
***************************************************************************
***************************************************************************
Parameters to compute the Sentinel 1 StripMap Interferogram
***************************************************************************
***************************************************************************
By Alexis Hrysiewicz (LMV / OPGC / UCA / INSU)

***************************************************************************
Doris Parameters
***************************************************************************
Path_of_Doris_processor:            /usr/local/Doris_5 #Directory of Doris core
Path_of_Doris_function:             /usr/local/Doris_utility #Directory of Cpxfiddle function
Path_of_MATLAB_Sentinel_Toolbox:    ???? #Directory of Sentinel SM Toolbox
***************************************************************************
Global parameters
***************************************************************************
Path_for_the_excecution:    ???? #Where do you can compute the interferogram
Data_of_the_Master:         ???? #Full path of the Master (.zip or .SAFE) 
Data_of_the_Slave:          ???? #Full path of the Slave (.zip or .SAFE) 
Choice_of_the_polarition:   vv #vv, hh, vh, or hv 

Path_of_orbits:             ???? #Directory of the orbit file 

Name_of_the_input:          date #date or #orbits
Colormap_of_the_display:    /home/alexis/Documents/MATLAB/Sentinel_SM_toolbox/cmap_sar.csv
RAM_Memory_(MB):            4000

***************************************************************************
Parameters of the DEM
***************************************************************************
Path_of_the_DEM:            ???? #Full path of the DEM file 
Number_of_lines:            1000
Number_of_pixels:           1000
Latitude_DEM:               0
Longitude_DEM:              0
Delta_Latitude:             0.1
Delta_Longitude:            0.1
Value_of_NO_data:           0

***************************************************************************
Processing Parameters #See Doris manuel for these parameters
***************************************************************************
Extraction ------------------------------- YES
First_line_master:     1
Last_line_master:      19680
First_pixel_master:    1
Last_pixel_master:     10120
First_line_slave:      1
Last_line_slave:       19680
First_pixel_slave:     1
Last_pixel_slave:      10120

Multilooking ----------------------------- YES
In_Azimuth:     2
In_Range:       2

Master_Timing ---------------------------- YES
MTE_NWIN        30
MTE_INITOFF     0 0
MTE_WINSIZE     1024 512

Coarse_Correlation ----------------------- YES
CC_NWIN         21
CC_WINSIZE      1024 512
CC_INITOFF      orbit

Fine_Correlation ------------------------- YES
The default parameters are fine. 

DEM_Correlation -------------------------- NO
For the next version

Resample_of_the_slave -------------------- YES
RS_METHOD       rc12p

Interferogram_formation ------------------ YES
FE_DEGREE       5
FE_NPOINTS      501
SRP_METHOD      polynomial

Computation_of_topo_Phase ---------------- YES
The default parameters are fine.

Substraction_of_the_topo_Phase ----------- YES
The default parameters are fine.

Filter_of_the_phase  --------------------- YES
PF_METHOD       goldstein
PF_ALPHA        0.125
PF_BLOCKSIZE    64
PF_OVERLAP      16

Coherence_formation ---------------------- YES
COH_WINSIZE     3 3

Geolocalisation  ------------------------- YES
On the DEM grid

***************************************************************************
Finalisation of the processing 
***************************************************************************
Removing_the_unused_files: YES

GEOTIFF_creating  ------------------------ YES
Latitude_extention:     0   0                              
Longitude_extention:    1   1  
Step_latitude:          0.1         
Step_longitude:         0.1
Method_Geotiff:         natural
Mode_O2I:               YES        dem.hdr #For the SNOV OI2
```

# Authors

The authors are:
  - Alexis Hrysiewicz (Laboratoire Magmas et Volcans)
  - Delphine Smittarello (IsTerre)
  
The authors thank Louis-Philippe Rousseau for the reading functions of the Snetinel-1 SLC. 
