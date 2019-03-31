function computation_Sentinel_SM(varargin)
% Major function to compute the interferogramme with Doris: 
%  The function permit to link Matlab to Doris. All interferometric
%  computation will be made with Doris and the managmenent of the
%  computation the reading of the data (SLC) and the DEM geometry products
%  using Matlab language. 
%
% The objectif of this function is the simplicity of the use. The inputs
% must be: 
% - "path_input_cart" with the parameters of the computation. 
% - "path_input_cart without the paths", "path of the executation directory",
% "path of the master", "path of the slave", "polarisation". This solution
% can be use for automatic processes. 
%
% The input card must be just a txt file. 
%
% EXAMPLES: 
% computation_Sentinel_SM('path_input_card.input')
% OR 
% computation_Sentinel_SM('path_input_card.input',path_EXEC,path_MASTER,path_SLAVE,pol)
%
% THE PATH of the INPUT_CARD must be complet ! ! 
%
% Required functions / programs: 
%   - detection_orbits_files.m
%   - getTIFFinfo.m by Louis-Philippe Rousseau
%   - read_data_Sentinel_SM.m
%   - readSent1Data.m by Louis-Philippe Rousseau
%   - xml2struct.m 
%   - ToolBox Envi for Matlab to create the UTM products.
%   - Doris in version > 4.6 of course. 
%   - ToolBox Mapping of Matlab to create the DEM geometry products
%
% System requirements: 
%   - Linux / Mac (using of the bash language)
%
% Information about the developpement: 
%Author: Alexis Hrysiewicz (Laboratoire Magmas et Volcans / OPGC / OI2)
%Updates: 
%       - Version 0.9 Alpha: creating the function (??)
%       - Version 1.0 Beta: Fixed some problems (??)
%       - Version 1.5 Beta: Implementation of the precise/restitued orbits
%       (October 2018)
%       - Version 1.5.5 Beta: Fixed the bugs with multilooking windows
%       (October 2018)
%       - Version 1.5.6 Beta: Fixed the bugs with the detection of the
%       orbit files (November 2018). 
%
% PLEASE READ the README for more information. 

%% Display informations:
current_directory = cd;
fprintf(1,'***********************************************************\n');
fprintf(1,'Computation of the SENTINEL 1 StripMap inferterogram\n');
fprintf(1,'***********************************************************\n');
fprintf(1,'For Doris 5.03 Beta (TU DELFT)\n');
fprintf(1,'Version 1.5.6 Beta: November 2018\n');
fprintf(1,'\n');

%% Read of the input file
if nargin == 1
    file_input = varargin{1};
    
    [a,path_Doris] = system(['grep ''Path_of_Doris_processor:'' ',file_input,' | awk ''END {print $2}''']); path_Doris = strtrim(path_Doris);
    [a,path_Doris_function] = system(['grep ''Path_of_Doris_function:'' ',file_input,' | awk ''END {print $2}''']); path_Doris_function = strtrim(path_Doris_function);
    
    [a,path_Tool] = system(['grep ''Path_of_MATLAB_Sentinel_Toolbox:'' ',file_input,' | awk ''END {print $2}''']); path_Tool = strtrim(path_Tool);
    addpath(path_Tool);
    
    [a,path_Exec] = system(['grep ''Path_for_the_excecution:'' ',file_input,' | awk ''END {print $2}''']); path_Exec = strtrim(path_Exec);
    [a,path_master] = system(['grep ''Data_of_the_Master:'' ',file_input,' | awk ''END {print $2}''']); path_master = strtrim(path_master);
    [a,path_slave] = system(['grep ''Data_of_the_Slave:'' ',file_input,' | awk ''END {print $2}''']); path_slave = strtrim(path_slave);
    
    [a,pol] = system(['grep ''Choice_of_the_polarition:'' ',file_input,' | awk ''END {print $2}''']); pol = strtrim(pol);
    
elseif nargin == 5
    
    file_input = varargin{1};
    
    [a,path_Doris] = system(['grep ''Path_of_Doris_processor:'' ',file_input,' | awk ''END {print $2}''']); path_Doris = strtrim(path_Doris);
    [a,path_Doris_function] = system(['grep ''Path_of_Doris_function:'' ',file_input,' | awk ''END {print $2}''']); path_Doris_function = strtrim(path_Doris_function);
    
    [a,path_Tool] = system(['grep ''Path_of_MATLAB_Sentinel_Toolbox:'' ',file_input,' | awk ''END {print $2}''']); path_Tool = strtrim(path_Tool);
    addpath(path_Tool);
    
    path_Exec = varargin{2};
    path_master = varargin{3};
    path_slave = varargin{4};
    pol = varargin{5};
    
    fprintf(1,'Excecution in:\n %s \n',path_Exec);
    fprintf(1,'The Master is:\n %s \n',path_master);
    fprintf(1,'The Slave is:\n %s \n',path_slave);
    fprintf(1,'The Polarisation is:\n %s \n',pol);
else
    error('Problem with the Input number');
end

[a,path_orb] = system(['grep ''Path_of_orbits:'' ',file_input,' | awk ''END {print $2}''']); path_orb = strtrim(path_orb);

[a,mod_save] = system(['grep ''Name_of_the_input:'' ',file_input,' | awk ''END {print $2}''']); mod_save = strtrim(mod_save);
[a,path_color] = system(['grep ''Colormap_of_the_display:'' ',file_input,' | awk ''END {print $2}''']); path_color = strtrim(path_color);
[a,val_ram] = system(['grep ''RAM_Memory_(MB):'' ',file_input,' | awk ''END {print $2}''']); val_ram = strtrim(val_ram);

fprintf(1,'Excecution in:\n %s \n',path_Exec);
fprintf(1,'The Master is:\n %s \n',path_master);
fprintf(1,'The Slave is:\n %s \n',path_slave);
fprintf(1,'The Polarisation is:\n %s \n',pol);

fprintf(1,'\nFor the Orbits:\n');
if strcmp(path_orb,'NONE') == 1
    fprintf(1,' Using of the Orbits in the SAFE directory\n');
    path_orb_master = 'NONE';
    path_orb_slave = 'NONE';
else
    fprintf(1,' Precise and Restitued Orbits: used\n');
    path_orb_master = path_orb;
    path_orb_slave = path_orb;
end

[a,path_dem] = system(['grep ''Path_of_the_DEM:'' ',file_input,' | awk ''END {print $2}''']); path_dem = strtrim(path_dem);
[a,nbl_dem] = system(['grep ''Number_of_lines:'' ',file_input,' | awk ''END {print $2}''']); nbl_dem = strtrim(nbl_dem);
[a,nbc_dem] = system(['grep ''Number_of_pixels:'' ',file_input,' | awk ''END {print $2}''']); nbc_dem = strtrim(nbc_dem);

[a,lat_dem] = system(['grep ''Latitude_DEM:'' ',file_input,' | awk ''END {print $2}''']); lat_dem = strtrim(lat_dem);
[a,lon_dem] = system(['grep ''Longitude_DEM:'' ',file_input,' | awk ''END {print $2}''']); lon_dem = strtrim(lon_dem);
[a,dlat_dem] = system(['grep ''Delta_Latitude:'' ',file_input,' | awk ''END {print $2}''']); dlat_dem = strtrim(dlat_dem);
[a,dlon_dem] = system(['grep ''Delta_Longitude:'' ',file_input,' | awk ''END {print $2}''']); dlon_dem = strtrim(dlon_dem);
[a,no_data_dem] = system(['grep ''Value_of_NO_data:'' ',file_input,' | awk ''END {print $2}''']); no_data_dem = strtrim(no_data_dem);

fprintf(1,'\nFor the DEM:\n');
fprintf(1,' The DEM is %s\n %s of lines and %s of pixels\n',path_dem,nbl_dem,nbc_dem);

fprintf(1,'\nRead of the Processing parameters:\n');

fprintf(1,'\tFor the extraction of the SLC:\n');
[a,val_extraction] = system(['grep ''Extraction'' ',file_input,' | awk ''END {print $3}''']); val_extraction = strtrim(val_extraction);
[a,first_l_m] = system(['grep ''First_line_master:'' ',file_input,' | awk ''END {print $2}''']); first_l_m = strtrim(first_l_m);
[a,last_l_m] = system(['grep ''Last_line_master:'' ',file_input,' | awk ''END {print $2}''']); last_l_m = strtrim(last_l_m);
[a,first_p_m] = system(['grep ''First_pixel_master:'' ',file_input,' | awk ''END {print $2}''']); first_p_m = strtrim(first_p_m);
[a,last_p_m] = system(['grep ''Last_pixel_master:'' ',file_input,' | awk ''END {print $2}''']); last_p_m = strtrim(last_p_m);

[a,first_l_s] = system(['grep ''First_line_slave:'' ',file_input,' | awk ''END {print $2}''']); first_l_s = strtrim(first_l_s);
[a,last_l_s] = system(['grep ''Last_line_slave:'' ',file_input,' | awk ''END {print $2}''']); last_l_s = strtrim(last_l_s);
[a,first_p_s] = system(['grep ''First_pixel_slave:'' ',file_input,' | awk ''END {print $2}''']); first_p_s = strtrim(first_p_s);
[a,last_p_s] = system(['grep ''Last_pixel_slave:'' ',file_input,' | awk ''END {print $2}''']); last_p_s = strtrim(last_p_s);

if strcmp(lower(val_extraction),'yes')==1
    fprintf(1,'\t\tExtraction: %s\n','OK');
    fprintf(1,'\t\tFirst line MASTER: %s\n',first_l_m);
    fprintf(1,'\t\tLast line MASTER: %s\n',last_l_m);
    fprintf(1,'\t\tFirst pixel MASTER: %s\n',first_p_m);
    fprintf(1,'\t\tLast pixel MASTER: %s\n',last_p_m);
    fprintf(1,'\t\tFirst line SLAVE: %s\n',first_l_s);
    fprintf(1,'\t\tLast line SLAVE: %s\n',last_l_s);
    fprintf(1,'\t\tFirst pixel SLAVE: %s\n',first_p_s);
    fprintf(1,'\t\tLast pixel SLAVE: %s\n',last_p_s);
else
    warning('Extraction: NO')
end

fprintf(1,'\tFor multilooking:\n');
[a,val_multilooking] = system(['grep ''Multilooking'' ',file_input,' | awk ''END {print $3}''']); val_multilooking = strtrim(val_multilooking);
[a,b] = system(['grep ''In_Azimuth:'' ',file_input,' | awk ''END {print $2}''']); ml_kernel(1)=str2num(b);
[a,b] = system(['grep ''In_Range:'' ',file_input,' | awk ''END {print $2}''']); ml_kernel(2)=str2num(b);
if strcmp(lower(val_multilooking),'yes')==1
    fprintf(1,'\t\tMultilooking: %s\n','OK');
    fprintf(1,'\t\tThe kernel is %d %d.\n',ml_kernel(1),ml_kernel(2));
else
    fprintf(1,'\t\tMultilooking: %s\n','NO');
end

fprintf(1,'\tFor Master Timing:\n');
[a,val_Master_Timing] = system(['grep ''Master_Timing'' ',file_input,' | awk ''END {print $3}''']); val_Master_Timing = strtrim(val_Master_Timing);
[a,MTE_NWIN] = system(['grep ''MTE_NWIN'' ',file_input,' | awk ''END {print $2}''']); MTE_NWIN = strtrim(MTE_NWIN);
[a,MTE_INITOFF] = system(['grep ''MTE_INITOFF'' ',file_input,' | awk ''END {print $2}''']); MTE_INITOFF = strtrim(MTE_INITOFF);
[a,b] = system(['grep ''MTE_INITOFF'' ',file_input,' | awk ''END {print $3}''']); MTE_INITOFF = [MTE_INITOFF,' ',strtrim(b)];
[a,MTE_WINSIZE] = system(['grep ''MTE_WINSIZE'' ',file_input,' | awk ''END {print $2}''']); MTE_WINSIZE = strtrim(MTE_WINSIZE);
[a,b] = system(['grep ''MTE_WINSIZE'' ',file_input,' | awk ''END {print $3}''']); MTE_WINSIZE = [MTE_WINSIZE,' ',strtrim(b)];
if strcmp(lower(val_Master_Timing),'yes')==1
    fprintf(1,'\t\tMaster Timing: %s\n','OK');
else
    fprintf(1,'\t\tMaster Timing: %s\n','NO');
end

fprintf(1,'\tFor Coarse Correlation:\n');
[a,val_Coarse_Correlation] = system(['grep ''Coarse_Correlation'' ',file_input,' | awk ''END {print $3}''']); val_Coarse_Correlation = strtrim(val_Coarse_Correlation);
[a,CC_NWIN] = system(['grep ''CC_NWIN'' ',file_input,' | awk ''END {print $2}''']); CC_NWIN = strtrim(CC_NWIN);
[a,CC_WINSIZE] = system(['grep ''CC_WINSIZE'' ',file_input,' | awk ''END {print $2}''']); CC_WINSIZE = strtrim(CC_WINSIZE);
[a,b] = system(['grep ''CC_WINSIZE'' ',file_input,' | awk ''END {print $3}''']); CC_WINSIZE = [CC_WINSIZE,' ',strtrim(b)];
[a,CC_INITOFF] = system(['grep ''CC_INITOFF'' ',file_input,' | awk ''END {print $2}''']); CC_INITOFF = strtrim(CC_INITOFF);
if strcmp(lower(val_Coarse_Correlation),'yes')==1
    fprintf(1,'\t\tCoarse Correlation: %s\n','OK');
else
    error('The Coarse Correlation is the necessary step !');
end

fprintf(1,'\tFor Fine Correlation:\n');
[a,val_Fine_Correlation] = system(['grep ''Fine_Correlation'' ',file_input,' | awk ''END {print $3}''']); val_Fine_Correlation = strtrim(val_Fine_Correlation);
if strcmp(lower(val_Fine_Correlation),'yes')==1
    fprintf(1,'\t\tFine Correlation: %s\n','OK');
else
    %error('The Fine Correlation is the necessary step !');
end

fprintf(1,'\tFor DEM Correlation:\n');
[a,val_DEM_Correlation] = system(['grep ''DEM_Correlation'' ',file_input,' | awk ''END {print $3}''']); val_DEM_Correlation = strtrim(val_DEM_Correlation);
if strcmp(lower(val_DEM_Correlation),'yes')==1
    fprintf(1,'\t\tDEM Correlation unused\n');
else
    fprintf(1,'\t\tDEM Correlation unused\n');
end

fprintf(1,'\tFor Resampling:\n');
[a,val_Resample] = system(['grep ''Resample_of_the_slave'' ',file_input,' | awk ''END {print $3}''']); val_Resample = strtrim(val_Resample);
[a,RS_METHOD] = system(['grep ''RS_METHOD'' ',file_input,' | awk ''END {print $2}''']); RS_METHOD = strtrim(RS_METHOD);
if strcmp(lower(val_Resample),'yes')==1
    fprintf(1,'\t\tResampling: %s\n','OK');
    fprintf(1,'\t\tThe Resampling method is %s\n',RS_METHOD);
else
    %error('The Resampling step is the necessary step !');
end

fprintf(1,'\tFor Interferogram Formation:\n');
[a,val_IFG] = system(['grep ''Interferogram_formation'' ',file_input,' | awk ''END {print $3}''']); val_IFG = strtrim(val_IFG);
[a,FE_DEGREE] = system(['grep ''FE_DEGREE'' ',file_input,' | awk ''END {print $2}''']); FE_DEGREE = strtrim(FE_DEGREE);
[a,FE_NPOINTS] = system(['grep ''FE_NPOINTS'' ',file_input,' | awk ''END {print $2}''']); FE_NPOINTS = strtrim(FE_NPOINTS);
[a,SRP_METHOD] = system(['grep ''SRP_METHOD'' ',file_input,' | awk ''END {print $2}''']); SRP_METHOD = strtrim(SRP_METHOD);
if strcmp(lower(val_IFG),'yes')==1
    fprintf(1,'\t\tInterferogram Formation: %s\n','OK');
else
    %error('The Interferogram Formation step is the necessary step !');
end

fprintf(1,'\tFor Computation of the topo phase:\n');
[a,val_COMPREFDEM] = system(['grep ''Computation_of_topo_Phase'' ',file_input,' | awk ''END {print $3}''']); val_COMPREFDEM = strtrim(val_COMPREFDEM);
if strcmp(lower(val_COMPREFDEM),'yes')==1
    fprintf(1,'\t\tComputation of the topo phase: %s\n','OK');
else
    fprintf(1,'\t\tComputation of the topo phase: %s\n','NO');
end

fprintf(1,'\tFor Substraction of the topo phase:\n');
[a,val_SUBPREFDEM] = system(['grep ''Substraction_of_the_topo_Phase'' ',file_input,' | awk ''END {print $3}''']); val_SUBPREFDEM = strtrim(val_SUBPREFDEM);
if strcmp(lower(val_SUBPREFDEM),'yes')==1
    fprintf(1,'\t\tSubstraction of the topo phase: %s\n','OK');
else
    fprintf(1,'\t\tSubstraction of the topo phase: %s\n','NO');
end

fprintf(1,'\tFor Filter of the phase:\n');
[a,val_Filtphase] = system(['grep ''Filter_of_the_phase'' ',file_input,' | awk ''END {print $3}''']); val_Filtphase = strtrim(val_Filtphase);
[a,PF_METHOD] = system(['grep ''PF_METHOD'' ',file_input,' | awk ''END {print $2}''']); PF_METHOD = strtrim(PF_METHOD);
[a,PF_ALPHA] = system(['grep ''PF_ALPHA'' ',file_input,' | awk ''END {print $2}''']); PF_ALPHA = strtrim(PF_ALPHA);
[a,PF_BLOCKSIZE] = system(['grep ''PF_BLOCKSIZE'' ',file_input,' | awk ''END {print $2}''']); PF_BLOCKSIZE = strtrim(PF_BLOCKSIZE);
[a,PF_OVERLAP] = system(['grep ''PF_OVERLAP'' ',file_input,' | awk ''END {print $2}''']); PF_OVERLAP = strtrim(PF_OVERLAP);
if strcmp(lower(val_Filtphase),'yes')==1
    fprintf(1,'\t\tFilter of the phase: %s\n','OK');
    fprintf(1,'\t\tThe method is %s, with a kernel %s.\n',PF_METHOD,PF_BLOCKSIZE);
    fprintf(1,'\t\tThe power is %s, with an overlap %s.\n',PF_ALPHA,PF_OVERLAP);
else
    fprintf(1,'\t\tFilter of the topo phase: %s\n','NO');
end

fprintf(1,'\tFor Coherence Formation:\n');
[a,val_coh] = system(['grep ''Coherence_formation'' ',file_input,' | awk ''END {print $3}''']); val_coh = strtrim(val_coh);
[a,COH_WINSIZE] = system(['grep ''COH_WINSIZE'' ',file_input,' | awk ''END {print $2}''']); COH_WINSIZE = strtrim(COH_WINSIZE);
[a,b] = system(['grep ''COH_WINSIZE'' ',file_input,' | awk ''END {print $3}''']); COH_WINSIZE = [COH_WINSIZE,' ',strtrim(b)];
if strcmp(lower(val_coh),'yes')==1
    fprintf(1,'\t\tCoherence formation: %s\n','OK');
    fprintf(1,'\t\tThe windows is %s.\n',COH_WINSIZE);
else
    fprintf(1,'\t\tCoherence formation: %s\n','NO');
end

fprintf(1,'\tFor Geolocalisation:\n');
[a,val_geo] = system(['grep ''Geolocalisation'' ',file_input,' | awk ''END {print $3}''']); val_geo = strtrim(val_geo);
if strcmp(lower(val_geo),'yes')==1
    if strcmp(lower(val_COMPREFDEM),'yes')==1
        fprintf(1,'\t\tGeolocalisation on the DEM: %s\n','OK');
    else
        error('Computation of the DEM Radar is NO');
    end
else
    fprintf(1,'\t\tGeolocalisation on the DEM: %s\n','NO');
end

[a,val_remove] = system(['grep ''Removing_the_unused_files:'' ',file_input,' | awk ''END {print $2}''']); val_remove = strtrim(val_remove);
if strcmp(lower(val_remove),'yes')==1
    fprintf(1,'Removing of the unused files: %s\n','OK');
else
    fprintf(1,'Removing of the unused files: %s\n','NO');
end

[a,val_geotiff] = system(['grep ''GEOTIFF_creating'' ',file_input,' | awk ''END {print $3}''']); val_geotiff = strtrim(val_geotiff);
if strcmp(lower(val_geotiff),'yes')==1
    fprintf(1,'Geotiff formation: %s\n','OK');
    v = ver;
    if any(strcmp(cellstr(char(v.Name)),'Mapping Toolbox')) == 1
        fprintf(1,'\tMapping Toolbox: OK\n');
        [status,errmsg] = license('checkout','Map_Toolbox');
        if status == 0
            error('Need to Mapping Toolbox');
        end
    else
        error('Need to Mapping Toolbox');
    end
else
    fprintf(1,'Geotiff formation: %s\n','NO');
end

%fprintf(1,'\nChecking of the parameters...\n'); I must to do it...

%% Creating of the Exec directory
if exist(path_Exec) == 0
    system(['mkdir ',path_Exec]);
    cd(path_Exec);
else
    cd(path_Exec);
end

%% Extraction of the SLC files
%For the Master
pos_pt = find(path_master=='.');
if strcmp(path_master(pos_pt(end)+1:end),'SAFE') == 1
    fprintf(1,'\nThe Master Format is .SAFE.\n');
    
    cd([path_master,'/measurement']);
    d = dir;
    for i1 = 1 : size(d,1)
        namei = d(i1).name;
        if isempty(strfind(namei,lower(pol)))==0
            slc_master = [path_master,'/measurement/',namei];
        end
    end
    cd([path_master,'/annotation']);
    d = dir;
    for i1 = 1 : size(d,1)
        namei = d(i1).name;
        if isempty(strfind(namei,lower(pol)))==0
            header_master = [path_master,'/annotation/',namei];
        end
    end
    
    fprintf(1,'Master SLC: \n%s\n',slc_master);
    fprintf(1,'Master Header: \n%s\n',header_master);
    
elseif strcmp(path_master(pos_pt(end)+1:end),'zip') == 1
    fprintf(1,'\nThe Master Format is .zip.\n');
    
    system(['unzip ',path_master,' -d ',path_Exec,'/slc_master/']);
    cd([path_Exec,'/slc_master/']);
    [a,b] = system('ls');
    b = strtrim(b);
    cd([path_Exec,'/slc_master/',b,'/measurement']);
    d = dir;
    for i1 = 1 : size(d,1)
        namei = d(i1).name;
        if isempty(strfind(namei,lower(pol)))==0
            slc_master = [path_Exec,'/slc_master/',b,'/measurement/',namei];
        end
    end
    cd([path_Exec,'/slc_master/',b,'/annotation']);
    d = dir;
    for i1 = 1 : size(d,1)
        namei = d(i1).name;
        if isempty(strfind(namei,lower(pol)))==0
            header_master = [path_Exec,'/slc_master/',b,'/annotation/',namei];
        end
    end
    
    fprintf(1,'Master SLC: \n%s\n',slc_master);
    fprintf(1,'Master Header: \n%s\n',header_master);
    
else
    error('Unknowed format for the Master');
end

%For the Slave
pos_pt = find(path_slave=='.');
if strcmp(path_slave(pos_pt(end)+1:end),'SAFE') == 1
    fprintf(1,'\nThe Slave Format is .SAFE.\n');
    
    cd([path_slave,'/measurement']);
    d = dir;
    for i1 = 1 : size(d,1)
        namei = d(i1).name;
        if isempty(strfind(namei,lower(pol)))==0
            slc_slave = [path_slave,'/measurement/',namei];
        end
    end
    cd([path_slave,'/annotation']);
    d = dir;
    for i1 = 1 : size(d,1)
        namei = d(i1).name;
        if isempty(strfind(namei,lower(pol)))==0
            header_slave = [path_slave,'/annotation/',namei];
        end
    end
    
    fprintf(1,'Slave SLC: \n%s\n',slc_slave);
    fprintf(1,'Slave Header: \n%s\n',header_slave);
    
elseif strcmp(path_slave(pos_pt(end)+1:end),'zip') == 1
    fprintf(1,'\nThe Slave Format is .zip.\n');
    
    system(['unzip ',path_slave,' -d ',path_Exec,'/slc_slave/']);
    cd([path_Exec,'/slc_slave/']);
    [a,b] = system('ls');
    b = strtrim(b);
    cd([path_Exec,'/slc_slave/',b,'/measurement']);
    d = dir;
    for i1 = 1 : size(d,1)
        namei = d(i1).name;
        if isempty(strfind(namei,lower(pol)))==0
            slc_slave = [path_Exec,'/slc_slave/',b,'/measurement/',namei];
        end
    end
    cd([path_Exec,'/slc_slave/',b,'/annotation']);
    d = dir;
    for i1 = 1 : size(d,1)
        namei = d(i1).name;
        if isempty(strfind(namei,lower(pol)))==0
            header_slave = [path_Exec,'/slc_slave/',b,'/annotation/',namei];
        end
    end
    
    fprintf(1,'Slave SLC: \n%s\n',slc_slave);
    fprintf(1,'Slave Header: \n%s\n',header_slave);
    
else
    error('Unknowed format for the Master');
end

if isempty(slc_master)==1 | isempty(header_master)==1 | isempty(slc_slave)==1 | isempty(header_slave)==1
    error('Problem with the SLC or Header path...');
end

cd(current_directory);

%% READING OF THE MASTER
fprintf(1,'\nReading of the Master DATA...\n');
read_data_Sentinel_SM(slc_master,header_master,path_Exec,'Master',val_extraction,first_l_m,last_l_m,first_p_m,last_p_m,path_orb_master);
fprintf(1,'Reading of the Master DATA: OK\n');

pos_pt = find(path_master=='.');
if strcmp(path_master(pos_pt(end)+1:end),'zip') == 1
    cd(path_Exec);
    !rm -R slc_master
    cd(current_directory);
end

%% READING OF THE SLAVE
fprintf(1,'\nReading of the Slave DATA...\n');
read_data_Sentinel_SM(slc_slave,header_slave,path_Exec,'Slave',val_extraction,first_l_s,last_l_s,first_p_s,last_p_s,path_orb_slave);
fprintf(1,'Reading of the Slave DATA: OK\n');

pos_pt = find(path_slave=='.');
if strcmp(path_slave(pos_pt(end)+1:end),'zip') == 1
    cd(path_Exec);
    !rm -R slc_slave
    cd(current_directory);
end

%% Computation of the file name
[a,date_master] = system(['grep ''First_pixel_azimuth_time (UTC):'' ',path_Exec,'/master.res | awk ''END {print $3}''']); date_master = strtrim(date_master);
date_master = datestr(date_master,'yyyymmdd');
[a,orb_abs_master] = system(['grep ''Scene identification:'' ',path_Exec,'/master.res | awk ''END {print $4}''']); orb_abs_master = strtrim(orb_abs_master);
[~,swath] = system(['grep ''SWATH:'' ',path_Exec,'/master.res | awk ''END {print $2}''']); swath = strtrim(swath);
[a,plt_master] = system(['grep ''Sensor platform mission identifer:'' ',path_Exec,'/master.res | awk ''END {print $5}''']); plt_master = strtrim(plt_master);

[a,date_slave] = system(['grep ''First_pixel_azimuth_time (UTC):'' ',path_Exec,'/slave.res | awk ''END {print $3}''']); date_slave = strtrim(date_slave);
date_slave = datestr(date_slave,'yyyymmdd');
[a,orb_abs_slave] = system(['grep ''Scene identification:'' ',path_Exec,'/slave.res | awk ''END {print $4}''']); orb_abs_slave = strtrim(orb_abs_slave);
[a,swath] = system(['grep ''SWATH:'' ',path_Exec,'/slave.res | awk ''END {print $2}''']); swath = strtrim(swath);
[a,plt_slave] = system(['grep ''Sensor platform mission identifer:'' ',path_Exec,'/slave.res | awk ''END {print $5}''']); plt_slave = strtrim(plt_slave);

%ifg
if strcmp(mod_save,'date')==1
    if strcmp(lower(val_multilooking),'yes')==1
        name_ifg = ['ifg_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',date_master,'_',date_slave,'_ml_',num2str(ml_kernel(1)),num2str(ml_kernel(2))];
        name_ifg_filt = ['ifg_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',date_master,'_',date_slave,'_ml_',num2str(ml_kernel(1)),num2str(ml_kernel(2)),'_filt'];
        name_coh = ['coh_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',date_master,'_',date_slave,'_ml_',num2str(ml_kernel(1)),num2str(ml_kernel(2))];
        name_lat = ['lat_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',date_master,'_',date_slave,'_ml_',num2str(ml_kernel(1)),num2str(ml_kernel(2))];
        name_lon = ['lon_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',date_master,'_',date_slave,'_ml_',num2str(ml_kernel(1)),num2str(ml_kernel(2))];
    else
        name_ifg = ['ifg_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',date_master,'_',date_slave];
        name_ifg_filt = ['ifg_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',date_master,'_',date_slave,'_filt'];
        name_coh = ['coh_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',date_master,'_',date_slave];
        name_lat = ['lat_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',date_master,'_',date_slave];
        name_lon = ['lon_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',date_master,'_',date_slave];
    end
else
    if strcmp(lower(val_multilooking),'yes')==1
        name_ifg = ['ifg_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',orb_abs_master,'_',orb_abs_slave,'_ml_',num2str(ml_kernel(1)),num2str(ml_kernel(2))];
        name_ifg_filt = ['ifg_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',orb_abs_master,'_',orb_abs_slave,'_ml_',num2str(ml_kernel(1)),num2str(ml_kernel(2)),'_filt'];
        name_coh = ['coh_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',orb_abs_master,'_',orb_abs_slave,'_ml_',num2str(ml_kernel(1)),num2str(ml_kernel(2))];
        name_lat = ['lat_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',orb_abs_master,'_',orb_abs_slave,'_ml_',num2str(ml_kernel(1)),num2str(ml_kernel(2))];
        name_lon = ['lon_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',orb_abs_master,'_',orb_abs_slave,'_ml_',num2str(ml_kernel(1)),num2str(ml_kernel(2))];
    else
        name_ifg = ['ifg_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',orb_abs_master,'_',orb_abs_slave];
        name_ifg_filt = ['ifg_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',orb_abs_master,'_',orb_abs_slave,'_filt'];
        name_coh = ['coh_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',orb_abs_master,'_',orb_abs_slave];
        name_lat = ['lat_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',orb_abs_master,'_',orb_abs_slave];
        name_lon = ['lon_',plt_master,'_',plt_slave,'_SM_',swath,'_',pol,'_',orb_abs_master,'_',orb_abs_slave];
    end
end

%% Master Timing STEP
if strcmp(lower(val_Master_Timing),'yes')==1
    f_temp = fopen([path_Exec,'/input_temp.dorisin'],'w');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,'c * Master Timing\n');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment  ___general options___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'c SCREEN          debug                         // level of output to standard out\n');
    fprintf(f_temp,'SCREEN          info                          // level of output to standard out\n');
    fprintf(f_temp,'MEMORY          %s                             // MB\n',val_ram);
    fprintf(f_temp,'BEEP            error                		 // level of beeping\n');
    fprintf(f_temp,'OVERWRITE                                       // overwrite existing files\n');
    fprintf(f_temp,'BATCH                                           // non-interactive\n');
    fprintf(f_temp,'PREVIEW ON                                 // prevents copy of this file to log\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'PROCESS          M_SIMAMP\n');
    fprintf(f_temp,'PROCESS          M_TIMING\n');
    fprintf(f_temp,'c                                               //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,' comment  ___the general io files___            //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,'LOGFILE         log.out                         // log file\n');
    fprintf(f_temp,'M_RESFILE       master.res                      // parameter file\n');
    fprintf(f_temp,'S_RESFILE       slave.res                       // parameter file\n');
    fprintf(f_temp,'I_RESFILE       timing.out                         // parameter file\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'comment ___SIMAMP___\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'SAM_IN_FORMAT   real4\n');
    fprintf(f_temp,'SAM_IN_DEM      %s\n',path_dem);
    fprintf(f_temp,'SAM_IN_SIZE     %s %s              // rows cols\n',nbl_dem,nbc_dem);
    fprintf(f_temp,'SAM_IN_DELTA    %s %s    // in degrees       \n',dlat_dem,dlon_dem);
    fprintf(f_temp,'SAM_IN_UL       %s %s     // lat and lon of upper left\n',lat_dem,lon_dem);
    fprintf(f_temp,'SAM_IN_NODATA   %s\n',no_data_dem);
    fprintf(f_temp,'SAM_OUT_FILE    master_sam.raw // synthetic amplitude\n');
    fprintf(f_temp,'SAM_OUT_DEM     dem_sam.raw    // cropped DEM for\n');
    fprintf(f_temp,'c       ___                           ___\n');
    fprintf(f_temp,'comment ___COMPUTE MASTER TIMING ERROR___\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'MTE_METHOD      magfft         // computes faster than magspace\n');
    fprintf(f_temp,'MTE_NWIN        %s             // number of large windows\n',MTE_NWIN);
    fprintf(f_temp,'MTE_INITOFF     %s            // initial offset\n',MTE_INITOFF);
    fprintf(f_temp,'MTE_WINSIZE     %s      // rectangular window\n',MTE_WINSIZE);
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'STOP\n');
    fclose(f_temp);
    
    cd(path_Exec);
    system([path_Doris,'/doris input_temp.dorisin']);
    !rm input_temp.dorisin
    cd(current_directory);
end

%% Coarse Correlation STEP
if strcmp(lower(val_Coarse_Correlation),'yes')==1
    f_temp = fopen([path_Exec,'/input_temp.dorisin'],'w');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,'c ***\n');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment  ___general options___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'c SCREEN          debug                           // level of output to standard out\n');
    fprintf(f_temp,'SCREEN          info                           // level of output to standard out\n');
    fprintf(f_temp,'MEMORY          %s                             // MB\n',val_ram);
    fprintf(f_temp,'BEEP            error                            // level of beeping\n');
    fprintf(f_temp,'OVERWRITE                                       // overwrite existing files\n');
    fprintf(f_temp,'BATCH                                           // non-interactive\n');
    fprintf(f_temp,'c LISTINPUT OFF                                 // prevents copy of this file to log\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'PROCESS          COARSEORB\n');
    fprintf(f_temp,'PROCESS          COARSECORR\n');
    fprintf(f_temp,'c                                              //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,' comment  ___the general io files___            //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,'LOGFILE         log.out                         // log file\n');
    fprintf(f_temp,'M_RESFILE       master.res  // parameter file\n');
    fprintf(f_temp,'S_RESFILE       slave.res                       // parameter file\n');
    fprintf(f_temp,'I_RESFILE       ifg.out              // parameter file\n');
    fprintf(f_temp,'DUMPBASELINE    50 50\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment ___COARSE CORR (COREGISTRATION)___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'CC_METHOD       magfft                          // default\n');
    fprintf(f_temp,'CC_NWIN         %s                             // number of windows\n',CC_NWIN);
    fprintf(f_temp,'CC_WINSIZE      %s                       // size of windows\n',CC_WINSIZE);
    fprintf(f_temp,'CC_INITOFF      %s                           // use result of orbits for initial offset\n',CC_INITOFF);
    fprintf(f_temp,'STOP\n');
    fclose(f_temp);
    
    cd(path_Exec);
    system([path_Doris,'/doris input_temp.dorisin']);
    !rm input_temp.dorisin
    cd(current_directory);
end

%% Fine Correlation STEP
if strcmp(lower(val_Fine_Correlation),'yes')==1
    f_temp = fopen([path_Exec,'/input_temp.dorisin'],'w');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,'c ***\n');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment  ___general options___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'c SCREEN          debug                           // level of output to standard out\n');
    fprintf(f_temp,'SCREEN          info                           // level of output to standard out\n');
    fprintf(f_temp,'c SCREEN          error                           // level of output to standard out\n');
    fprintf(f_temp,'MEMORY          %s                             // MB\n',val_ram);
    fprintf(f_temp,'BEEP            error                            // level of beeping\n');
    fprintf(f_temp,'OVERWRITE                                       // overwrite existing files\n');
    fprintf(f_temp,'BATCH                                           // non-interactive\n');
    fprintf(f_temp,'c LISTINPUT OFF                                 // prevents copy of this file to log\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'PROCESS          FINE\n');
    fprintf(f_temp,'PROCESS          COREGPM\n');
    fprintf(f_temp,'c                                              //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,' comment  ___the general io files___            //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,'LOGFILE         log.out                         // log file\n');
    fprintf(f_temp,'M_RESFILE       master.res  // parameter file\n');
    fprintf(f_temp,'S_RESFILE       slave.res                       // parameter file\n');
    fprintf(f_temp,'I_RESFILE       ifg.out               // parameter file\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment ___FINE COREGISTRATION___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'FC_METHOD       oversample                          //\n');
    fprintf(f_temp,'c FC_METHOD       magfft                          //\n');
    fprintf(f_temp,'c FC_METHOD     magspace                        //\n');
    fprintf(f_temp,'c FC_NWIN         8000                             // number of windows\n');
    fprintf(f_temp,'c FC_IN_POS       fc_pos.in                // file containing position of windows\n');
    fprintf(f_temp,'FC_WINSIZE      64 64                           // size of windows\n');
    fprintf(f_temp,'FC_ACC          8 8                             \n');
    fprintf(f_temp,'FC_INITOFF      coarsecorr                      // use result of coarse to compute first\n');
    fprintf(f_temp,'FC_OSFACTOR     32                              // oversampling factor\n');
    fprintf(f_temp,'c FC_PLOT         0.65 BG\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment ___COMPUTE COREGISTRATION PARAMETERS___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'CPM_THRESHOLD   0.30\n');
    fprintf(f_temp,'CPM_DEGREE      2\n');
    fprintf(f_temp,'c CPM_WEIGHT      bamler                          // none\n');
    fprintf(f_temp,'c CPM_WEIGHT      linear                          // none\n');
    fprintf(f_temp,'CPM_WEIGHT      quadratic                          // none\n');
    fprintf(f_temp,'CPM_MAXITER     8000\n');
    fprintf(f_temp,'c CPM_PLOT        BG\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'STOP\n');
    fclose(f_temp);
    
    cd(path_Exec);
    system([path_Doris,'/doris input_temp.dorisin']);
    !rm input_temp.dorisin
    cd(current_directory);
end

%% Resample STEP
if strcmp(lower(val_Resample),'yes')==1
    f_temp = fopen([path_Exec,'/input_temp.dorisin'],'w');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,'c ***  Doris inputfile generated by: run at: Nov 27, 2000 (Monday) *****\n');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,'c ***\n');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment  ___general options___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'c SCREEN          debug                           // level of output to standard out\n');
    fprintf(f_temp,'SCREEN          info                           // level of output to standard out\n');
    fprintf(f_temp,'MEMORY          %s                             // MB\n',val_ram);
    fprintf(f_temp,'BEEP            error                            // level of beeping\n');
    fprintf(f_temp,'OVERWRITE                                       // overwrite existing files\n');
    fprintf(f_temp,'BATCH                                           // non-interactive\n');
    fprintf(f_temp,'c LISTINPUT OFF                                 // prevents copy of this file to log\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'PROCESS          RESAMPLE\n');
    fprintf(f_temp,'c                                              //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,' comment  ___the general io files___            //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,'LOGFILE         log.out                         // log file\n');
    fprintf(f_temp,'M_RESFILE       master.res  // parameter file\n');
    fprintf(f_temp,'S_RESFILE       slave.res                       // parameter file\n');
    fprintf(f_temp,'I_RESFILE       ifg.out               // parameter file\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment ___RESAMPLING SLAVE___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'RS_METHOD       %s\n',RS_METHOD);
    fprintf(f_temp,'RS_OUT_FILE     slave_res.slc\n');
    fprintf(f_temp,'RS_OUT_FORMAT   cr4\n');
    fprintf(f_temp,'STOP\n');
    fclose(f_temp);
    
    cd(path_Exec);
    system([path_Doris,'/doris input_temp.dorisin']);
    !rm input_temp.dorisin
    cd(current_directory);
end

%% Interferogram STEP
if strcmp(lower(val_IFG),'yes')==1
    f_temp = fopen([path_Exec,'/input_temp.dorisin'],'w');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,'c ***\n');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment  ___general options___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'c SCREEN          debug                           // level of output to standard out\n');
    fprintf(f_temp,'SCREEN          info                           // level of output to standard out\n');
    fprintf(f_temp,'MEMORY          %s                             // MB\n',val_ram);
    fprintf(f_temp,'OVERWRITE                                       // overwrite existing files\n');
    fprintf(f_temp,'BATCH                                           // non-interactive\n');
    fprintf(f_temp,'c LISTINPUT OFF                                 // prevents copy of this file to log\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'PROCESS          INTERFERO\n');
    fprintf(f_temp,'PROCESS          COMPREFPHA\n');
    fprintf(f_temp,'PROCESS          SUBTRREFPHA\n');
    fprintf(f_temp,'c                                              //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,' comment  ___the general io files___            //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,'LOGFILE         log.out                         // log file\n');
    fprintf(f_temp,'M_RESFILE       master.res  // parameter file\n');
    fprintf(f_temp,'S_RESFILE       slave.res                       // parameter file\n');
    fprintf(f_temp,'I_RESFILE       ifg.out               // parameter file\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment ___interferogram generation___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'INT_OUT_CINT    cint.raw                 // optional\n');
    fprintf(f_temp,'INT_MULTILOOK   1 1                            // line, pixel\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment ___ COMPREFPHA ___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'FE_METHOD       porbits\n');
    fprintf(f_temp,'FE_DEGREE       %s\n',FE_DEGREE);
    fprintf(f_temp,'FE_NPOINTS      %s\n',FE_NPOINTS);
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'c ___ step subtrrefpha ___\n');
    fprintf(f_temp,'c \n');
    fprintf(f_temp,'SRP_METHOD      %s\n',SRP_METHOD);
    fprintf(f_temp,'SRP_OUT_CINT    cint.minrefpha.raw\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'STOP\n');
    fclose(f_temp);
    
    cd(path_Exec);
    system([path_Doris,'/doris input_temp.dorisin']);
    !rm input_temp.dorisin
    cd(current_directory);
end

%% DEM phase STEP
if strcmp(lower(val_COMPREFDEM),'yes')==1
    f_temp = fopen([path_Exec,'/input_temp.dorisin'],'w');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,'c ***\n');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment  ___general options___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'c SCREEN          debug                         // level of output to standard out\n');
    fprintf(f_temp,'SCREEN          info                          // level of output to standard out\n');
    fprintf(f_temp,'c SCREEN          error                           // level of output to standard out\n');
    fprintf(f_temp,'MEMORY          %s                            // MB\n',val_ram);
    fprintf(f_temp,'BEEP            error                            // level of beeping\n');
    fprintf(f_temp,'OVERWRITE                                       // overwrite existing files\n');
    fprintf(f_temp,'BATCH                                           // non-interactive\n');
    fprintf(f_temp,'c LISTINPUT OFF                                 // prevents copy of this file to log\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'PROCESS          COMPREFDEM\n');
    fprintf(f_temp,'c                                               //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,' comment  ___the general io files___            //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,'LOGFILE         log.out                         // log file\n');
    fprintf(f_temp,'M_RESFILE       master.res                      // parameter file\n');
    fprintf(f_temp,'S_RESFILE       slave.res                       // parameter file\n');
    fprintf(f_temp,'I_RESFILE       ifg.out                        // parameter file\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'c ___ step comprefdem ___\n');
    fprintf(f_temp,'CRD_INCLUDE_FE  OFF                     // phase w.r.t. ellipsoid\n');
    fprintf(f_temp,'CRD_OUT_FILE    refdem_1l.raw           //\n');
    fprintf(f_temp,'CRD_OUT_DEM_LP  dem_radar.raw\n');
    fprintf(f_temp,'CRD_IN_FORMAT   real4\n');
    fprintf(f_temp,'CRD_IN_DEM      %s\n',path_dem);
    fprintf(f_temp,'CRD_IN_SIZE     %s %s               // rows cols\n',nbl_dem,nbc_dem);
    fprintf(f_temp,'CRD_IN_DELTA    %s %s    // in degrees       \n',dlat_dem,dlon_dem);
    fprintf(f_temp,'CRD_IN_UL       %s %s     // lat and lon of upper left\n',lat_dem,lon_dem);
    fprintf(f_temp,'CRD_IN_NODATA   %s\n',no_data_dem);
    fprintf(f_temp,'CRD_OUT_FILE    master_sam.raw // synthetic amplitude\n');
    fprintf(f_temp,'CRD_OUT_DEM     dem_sam.raw    // cropped DEM for\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'STOP\n');
    fclose(f_temp);
    
    cd(path_Exec);
    system([path_Doris,'/doris input_temp.dorisin']);
    !rm input_temp.dorisin
    cd(current_directory);
end

%% Substraction of DEM phase STEP
if strcmp(lower(val_SUBPREFDEM),'yes')==1
    f_temp = fopen([path_Exec,'/input_temp.dorisin'],'w');
    fprintf(f_temp,'c ***\n');
    fprintf(f_temp,'c         ___           ___ \n');
    fprintf(f_temp,'comment   ___SUBTRREFDEM___ \n');
    fprintf(f_temp,'c                             \n');
    fprintf(f_temp,'MEMORY              	%s                	 // MB\n',val_ram);
    fprintf(f_temp,'BEEP                	error               	 // level of beeping\n');
    fprintf(f_temp,'OVERWRITE           	on                  	 // overwrite existing files\n');
    fprintf(f_temp,'PREVIEW             	off                  	 // on\n');
    fprintf(f_temp,'BATCH               	on                	 // non-interactive\n');
    fprintf(f_temp,'LISTINPUT           	on                  	 // prevents copy of this file to log\n');
    fprintf(f_temp,'SCREEN              	info                	 // level of output to standard out\n');
    fprintf(f_temp,'c \n');
    fprintf(f_temp,'PROCESS        SUBTRREFDEM \n');
    fprintf(f_temp,'c \n');
    fprintf(f_temp,'LOGFILE             	log.out             	 // log file\n');
    fprintf(f_temp,'I_RESFILE           	ifg.out            	 // interferogram parameter file\n');
    fprintf(f_temp,'M_RESFILE           	master.res          	 // master parameter file\n');
    fprintf(f_temp,'S_RESFILE           	slave.res           	 // slave parameter file\n');
    fprintf(f_temp,'HEIGHT              	0.0                 	 // average WGS84 height\n');
    fprintf(f_temp,'ORB_INTERP          	POLYFIT             	 // orbit interpolation method\n');
    fprintf(f_temp,'ELLIPSOID           	WGS84               	 // WGS84, GRS80, BESSEL or define major and minor axis\n');
    fprintf(f_temp,'c                             \n');
    fprintf(f_temp,'SRD_OUT_CINT        	cint.minrefdem.raw        \n');
    fprintf(f_temp,'SRD_OFFSET          	0 0                 \n');
    fprintf(f_temp,'STOP                          \n');
    fclose(f_temp);
    
    cd(path_Exec);
    system([path_Doris,'/doris input_temp.dorisin']);
    !rm input_temp.dorisin
    cd(current_directory);
end

%% Filter Phase STEP
if strcmp(lower(val_Filtphase),'yes')==1
    f_temp = fopen([path_Exec,'/input_temp.dorisin'],'w');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,'c ***\n');
    fprintf(f_temp,'c **********************************************************************\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,' comment  ___general options___\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'c SCREEN          debug                           // level of output to standard out\n');
    fprintf(f_temp,'c SCREEN          info                           // level of output to standard out\n');
    fprintf(f_temp,'SCREEN          error                           // level of output to standard out\n');
    fprintf(f_temp,'MEMORY          %s                            // MB\n',val_ram);
    fprintf(f_temp,'OVERWRITE         on                              // overwrite existing files\n');
    fprintf(f_temp,'BATCH                                           // non-interactive\n');
    fprintf(f_temp,'c LISTINPUT OFF                                 // prevents copy of this file to log\n');
    fprintf(f_temp,'c\n');
    fprintf(f_temp,'PROCESS         filtphase\n');
    fprintf(f_temp,'c                                              //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,' comment  ___the general io files___            //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,'LOGFILE         log.out                         // log file\n');
    fprintf(f_temp,'M_RESFILE       master.res  // parameter file\n');
    fprintf(f_temp,'S_RESFILE       slave.res                       // parameter file\n');
    fprintf(f_temp,'I_RESFILE       ifg.out               // parameter file\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,' c                                              //\n');
    fprintf(f_temp,' comment PHASEFILT\n');
    fprintf(f_temp,' c\n');
    fprintf(f_temp,'PF_METHOD      %s\n',PF_METHOD);
    fprintf(f_temp,'PF_OUT_FILE    cint.filt.raw\n');
    fprintf(f_temp,'PF_ALPHA       %s\n',PF_ALPHA);
    fprintf(f_temp,'PF_BLOCKSIZE   %s\n',PF_BLOCKSIZE);
    fprintf(f_temp,'PF_OVERLAP     %s\n',PF_OVERLAP);
    fprintf(f_temp,'STOP\n');
    fclose(f_temp);
    
    cd(path_Exec);
    system([path_Doris,'/doris input_temp.dorisin']);
    !rm input_temp.dorisin
    cd(current_directory);
end

%% Coherence STEP
if strcmp(lower(val_coh),'yes')==1
    
    if strcmp(lower(val_multilooking),'yes')==1
        ml_kernel_bis = ml_kernel;
    else
        ml_kernel_bis = [1 1];
    end
    
    f_temp = fopen([path_Exec,'/input_temp.dorisin'],'w');
    fprintf(f_temp,'c ********\n');
    fprintf(f_temp,'c         ___         ___ \n');
    fprintf(f_temp,'comment   ___COHERENCE___ \n');
    fprintf(f_temp,'c                             \n');
    fprintf(f_temp,'MEMORY              	%s                	 // MB\n',val_ram);
    fprintf(f_temp,'BEEP                	error               	 // level of beeping\n')
    fprintf(f_temp,'OVERWRITE           	on                  	 // overwrite existing files\n');
    fprintf(f_temp,'PREVIEW             	on                  	 // on\n');
    fprintf(f_temp,'BATCH               	                   	 // non-interactive\n');
    fprintf(f_temp,'SCREEN              	info                	 // level of output to standard out\n');
    fprintf(f_temp,'c \n');
    fprintf(f_temp,'PROCESS        COHERENCE \n');
    fprintf(f_temp,'c \n');
    fprintf(f_temp,'LOGFILE             	log.out             	 // log file\n');
    fprintf(f_temp,'I_RESFILE           	ifg.out            	 // interferogram parameter file\n');
    fprintf(f_temp,'M_RESFILE           	master.res          	 // master parameter file\n');
    fprintf(f_temp,'S_RESFILE           	slave.res           	 // slave parameter file\n');
    fprintf(f_temp,'HEIGHT              	0.0                 	 // average WGS84 height\n');
    fprintf(f_temp,'ORB_INTERP          	POLYFIT             	 // orbit interpolation method\n');
    fprintf(f_temp,'ELLIPSOID           	WGS84               	 // WGS84, GRS80, BESSEL or define major and minor axis\n');
    fprintf(f_temp,'c                             \n');
    if strcmp(lower(val_SUBPREFDEM),'yes')==1
        fprintf(f_temp,'COH_METHOD          	INCLUDE_REFDEM      \n');
    else
        fprintf(f_temp,'COH_METHOD          	REFPHASE_ONLY      \n');
    end
    fprintf(f_temp,'COH_OUT_COH         	%s      	 // output image\n',[name_coh,'.raw']);
    fprintf(f_temp,'COH_MULTILOOK       	%d %d                 \n',ml_kernel_bis(1),ml_kernel_bis(2));
    fprintf(f_temp,'COH_WINSIZE         	%s               \n',COH_WINSIZE);
    fprintf(f_temp,'STOP                          \n');
    fclose(f_temp);
    
    cd(path_Exec);
    system([path_Doris,'/doris input_temp.dorisin']);
    !rm input_temp.dorisin
    cd(current_directory);
end

%% Geolocalisation STEP
if strcmp(lower(val_geo),'yes')==1
    
    %Fake Slant2h
    f_temp = fopen([path_Exec,'/fake_slant2h.out'],'w');
    fprintf(f_temp,' \n');
    fprintf(f_temp,'*****************************************************\n');
    fprintf(f_temp,'*_Start_slant2h:\n');
    fprintf(f_temp,'*****************************************************\n');
    fprintf(f_temp,'Method:                                 schwabisch\n');
    fprintf(f_temp,'Data_output_file:                       dem_radar.raw\n');
    fprintf(f_temp,'Data_output_format:                     real4\n');
    fprintf(f_temp,'First_line (w.r.t. original_master):   %s\n',first_l_m);
    fprintf(f_temp,'Last_line (w.r.t. original_master):    %s\n',last_l_m);
    fprintf(f_temp,'First_pixel (w.r.t. original_master):  %s\n',first_p_m);
    fprintf(f_temp,'Last_pixel (w.r.t. original_master):   %s\n',last_p_m);
    fprintf(f_temp,'Multilookfactor_azimuth_direction:      1\n');
    fprintf(f_temp,'Multilookfactor_range_direction:        1\n');
    fprintf(f_temp,'Ellipsoid (name,a,b):                   WGS84 6.37814e+06 6.35675e+06\n');
    fprintf(f_temp,'*****************************************************\n');
    fprintf(f_temp,'* End_slant2h:_NORMAL\n');
    fprintf(f_temp,'*****************************************************\n');
    fclose(f_temp);
    
    cd(path_Exec);
    system('cat ifg.out fake_slant2h.out > ifg_temp.out');
    !rm ifg.out
    !rm fake_slant2h.out
    !mv ifg_temp.out ifg.out
    cd(current_directory);
    
    %Run geo step
    f_temp = fopen([path_Exec,'/input_temp.dorisin'],'w');
    fprintf(f_temp,'c *****\n');
    fprintf(f_temp,'c         ___       ___ \n');
    fprintf(f_temp,'comment   ___GEOCODE___ \n');
    fprintf(f_temp,'c                             \n');
    fprintf(f_temp,'MEMORY              	%s                	 // MB\n',val_ram);
    fprintf(f_temp,'BEEP                	error               	 // level of beeping\n');
    fprintf(f_temp,'OVERWRITE           	on                  	 // overwrite existing files\n');
    fprintf(f_temp,'PREVIEW             	on                  	 // on\n');
    fprintf(f_temp,'BATCH               	on                  	 // non-interactive\n');
    fprintf(f_temp,'LISTINPUT           	on                  	 // prevents copy of this file to log\n');
    fprintf(f_temp,'SCREEN              	info                	 // level of output to standard out\n');
    fprintf(f_temp,'c \n');
    fprintf(f_temp,'PROCESS        GEOCODE \n');
    fprintf(f_temp,'c \n');
    fprintf(f_temp,'LOGFILE             	log.out             	 // log file\n');
    fprintf(f_temp,'I_RESFILE           	ifg.out            	 // interferogram parameter file\n');
    fprintf(f_temp,'M_RESFILE           	master.res          	 // master parameter file\n');
    fprintf(f_temp,'S_RESFILE           	slave.res           	 // slave parameter file\n');
    fprintf(f_temp,'HEIGHT              	0.0                 	 // average WGS84 height\n');
    fprintf(f_temp,'ORB_INTERP          	POLYFIT             	 // orbit interpolation method\n');
    fprintf(f_temp,'ELLIPSOID           	WGS84               	 // WGS84, GRS80, BESSEL or define major and minor axis\n');
    fprintf(f_temp,'c                             \n');
    fprintf(f_temp,'GEO_OUT_LAM         	lon.raw             	 // longitude coordinates\n');
    fprintf(f_temp,'GEO_OUT_PHI         	lat.raw             	 // latitude coordinates\n');
    fprintf(f_temp,'STOP                          \n');
    fclose(f_temp);
    
    cd(path_Exec);
    system([path_Doris,'/doris input_temp.dorisin']);
    !rm input_temp.dorisin
    cd(current_directory);
end

%% Multilooking and visualisation creating
fprintf(1,'\nMultilooking and visulation creating\n\n');
cd(path_Exec);
nc = str2num(last_p_m)-str2num(first_p_m)+1;
if strcmp(lower(val_multilooking),'yes')==1
    if strcmp(lower(val_Filtphase),'yes')==1
        if isempty(path_color)==0
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M',num2str(ml_kernel(2)),'/',num2str(ml_kernel(1)),' -osunraster -c',path_color,' -e0.2 -s1.8 cint.filt.raw > ',name_ifg_filt,'.ras']);
        else
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M',num2str(ml_kernel(2)),'/',num2str(ml_kernel(1)),' -osunraster -e0.2 -s1.8 cint.filt.raw > ',name_ifg_filt,'.ras']);
        end
        system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qnormal -M',num2str(ml_kernel(2)),'/',num2str(ml_kernel(1)),' -ofloat cint.filt.raw > ',name_ifg_filt,'.raw']);
    end
    if strcmp(lower(val_SUBPREFDEM),'yes')==1
        if isempty(path_color)==0
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M',num2str(ml_kernel(2)),'/',num2str(ml_kernel(1)),' -osunraster -c',path_color,' -e0.2 -s1.8 cint.minrefdem.raw > ',name_ifg,'.ras']);
        else
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M',num2str(ml_kernel(2)),'/',num2str(ml_kernel(1)),' -osunraster -e0.2 -s1.8 cint.minrefdem.raw > ',name_ifg,'.ras']);
        end
        system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qnormal -M',num2str(ml_kernel(2)),'/',num2str(ml_kernel(1)),' -ofloat cint.minrefdem.raw > ',name_ifg,'.raw']);
    else
        if isempty(path_color)==0
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M',num2str(ml_kernel(2)),'/',num2str(ml_kernel(1)),' -osunraster -c',path_color,' -e0.2 -s1.8 cint.minrefpha.raw > ',name_ifg,'.ras']);
        else
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M',num2str(ml_kernel(2)),'/',num2str(ml_kernel(1)),' -osunraster -e0.2 -s1.8 cint.minrefpha.raw > ',name_ifg,'.ras']);
        end
        system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qnormal -M',num2str(ml_kernel(2)),'/',num2str(ml_kernel(1)),' -ofloat cint.minrefpha.raw > ',name_ifg,'.raw']);
    end
    if strcmp(lower(val_geo),'yes')==1
        system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fr4 -qnormal -M',num2str(ml_kernel(2)),'/',num2str(ml_kernel(1)),' -ofloat lat.raw > ',name_lat,'.raw']);
        system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fr4 -qnormal -M',num2str(ml_kernel(2)),'/',num2str(ml_kernel(1)),' -ofloat lon.raw > ',name_lon,'.raw']);
    end
else
    if strcmp(lower(val_Filtphase),'yes')==1
        if isempty(path_color)==0
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M1/1 -osunraster -c',path_color,' -e0.2 -s1.8 cint.filt.raw > ',name_ifg_filt,'.ras']);
        else
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M1/1 -osunraster -e0.2 -s1.8 cint.filt.raw > ',name_ifg_filt,'.ras']);
        end
        system(['mv cint.filt.raw ',name_ifg_filt,'.raw']);
    end
    if strcmp(lower(val_SUBPREFDEM),'yes')==1
        if isempty(path_color)==0
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M1/1 -osunraster -c',path_color,' -e0.2 -s1.8 cint.minrefdem.raw > ',name_ifg,'.ras']);
        else
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M1/1 -osunraster -e0.2 -s1.8 cint.minrefdem.raw > ',name_ifg,'.ras']);
        end
        system(['mv cint.minrefdem.raw ',name_ifg,'.raw']);
    else
        if isempty(path_color)==0
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M1/1 -osunraster -c',path_color,' -e0.2 -s1.8 cint.minrefpha.raw > ',name_ifg,'.ras']);
        else
            system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -fcr4 -qmixed -M1/1 -osunraster -e0.2 -s1.8 cint.minrefpha.raw > ',name_ifg,'.ras']);
        end
        system(['mv cint.minrefpha.raw ',name_ifg,'.raw']);
    end
    if strcmp(lower(val_geo),'yes')==1
        system(['mv lat.raw ',name_lat,'.raw']);
        system(['mv lon.raw ',name_lon,'.raw']);
    end
end

if strcmp(lower(val_coh),'yes')==1
    if strcmp(lower(val_multilooking),'yes')==1
        %system([path_Doris_function,'/cpxfiddle -w',num2str(fix(nc./ml_kernel(2))),' -qnormal -osunraster -b -c gray -r 0.0/1.0  -f r4 -l1 -p1 -P1427 ',name_coh,'.raw > ',name_coh,'.ras']);
        system([path_Doris_function,'/cpxfiddle -w',num2str(fix(nc./ml_kernel(2))),' -qnormal -osunraster -b -c gray -r 0.0/1.0  -f r4 -l1 -p1 ',name_coh,'.raw > ',name_coh,'.ras']);
    else
        %system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -qnormal -osunraster -b -c gray -r 0.0/1.0  -f r4 -l1 -p1 -P1427 ',name_coh,'.raw > ',name_coh,'.ras']);
        system([path_Doris_function,'/cpxfiddle -w',num2str(nc),' -qnormal -osunraster -b -c gray -r 0.0/1.0  -f r4 -l1 -p1 ',name_coh,'.raw > ',name_coh,'.ras']);
    end
end

cd(current_directory);

%% Removing of the files
cd(path_Exec);
if strcmp(lower(val_remove),'yes')==1
    fprintf(1,'\nRemoving of the files\n');
    if exist('cint.minrefdem.raw')==2;
        !rm cint.minrefdem.raw
    end
    if exist('cint.minrefpha.raw')==2;
        !rm cint.minrefpha.raw
    end
    if exist('cint.filt.raw')==2;
        !rm cint.filt.raw
    end
    if exist('cint.raw')==2;
        !rm cint.raw
    end
    if exist('coherence.ras')==2;
        !rm coherence.ras
    end
    if exist('coherence.ras.sh')==2;
        !rm coherence.ras.sh
    end
    if exist('CPM_Data')==2;
        !rm CPM_Data
    end
    if exist('dem_radar.raw')==2;
        !rm dem_radar.raw
    end
    if exist('dem_sam.raw')==2;
        !rm dem_sam.raw
    end
    if exist('image_master.slc')==2;
        !rm image_master.slc
    end
    if exist('image_slave.slc')==2;
        !rm image_slave.slc
    end
    if exist('lat.raw')==2;
        !rm lat.raw
    end
    if exist('lon.raw')==2;
        !rm lon.raw
    end
    if exist('master_sam.raw')==2;
        !rm master_sam.raw
    end
    if exist('master_simamp.ras')==2;
        !rm master_simamp.ras
    end
    if exist('master_simamp.ras.sh')==2;
        !rm master_simamp.ras.sh
    end
    if exist('resample.kml')==2;
        !rm resample.kml
    end
    if exist('sam_m_demline.temp')==2;
        !rm sam_m_demline.temp
    end
    if exist('sam_m_dempixel.tem')==2;
        !rm sam_m_dempixel.tem
    end
    if exist('sam_m_theta.success')==2;
        !rm sam_m_theta.success
    end
    if exist('sam_m_theta.temp')==2;
        !rm sam_m_theta.temp
    end
    if exist('slave_res.slc')==2;
        !rm slave_res.slc
    end
    if exist('rsmp_orig_slave_pixel.raw')==2;
        !rm rsmp_orig_slave_pixel.raw
    end
    if exist('rsmp_orig_slave_line.raw')==2;
        !rm rsmp_orig_slave_line.raw
    end
    if exist('sam_m_dempixel.temp')==2;
        !rm sam_m_dempixel.temp
    end
end
cd(current_directory);

%% Writting of the ifg.out
cd(path_Exec);
fifg = fopen('ifg.out','a');
fprintf(fifg,'\n\n\n  ------------------------------------------------------- \n');
fprintf(fifg,'Output of the Sentinel S1 StripMap computation:\n');
fprintf(fifg,' ------------------------------------------------------- \n');
fprintf(fifg,'Master_information:\t\t\t\t%s\n','master.res');
fprintf(fifg,'Slave_information:\t\t\t\t%s\n','slave.res');
fprintf(fifg,'IFG_information:\t\t\t\t%s\n','ifg.out.res');
if strcmp(lower(val_multilooking),'yes')==1
    fprintf(fifg,'Multilooking:\t\t\t\t%d %d\n',ml_kernel(1),ml_kernel(2));
    fprintf(fifg,'Size:\t\t\t\t%d %d\n',fix((str2num(last_l_m)-str2num(first_l_m)+1)./ml_kernel(1)),fix((str2num(last_p_m)-str2num(first_p_m)+1)./ml_kernel(2)));
else
    fprintf(fifg,'Multilooking:\t\t\t\t%s\n','NONE');
    fprintf(fifg,'Size:\t\t\t\t%d %d\n',fix((str2num(last_l_m)-str2num(first_l_m)+1)),fix((str2num(last_p_m)-str2num(first_p_m)+1)));
end
if strcmp(lower(val_Filtphase),'yes')==1
    fprintf(fifg,'Complex_Filter_IFG:\t\t\t\t%s\n',[name_ifg_filt,'.raw']);
    fprintf(fifg,'Ras_Filter_IFG:\t\t\t\t%s\n',[name_ifg_filt,'.ras']);
else
    fprintf(fifg,'Complex_Filter_IFG:\t\t\t\t%s\n','NONE');
    fprintf(fifg,'Ras_Filter_IFG:\t\t\t\t%s\n','NONE');
end
fprintf(fifg,'Complex_IFG:\t\t\t\t%s\n',[name_ifg,'.raw']);
fprintf(fifg,'Ras_IFG:\t\t\t\t%s\n',[name_ifg,'.ras']);
if strcmp(lower(val_coh),'yes')==1
    fprintf(fifg,'Real4_Coherence:\t\t\t\t%s\n',[name_coh,'.raw']);
    fprintf(fifg,'Ras_Coherence:\t\t\t\t%s\n',[name_coh,'.ras']);
else
    fprintf(fifg,'Real4_Coherence:\t\t\t\t%s\n','NONE');
    fprintf(fifg,'Ras_Coherence:\t\t\t\t%s\n','NONE');
end
if strcmp(lower(val_geo),'yes')==1
    fprintf(fifg,'Lat_file:\t\t\t\t%s\n',[name_lat,'.raw']);
    fprintf(fifg,'Lon_file:\t\t\t\t%s\n',[name_lon,'.raw']);
else
    fprintf(fifg,'Lat_file:\t\t\t\t%s\n','NONE');
    fprintf(fifg,'Lon_file:\t\t\t\t%s\n','NONE');
end
fprintf(fifg,' ------------------------------------------------------- \n');
fprintf(fifg,' Finish at: %s \n',datestr(now));
fprintf(fifg,' ------------------------------------------------------- \n');
fclose(fifg);

cd(current_directory);

%% Geotiff creating
if strcmp(lower(val_geotiff),'yes')==1
    
    [a,mode_OI2] = system(['grep ''Mode_O2I:'' ',file_input,' | awk ''END {print $2}''']); mode_OI2 = strtrim(mode_OI2);
    [a,interp_method] = system(['grep ''Method_Geotiff:'' ',file_input,' | awk ''END {print $2}''']); interp_method = strtrim(interp_method);
    
    if strcmp(lower(mode_OI2),'yes')==1
        fprintf(1,'\nCreating of the data for O2I...\n');
        
        if strcmp(lower(val_Filtphase),'yes')==1
            f_cplx_ifg = fopen([path_Exec,'/',name_ifg_filt,'.raw'],'r');
        else
            f_cplx_ifg = fopen([path_Exec,'/',name_ifg,'.raw'],'r');
        end
        if strcmp(lower(val_multilooking),'yes')==1
            nbl = fix((str2num(last_l_m)-str2num(first_l_m)+1)./ml_kernel(1));
            nbc = fix((str2num(last_p_m)-str2num(first_p_m)+1)./ml_kernel(2));
        else
            nbl = fix((str2num(last_l_m)-str2num(first_l_m)+1));
            nbc = fix((str2num(last_p_m)-str2num(first_p_m)+1));
        end
        cplx_ifg = fread(f_cplx_ifg,[nbc.*2 nbl],'float32')'; fclose(f_cplx_ifg);
        realpart = cplx_ifg(:,1:2:nbc.*2);
        cplxpart = cplx_ifg(:,2:2:nbc.*2);
        cplx_ifg = complex(realpart,cplxpart);
        
        if strcmp(lower(val_coh),'yes')==1
            f_coh = fopen([path_Exec,'/',name_coh,'.raw'],'r'); coh = fread(f_coh,[nbc nbl],'float32')'; fclose(f_coh);
        end
        
        [a,path_hdr_dem] = system(['grep ''Mode_O2I:'' ',file_input,' | awk ''END {print $3}''']); path_hdr_dem = strtrim(path_hdr_dem);
        info_dem = envihdrread(path_hdr_dem);
        
        f_lat = fopen([path_Exec,'/',name_lat,'.raw'],'r'); lat = fread(f_lat,[nbc nbl],'float32')'; fclose(f_lat);
        f_lon = fopen([path_Exec,'/',name_lon,'.raw'],'r'); lon = fread(f_lon,[nbc nbl],'float32')'; fclose(f_lon);
        
        dczone = utmzone(mean(lat(:),'omitnan'),mean(lon(:),'omitnan'));
        utmstruct = defaultm('utm');
        utmstruct.zone = dczone;
        utmstruct.geoid = wgs84Ellipsoid;
        utmstruct = defaultm(utmstruct);
        [x_utm,y_utm] = mfwdtran(utmstruct,lat,lon);
        
        [X,Y] = meshgrid(info_dem.x,info_dem.y);
        
        F = scatteredInterpolant(x_utm(:),y_utm(:),cplx_ifg(:),interp_method,'none');
        cplx_ifg_interp = F(X,Y);
        
        if strcmp(lower(val_coh),'yes')==1
            F = scatteredInterpolant(x_utm(:),y_utm(:),coh(:),interp_method,'none');
            coh_int = F(X,Y);
        end
        ifg_inter = angle(cplx_ifg_interp);
        
        I = uint8(round(((ifg_inter+pi)./(2.*pi)).*256));
        C = uint8(round(((coh_int)).*255));
        
        if strcmp(lower(val_Filtphase),'yes')==1
            f1 = fopen([path_Exec,'/',name_ifg_filt,'.pha'],'w');
        else
            f1 = fopen([path_Exec,'/',name_ifg,'.pha'],'w');
        end
        fwrite(f1,I','uint8'); fclose(f1);
        
        if strcmp(lower(val_coh),'yes')==1
            f2 = fopen([path_Exec,'/',name_coh,'.byt'],'w');
            fwrite(f2,C','uint8'); fclose(f2);
        end
        
        fprintf(1,'\nCreating of the data for O2I: OK\n');
        
    end
    
    fprintf(1,'\nCreating of the GEOTIFF...\n');
    
    if strcmp(lower(val_Filtphase),'yes')==1
        f_cplx_ifg = fopen([path_Exec,'/',name_ifg_filt,'.raw'],'r');
    else
        f_cplx_ifg = fopen([path_Exec,'/',name_ifg,'.raw'],'r');
    end
    if strcmp(lower(val_multilooking),'yes')==1
        nbl = fix((str2num(last_l_m)-str2num(first_l_m)+1)./ml_kernel(1));
        nbc = fix((str2num(last_p_m)-str2num(first_p_m)+1)./ml_kernel(2));
    else
        nbl = fix((str2num(last_l_m)-str2num(first_l_m)+1));
        nbc = fix((str2num(last_p_m)-str2num(first_p_m)+1));
    end
    cplx_ifg = fread(f_cplx_ifg,[nbc.*2 nbl],'float32')'; fclose(f_cplx_ifg);
    realpart = cplx_ifg(:,1:2:nbc.*2);
    cplxpart = cplx_ifg(:,2:2:nbc.*2);
    cplx_ifg = complex(realpart,cplxpart);
    
    if strcmp(lower(val_coh),'yes')==1
        f_coh = fopen([path_Exec,'/',name_coh,'.raw'],'r'); coh = fread(f_coh,[nbc nbl],'float32')'; fclose(f_coh);
    end
    
    f_lat = fopen([path_Exec,'/',name_lat,'.raw'],'r'); lat = fread(f_lat,[nbc nbl],'float32')'; fclose(f_lat);
    f_lon = fopen([path_Exec,'/',name_lon,'.raw'],'r'); lon = fread(f_lon,[nbc nbl],'float32')'; fclose(f_lon);
    
    [a,lat_1] = system(['grep ''Latitude_extention:'' ',file_input,' | awk ''END {print $2}''']); lat_1 = strtrim(lat_1); lat_1 = str2num(lat_1);
    [a,lat_2] = system(['grep ''Latitude_extention:'' ',file_input,' | awk ''END {print $3}''']); lat_2 = strtrim(lat_2); lat_2 = str2num(lat_2);
    
    [a,lon_1] = system(['grep ''Longitude_extention:'' ',file_input,' | awk ''END {print $2}''']); lon_1 = strtrim(lon_1); lon_1 = str2num(lon_1);
    [a,lon_2] = system(['grep ''Longitude_extention:'' ',file_input,' | awk ''END {print $3}''']); lon_2 = strtrim(lon_2); lon_2 = str2num(lon_2);
    
    [a,dlat_geo] = system(['grep ''Step_latitude:'' ',file_input,' | awk ''END {print $2}''']); dlat_geo = strtrim(dlat_geo); dlat_geo = str2num(dlat_geo);
    [a,dlon_geo] = system(['grep ''Step_longitude:'' ',file_input,' | awk ''END {print $2}''']); dlon_geo = strtrim(dlon_geo); dlon_geo = str2num(dlon_geo);
    
    [a,interp_method] = system(['grep ''Method_Geotiff:'' ',file_input,' | awk ''END {print $2}''']); interp_method = strtrim(interp_method);
    
    [LON,LAT] = meshgrid(lon_1:dlon_geo:lon_2,lat_1:dlat_geo:lat_2);
    
    F = scatteredInterpolant(lon(:),lat(:),cplx_ifg(:),interp_method,'none');
    cplx_ifg_interp = F(LON,LAT);
    
    if strcmp(lower(val_coh),'yes')==1
        F = scatteredInterpolant(lon(:),lat(:),coh(:),interp_method,'none');
        coh_int = F(LON,LAT);
    end
    ifg_inter = angle(cplx_ifg_interp);
    
    R = georefcells();
    R.LatitudeLimits = [min(LAT(:)) max(LAT(:))];
    R.LongitudeLimits = [min(LON(:)) max(LON(:))];
    R.RasterSize = [size(ifg_inter)];
    R.ColumnsStartFrom = 'south';
    R.RowsStartFrom =  'west';
    R.CellExtentInLatitude = dlat_geo;
    R.CellExtentInLongitude = dlon_geo;
    
    I = uint8(round(((ifg_inter+pi)./(2.*pi)).*256));
    C = uint8(round(((coh_int)).*255));
    
    if exist('sar')==2
        if strcmp(lower(val_Filtphase),'yes')==1
            geotiffwrite([path_Exec,'/',name_ifg_filt,'.tif'],I,sar,R);
        else
            geotiffwrite([path_Exec,'/',name_ifg,'.tif'],I,sar,R);
        end
    else
        if strcmp(lower(val_Filtphase),'yes')==1
            geotiffwrite([path_Exec,'/',name_ifg_filt,'.tif'],I,jet,R);
        else
            geotiffwrite([path_Exec,'/',name_ifg,'.tif'],I,jet,R);
        end
    end
    
    if strcmp(lower(val_coh),'yes')==1
        geotiffwrite([path_Exec,'/',name_coh,'.tif'],C,R);
    end
    
    fprintf(1,'\t\tGEOTIFF: OK\n');
end
end

