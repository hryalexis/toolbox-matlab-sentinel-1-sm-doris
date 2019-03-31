function read_data_Sentinel_SM(slc_path,header_path,path_Exec,Mode,val_extraction,first_line,last_line,first_pixel,last_pixel,orb_path)
% Function to read the Sentinel StripMap data and to create the files (.slc and .res) for
% Doris. The extraction of ROI and orbit reading are done here. The
% function is the same for the master and the slave (input "Mode"). The
% user must not manually run this function.
%
% This function is ran by computation_Sentinel_SM
%
% Required function / programs: 
%   - detection_orbits_files.m
%   - getTIFFinfo.m by Louis-Philippe Rousseau
%   - readSent1Data.m by Louis-Philippe Rousseau
%   - xml2struct.m 
%
% System requirements: 
%   - Linux / Mac (using of the bash language)
%
% Information about the developpement: 
%Author: Alexis Hrysiewicz (Laboratoire Magmas et Volcans / OPGC / OI2)
%Updates: 
%   - Created : (??)
%   - Implementation for the orbit files (October 2018)
%
% PLEASE READ the README for more information

file_SLC_SM = xml2struct(header_path);
% Initialisation
if strcmp(Mode,'Master') == 1
    f_res = fopen([path_Exec,'/master.res'],'w');
    fprintf(f_res,'===============================================\n');
    fprintf(f_res,'MASTER RESULTFILE:                      master.res\n');
else
    f_res = fopen([path_Exec,'/slave.res'],'w');
    fprintf(f_res,'===============================================\n');
    fprintf(f_res,'SLAVE RESULTFILE:                      slave.res\n');
end
fprintf(f_res,'Created by:                             Alexis Hrysiewicz\n');
fprintf(f_res,'DVersion:                               Version (2015)\n');
fprintf(f_res,'FFTW library:                           used\n');
fprintf(f_res,'VECLIB library:                         not used\n');
fprintf(f_res,'LAPACK library:                         used\n');
fprintf(f_res,'Compiled at:                            XXXXXXXX\n');
fprintf(f_res,'By GUN gcc:                             XXXXXXXX\n');
fprintf(f_res,'===============================================\n');
fprintf(f_res,'File creation at:       %s\n\n',datestr(now));
fprintf(f_res,' -------------------------------------------------------\n');
fprintf(f_res,'| Delft Institute of Earth Observation & Space Systems  |\n');
fprintf(f_res,'|          Delft University of Technology               |\n');
fprintf(f_res,'|              http://doris.tudelft.nl                  |\n');
fprintf(f_res,'|                                                       |\n');
fprintf(f_res,'| Author: (c) TUDelft - DEOS Radar Group                |\n');
fprintf(f_res,' -------------------------------------------------------\n\n\n');
fprintf(f_res,'Start_process_control\n');
fprintf(f_res,'readfiles:\t\t\t1\n');
fprintf(f_res,'precise_orbits:\t\t1\n');
fprintf(f_res,'crop:\t\t\t\t1\n');
fprintf(f_res,'sim_amplitude:\t\t0\n');
fprintf(f_res,'master_timing:\t\t0\n');
fprintf(f_res,'oversample:\t\t\t0\n');
fprintf(f_res,'resample:\t\t\t0\n');
fprintf(f_res,'filt_azi:\t\t\t0\n');
fprintf(f_res,'filt_range:\t\t\t0\n');
fprintf(f_res,'NOT_USED:\t\t\t0\n');
fprintf(f_res,'End_process_control\n\n\n');

% Read of xml file
fprintf(f_res,'*******************************************************************\n');
fprintf(f_res,'*_Start_readfiles:\n');
fprintf(f_res,'*******************************************************************\n');
fprintf(f_res,'Volume_file: \t\t\t%s\n','dummy');
fprintf(f_res,'Volume_ID: \t\t\t%s\n',file_SLC_SM.product.adsHeader.missionDataTakeId.Text);
fprintf(f_res,'Volume_identifier: \t\t\t%s\n','dummy');
fprintf(f_res,'Volume_set_identifier: \t\t\t%s\n','dummy');
fprintf(f_res,'Number of records in ref. file: \t\t\t%s\n','dummy');
fprintf(f_res,'SAR_PROCESSOR: \t\t\t%s\n','dummy');
fprintf(f_res,'SWATH: \t\t\t%s\n',file_SLC_SM.product.adsHeader.swath.Text);
fprintf(f_res,'PASS: \t\t\t%s\n',file_SLC_SM.product.generalAnnotation.productInformation.pass.Text);
fprintf(f_res,'IMAGE_MODE: \t\t\t%s\n',file_SLC_SM.product.adsHeader.mode.Text);
fprintf(f_res,'polarisation: \t\t\t%s\n',file_SLC_SM.product.adsHeader.polarisation.Text);
fprintf(f_res,'Product type specifier: \t\t\t%s\n',file_SLC_SM.product.adsHeader.missionId.Text);
fprintf(f_res,'Logical volume generating facility: \t\t\t%s\n','dummy');
fprintf(f_res,'Location and date/time of product creation: \t\t\t%s\n','dummy');
fprintf(f_res,'RADAR_FREQUENCY (HZ): \t\t\t%s\n',file_SLC_SM.product.generalAnnotation.productInformation.radarFrequency.Text);
fprintf(f_res,'Scene identification: \t\t\t%s %s\n','Orbit',file_SLC_SM.product.adsHeader.absoluteOrbitNumber.Text);

%Extraction of the coordinate grid
azimuthTime=cell(1);
slantRangeTime = cell(1);
line = [];
pixel = [];
latitude = [];
longitude = [];
height = [];
incidenceAngle = [];
elevationAngle = [];
for i1 = 1 : length(file_SLC_SM.product.geolocationGrid.geolocationGridPointList.geolocationGridPoint)
    azimuthTime{i1,1} = file_SLC_SM.product.geolocationGrid.geolocationGridPointList.geolocationGridPoint{i1}.azimuthTime.Text;
    slantRangeTime{i1,1} = file_SLC_SM.product.geolocationGrid.geolocationGridPointList.geolocationGridPoint{i1}.slantRangeTime.Text;
    line(i1,1) = str2num(file_SLC_SM.product.geolocationGrid.geolocationGridPointList.geolocationGridPoint{i1}.line.Text);
    pixel(i1,1) = str2num(file_SLC_SM.product.geolocationGrid.geolocationGridPointList.geolocationGridPoint{i1}.pixel.Text);
    latitude(i1,1) = str2num(file_SLC_SM.product.geolocationGrid.geolocationGridPointList.geolocationGridPoint{i1}.latitude.Text);
    longitude(i1,1) = str2num(file_SLC_SM.product.geolocationGrid.geolocationGridPointList.geolocationGridPoint{i1}.longitude.Text);
    height(i1,1) = str2num(file_SLC_SM.product.geolocationGrid.geolocationGridPointList.geolocationGridPoint{i1}.height.Text);
    incidenceAngle(i1,1) = str2num(file_SLC_SM.product.geolocationGrid.geolocationGridPointList.geolocationGridPoint{i1}.incidenceAngle.Text);
    elevationAngle(i1,1) = str2num(file_SLC_SM.product.geolocationGrid.geolocationGridPointList.geolocationGridPoint{i1}.elevationAngle.Text);
end

nbl = str2num(file_SLC_SM.product.imageAnnotation.imageInformation.numberOfLines.Text);
nbc = str2num(file_SLC_SM.product.imageAnnotation.imageInformation.numberOfSamples.Text);

%Write of the parameters
fprintf(f_res,'Scene location: \t\t\t lat: %f lon: %f\n',latitude(1),longitude(1));
fprintf(f_res,'Sensor platform mission identifer: \t\t\t%s\n',file_SLC_SM.product.adsHeader.missionId.Text);
fprintf(f_res,'Scene_center_heading: \t\t\t%s\n',file_SLC_SM.product.generalAnnotation.productInformation.platformHeading.Text);
fprintf(f_res,'Scene_centre_latitude: \t\t\t%10.9f\n',mean(latitude));
fprintf(f_res,'Scene_centre_longitude: \t\t\t%10.9f\n',mean(longitude));
fprintf(f_res,'Radar_wavelength (m): \t\t\t%10.9f\n',299792458.0/str2num(file_SLC_SM.product.generalAnnotation.productInformation.radarFrequency.Text));
a = azimuthTime{1}; a(a=='T')=' ';
posa = find(a=='.');
b = a(posa+1:end);
fprintf(f_res,'First_pixel_azimuth_time (UTC): \t\t\t%s\n',[datestr(a,'dd-mmm-yyyy HH:MM:ss'),'.',b]);

fprintf(f_res,'Pulse_Repetition_Frequency (computed, Hz): \t\t\t%f\n',str2num(file_SLC_SM.product.imageAnnotation.imageInformation.azimuthFrequency.Text));
fprintf(f_res,'Total_azimuth_band_width (Hz): \t\t\t%f\n',str2num(file_SLC_SM.product.imageAnnotation.processingInformation.swathProcParamsList.swathProcParams.azimuthProcessing.totalBandwidth.Text));
fprintf(f_res,'Weighting_azimuth: \t\t\t%s\n',file_SLC_SM.product.imageAnnotation.processingInformation.swathProcParamsList.swathProcParams.azimuthProcessing.windowType.Text);
fprintf(f_res,'Range_time_to_first_pixel (2way) (ms): \t\t\t%10.20f\n',str2num(file_SLC_SM.product.imageAnnotation.imageInformation.slantRangeTime.Text).*1000);
fprintf(f_res,'Range_sampling_rate (computed, MHz): \t\t\t%10.9f\n',str2num(file_SLC_SM.product.generalAnnotation.productInformation.rangeSamplingRate.Text)./1000000);
fprintf(f_res,'Total_range_band_width (MHz): \t\t\t%10.9f\n',str2num(file_SLC_SM.product.imageAnnotation.processingInformation.swathProcParamsList.swathProcParams.rangeProcessing.processingBandwidth.Text)./1000000);
fprintf(f_res,'Weighting_range: \t\t\t%s\n',file_SLC_SM.product.imageAnnotation.processingInformation.swathProcParamsList.swathProcParams.rangeProcessing.windowType.Text);

%For the Doppler centroid
a = file_SLC_SM.product.dopplerCentroid.dcEstimateList.dcEstimate{1,2}.dataDcPolynomial.Text;
posa = find(a==' ');
fprintf(f_res,'Xtrack_f_DC_constant (Hz, early edge): \t\t\t%s\n',a(1:posa(1)-1));
fprintf(f_res,'Xtrack_f_DC_linear (Hz/s, early edge): \t\t\t%s\n',a(posa(1)+1:posa(2)-1));
fprintf(f_res,'Xtrack_f_DC_quadratic (Hz/s/s, early edge): \t\t\t%s\n',a(posa(2)+1:end));

%Finalisation
fprintf(f_res,'\n*******************************************************************\n');
fprintf(f_res,'Datafile: \t\t\t%s\n','NaN');
fprintf(f_res,'Dataformat: \t\t\t%s\n','tiff');
fprintf(f_res,'Number_of_lines_original: \t\t\t%d\n',nbl);
fprintf(f_res,'Number_of_pixels_original: \t\t\t%d\n',nbc);
fprintf(f_res,'*******************************************************************\n');
fprintf(f_res,'* End_readfiles:_NORMAL\n');
fprintf(f_res,'*******************************************************************\n');

% Orbites
if strcmp(orb_path,'NONE') == 0
    %% Reading of the orbit directory
    a = azimuthTime{1}; a(a=='T')=' ';
    posa = find(a=='.');
    b = a(posa+1:end);
    if isempty(strfind(slc_path,'S1A'))==0
        [precised,restitued] = detection_orbit_files(orb_path,datetime([datestr(a,'dd-mmm-yyyy HH:MM:ss')],'InputFormat','dd-MMM-yyyy HH:mm:SS'),'S1A');
    else
        [precised,restitued] = detection_orbit_files(orb_path,datetime([datestr(a,'dd-mmm-yyyy HH:MM:ss')],'InputFormat','dd-MMM-yyyy HH:mm:SS'),'S1B');
    end
    
    if strcmp(precised,'NONE')==0
        mode_orbit = 'precised';
        path_orbit_file = precised;
    elseif strcmp(restitued,'NONE')==0
        mode_orbit = 'restitued';
        path_orbit_file = restitued;
    else
        mode_orbit = 'RAW';
        path_orbit_file = 'NONE';
    end
    if strcmp(path_orbit_file,'NONE')==0
        
        [a,abs_orb] = system(['grep ''<Absolute_Orbit>'' ',[path_orbit_file],' | awk '' {print $1}''']); abs_orb = strsplit(abs_orb);
        [a,UTC] = system(['grep ''<UTC>UTC='' ',[path_orbit_file],' | awk '' {print $1}''']); UTC = strsplit(UTC);
        [a,X] = system(['grep ''<X unit="m">'' ',[path_orbit_file],' | awk '' {print $2}''']); X = strsplit(X);
        [a,Y] = system(['grep ''<Y unit="m">'' ',[path_orbit_file],' | awk '' {print $2}''']); Y = strsplit(Y);
        [a,Z] = system(['grep ''<Z unit="m">'' ',[path_orbit_file],' | awk '' {print $2}''']); Z = strsplit(Z);
        [a,Vx] = system(['grep ''<VX unit="m/s">'' ',[path_orbit_file],' | awk '' {print $2}''']); Vx = strsplit(Vx);
        [a,Vy] = system(['grep ''<VY unit="m/s">'' ',[path_orbit_file],' | awk '' {print $2}''']); Vy = strsplit(Vy);
        [a,Vz] = system(['grep ''<VZ unit="m/s">'' ',[path_orbit_file],' | awk '' {print $2}''']); Vz = strsplit(Vz);
        
        %Creating the cell matrix
        orb_mat = cell(1);
        for oi1 = 1 : length(abs_orb)
            if isempty(abs_orb{oi1})==0
                posi1 = find(abs_orb{oi1}=='<');
                posi2 = find(abs_orb{oi1}=='>');
                orb_mat{oi1,1} = str2num(abs_orb{oi1}(posi2(1)+1:posi1(end)-1));
                
                posi1 = find(UTC{oi1}=='=');
                posi2 = find(UTC{oi1}=='T'); posi2 = posi2(3);
                yi = UTC{oi1}(posi1+1:posi1+4);
                moi = UTC{oi1}(posi1+6:posi1+7);
                dai = UTC{oi1}(posi1+9:posi1+10);
                hi = UTC{oi1}(posi2+1:posi2+2);
                mi = UTC{oi1}(posi2+4:posi2+5);
                si = UTC{oi1}(posi2+7:posi2+8);
                dsi = UTC{oi1}(posi2+10:posi2+15);
                date_orb(oi1) = datetime([yi,'-',moi,'-',dai,' ',hi,':',mi,':',si]);
                date_int(oi1) = str2num(hi).*3600 + str2num(mi).*60 + str2num([si,'.',dsi]);
                
                posi1 = find(X{oi1}=='<');
                posi2 = find(X{oi1}=='>');
                X_mat(oi1) = str2num(X{oi1}(posi2(1)+1:posi1(1)-1));
                
                posi1 = find(Y{oi1}=='<');
                posi2 = find(Y{oi1}=='>');
                Y_mat(oi1) = str2num(Y{oi1}(posi2(1)+1:posi1(1)-1));
                
                posi1 = find(Z{oi1}=='<');
                posi2 = find(Z{oi1}=='>');
                Z_mat(oi1) = str2num(Z{oi1}(posi2(1)+1:posi1(1)-1));
                
                posi1 = find(Vx{oi1}=='<');
                posi2 = find(Vx{oi1}=='>');
                Vx_mat(oi1) = str2num(Vx{oi1}(posi2(1)+1:posi1(1)-1));
                
                posi1 = find(Vy{oi1}=='<');
                posi2 = find(Vy{oi1}=='>');
                Vy_mat(oi1) = str2num(Vy{oi1}(posi2(1)+1:posi1(1)-1));
                
                posi1 = find(Vz{oi1}=='<');
                posi2 = find(Vz{oi1}=='>');
                Vz_mat(oi1) = str2num(Vz{oi1}(posi2(1)+1:posi1(1)-1));
                
            end
            if mod(oi1,100)==0; fprintf(1,'\t\t\tReading of orbit files: %d orbits on %d orbits\n',oi1,length(abs_orb)); end;
        end
        
        %Detection of the intersected orbits
        a = azimuthTime{1}; a(a=='T')=' ';
        posa = find(a=='.');
        b = a(posa+1:end);
        a = datetime([datestr(a,'dd-mmm-yyyy HH:MM:ss')],'InputFormat','dd-MMM-yyyy HH:mm:SS'); 
        pos_date = find(min(abs(datenum(a)-datenum(date_orb)))==abs(datenum(a)-datenum(date_orb)));
        if pos_date + 25 > length(date_orb)
            select = [pos_date-25 : length(date_orb)];
        elseif pos_date -25 < 1
            select = [1 : pos_date+25];
        else
            select = [pos_date-25 : pos_date+25];
        end
        
        %Writting of the .res file
        fprintf(f_res,'\n\n*******************************************************************\n');
        fprintf(f_res,'*_Start_precise_orbits:\n');
        fprintf(f_res,'*******************************************************************\n');
        fprintf(f_res,' t(s)            X(m)            Y(m)            Z(m)            X_V(m/s)            Y_V(m/s)            Z_V(m/s)\n');
        fprintf(f_res,'NUMBER_OF_DATAPOINTS:                    %d\n\n',length(select));
        for i3 = 1  : length(select)
            date_inti = date_int(select(i3));
            orb_X_int = X_mat(select(i3));
            orb_Y_int = Y_mat(select(i3));
            orb_Z_int = Z_mat(select(i3));
            
            orb_VX_int = Vx_mat(select(i3));
            orb_VY_int = Vy_mat(select(i3));
            orb_VZ_int = Vz_mat(select(i3));
            
            fprintf(f_res,'%f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n',[date_inti orb_X_int orb_Y_int orb_Z_int orb_VX_int orb_VY_int orb_VZ_int]');
        end
        fprintf(f_res,'\n*******************************************************************\n');
        fprintf(f_res,'* End_precise_orbits:_NORMAL\n');
        fprintf(f_res,'*******************************************************************\n');
        
    end
end
if strcmp(orb_path,'NONE') == 1 | strcmp(path_orbit_file,'NONE')==1
    mode_orbit = 'RAW';
    
    fprintf(f_res,'\n\n*******************************************************************\n');
    fprintf(f_res,'*_Start_precise_orbits:\n');
    fprintf(f_res,'*******************************************************************\n');
    fprintf(f_res,' t(s)            X(m)            Y(m)            Z(m)            X_V(m/s)            Y_V(m/s)            Z_V(m/s)\n');
    fprintf(f_res,'NUMBER_OF_DATAPOINTS:                    %d\n\n',length(file_SLC_SM.product.generalAnnotation.orbitList.orbit));
    for i3 = 1  : length(file_SLC_SM.product.generalAnnotation.orbitList.orbit)
        a = file_SLC_SM.product.generalAnnotation.orbitList.orbit{i3}.time.Text; pos = find(a=='T'); a = a(pos+1:end);
        date_int = str2num(a(1:2)).*3600 + str2num(a(4:5)).*60 + str2num(a(7:end));
        orb_X_int = str2num(file_SLC_SM.product.generalAnnotation.orbitList.orbit{i3}.position.x.Text);
        orb_Y_int = str2num(file_SLC_SM.product.generalAnnotation.orbitList.orbit{i3}.position.y.Text);
        orb_Z_int = str2num(file_SLC_SM.product.generalAnnotation.orbitList.orbit{i3}.position.z.Text);
        
        orb_VX_int = str2num(file_SLC_SM.product.generalAnnotation.orbitList.orbit{i3}.velocity.x.Text);
        orb_VY_int = str2num(file_SLC_SM.product.generalAnnotation.orbitList.orbit{i3}.velocity.y.Text);
        orb_VZ_int = str2num(file_SLC_SM.product.generalAnnotation.orbitList.orbit{i3}.velocity.z.Text);
        
        fprintf(f_res,'%f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n',[date_int orb_X_int orb_Y_int orb_Z_int orb_VX_int orb_VY_int orb_VZ_int]');
    end
    fprintf(f_res,'\n*******************************************************************\n');
    fprintf(f_res,'* End_precise_orbits:_NORMAL\n');
    fprintf(f_res,'*******************************************************************\n');
end

% Copy of the tiff SLC
if strcmp(lower(val_extraction),'yes')==1
    roi.firstSample = str2num(first_pixel);
    roi.numSample = str2num(last_pixel)-str2num(first_pixel)+1;
    roi.firstLine = str2num(first_line);
    roi.numLine = str2num(last_line)-str2num(first_line)+1;
    data = readSent1Data(slc_path,roi,true);
    real_part = real(data)';
    imag_part = imag(data)';
else
    data = readSent1Data(slc_path,true);
    real_part = real(data)';
    imag_part = imag(data)';
end

fprintf(f_res,'\n\n\n*******************************************************************\n');
if strcmp(Mode,'Master') == 1
    fprintf(f_res,'*_Start_crop:			master step01\n');
    fprintf(f_res,'*******************************************************************\n');
    fprintf(f_res,'Data_output_file: 				%s\n',[path_Exec,'/image_master.slc']);
else
    fprintf(f_res,'*_Start_crop:			slave step01\n');
    fprintf(f_res,'*******************************************************************\n');
    fprintf(f_res,'Data_output_file: 				%s\n',[path_Exec,'/image_slave.slc']);
end
fprintf(f_res,'Data_output_format: 				complex_short\n');
if strcmp(lower(val_extraction),'yes')==1
    fprintf(f_res,'First_line (w.r.t. original_image): 		%s\n',first_line);
    fprintf(f_res,'Last_line (w.r.t. original_image): 		%s\n',last_line);
    fprintf(f_res,'First_pixel (w.r.t. original_image): 	%s\n',first_pixel);
    fprintf(f_res,'Last_pixel (w.r.t. original_image): 		%s\n',last_pixel);
    fprintf(f_res,'Number of lines (non-multilooked): 		%d\n',str2num(last_line)-str2num(first_line)+1);
    fprintf(f_res,'Number of pixels (non-multilooked): 		%d\n',str2num(last_pixel)-str2num(first_pixel)+1);
else
    fprintf(f_res,'First_line (w.r.t. original_image): 		%s\n',1);
    fprintf(f_res,'Last_line (w.r.t. original_image): 		%s\n',size(real_part,1));
    fprintf(f_res,'First_pixel (w.r.t. original_image): 	%s\n',1);
    fprintf(f_res,'Last_pixel (w.r.t. original_image): 		%s\n',size(real_part,2));
    fprintf(f_res,'Number of lines (non-multilooked): 		%d\n',size(real_part,1));
    fprintf(f_res,'Number of pixels (non-multilooked): 		%d\n',size(real_part,2));
end
fprintf(f_res,'*******************************************************************\n');
fprintf(f_res,'* End_crop:_NORMAL\n');
fprintf(f_res,'*******************************************************************\n');
mat = zeros(size(real_part,1),size(real_part,2).*2);
mat(:,1:2:size(real_part,2).*2)=real_part;
mat(:,2:2:size(real_part,2).*2)=imag_part;
clearvars real_part;
clearvars imag_part;
if strcmp(Mode,'Master') == 1
    fid = fopen([path_Exec,'/image_master.slc'],'w');
else
    fid = fopen([path_Exec,'/image_slave.slc'],'w');
end
fwrite(fid,mat','int16');
fclose(fid);
fclose(f_res);
clearvars mat;


end
