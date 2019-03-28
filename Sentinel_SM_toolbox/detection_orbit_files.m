%% Function to detecte the orbit file inside the orbit directory
%By Alexis HRYSIEWICZ Laboratoire Magmas et Volcans

function [precised,restitued] = detection_orbit_files(path_orbit,date,sat)
% Function to detect the good orbit file by acquisition. The function
% creates the list of the orbits files stocked in the same directory
% (precise and restitued) then detects the files using the dates (with the
% names of files).
%
% This function is ran by read_data_Sentinel_SM.m
%
% The order: 
%   - To search a possible precise orbit file 
%       if not
%           - To search a possible restitued orbit file 
%               if not 
%                   - To use the orbit in the .SAFE directory. 
%
% System requirements: 
%   - Linux / Mac (using of the bash language)
%
% Information about the developpement: 
%Author: Alexis Hrysiewicz (Laboratoire Magmas et Volcans / OPGC / OI2)
%Updates: 
%   - Created : October 2018
%   - Fixed the bugs with the name of satellite (November 2018)
%   - Fixed for the bugs with the matrix initialisation (March 2018) by Delphine Smittarello
%
% PLEASE READ the README for more information.

cur_no_orbit = cd;
fprintf(1,'####################\n'); 
fprintf(1,'Detection of the orbit files in \n\t%s\n',path_orbit);

%% Creating of the orbit lists; 
cd(path_orbit);
[a,b] = system('ls *.EOF');
list = strsplit(b);
cd(cur_no_orbit);

list_restitued_S1A = cell(1); h1 = 1; 
list_restitued_S1B = cell(1); h2 = 1; 
list_precise_S1A = cell(1); h3 = 1; 
list_precise_S1B = cell(1); h4 = 1; 
for i1 = 1 : length(list)
    if isempty(list{i1}) == 0
        if isempty(strfind(list{i1},'S1A'))==0
            if isempty(strfind(list{i1},'RESORB'))==0
                list_restitued_S1A{h1} = strtrim(list{i1}); h1 = h1 + 1;
            else
                list_precise_S1A{h3} = strtrim(list{i1}); h3 = h3 + 1;
            end
        elseif isempty(strfind(list{i1},'S1B'))==0
            if isempty(strfind(list{i1},'RESORB'))==0
                list_restitued_S1B{h2} = strtrim(list{i1}); h2 = h2 + 1;
            else
                list_precise_S1B{h4} = strtrim(list{i1}); h4 = h4 + 1;
            end
        end
    end
end
if strcmp(sat,'S1A')
    list_precise = list_precise_S1A;
    list_restitued = list_restitued_S1A;
else
    list_precise = list_precise_S1B;
    list_restitued = list_restitued_S1B;
end

%% Detection of the possible precised orbit
res_precise = [];
for i1 = 1 : length(list_precise)
    pos = find(list_precise{i1}=='_');
    if pos
    di1 = list_precise{i1}(pos(end-1)+2:pos(end)-1);
    di2 = list_precise{i1}(pos(end)+1:end-4);
    d1 = datetime([di1(1:4),'-',di1(5:6),'-',di1(7:8),' ',di1(10:11),':',di1(12:13),':',di1(14:15)],'InputFormat','yyyy-MM-dd HH:mm:SS');
    d2 = datetime([di2(1:4),'-',di2(5:6),'-',di2(7:8),' ',di2(10:11),':',di2(12:13),':',di2(14:15)],'InputFormat','yyyy-MM-dd HH:mm:SS');
    % Test of validity
    if date > d1 & date < d2
        res_precise(i1,1) = 1;
        res_precise(i1,2) = datenum(d1 - date);
        res_precise(i1,3) = datenum(date - d2);
        res_precise(i1,4) = datenum(d1);
        res_precise(i1,5) = datenum(d2);
    else
        res_precise(i1,1) = 0;
        res_precise(i1,2) = datenum(d1 - date);
        res_precise(i1,3) = datenum(date - d2);
        res_precise(i1,4) = datenum(d1);
        res_precise(i1,5) = datenum(d2);
    end
    end
end
if res_precise
if isempty(find(res_precise(:,1)==1))==0
    diff = abs(datenum(date) - (res_precise(:,4) + (res_precise(:,5) - res_precise(:,4))./2));
    diff(res_precise(:,1)==0)= NaN;
    precised = find(min(diff)==diff);
    precised = [path_orbit,'/',list_precise{precised(end)}];
    fprintf(1,'\t\tPRECISED ORBIT DETECTED in\n\t%s\n',precised);
    restitued = 'NONE';
end
else
    fprintf(1,'\t\tNO PRECISED ORBIT DETECTED\n');
    precised = 'NONE';
    %% Detection of the possible restitued orbit
    res_restitued = [];
    for i1 = 1 : length(list_restitued)
        pos = find(list_restitued{i1}=='_');
        di1 = list_restitued{i1}(pos(end-1)+2:pos(end)-1);
        di2 = list_restitued{i1}(pos(end)+1:end-4);
        d1 = datetime([di1(1:4),'-',di1(5:6),'-',di1(7:8),' ',di1(10:11),':',di1(12:13),':',di1(14:15)],'InputFormat','yyyy-MM-dd HH:mm:SS');
        d2 = datetime([di2(1:4),'-',di2(5:6),'-',di2(7:8),' ',di2(10:11),':',di2(12:13),':',di2(14:15)],'InputFormat','yyyy-MM-dd HH:mm:SS');
        % Test of validity
        if date > d1 & date < d2
            res_restitued(i1,1) = 1;
            res_restitued(i1,2) = datenum(d1 - date);
            res_restitued(i1,3) = datenum(date - d2);
            res_restitued(i1,4) = datenum(d1);
            res_restitued(i1,5) = datenum(d2);
        else
            res_restitued(i1,1) = 0;
            res_restitued(i1,2) = datenum(d1 - date);
            res_restitued(i1,3) = datenum(date - d2);
            res_restitued(i1,4) = datenum(d1);
            res_restitued(i1,5) = datenum(d2);
        end
    end
    if isempty(find(res_restitued(:,1)==1))==0
        diff = abs(datenum(date) - (res_restitued(:,4) + (res_restitued(:,5) - res_restitued(:,4))./2));
        diff(res_restitued(:,1)==0)= NaN;
        restitued = find(min(diff)==diff);
        restitued = [path_orbit,'/',list_restitued{restitued(end)}];
        fprintf(1,'\t\tRESTITUED ORBIT DETECTED in\n\t%s\n',restitued);
    else
        restitued = 'NONE';
        fprintf(1,'\t\tNO RESTITUED ORBIT DETECTED\n');
    end  
end
fprintf(1,'Detection of the orbit files: OVER\n');
fprintf(1,'####################\n'); 

