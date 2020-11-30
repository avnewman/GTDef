function [] = Abaqus2GTdef(fin_name,lon0,lat0,rot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              		         Abaqus2GTdef.m  			       %
%									       %
% Convert Abaqus .rpt format to GTdef format			  	       %
% Currently only work for nodes						       %
%									       %
% INPUT									       %
%   fin_name - *.rpt input file name                                           %
%   lon0,lat0 - geographyic coordinates of origion in Abaqus cartesian system  %
%   rot - the rotation angle from Cartesian to Geographic coordinate	       %
%   (E=x+; N=y+; CCW=+; CW=-) [degree]					       %
% Note: 								       %
% (1) set lon0,lat0 or rot as nan if no conversion needed		       %	
% (2) Usually, km is used in Abaqus, so need to convert it to m in GTdef       %
%									       %
% e.g. lon0=21.490;lat0=37.925;rot=30 for 2008 Patras EQ		       %
%      									       %
% first created by lfeng Tue May  4 15:42:34 EDT 2010			       %
% use cell array of strings for names lfeng Wed Dec  1 17:58:23 EST 2010       %
% last modified by lfeng Wed Dec  1 17:58:33 EST 2010			       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read in Abaqus Field Output Report %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fin = fopen(fin_name,'r');
nod_num = 0; nod_data = [];
while(1)   
    % read in one line
    tline = fgetl(fin);
    % test if it is the end of file; exit if yes
    if ischar(tline)~=1, break; end
    % read the 1st term that is a flag for data type; could be some white-spaces before the 1st term
    [flag,remain] = strtok(tline);		% the default delimiter is white-space
    % omit '#' comment lines
    if strncmp(flag,'#',1), continue; end

    %%%%% node %%%%%
    if strcmp(flag,'Node')
	% read in the variable description
        var_num = 1; var_name = {'Node.Label'};
	while(1)
            [var,remain] = strtok(remain);
    	    if isempty(var)||strncmp(var,'#',1), break; end
            var_num = var_num+1; 
	    var_name = [ var_name; var ];
	end
        % skip two lines after node
        tline = fgetl(fin); tline = fgetl(fin);
	% read in the node data
	while(1)
    	    % read in one line
    	    tline = fgetl(fin);
    	    % test if this line still includes the data
            [str1,remain] = strtok(tline); var = str2double(str1);
	    if isnan(var), break; end
	    nod_num = nod_num+1;
	    nod_data(nod_num,1) = var;
            for ii = 2:var_num
 	        [ nod_data(nod_num,ii),remain ] = GTdef_read1double(remain);
            end
	end
	continue
    end
end
fclose(fin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Some conversions if necessary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create an index list for (1)Node.Lable (2)COORD.COOR1 (3)COORD.COOR2 (4)U.U1 (5)U.U2 (6)U.U3
ind = [];
for ii = 1:var_num
    if strncmp(var_name{ii},'Node.Label',10),  ind(1) = ii; end
    if strncmp(var_name{ii},'COORD.COOR1',11), ind(2) = ii; end
    if strncmp(var_name{ii},'COORD.COOR2',11), ind(3) = ii; end
    if strncmp(var_name{ii},'U.U1',4), ind(4) = ii; end
    if strncmp(var_name{ii},'U.U2',4), ind(5) = ii; end
    if strncmp(var_name{ii},'U.U3',4), ind(6) = ii; end
end
% convert from Abaqus [km] to GTdef [m]
nod_data(:,ind(2:6)) = nod_data(:,ind(2:6))*1000;
% convert from Abaqus Cartesian to Geographic coordinate
if ~isnan(lon0)&&~isnan(lat0)&&~isnan(rot)
    [nod_data(:,ind(2)),nod_data(:,ind(3))] = ckm2LLd(nod_data(:,ind(2)),nod_data(:,ind(3)),lon0,lat0,rot);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Write out GTdef file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basename = strtok(fin_name,'.');
fout_name = strcat(basename,'.out');
fout = fopen(fout_name,'w');
%%%%%%%%%%  point  %%%%%%%%%%
for ii =1:nod_num
    fprintf(fout,'point 3 %d     %-14.8f %-12.8f nan  %10.5f %-8.5f %-8.5f  nan	 nan  nan  1.0\n', nod_data(ii,ind));
end
fclose(fout);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_read1double.m				%
% 	      function to read in one double number from a string		%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ value,remain ] = GTdef_read1double(str)

[str1,remain] = strtok(str);
if isempty(str1)||strncmp(str1,'#',1)  
    error('The input file is wrong!');
end
value = str2double(str1);
