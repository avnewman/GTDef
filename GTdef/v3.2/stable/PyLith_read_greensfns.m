function [ siteList,siteloc,vertices,grnList,grnfns ] = PyLith_read_greensfns(fileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             PyLith_read_greensfns                             %
% read PyLith greensfns after reordering vertices and stations                  %
%                                                                               %
% INPUT:                                                                        %
% a green function database file listing                                        %
% (1) location of gps sites                                                     %
%     1    2    3   4   5                                                       %
%     Site LOC  X   Y   Z                                                       %
% (2) location of vertices                                                      %
%     1  2    3    4     5   6    7                                             %
%     id dnum snum LOC   X   Y    Z                                             %
% (3) greens functions of vertex-site pairs                                     %
%     1  2    3    4     5     6     7     8     9     10    11    12    13     %
%     id dnum snum Site  ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz  %
%                                                                               %
% OUTPUT:                                                                       %
% siteList - site names stored as a cell          (pntNum*1)                    %
% siteloc  - [ xx yy zz ]                         (pntNum*3)                    %
% vertices - [ id dnum snum xx yy zz ]            (vertNum*6)                   %
% grnList  - site names corresponding to grnfns   (vertexNum*pntNum*1)          %
% grnfns   - array of size                        (vertexNum*pntNum*9)          %
%            for each patch-site pair                                           % 
%          = [ ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz ]          %
%                                                                               %
% first created by lfeng Fri Nov 30 15:56:07 SGT 2012                           %
% last modified by lfeng Fri Nov 30 17:59:15 SGT 2012                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% check if this file exists %%%%%%%%%%
if ~exist(fileName,'file'), error('PyLith_read_greensfns ERROR: %s does not exist!',fileName); end

fin = fopen(fileName,'r');
sitestd  = '^\w{4}\s+LOC'; 			% \w = [a-zA-Z_0-9]
vertstd  = '^\d+\s+\d+\s+\d+\s+LOC';            % \d = [0-9]
grnstd   = '^\d+\s+\d+\s+\d+\s+\w{4}\s+';
siteList = {};
siteloc  = [];
vertices = [];
grnList  = {};
grnfns   = [];

%%%%%%%%%%%%%%%%%% start reading %%%%%%%%%%%%%%
while(1)   
    % read in one line
    tline = fgetl(fin);
    % test if it is the end of file; exit if yes
    if ischar(tline)~=1, break; end
    % read site location
    isite = regexp(tline,sitestd,'match');
    if ~isempty(isite)
       dataCell = regexp(tline,'\s+','split');
       dataStr  = char(dataCell);
       siteList = [ siteList; dataStr(1,1:4) ];
       data     = str2num(dataStr(3:5,:));
       siteloc  = [ siteloc; data' ];
    end

    % read vertices location
    isvert = regexp(tline,vertstd,'match');
    if ~isempty(isvert)
       dataCell = regexp(tline,'\s+','split');
       dataStr  = char(dataCell);
       data     = str2num(dataStr([1:3 5:7],:));
       vertices = [ vertices; data' ];
    end

    % read greens function pairs
    isgrn = regexp(tline,grnstd,'match');
    if ~isempty(isgrn)
       dataCell = regexp(tline,'\s+','split');
       dataStr  = char(dataCell);
       grnList  = [ grnList; dataStr(4,1:4) ];
       data     = str2num(dataStr(5:end,:));
       grnfns   = [ grnfns; data' ];
    end
end

fclose(fin);

% reshape greens functions  
vertNum = size(vertices,1);
pntNum  = size(siteloc,1);
grnfns  = reshape(grnfns',9,pntNum,vertNum);
grnfns  = permute(grnfns,[ 3 2 1 ]);
