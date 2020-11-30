function [ siteList,siteloc,patchloc,grnList,grnfns ] = GTdef_read_greensfns(fileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_read_greensfns                              %
% read GTdef greensfns after reordering vertices and stations                   %
%                                                                               %
% INPUT:                                                                        %
% GTdef green function file                                                     %
% (1) location of gps sites                                                     %
%     1     2    3    4   5   6                                                 %
%     point num  SITE lon lat z                                                 %
% (2) location of patches                                                       %
%     1     2  3    4    5     6     7     8     9     10    11    12    13     %
%     patch id dnum snum xtop1 ytop1 ztop1 xbot1 ybot1 zbot1 xbot2 ybot2 zbot2  %
%     14    15    16    17   18   19                                            %
%     xtop2 ytop2 ztop2 xctr yctr zcrt                                          %
% (3) greens functions of patch-site pairs                                      %
%     1  2    3    4     5     6     7     8     9     10    11    12    13     %
%     id dnum snum Site  ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz  %
%                                                                               %
% OUTPUT:                                                                       %
% siteList - site names stored as a cell                                        %
% siteloc  - site location                                                      %
%          = [ lon lat height ]                  (siteNumx3)                    %
% patchloc - patch location                      (patchNum*18)                  %
%          = [ id dnum snum xtop1 ytop1 ztop1 xbot1 ybot1 zbot1                 %
%              xbot2 ybot2 zbot2 xtop2 ytop2 ztop2 xctr yctr zcrt ]             %
% grnList - site names stored as a cell          (pntNum*1)                     %
% grnfns  - array of size                        (patchNum*pntNum*12)           %
%            for each patch-site pair                                           % 
%   = [ id dnum snum ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz ]    %
%                                                                               %
% first created by Lujia Feng Tue Aug  6 10:25:58 SGT 2013                      %
% modified based on PyLith_read_greensfns.m lfeng Tue Aug  6 10:26:36 SGT 2013  %
% scanning line is too slow, so use textscan lfeng Tue Aug  6 10:48:07 SGT 2013 %
% added read in point & patch too lfeng Fri Nov 14 16:52:19 SGT 2014            %
% last modified by Lujia Feng Fri Nov 14 17:30:38 SGT 2014                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fileName,'file'), error('GTdef_read_greensfns ERROR: %s does not exist!',fileName); end

%%                        1  2  3  4   5  6  7  8  9  10 11 12 13
%dataCell = textscan(fin,'%f %f %f %*s %f %f %f %f %f %f %f %f %f','CommentStyle',{'#','p'});
%grnfns   = cell2mat(dataCell);
%frewind(fin);
%%                        1   2   3   4  5   6   7   8   9   10  11  12  13
%nameCell = textscan(fin,'%*f %*f %*f %s %*f %*f %*f %*f %*f %*f %*f %*f %*f','CommentStyle',{'#','p'});
%grnList  = nameCell{1};

fin = fopen(fileName,'r');
sitestd  = 'point';
patchstd = 'patch';
grnstd   = '^\d+\s+\d+\s+\d+\s+\w{4}\s+';
siteList = {};
siteloc  = [];
patchloc = [];
grnList  = {};
grnfns   = [];

while(1)   
    % read in one line
    tline = fgetl(fin);
    % test if it is the end of file; exit if yes
    if ischar(tline)~=1, break; end
    % exclude comment lines
    if strncmp(tline,'#',1), continue; end

    % read site location
    isite = regexp(tline,sitestd,'match');
    if ~isempty(isite)
       dataCell = regexp(tline,'\s+','split');
       dataStr  = char(dataCell);
       siteList = [ siteList; dataStr(3,1:4) ];
       data     = str2num(dataStr(4:6,:));
       siteloc  = [ siteloc; data' ];
    end

    % read patch location
    ispatch = regexp(tline,patchstd,'match');
    if ~isempty(ispatch)
       dataCell = regexp(tline,'\s+','split');
       dataStr  = char(dataCell);
       data     = str2num(dataStr(2:end,:));
       patchloc = [ patchloc; data' ];
    end

    % read greens function pairs
    isgrn = regexp(tline,grnstd,'match');
    if ~isempty(isgrn)
       dataCell = regexp(tline,'\s+','split');
       dataStr  = char(dataCell);
       grnList  = [ grnList; dataStr(4,1:4) ];
       data     = str2num(dataStr([1:3 5:end],:));
       grnfns   = [ grnfns; data' ];
    end
end

fclose(fin);
