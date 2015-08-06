function [ grnList,grnfns ] = GTdef_read_greensfns(fileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_read_greensfns                              %
% read GTdef greensfns after reordering vertices and stations                   %
%                                                                               %
% INPUT:                                                                        %
% GTdef green function file                                                     %
% (*1) location of gps sites                                                    %
%     1    2    3   4   5                                                       %
%     Site LOC  X   Y   Z                                                       %
% (*2) location of vertices                                                     %
%     1  2    3    4     5   6    7                                             %
%     id dnum snum LOC   X   Y    Z                                             %
% (3) greens functions of vertex-site pairs                                     %
%     1  2    3    4     5     6     7     8     9     10    11    12    13     %
%     id dnum snum Site  ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz  %
%                                                                               %
% OUTPUT:                                                                       %
% grnList - site names stored as a cell          (pntNum*1)                     %
% grnfns  - array of size                        (patchNum*pntNum*12)           %
%            for each patch-site pair                                           % 
%   = [ id dnum snum ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz ]    %
%                                                                               %
% first created by Lujia Feng Tue Aug  6 10:25:58 SGT 2013                      %
% modified based on PyLith_read_greensfns.m lfeng Tue Aug  6 10:26:36 SGT 2013  %
% scanning line is too slow, so use textscan lfeng Tue Aug  6 10:48:07 SGT 2013 %
% last modified by Lujia Feng Tue Aug  6 11:11:48 SGT 2013                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% check if this file exists %%%%%%%%%%
if ~exist(fileName,'file'), error('GTdef_read_greensfns ERROR: %s does not exist!',fileName); end

fin = fopen(fileName,'r');
grnList = {};
grnfns   = [];

%%%%%%%%%%%%%%%%%% start reading %%%%%%%%%%%%%%
dataCell = textscan(fin,'%f %f %f %*s %f %f %f %f %f %f %f %f %f','CommentStyle','#');
grnfns   = cell2mat(dataCell);
frewind(fin);
nameCell = textscan(fin,'%*f %*f %*f %s %*f %*f %*f %*f %*f %*f %*f %*f %*f','CommentStyle','#');
grnList  = nameCell{1};

fclose(fin);
