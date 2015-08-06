function [ limits ] = GTdef_detect_greensfns(fgrnName,threshold,area)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_detect_greensfns                             %
% read in *.grnfns file and calculate the minimum slip needed for each patch    %
% to be detected by at least one station                                        %
%                                                                               %
% INPUT:                                                                        %
% fgrnName  - GTdef greens function file                                        %
% (1) location of gps sites                                                     %
%     1    2    3   4   5                                                       %
%     Site LOC  lon lat z                                                       %
% (2) location of vertices                                                      %
%     1  2    3    4    5   6    7                                              %
%     id dnum snum LOC  lon lat  z                                              %
% (3) greens functions of patch-site pairs                                      %
%     1  2    3    4     5     6     7     8     9     10    11    12    13     %
%     id dnum snum Site  ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz  %
% threshold - minimul level of detection for each component [m]                 %
% area      - area of each patch [m^2]                                          %
%    if area = [], M0 & Mw are not calculated                                   %
%                                                                               %
% grnList - site names stored as a cell          (pntNum*1)                     %
% grnfns  - array of size                        (patchNum*pntNum*12)           %
%            for each patch-site pair                                           % 
%       1  2    3    4     5     6     7     8     9     10    11    12         %
%   = [ id dnum snum ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz ]    %
%                                                                               %
% OUTPUT:                                                                       %
% *.min records minimum slip & magnitdue needed for detection of threshold      %
%            1  2    3    4     5     6     7       8       9                   %
% limits = [ id dnum snum ssmin dsmin tsmin smagmin dmagmin tmagmin ]           %
%                                                                               %
% first created by Lujia Feng Fri Dec  6 09:05:58 SGT 2013                      %
% last modified by Lujia Feng Fri Dec  6 10:39:01 SGT 2013                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read green's functions
fprintf(1,'\n.......... reading %s ...........\t',fgrnName);
tic
[ ~,grnfns ] = GTdef_read_greensfns(fgrnName);
toc

% read in parameters
patchNum  = grnfns(end,1);
patchList = grnfns(:,1); 

% loop through all patches
fprintf(1,'\n.......... processing %d patches ...........\t',patchNum);
tic
ssM0 = nan; ddM0 = nan; ttM0 = nan;
ssMw = nan; ddMw = nan; ttMw = nan;
limits = zeros(patchNum,9);
for ii=1:patchNum
    ind = patchList==ii;
    grnpatch = grnfns(ind,:);
    dd       = grnpatch(1,2);
    ss       = grnpatch(1,3);
    disps    = grnpatch(:,4:end);
    slips    = abs(threshold./disps); % absolute values!!
    slipmin  = min(slips); % detected by at least one site
    % minimum slip
    ssmin    = min(slipmin(1:3));
    ddmin    = min(slipmin(4:6));
    ttmin    = min(slipmin(7:9));
    if ~isempty(area)
        % minimum M0
        [ ssM0,~,~ ] = EQS_M0AreaSlip([],[],area,ssmin);
        [ ddM0,~,~ ] = EQS_M0AreaSlip([],[],area,ddmin);
        [ ttM0,~,~ ] = EQS_M0AreaSlip([],[],area,ttmin);
        % minimum Mw
        [ ssMw,~ ]   = EQS_MwM([],ssM0);
        [ ddMw,~ ]   = EQS_MwM([],ddM0);
        [ ttMw,~ ]   = EQS_MwM([],ttM0);
    end
    limits(ii,:) = [ ii dd ss ssmin ddmin ttmin ssMw ddMw ttMw ];
end
toc

% output limits file
[ ~,basename,~ ] = fileparts(fgrnName);
foutName = [ basename '_disp' num2str(threshold*1e3,'%04.0f') 'mm.min' ];
fout = fopen(foutName,'w');
fprintf(fout,'# Output from GTdef_detect_greensfns.m\n');
fprintf(fout,'# 1   2     3     4      5      6      7        8        9\n');
fprintf(fout,'# ID  dnum  snum  ssmin  dsmin  tsmin  smagmin  dmagmin  tmagmin\n');
fprintf(fout,'%-6d %-4d %-4d %14.5e %14.5e %14.5e %10.2f %8.2f %8.2f\n',limits');
fclose(fout);
