function [ limits ] = GTdef_detect_greensfns(fgrnName,threshold,sitenum,area,mu)

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
%     id dnum snum Site  ss_E  ss_N  ss_U  ds_E  ds_N  ds_U  ts_E  ts_N  ts_u   %
% threshold - minimul level of detection for each component [m]                 %
%           = [ ee nn uu ]       (1*3)                                          %
% sitenum   - number of sites that detect the signal                            %
% area      - area of each patch [m^2]           (scalar or vector)             %
%    if area = [], M0 & Mw are not calculated                                   %
%                                                                               %
% grnList - site names stored as a cell          (pntNum*1)                     %
% grnfns  - array of size                        (patchNum*pntNum*12)           %
%            for each patch-site pair                                           % 
%       1  2    3    4    5    6    7    8    9    10   11   12                 %
%   = [ id dnum snum ss_E ss_N ss_U ds_E ds_N ds_U ts_E ts_N ts_U ]             %
%                                                                               %
% OUTPUT:                                                                       %
% *.min records minimum slip & magnitdue needed for detection of threshold      %
%            1  2    3    4     5     6     7       8       9                   %
% limits = [ id dnum snum ssmin dsmin tsmin smagmin dmagmin tmagmin ]           %
%                                                                               %
% first created by Lujia Feng Fri Dec  6 09:05:58 SGT 2013                      %
% output patch location lfeng Fri Nov 14 17:39:39 SGT 2014                      %
% added thres3 lfeng Fri Jun 12 14:26:30 SGT 2015                               %
% added and removed fsiteName lfeng Fri Jun 12 19:25:35 SGT 2015                %
% area can be scalar or vector lfeng Wed Aug  5 17:47:47 SGT 2015               %
% added threshold for number of stations lfeng Wed Aug 12 13:18:49 SGT 2015     %
% added rigidity mu lfeng Thu Aug 27 15:19:08 SGT 2015                          %
% last modified by Lujia Feng Thu Aug 27 15:20:44 SGT 2015                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(threshold)~=[1 3], error('GTdef_detect_greensfns ERROR: need a 1*3 vector for threshold as input!'); end

% read green's functions
fprintf(1,'\n.......... reading %s ...........\t',fgrnName);
tic
[ ~,~,patchloc,~,grnfns ] = GTdef_read_greensfns(fgrnName);
toc

% read in parameters
patchNum  = grnfns(end,1);
patchList = grnfns(:,1); 

thres3 = [ threshold threshold threshold ];

% check area
if ~isscalar(area)
   areaNum = length(area);
   if areaNum~=patchNum, error('GTdef_detect_greensfns ERROR: the number of areas does not match that of patches!'); end
end

% loop through all patches
fprintf(1,'\n.......... processing %d patches ...........\t',patchNum);
tic
ssM0 = nan; ddM0 = nan; ttM0 = nan;
ssMw = nan; ddMw = nan; ttMw = nan;
limits = zeros(patchNum,9);
for ii=1:patchNum
    ind = patchList==ii;
    grnpatch = grnfns(ind,:);
    id       = grnpatch(1,1);
    dd       = grnpatch(1,2);
    ss       = grnpatch(1,3);
    disps    = grnpatch(:,4:end);
    %slips    = abs(threshold./disps);
    slips    = abs(bsxfun(@rdivide,thres3,disps)); % absolute values!! slips = matrix [siteNum 9]
    sslist   = min(slips(:,1:3),[],2);
    ddlist   = min(slips(:,4:6),[],2);
    ttlist   = min(slips(:,7:9),[],2);
    sslist   = sort(sslist);
    ddlist   = sort(ddlist);
    ttlist   = sort(ttlist);
    ssmin    = sslist(sitenum);
    ddmin    = ddlist(sitenum);
    ttmin    = ttlist(sitenum);
    %%-------- minimum --------
    %slipmin  = min(slips); % detected by at least one site
    %% minimum slip
    %ssmin    = min(slipmin(1:3));
    %ddmin    = min(slipmin(4:6));
    %ttmin    = min(slipmin(7:9));
    if ~isempty(area)
        if ~isscalar(area)
           parea = area(ii);
	else
	   parea = area;
	end
        % minimum M0
        [ ssM0,~,~ ] = EQS_M0AreaSlip(mu,[],parea,ssmin);
        [ ddM0,~,~ ] = EQS_M0AreaSlip(mu,[],parea,ddmin);
        [ ttM0,~,~ ] = EQS_M0AreaSlip(mu,[],parea,ttmin);
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
foutName = [ basename '_disp' num2str(threshold*1e3,'%03.0f') 'mm_' num2str(sitenum,'%d') 'sites.min' ];
fout = fopen(foutName,'w');
fprintf(fout,'# Output from GTdef_detect_greensfns.m\n');
fprintf(fout,'# (1)ID (2)dnum (3)snum (4)xtop1 (5)ytop1 (6)ztop1[m] (7)xbot1 (8)ybot1 (9)zbot1[m] (10)xbot2 (11)ybot2 (12)zbot2[m] (13)xtop2 (14)ytop2 (15)ztop2[m] (16)xctr (17)yctr (18)zcrt[m]\n');
fprintf(fout,'# (19)l (20)ssmin[m] (21)dsmin[m] (22)tsmin[m]  (23)smagmin  (24)dmagmin  (25)tmagmin\n');
out = [ patchloc limits(:,4:end) ];
fprintf(fout,'%-8d %-4d %-4d %10.4f %9.4f %10.1e %10.4f %9.4f %10.1e %10.4f %9.4f %10.1e %10.4f %9.4f %10.1e %10.4f %9.4f %10.1e l %14.5e %14.5e %14.5e %7.2f %7.2f %7.2f \n',out');
fclose(fout);
