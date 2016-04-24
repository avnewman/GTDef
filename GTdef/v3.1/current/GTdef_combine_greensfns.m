function [ ] = GTdef_combine_greensfns(fsiteName,fgrnName,fsurfName,parsMat)

% not updated after GTdef_save_greensfns.m & GTdef_read_greensfns.m changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              GTdef_combine_greensfns			             %
% combine green's functions to generate deformation at sites                         %
%                                                                                    %
% INPUT:                                                                             %
% fsiteName - *.site location file                                                   %
% 1      2              3              4    			                     %
% Site	 Lon		Lat	       Height [m]			             %
% ABGS   99.387520914   0.220824642    236.2533			                     %
%                                                                                    %
% fgrnName  - GTdef greens function file                                             %
% (1) location of gps sites                                                          %
%     1    2    3   4   5                                                            %
%     Site LOC  lon lat z                                                            %
% (2) location of vertices                                                           %
%     1  2    3    4    5   6    7                                                   %
%     id dnum snum LOC  lon lat  z                                                   %
% (3) greens functions of patch-site pairs                                           %
%     1  2    3    4     5     6     7     8     9     10    11    12    13          %
%     id dnum snum Site  ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz       %
%                                                                                    %
% grnList - site names stored as a cell          (pntNum*1)                          %
% grnfns  - array of size                        (patchNum*pntNum*12)                %
%            for each patch-site pair                                                %
%       1  2    3    4     5     6     7     8     9     10    11    12              %
%   = [ id dnum snum ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz ]         %
%                                                                                    %
% fsurfName - GTdef surface projection file                                          %
% 1    2    3    4     5     6     7     8     9     10    11    12   13   14 15 16  %
% name dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2 xctr yctr ss ds ts  %
%                                                                                    %
% parsMat - parameter matrix                                                         %
%         = [ ssize dsize sstep dstep sslip dslip tslip ]                            % 
%   combine dsize patches along dip & ssize patches along strike                     %
%   and multiply values in grnfns with dslip, sslip & tslip                          %
%                                                                                    %
% OUTPUT:                                                                            %
% a new greens function file and a corresponding surface projection file             %
%                                                                                    %
% first created by Lujia Feng Tue Aug  6 09:15:07 SGT 2013                           %
% added dslip sslip tslip & coded odd-even mix by lfeng Thu Dec  5 20:16:41 SGT 2013 %
% last modified by Lujia Feng Fri Dec  6 01:41:20 SGT 2013                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read sites
fprintf(1,'\n.......... reading %s ...........\t',fsiteName);
tic
[ siteList,locList ] = GPS_readsites(fsiteName);
toc

% read green's functions
fprintf(1,'\n.......... reading %s ...........\t',fgrnName);
tic
[ grnList,grnfns ] = GTdef_read_greensfns(fgrnName);
toc

% read surface projection
fprintf(1,'\n.......... reading %s ...........\t',fsurfName);
tic
fin = fopen(fsurfName,'r');
%                        1   2  3  4  5  6  7  8  9  10 11 12 13 14 15 16
surfCell = textscan(fin,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
surfMat  = cell2mat(surfCell);
fclose(fin);
toc

% read in parameters
patchNum = grnfns(end,1);
dnum     = grnfns(end,2);
snum     = grnfns(end,3);
% check file consistency
if dnum ~= surfMat(end,1) || snum ~= surfMat(end,2) 
    error('GTdef_combine_greensfns ERROR: %s and %s do not match!',fgrnName,fsurfName);
end
dlist    = grnfns(:,2);
slist    = grnfns(:,3);
ind      = dlist==1 & slist==1;
siteNum  = sum(ind);
siteCell = grnList(ind);

fprintf(1,'\n total patch number along fault  = %.0f\n',patchNum);
fprintf(1,' total patch number along dip    = %.0f\n',dnum);
fprintf(1,' total patch number along strike = %.0f\n',snum);
fprintf(1,' total site number               = %d\n',siteNum);

ssize = parsMat(1); 
dsize = parsMat(2); 
sstep = parsMat(3); 
dstep = parsMat(4); 
sslip = parsMat(5);
dslip = parsMat(6);
tslip = parsMat(7);

% position of center point in the combined patch
if mod(dsize,2)==1 && mod(ssize,2)==1  % odd numbers
    dcenter = 0.5*(dsize-1);
    scenter = 0.5*(ssize-1);
    cposLon = 11; % index of xctr in one line of surface file
    cposLat = 12; % index of yctr in one line of surface file
elseif mod(dsize,2)==0 && mod(ssize,2)==0  % even numbers
    dcenter = 0.5*dsize-1;
    scenter = 0.5*ssize-1;
    cposLon = 7; % index of xbot2 in one line of surface file
    cposLat = 8; % index of ybot2 in one line of surface file
elseif mod(dsize,2)==1 && mod(ssize,2)==0  % dsize odd; ssize even 
    dcenter = 0.5*(dsize-1);
    scenter = 0.5*ssize-1;
    cposLon = [7 9];  % index of xbot2 xtop2 in one line of surface file
    cposLat = [8 10]; % index of ybot2 ytop2 in one line of surface file
elseif mod(dsize,2)==0 && mod(ssize,2)==1  % dsize even; ssize odd
    dcenter = 0.5*dsize-1;
    scenter = 0.5*(ssize-1);
    cposLon = [5 7]; % index of xbot1 xbot2 in one line of surface file
    cposLat = [6 8]; % index of ybot1 ybot2 in one line of surface file
end

% output grnfns file
[ ~,basename,~ ] = fileparts(fgrnName);
fout1Name = [ basename '_' num2str(dsize,'%02d') 'x' num2str(ssize,'%02d') ...
              '_ss' num2str(sslip*1e3,'%04d') 'mm_ds' num2str(dslip*1e3,'%04d') 'mm.grnfns' ];
fout1 = fopen(fout1Name,'w');

%% output station position (optional)
fprintf(fout1,'# Output from GTdef_combine_greensfns.m\n');

%fprintf(fout1,'# 1     2    3    4    5\n');
%fprintf(fout1,'# Site  LOC  Lon  Lat  Depth (m)\n');
%for ii=1:siteNum
%    site = siteCell{ii};
%    ind  = strcmpi(site,siteList);
%    loc  = locList(ind,:);
%    fprintf(fout1,'%4s LOC %12.5f %12.5f %14.5f\n',site,loc);
%end

fprintf(fout1,'# 1   2     3     4     5    6    7     8    9    10    11   12   13\n');
fprintf(fout1,'# ID  dnum  snum  Site  ss_E ss_N ss_U  ds_E ds_N ds_U  ts_E ts_N ts_U (m)\n');

% output surface file
[ ~,basename,~ ] = fileparts(fout1Name);
fout2Name = [ basename '_surface.out' ];
fout2 = fopen(fout2Name,'w');
fprintf(fout2,'# Output from GTdef_combine_greensfns.m\n');
fprintf(fout2,'#(1)name (2)dnum (3)snum (4)xtop1 (5)ytop1 (6)xbot1 (7)ybot1 (8)xbot2 (9)ybot2 (10)xtop2 (11)ytop2 (12)xctr (13)yctr (14)ss[m] (15)ds[m] (16)ts[m]\n'); 

% loop along strike
fprintf(1,'\n.......... calculating %d by %d patches with %d by %d offsets ...........\t',dsize,ssize,dstep,sstep);
tic
ss = 0;
tt = 0;   % total
sstart = 1;
while(1)
    ss   = ss + 1;
    send = sstart + ssize - 1;
    if(send>snum), break; end
    % loop along dip
    dd = 0;
    dstart = 1; 
    while(1)
        dd = dd + 1;
        dend = dstart + dsize - 1;
        if(dend>dnum), break; end
        tt = tt + 1;
        % 1  2    3    4     5     6     7     8     9     10    11    12    13
        % id dnum snum Site  ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz
        ind = dlist>=dstart & dlist<=dend & slist>=sstart & slist<=send;
        grnpatch  = grnfns(ind,:); % deformation from combined patch
        sitepatch = grnList(ind);  % corresponding sites
        % loop through stations
        for ii = 1:siteNum 
            site = siteCell{ii};
            ind  = strcmpi(site,sitepatch);
            grnsite = grnpatch(ind,:); % deformation for one site from combined patch
            if (dsize*ssize)~=size(grnsite,1) % check if number of patches is correct
                error('GTdef_combine_greensfns ERROR: patch size does not match!');
            end
            % sum up all rows into one row
            sgrn = sslip*sum(grnsite(:,4:6),1); % strike slip
            dgrn = dslip*sum(grnsite(:,7:9),1); % dip sip
	    tgrn = tslip*sum(grnsite(:,10:end),1); % tensile slip
            fprintf(fout1,'%-6d %-4d %-4d %4s  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',...
                    tt,dd,ss,site,sgrn,dgrn,tgrn);
        end
        % top left point of the combined patch
        ind   = surfMat(:,1)==dstart & surfMat(:,2)==sstart; 
        xtop1 = surfMat(ind,3);
        ytop1 = surfMat(ind,4); 
        % bottom left point of the combined patch
        ind   = surfMat(:,1)==dend & surfMat(:,2)==sstart; 
        xbot1 = surfMat(ind,5);
        ybot1 = surfMat(ind,6); 
        % bottom right point of the combined patch
        ind   = surfMat(:,1)==dend & surfMat(:,2)==send; 
        xbot2 = surfMat(ind,7);
        ybot2 = surfMat(ind,8); 
        % top right point of the combined patch
        ind   = surfMat(:,1)==dstart & surfMat(:,2)==send; 
        xtop2 = surfMat(ind,9);
        ytop2 = surfMat(ind,10); 
        % center point of the combined patch
        ind   = surfMat(:,1)==(dstart+dcenter) & surfMat(:,2)==(sstart+scenter); 
	if isscalar(cposLon)
            xctr = surfMat(ind,cposLon);
            yctr = surfMat(ind,cposLat); 
        else
            lon = surfMat(ind,cposLon);
            lat = surfMat(ind,cposLat); 
            [ xx,yy ]     = LL2ckmd(lon(2),lat(2),lon(1),lat(1),0);
            [ xctr,yctr ] = ckm2LLd(0.5*xx,0.5*yy,lon(1),lat(1),0);
        end
        patch = [ tt dd ss xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2 xctr yctr sslip dslip tslip ];
        fprintf(fout2,'%-6d %-3d %-3d %12.5f %11.5f %12.5f %11.5f %12.5f %11.5f %12.5f %11.5f %12.5f %11.5f %10.5f %10.5f %10.5f\n',patch);
        dstart = dstart + dstep;
    end
    sstart = sstart + sstep;
end
toc

fclose(fout1); fclose(fout2);
