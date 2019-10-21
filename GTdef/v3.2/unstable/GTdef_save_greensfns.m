function [  ] = GTdef_save_greensfns(fltName,pnt,prjflt,Xgrn,coord,lon0,lat0,rot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_save_greensfns                              %
% save greensfns for one fault                                                  %
%                                                                               %
% INPUT:                                                                        %
% fltName - fault name                                                          %
% pnt     - point structure                                                     %
%           pnt.num  - number of point data             (scalar)                %
%           pnt.name - names of point data              (cell array)            %
%           pnt.loc  - [lon lat z]                      (nn*3)                  % 
%             1    2    3     4     5     6     7     8                         %
% prjflt  = [ dnum snum xtop1 ytop1 ztop1 xbot1 ybot1 zbot1                     %
%            9     10    11    12    13    14    15   16   17                   %
%            xbot2 ybot2 zbot2 xtop2 ytop2 ztop2 xctr yctr zctr                 %
%            18 19 20 21   22                                                   %
%            ss ds ts rake rs ]                                                 %
%                                                                               %
% Xgrn    - green's function array                   (3*siteNum)*slipNum        %
%         = [east;north;vertical] for different sites                           % 
%           from unit slips [strike dip tensile] slipNum = 3*patchNum           %
% coord,lon0,lat0,rot - parameters for coordinate conversion                    %
% ----------------------------------------------------------------------------- %
%                                                                               %
% OUTPUT:                                                                       %
% a green function database file listing                                        %
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
% first created by Lujia Feng Mon Aug  5 14:07:48 SGT 2013                      %
% modified based on PyLith_save_greensfns.m lfeng Mon Aug  5 14:08:24 SGT 2013  %
% output point & patch lfeng Fri Nov 14 15:47:17 SGT 2014                       %
% added coord, lon0, lat0, and rot lfeng Fri Nov 14 16:29:01 SGT 2014           %
% modified prjflt to include [ ss ds ts rake rs ] lfeng Wed Apr 27 SGT 2016     %
% last modified by Lujia Feng Wed Apr 27 23:07:41 SGT 2016                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

foutName = [ fltName '.grnfns' ];
fout     = fopen(foutName,'w');
dnum     = prjflt(end,1);
snum     = prjflt(end,2);
patchNum = dnum*snum;

fprintf(fout,'# Output from GTdef_save_greensfns.m\n');

% site location
fprintf(fout,'# (1)point (2)num (3)site (4)lon (5)lat (6)height[m]\n');
pntNum = pnt.num;
for ii=1:pntNum
   site = pnt.name{ii};
   loc  = pnt.loc(ii,:);
   fprintf(fout,'point  %5d %8s %12.5f %11.5f %12.3e\n',ii,site,loc);
end

% patch location
switch coord
   case 'geo'
      xx = prjflt(:,[3 6 9 12 15]);  yy = prjflt(:,[4 7 10 13 16]);
      [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,rot);
      newprjflt = [ prjflt(:,1:2) lon(:,1) lat(:,1) prjflt(:,5) lon(:,2) lat(:,2) prjflt(:,8) lon(:,3) lat(:,3) prjflt(:,11) ...
                       lon(:,4) lat(:,4) prjflt(:,14) lon(:,5) lat(:,5) prjflt(:,17) prjflt(:,18:22) ];
   case 'geo_polyconic'
      xx = prjflt(:,[3 6 9 12 15]);  yy = prjflt(:,[4 7 10 13 16]);
      lon = zeros(size(xx));         lat = zeros(size(yy));
      for ii=1:5
         [lon(:,ii),lat(:,ii)] = xy_to_latlon(xx(:,ii),yy(:,ii),lon0,lat0,rot);
      end        
      newprjflt = [ prjflt(:,1:2) lon(:,1) lat(:,1) prjflt(:,5) lon(:,2) lat(:,2) prjflt(:,8) lon(:,3) lat(:,3) prjflt(:,11) ...
                    lon(:,4) lat(:,4) prjflt(:,14) lon(:,5) lat(:,5) prjflt(:,17) prjflt(:,18:22) ];
   case 'local'
      newprjflt = prjflt;
end

fprintf(fout,'# (1)patch (2)ID (3)dnum (4)snum (5)xtop1 (6)ytop1 (7)ztop1[m] (8)xbot1 (9)ybot1 (10)zbot1[m] (11)xbot2 (12)ybot2 (13)zbot2[m] (14)xtop2 (15)ytop2 (16)ztop2[m] (17)xctr (18)yctr (19)zcrt[m]\n');
id  = [ 1:patchNum ]';
out = [ id newprjflt(:,1:17) ];
%             1     2    3   4   5      6      7      8      9      10     11     12     13     14     15     16     17     18     19
fprintf(fout,'patch %10d %4d %4d %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e \n',out');

% green fucntions for subpatch-site pairs
fprintf(fout,'# 1   2     3     4     5    6    7     8    9    10    11   12   13\n');
fprintf(fout,'# ID  dnum  snum  Site  ss_E ss_N ss_U  ds_E ds_N ds_U  ts_E ts_N ts_U [m]\n');
patchDSNum = [ dnum snum ];
% the number of subpatches
for ii=1:patchNum
    [ dd,ss ] = ind2sub(patchDSNum,ii);
    % the number of points
    for jj=1:pntNum
        site = pnt.name{jj};
        % individual components
        ss_de = Xgrn(jj,ii);
        ss_dn = Xgrn(jj+pntNum,ii);
        ss_du = Xgrn(jj+pntNum*2,ii);
        ds_de = Xgrn(jj,ii+patchNum);
        ds_dn = Xgrn(jj+pntNum,ii+patchNum);
        ds_du = Xgrn(jj+pntNum*2,ii+patchNum);
        ts_de = Xgrn(jj,ii+patchNum*2);
        ts_dn = Xgrn(jj+pntNum,ii+patchNum*2);
        ts_du = Xgrn(jj+pntNum*2,ii+patchNum*2);
        grn = [ ss_de ss_dn ss_du ds_de ds_dn ds_du ts_de ts_dn ts_du ];
        fprintf(fout,'%-6d %-4d %-4d %4s  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',...
                ii,dd,ss,site,grn);
    end
end

fclose(fout);
