function [  ] = GTdef_save_greensfns(fltName,pnt,flt,Xgrn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_save_greensfns                              %
% save greensfns for one fault                                                  %
%                                                                               %
% INPUT:                                                                        %
% fltName - fault name                                                          %
% pnt - point structure                                                         %
%   read in from input                                                          %
%     pnt.num  - number of point data             (scalar)                      %
%     pnt.name - names of point data              (cell array)                  %
%     pnt.loc  - [lon lat z]                      (nn*3)                        % 
%     pnt.disp - [east north vert]                (nn*3)                        %
%     pnt.err  - [east north vert]                (nn*3)                        %
%     pnt.wgt  - [weight]                         (nn*1)                        % 
%   added later                                                                 %
%     pnt.crt  -                                                                %   
%     pnt.obs  -                                                                %
%     pnt.obs_err -                                                             %
%     pnt.obs_wgt -                                                             %
%     pnt.coef -                                                                %
% flt - flt1, flt2, flt3, flt4, and flt5 have different formats                 %
%       but they all have last two terms as Nd & Ns                             %
% Xgrn - green's function array                   (3*siteNum)*slipNum           %
%      = [east;north;vertical] for different sites                              % 
%        from unit slips [strike dip tensile] slipNum = 3*patchNum              %
% ----------------------------------------------------------------------------- %
%                                                                               %
% OUTPUT:                                                                       %
% a green function database file listing                                        %
% (1) location of gps sites                                                     %
%     1    2    3   4   5                                                       %
%     Site LOC  lon lat z                                                       %
% (2) location of vertices                                                      %
%     1  2    3    4    5   6    7                                              %
%     id dnum snum LOC  lon lat  z                                              %
% (3) greens functions of patch-site pairs                                      %
%     1  2    3    4     5     6     7     8     9     10    11    12    13     %
%     id dnum snum Site  ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz  %
%                                                                               %
% first created by Lujia Feng Mon Aug  5 14:07:48 SGT 2013                      %
% modified based on PyLith_save_greensfns.m lfeng Mon Aug  5 14:08:24 SGT 2013  %
% last modified by Lujia Feng Mon Aug  5 17:06:55 SGT 2013                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

foutName = [ fltName '.grnfns' ];
fout     = fopen(foutName,'w');
pntNum   = pnt.num;

fprintf(fout,'# Output from GTdef_save_greensfns.m\n');
%fprintf(fout,'# 1     2    3   4   5\n');
%fprintf(fout,'# Site  LOC  Lon Lat ZZ [m]\n');

%% site location
%for ii=1:pntNum
%    site = pnt.name{ii};
%    loc  = pnt.loc(ii,:);
%    fprintf(fout,'%4s LOC %14.5e %14.5e %14.5e\n',site,loc);
%end

% haven't corrected for GTdef yet
%% patch location
%fprintf(fout,'# 1   2     3     4    5   6   7\n');
%fprintf(fout,'# ID  dnum  snum  LOC  XX  YY  ZZ (m)\n');
%vertSize = [ vertices.Nd vertices.Ns ];
%for ii=1:vertices.num
%    [ dd,ss ] = ind2sub(vertSize,ii);
%    loc = vertices.loc(ii,:);
%    fprintf(fout,'%-6d %-4d %-4d LOC %14.5e %14.5e %14.5e\n',ii,dd,ss,loc);
%end

% green fucntions for subpatch-site pairs
fprintf(fout,'# 1   2     3     4     5    6    7     8    9    10    11   12   13\n');
fprintf(fout,'# ID  dnum  snum  Site  ss_E ss_N ss_U  ds_E ds_N ds_U  ts_E ts_N ts_U [m]\n');
dnum      = flt(end-1);
snum      = flt(end);
patchNum  = dnum*snum;
patchSize = [ dnum snum ];
% the number of subpatches
for ii=1:patchNum
    [ dd,ss ] = ind2sub(patchSize,ii);
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
