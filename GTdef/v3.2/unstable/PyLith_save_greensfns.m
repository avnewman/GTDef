function [  ] = PyLith_save_greensfns(vertices,gpsites,siteList,grnfns)

%%%%%% need to modify according to GTdef_save_greensfns.m!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             PyLith_save_greensfns                             %
% save PyLith greensfns after reordering vertices and stations                  %
%                                                                               %
% INPUT:                                                                        %
% vertices.num  - number of vertices               (scalar)                     %
% vertices.step - number of time steps             (scalar)                     %
% vertices.id   - 1-based id number for vertices   (vertexNum*1)                %
% vertices.loc  - location of vertices in global coordinate system              %
%               = [ x y z ]                        (vertexNum*3)                %
% vertices.slip - slip of vertices in fault coordinate                          %
%               = [ leftlateral reverse opening ]  (timeStep*vertexNum*3)       %
% vertices.dT   - traction change                  (timeStep*vertexNum*3)       %
%               = [ shear-leftlateral shear-updip normal ]                      %
% cells.num     - number of cells                  (scalar)                     %
% cells.id      - 1-based id number for cells      (vertexNum*1)                %
% cells.topo    - topology of cells (4 vertices)   (cellNum*4)                  %
% gpsites.num   - number of gps sites              (scalar)                     %
% gpsites.step  - number of time steps             (scalar)                     %
% gpsites.id    - 1-based id number for gps sites (PyLith randomizes site order)%
% gpsites.loc   - location of gps sites in global coordinate system             %
%               = [ x y z ]                        (pntNum*3)                   %
% gpsites.disp  - greens function for gps sites                                 %
%               = [ dx dy dz ]                     (timeStep*pntNum*3)          %
% siteList      - site names stored as a cell      (pntNum*1)                   %
% grnfns        - array of size                    (vertexNum*pntNum*9)         %
%               for each patch-site pair                                        % 
%               = [ ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz ]     %
% ----------------------------------------------------------------------------- %
%                                                                               %
% OUTPUT:                                                                       %
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
% first created by lfeng Fri Nov 30 12:38:16 SGT 2012                           %
% last modified by lfeng Fri Nov 30 13:13:55 SGT 2012                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout = fopen('pylith.grnfns','w');

fprintf(fout,'# Output from PyLith_save_greensfns.m\n');
fprintf(fout,'# 1     2    3   4   5\n');
fprintf(fout,'# Site  LOC  XX  YY  ZZ (m)\n');
% site location in pylith coordinate system
for ii=1:gpsites.num
    site = siteList{ii};
    loc  = gpsites.loc(ii,:);
    fprintf(fout,'%4s LOC %14.5e %14.5e %14.5e\n',site,loc);
end

% vertex location in pylith coordinate system
fprintf(fout,'# 1   2     3     4    5   6   7\n');
fprintf(fout,'# ID  dnum  snum  LOC  XX  YY  ZZ (m)\n');
vertSize = [ vertices.Nd vertices.Ns ];
for ii=1:vertices.num
    [ dd,ss ] = ind2sub(vertSize,ii);
    loc = vertices.loc(ii,:);
    fprintf(fout,'%-6d %-4d %-4d LOC %14.5e %14.5e %14.5e\n',ii,dd,ss,loc);
end

% green fucntions for vertex-site pairs
fprintf(fout,'# 1   2     3     4     5    6    7     8    9    10    11   12   13\n');
fprintf(fout,'# ID  dnum  snum  Site  ss_X ss_Y ss_Z  ds_X ds_Y ds_Z  ts_X ts_Y ts_Z (m)\n');
for ii=1:vertices.num
    [ dd,ss ] = ind2sub(vertSize,ii);
    for jj=1:gpsites.num
        site = siteList{jj};
        grn = grnfns(ii,jj,:);
        fprintf(fout,'%-6d %-4d %-4d %4s  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',...
	        ii,dd,ss,site,grn);
    end
end

fclose(fout);
