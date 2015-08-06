function [x_rot,y_rot] = latlon_to_xy (lon,lat,lon0,lat0,rot)

% Convert lat-lon-height to local xy via a given origin
% Based on code from PCAIM (Kositsky/Perfettini)
% Edited by Emma, 19 Mar 2011
%
% Inputs list of lat/lon coords and origin coords

% modified to be compatible with GTdef Lujia Feng Thu Jun  7 12:42:22 SGT 2012
% input lon & lat can be row or column vectors
% output x & y would be row or column vectors accordingly
% no need for unit conversion

% added rotation to be compatible with LL2ckmd.m Lujia Feng Tue Sep 10 14:31:45 SGT 2013
% rot [deg] is the rotation angle from the old coord system (east is x+ ;north is y+) to the new one.
% Counterclockwise rotation is positive; Clockwise rotation is negative
% use rot for latlon_to_xy(); use -rot for xy_to_latlon()

Nstn = length(lat);

%% Convert to decimal seconds
lat = 3600 * lat;
lon = 3600 * lon;
lat0 = 3600 * lat0;

diffLon = 3600 * lon0*ones(size(lon)) - lon;

%% Make lat or lon values that are exactly zero a very small number, otherwise it gives NaN
lat(find(lat==0)) = 0.000000001;
lon(find(lon==0)) = 0.000000001;

%% Initialize output
xy = zeros(Nstn,2);

%% Loop through and do projection
for k = 1:Nstn
     xy(k,:) = polyconic(lat(k), diffLon(k), lat0);
end

%% flip x-axis x & y are column vectors
xx = -xy(:,1);
yy =  xy(:,2);

[mm,nn] = size(lat);
% convert x & y to row vectors if inputs are row vectors
if mm<nn
   xx = xx';  
   yy = yy';
end

% rotate the coordinate system by rot [degree]
cos_rot = cosd(rot);        % cos_rot - cos of rotation angle
sin_rot = sind(rot);        % sin_rot - sin of rotation angle
x_rot = xx*cos_rot + yy*sin_rot;
y_rot =-xx*sin_rot + yy*cos_rot;
