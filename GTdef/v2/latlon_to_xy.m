function [x, y] = latlon_to_xy (lon,lat,lon0,lat0)

%% Convert lat-lon-height to local xy via a given origin
%% Based on code from PCAIM (Kositsky/Perfettini)
%% Edited by Emma, 19 Mar 2011
%%
%% Inputs list of lat/lon coords and origin coords

% modified to be compatible with GTdef lfeng Thu Jun  7 12:42:22 SGT 2012
% input lon & lat can be row or column vectors
% output x & y would be row or column vectors accordingly
% no need for unit conversion

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
x = -xy(:,1);
y =  xy(:,2);

[mm,nn] = size(lat);
% convert x & y to row vectors if inputs are row vectors
if mm<nn
   x = x';  y = y';
end
