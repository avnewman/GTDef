function [ lon,lat ] = Haversine( lon0,lat0,bearing,distance )
%lon0, lat0, and bearing should be in degrees, distance in km0
% copied from awilliamson's tectonic_toolbox (on masaya)
% by AVN 9/25/2020

RR = 6378137;     % m   raidus at equator
ff=1/298.257 ;     %  flattening factor
R=RR*(1-ff*sind(lat0)^2);  %radius at point

d=distance*1000;
lon0=deg2rad(lon0);
lat0=deg2rad(lat0);
b=deg2rad(bearing);




lat = asin(sin(lat0)*cos(d/R) + cos(lat0)*sin(d/R)*cos(b));


lon = lon0 + atan2(sin(b)*sin(d/R)*cos(lat0),(cos(d/R)-sin(lat0)*sin(lat)));

lat=rad2deg(lat);  lon=rad2deg(lon);


end
