function [lon,lat] = xy_to_latlon(xx,yy,lon0,lat0,rot)

% From PCAIM
%
%Converts from local coorindates to longitude and latitude 
%given the [lon, lat] of an origin. 'origin' should be in 
%decimal degrees. Note that heights are ignored and that 
%xy is in km.  llh is [lon, lat, height] in decimal 
%degrees. This is an iterative solution for the inverse of 
%a polyconic projection.

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Aug 23, 2001  Jessica Murray        Clarification to help.
%
%   Apr 4, 2001   Peter Cervelli        Added failsafe to avoid
%                                       infinite loop because of
%                                       covergence failure.
%   Sep 7, 2000   Peter Cervelli		Original Code
%
%-------------------------------------------------------------

% modified to be compatible with GTdef.m lfeng Thu Jun  7 12:41:44 SGT 2012
% input x & y can be row or column vectors
% output lat & lon would be row or column vectors accordingly
% no need for unit conversion

% added rotation to be compatible with ckm2LLd.m Lujia Feng Tue Sep 10 14:31:45 SGT 2013
% rot [deg] is the rotation angle from the old coord system (east is x+ ;north is y+) to the new one.
% Counterclockwise rotation is positive; Clockwise rotation is negative
% use rot for latlon_to_xy(); use -rot for xy_to_latlon()

cos_rot = cosd(rot);        % cos_rot - cos of rotation angle
sin_rot = sind(rot);        % sin_rot - sin of rotation angle
% rotate the coordinate system by rot [degree]
% back to the state of horizontal x-axis
x_rot = xx*cos_rot + yy*sin_rot;
y_rot =-xx*sin_rot + yy*cos_rot;

origin = [lon0 lat0];
[mm,nn] = size(x_rot);
% column vectors
if mm>=nn
    xy = [x_rot y_rot]';
% row vectors (need to be converted to row vectors)
else
    xy = [x_rot;y_rot];
end

%Set ellipsoid constants (WGS84)

   a=6378137.0;
   e=0.08209443794970;

%Convert to radian

   origin=origin*pi/180;

%Iterate to perform inverse projection

   M0=a*((1-e^2/4-3*e^4/64-5*e^6/256)*origin(2) - ...
        (3*e^2/8+3*e^4/32+45*e^6/1024)*sin(2*origin(2)) + ...
        (15*e^4/256 +45*e^6/1024)*sin(4*origin(2)) - ...
        (35*e^6/3072)*sin(6*origin(2)));

   z=xy(2,:)~=-M0;

   A=(M0+xy(2,z))/a;
   B=xy(1,z).^2./a^2+A.^2;

   llh(2,z)=A;

   delta=Inf;

   c=0;
   
   while max(abs(delta))>1e-8

      C=sqrt((1-e^2*sin(llh(2,z)).^2)).*tan(llh(2,z));

      M=a*((1-e^2/4-3*e^4/64-5*e^6/256)*llh(2,z) - ...
           (3*e^2/8+3*e^4/32+45*e^6/1024)*sin(2*llh(2,z)) + ...
           (15*e^4/256 +45*e^6/1024)*sin(4*llh(2,z)) - ...
           (35*e^6/3072)*sin(6*llh(2,z)));

      Mn=1-e^2/4-3*e^4/64-5*e^6/256 - ...
         -2*(3*e^2/8+3*e^4/32+45*e^6/1024)*cos(2*llh(2,z)) + ...
         4*(15*e^4/256 +45*e^6/1024)*cos(4*llh(2,z)) + ...
         -6*(35*e^6/3072)*cos(6*llh(2,z));

      Ma=M/a;
   
      delta=-(A.*(C.*Ma+1)-Ma-0.5*(Ma.^2+B).*C)./ ...
           (e^2*sin(2*llh(2,z)).*(Ma.^2+B-2*A.*Ma)./(4*C)+(A-Ma).*(C.*Mn-2./sin(2*llh(2,z)))-Mn);

      llh(2,z)=llh(2,z)+delta;

      c=c+1;
      if c>100
          error('xy_to_latlon ERROR: Convergence failure!');
      end
   end

   llh(1,z)=(asin(xy(1,z).*C/a))./sin(llh(2,z))+origin(1);

%Handle special case of latitude = 0

   llh(1,~z)=xy(1,~z)/a+origin(1);
   llh(2,~z)=0;

%Convert back to decimal degrees

   llh=llh*180/pi;
   
   
%Variables to return
% column vectors
if mm>=nn
    lon = llh(1,:)';
    lat = llh(2,:)';
% row vectors
else
    lon = llh(1,:);
    lat = llh(2,:);
end
