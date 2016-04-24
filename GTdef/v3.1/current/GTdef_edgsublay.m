function [ sublayer ] = GTdef_edgsublay(layer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   GTdef_edgsublay.m                                        %
%                                                                                            %
% Create sublayers that are thin enough to represent the change of properties                %
% This code is converted from fortran code edgsublay.F                                       %
%                                                                                            %
% INPUT											     %
% (earth.)layer = [ id depth vp vs ro ]	(nn*5)                                               %
%                                                                                            %
% OUTPUT                                                                                     %
% earth.sublayer structure (model sublayers)                                                 %
% corresponds to /model/ h,ro,vp,vs,n0 [h->hh,n0->nl]                                        %
% (earth.)sublayer.nl   - number of sublayers                                      (scalar)  %
% (earth.)sublayer.topz - top depth of each sublayer [m]               (nl*1 column vector)  % 
% (earth.)sublayer.botz - bottom depth of each sbulayer [m]            (nl*1 column vector)  % 
% (earth.)sublayer.hh   - thickness of each sublayer [m]               (nl*1 column vector)  %    
% (earth.)sublayer.vp   - P-wave velocity for each sublayer [m/s]      (nl*1 column vector)  %    
% (earth.)sublayer.vs   - S-wave velocity for each sublayer [m/s]      (nl*1 column vector)  %    
% (earth.)sublayer.ro   - density for each sublayer [kg/m^3]           (nl*1 column vector)  %    
%                                                                                            %
% REFERENCE  										     %
% Wang, R., Martin, F. L., & Roth, F. (2003)						     %
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5		     %
%                                                                                            %
% first created by Lujia Feng Tue Sep  8 19:11:20 SGT 2015                                   %
% changed topz,botz,hh,vp,vs,ro from rowwise to columnwise lfeng Wed Jan 20 SGT 2016         %
% changed halfspace thicknness from Inf to 0 lfeng Fri Jan 22 17:24:03 SGT 2016              %
% last modified by Lujia Feng Wed Jan 27 12:58:24 SGT 2016                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate parameters for layers
% common /model0/ z1,z2,ro1,ro2,vp1,vp2,vs1,vs2,l0
l0 = 1;
% upper bound for the 1st layer
z1(l0)  = layer(1,2);
vp1(l0) = layer(1,3);
vs1(l0) = layer(1,4);
ro1(l0) = layer(1,5);
for ii=2:size(layer,1)
   if layer(ii,2)>layer(ii-1,2)
      % lower bound
      z2(l0)  = layer(ii,2);
      vp2(l0) = layer(ii,3);
      vs2(l0) = layer(ii,4);
      ro2(l0) = layer(ii,5);
      l0 = l0 + 1;
   else
      % upper bound
      z1(l0)  = layer(ii,2);
      vp1(l0) = layer(ii,3);
      vs1(l0) = layer(ii,4);
      ro1(l0) = layer(ii,5);
   end
end

% resolution for (1) p-wave velocity, (2) s-wave velocity, and (3) density
resolut = 5e-2*ones(3,1);
% loop through layers except halfspace
nl = 0;
for ii=1:l0-1
   dz  = z2(ii) - z1(ii);
   % calculate the number of layers according to resolution
   dvp = 2.0*abs(vp2(ii)-vp1(ii))/(vp2(ii)+vp1(ii));
   dvs = 2.0*abs(vs2(ii)-vs1(ii))/(vs2(ii)+vs1(ii));
   dro = 2.0*abs(ro2(ii)-ro1(ii))/(ro2(ii)+ro1(ii));
   ll  = round(max([dvp/resolut(1) dvs/resolut(2) dro/resolut(3)])); % idnint in fortran = round in matlab
   ll  = max(1,ll); % in case ll = 0
   % calcuate the change of vp, vs, ro per meter
   dvp = (vp2(ii)-vp1(ii))/dz;
   dvs = (vs2(ii)-vs1(ii))/dz;
   dro = (ro2(ii)-ro1(ii))/dz;
   % thickness of each layer
   dh  = dz/ll;
   % loop through internal layers
   for jj=1:ll
      nl     = nl + 1;
      hh(nl,1) = dh; 
      zz     = (jj-0.5)*dh;
      topz(nl,1) = z1(ii)  + (jj-1)*dh;
      botz(nl,1) = z1(ii)  + jj*dh;
      % use properties at the middle depth of each layer to represent the layer
      ro(nl,1)   = ro1(ii) + dro*zz;
      vp(nl,1)   = vp1(ii) + dvp*zz;
      vs(nl,1)   = vs1(ii) + dvs*zz;
   end
end

% last layer is halfspace
nl = nl + 1;
topz(nl,1) = z1(l0);
hh(nl,1)   = 0; % halfspace thickness fortran code sets it 0.0
botz(nl,1) = Inf; % halfspace bottom depth is positive infinite
vp(nl,1)   = vp1(l0);
vs(nl,1)   = vs1(l0);
ro(nl,1)   = ro1(l0);

% integer n0: number of homogeneous layers
% common /model/ h,ro,vp,vs,n0
sublayer.nl   = nl;
sublayer.topz = topz;
sublayer.botz = botz;
sublayer.hh   = hh;
sublayer.vp   = vp;
sublayer.vs   = vs;
sublayer.ro   = ro;
