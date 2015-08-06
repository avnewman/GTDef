function [ hk ] = GTdef_edghask(kk,zz,vp,vs,ro,eps0,mm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 GTdef_edghask.m                                       %
%                                                                                       %
% Calculate Hask matrix based on EDGRN edghask.F                                        %
%                                                                                       %
% INPUT										        %
% (1) kk   - wave number                                                                %
% (2) zz   - layer depth [m]                                                            %
% (3) vp   - P-wave velocity [m/s]                                                      %
% (4) vs   - P-wave velocity [m/s]                                                      %
% (5) ro   - density [kg/m^3]                                                           %
% (6) eps0 - relative accuracy  eps0 = 1.0e-3                                           %
% (7) mm   - size of output hask matrix [ 2 or 4 ]                                      %
%										        %
% OUTPUT                                                                                %
% hk     - Hask matrix [2x2] or [4x4]                                                   %
%                                                                                       %
% REFERENCE  									        %
% Wang, R., Martin, F. L., & Roth, F. (2003)					        %
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5	        %
%		                                                                        %
% first created by Lujia Feng Thu Dec 11 04:42:13 SGT 2014                              %
% last modified by Lujia Feng Thu Dec 11 04:45:19 SGT 2014                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k2  = kk*kk;
xx  = kk*z;
x2  = xx*xx;
mu  = ro*vs^2;
la  = ro*vp^2-2.0*mu;
eta = mu/(la+mu);
ex  = exp(xx);

%	ch=(e^x+1/e^x)/2
%	sh=(e^x-1/e^x)/2x
ch = 0.5*(ex+1.0/ex);
if(abs(xx)>=eps0)
   sh = 1.0 + x2/6.0*(1.0+x2/20.0);
else
   sh = 0.5*(ex-1.0/ex)/xx;
end

% propagator matrix for SH waves
if mm==2
   hk(1,1) = ch;
   hk(1,2) = zz*sh/mu;
   hk(2,1) = kk*xx*sh*mu;
   hk(2,2) = ch;
% propagatior matrix for P-SV waves.
elseif mm==4
   hk(1,1) = ch-x2*sh/(1.0+eta);
   hk(1,2) = 0.50*zz*(-ch+(1.0+2.0*eta)*sh)/(1.0+eta)/mu;
   hk(1,3) = xx*(ch-eta*sh)/(1.0+eta);
   hk(1,4) = x2*sh/(1.0+eta);
   hk(2,1) = 2.0*mu*kk*xx*(-ch+sh)/(1.0+eta);
   hk(2,2) = hk(1,1);
   hk(2,3) = 2.0*mu*kk*hk(1,4);
   hk(2,4) = xx*(ch+eta*sh)/(1.0+eta);
   hk(3,1) = -hk(2,4);
   hk(3,2) = -0.50*xx*zz*sh/(1.0+eta)/mu;
   hk(3,3) = ch+x2*sh/(1.0+eta);
   hk(3,4) = 0.50*zz*(ch+(1.0+2.0*eta)*sh)/(1.0+eta)/mu;
   hk(4,1) = -hk(2,3);
   hk(4,2) = -hk(1,3);
   hk(4,3) = 2.0*mu*kk*xx*(ch+sh)/(1.0+eta);
   hk(4,4) = hk(3,3);
end
