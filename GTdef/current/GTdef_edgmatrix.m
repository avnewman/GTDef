function [ AA ] = GTdef_edgmatrix(kk,zz,vp,vs,ro)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 GTdef_edgmatrix.m                                     %
%                                                                                       %
% Calculate A matrix based on EDGRN edgmatrix.F                                         %
%                                                                                       %
% INPUT										        %
% (1) kk - wave number                                                                  %
% (2) zz - receiver depth [m]                                                           %
% (3) vp - P-wave velocity [m/s]                                                        %
% (4) vs - P-wave velocity [m/s]                                                        %
% (5) ro - density [kg/m^3]                                                             %
%										        %
% OUTPUT                                                                                %
% AA     - AA matrix [4x4]                                                              %
%                                                                                       %
% REFERENCE  									        %
% Wang, R., Martin, F. L., & Roth, F. (2003)					        %
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5	        %
%		                                                                        %
% first created by Lujia Feng Thu Dec 11 04:04:19 SGT 2014                              %
% last modified by Lujia Feng Thu Dec 11 04:10:10 SGT 2014                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx      =  kk*zz;
mu      =  ro*vs^2;
la      =  ro*vp^2 - 2.0*mu;
eta     =  mu/(la+mu);
alfa    =  2.0*mu*kk;
AA(1,1) =  1.0;
AA(1,2) = -1.0;
AA(1,3) = -xx + 1.0 + 2.0*eta;
AA(1,4) =  xx + 1.0 + 2.0*eta;
AA(2,1) =  alfa;
AA(2,2) =  alfa;
AA(2,3) =  alfa*(-x+1.0+eta);
AA(2,4) = -alfa*(x+1.0+eta);
AA(3,1) =  1.0;
AA(3,2) =  1.0;
AA(3,3) = -xx;
AA(3,4) = -xx;
AA(4,1) =  alfa;
AA(4,2) = -alfa;
AA(4,3) =  alfa*(-xx+eta);
AA(4,4) =  alfa*(xx+eta);
