function [ AAi ] = GTdef_edgmatinv(kk,hh,vp,vs,ro)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 GTdef_edgmatinv.m                                     %
%                                                                                       %
% Calculate inverse layer matrix based on fortran code EDGRN edgmatinv.F                %
% The matrix is a modification of L(z) in equation B.1 of Want et al. 2003              %
% The four columns of the poloidal layer matrix is the four eigenvectors of             %
% the coefficient A in equation 16                                                      %
%                                                                                       %
% INPUT										        %
% (1) kk - wave number                                                                  %
% (2) hh - layer thickness [m]                                                          %
% (3) vp - P-wave velocity [m/s]                                                        %
% (4) vs - P-wave velocity [m/s]                                                        %
% (5) ro - density [kg/m^3]                                                             %
%										        %
% OUTPUT                                                                                %
% AAi    - AAi matrix [4x4]                                                             %
%                                                                                       %
% REFERENCE  									        %
% Wang, R., Martin, F. L., & Roth, F. (2003)					        %
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5	        %
%		                                                                        %
% first created by Lujia Feng Thu Dec 11 04:52:53 SGT 2014                              %
% last modified by Lujia Feng Fri Jan 22 10:15:55 SGT 2016                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx       =  kk*hh;
mu       =  ro*vs^2;
la       =  ro*vp^2-2.0*mu;
eta      =  mu/(la+mu);
alfa     =  2.0*mu*kk;
psi      =  1.0+eta;
AAi(1,1) =  xx-eta;
AAi(1,2) =  xx/alfa;
AAi(1,3) = -xx+psi;
AAi(1,4) =  (psi+eta-xx)/alfa;
AAi(2,1) =  xx+eta;
AAi(2,2) = -xx/alfa;
AAi(2,3) =  xx+psi;
AAi(2,4) = -(psi+eta+xx)/alfa;
AAi(3,1) =  1.0;
AAi(3,2) =  1.0/alfa;
AAi(3,3) = -1.0;
AAi(3,4) = -1.0/alfa;
AAi(4,1) =  1.0;
AAi(4,2) = -1.0/alfa;
AAi(4,3) =  1.0;
AAi(4,4) = -1.0/alfa;

AAi=0.5.*AAi./psi;
