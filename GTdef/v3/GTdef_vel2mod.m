function [mu,lambda,nu] = GTdef_vel2mod(vp,vs,ro)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          GTdef_vel2mod			 %
% convert p-wave and s-wave velocities to shear moduli           %
%                                                                %
% INPUT:                                                         %
% vp - p-wave velocity [m/s]                                     %
% vs - s-wave velocity [m/s]                                     %
% ro - density         [kg/m^3]                                  %
%                                                                %
% OUTPUT:                                                        %
% mu - shear modulus/rigidity                                    %
% nu - poisson's ratio                                           %
% lambda - lame's constant                                       %
%                                                                %
% first created by Lujia Feng Sat May 19 14:49:58 SGT 2012       %
% last modified by lfeng Sat May 19 15:12:08 SGT 2012            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = vs*vs*ro;
lambda = vp*vp*ro-2*mu;
nu = 0.5*lambda/(lambda+mu);
