function [ yy ] = GTdef_edgkern(psflag,shflag,kk,eps,sublayer,depths,source) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  GTdef_edgkern.m                                      %
%                                                                                       %
% function to call GTdef_edgpsv.m & GTdef_edgsh.m                                       %
%                                                                                       %
% INPUT										        %
% psflag - is p-sv (poloidal) mode important?                                           %
% shflag - is sh (toroidal) mode important?                                             %
% kk     - wavenumber                                                                   %
% eps    - relative accuracy                                                            %
% sublayer structure corresponds to /model/  h,ro,vp,vs,n0 [h->hh,n0->nl]               %
% depths structure corresponds to /sublayer/ hp,lp,nno     [hp->hh,lp->nl]              %
%                                 /receiver/ zrec,lzrec    [zrec->recz,lzrec->lrec]     %
%                                  /source/  zs,ls         [zs->srcz,ls->lsrc]          %
% source structure corresponds to /source/ r0,ms,ics,sfct,kpower [ics->cs]              %
%										        %
% OUTPUT                                                                                %
% yy      - solution vector for P-SV waves                                              %
%         = (U,E,V,F,W,G)^T                                                             %
%   U,V,W - two displacement components in the frequency-wavenumber domain              %
%   P,S,G - two components of elastic surface force on a horizontal plane               %
% see Equation 15 & 20 of Wang et al. 2003                                              %
%										        %
% REFERENCE  									        %
% Wang, R., Martin, F. L., & Roth, F. (2003)					        %
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5.           %
%		                                                                        %
% first created by Lujia Feng Fri Sep  4 14:58:12 SGT 2015                              %
% last modified by Lujia Feng Wed Jan 20 18:40:04 SGT 2016                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yy = zeros(6,1);

if (psflag)
   [ yy ] = GTdef_edgpsv(yy,kk,eps,sublayer,depths,source);
end

if (shflag)
   [ yy ] = GTdef_edgsh(yy,kk,eps,sublayer,depths,source);
end
