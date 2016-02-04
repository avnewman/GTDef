function [ depths ] = GTdef_edgdepth(topz,botz,srcz,recz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  GTdef_edgdepth.m                                     %
%                                                                                       %
% Order depths of sublayers, a receiver and a source                                    %
% This code is converted from fortran code edglayer.F                                   %
%                                                                                       %
% INPUT											%
% srcz               - depth of source                                       (scalar)   %
% recz               - depth of receiver                                     (scalar)   %
% (earth.edgrn.)topz - top depth of each sublayer [m]            (nl*1 column vector)   %
% (earth.edgrn.)botz - bottom depth of each sublayer [m]         (nl*1 column vector)   %
%                                                                                       %
% OUTPUT                                                                                %
% depths structure (changing layers including source and receiver)                      %
% corresponds to sublayer(hp,l,nno) [hp->hh,l->nl]                                      %
% depths.nl        - the number of depths including source depth and receiver depth     %
% depths.hh        - thickness of each depth                                            %
% depths.nno       - layer number for each depth                                        %
% depths.srcz      - depth of source                                                    %
% depths.recz      - depth of receiver                                                  %
% depths.lsrc      - depth number for the source                                        %
% depths.lrec      - depth number for the receiver                                      %
%                                                                                       %
% REFERENCE  										%
% Wang, R., Martin, F. L., & Roth, F. (2003)						%
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5		%
%                                                                                       %
% first created by Lujia Feng Wed Sep  9 19:34:03 SGT 2015                              %
% added depths.srcz & depths.recz lfeng Fri Jan 15 18:18:43 SGT 2016                    %
% changed halfspace thickness from Inf to 0 lfeng Fri Jan 22 17:28:37 SGT 2016          %
% last modified by Lujia Feng Fri Jan 22 18:31:17 SGT 2016                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% depths including all layer boundaries and a source and a receiver
zz = [ topz; srcz; recz ];

% sort depths in an ascending order and delete duplicates
zz = sort(unique(zz));
hh = zz(2:end) - zz(1:end-1);
hh = [ hh; 0 ]; % add thickness of halfspace 0

% determine lsrc,lrec
lsrc = find(zz==srcz);
lrec = find(zz==recz);

% determine layer no of each depth
nl  = size(zz,1);
nno = zeros(nl,1);
nno(1) = 1;
for ii=2:nl
   curz    = zz(ii);
   ind     = find(curz>=topz & curz<botz);
   nno(ii) = ind;
end

% integer lp,nno(nzmax)
% common /sublayer/ hp,lp,nno
depths.nl   = nl;
depths.hh   = hh;
depths.nno  = nno;
depths.srcz = srcz;
depths.recz = recz;
depths.lsrc = lsrc;
depths.lrec = lrec;
