function [ newflt,prjflt,xyzflt ] = GTdef_prjflt6(modspace,geoname,colname,flt,subflt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_prjflt6				  %
%									  %
% Project distributed-slip fault6 geometry and slip information           %
% onto surface geographic coordinate                                      %
% Convert type-6 fault to type-3 fault                                    %
%									  %
% INPUT:					  		  	  %
% geoname - geometry file name                                            %
% colName - name of each column                                           %
% flt = [ rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns ]                  %
% subflt = [ dnum snum rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]           %
%   dnum - row number for subfaults					  %
%   snum - column number for subfaults				  	  %
%   rake - rake direction                                                 %
%  rs,ts - subfault slips						  %
%      rake0,rakeX - rake is usually fixed, currently dummy parameters    %
%  rs0,ts0,rsX,tsX - subfault slip bounds			          %
% ----------------------------------------------------------------------- %
% modspace structure                                                      %
% coord  - coordnate system                                               %
% origin = [ lon0 lat0 ]                                                  %
% ----------------------------------------------------------------------- %
%                                                                         %
% OUTPUT:                                                                 %
% newflt - all fault patches                                              %
%        = [ dnum snum x1 y1 z1 z2 len str dip slips ]                    %
% prjflt = [ dnum snum xtop1 ytop1 ztop1 xbot1 ybot1 zbot1                %
%            xbot2 ybot2 zbot2 xtop2 ytop2 ztop2 xctr yctr zctr           %
%            ss ds ts rake rs ]                                           %
%   [xtop1 ytop1], [xbot1 ybot1], [xbot2 ybot2], and [xtop2 ytop2]	  %
%   are the surface projection of four points that confine 		  %
%   the fault interface 						  %
%   They are in a counterclockwise sense looking from the RHS of endpoint %
%   [xtop1 ytop1] and [xbot1 ybot1] correspond [xx yy] at z1 and z2	  %
%   [xx yy] is surface projection of one fault endpoint 	  	  %
% ----------------------------------------------------------------------- %
% flt?.xyzflt stucture for fault patch info in cartesian coordinate       %
% xyzflt.xyzctr  - location of center points in cartesian coordinate      %
%                = [ xx yy zz ]                    (flt_num*3)            %
% xyzflt.xyztop1 - location of top upper corners in cartesian coordinate  %
%                = [ xx yy zz ]                    (flt_num*3)            %
% xyzflt.suv     - unit vector in the strike direction                    %
%                = [ xx yy zz ]                    (flt_num*3)            %
% xyzflt.duv     - unit vector in the dip direction                       %
%                = [ xx yy zz ]                    (flt_num*3)            %
% xyzflt.nuv     - normal unit vector                                     %
%                = [ xx yy zz ]                    (flt_num*3)            %
% Note: zz elevation positive upward                                      %
%                                                                         %
% first created based on GTdef_fault5.m by Lujia Feng Wed Jun 1 SGT 2016  %
% allowed rake to be not 0 in external geometry lfeng Thu Jun 2 SGT 2016  %
% last modified by Lujia Feng Thu Jun  2 16:24:27 SGT 2016                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 11], error('GTdef_prjflt6 ERROR: need a 1*11 fault vector as input!'); end

smooth = modspace.smooth;
surf   = modspace.surf;

Nd = flt(10);
Ns = flt(11);

% read geometry file
% newflt3 = [ x1 y1 z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]
if ~isnan(Nd) && ~isnan(Ns)
    [ ~,newflt3,Nd,Ns ] = GTdef_read_geometry(geoname,colname,modspace.origin,modspace.coord,[],Nd,Ns);
else
    [ ~,newflt3,Nd,Ns ] = GTdef_read_geometry(geoname,colname,modspace.origin,modspace.coord,[]);
    flt(10) = Nd;
    flt(11) = Ns;
end

% initialization
% Note: flt is a row vector for the master fault
x1    = newflt3(:,1); y1  = newflt3(:,2); 
z1    = newflt3(:,3); z2  = newflt3(:,4); 
len   = newflt3(:,5); str = newflt3(:,6);  
dip   = newflt3(:,7);
slips = newflt3(:,8:end);                      % slip block
subfltnum = Nd*Ns;                             % subfault num
unit = ones(subfltnum,1);

% only when slips from geometry file are zero (not specified), slips from GTdef input will be used
if ~any(slips(:,2:end)) % excluding rake, so allow rake to be non-zero
    % master fault
    mslips = flt(1:end-2);
    slips  = mslips(unit,:);

    % subflt = [ dnum snum rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]
    if ~isempty(subflt)
        num = size(subflt,1); mat = [Nd Ns];
        for ii = 1:num
            dnum = round(subflt(ii,1)); snum = round(subflt(ii,2));
            jj = sub2ind(mat,dnum,snum);
            slips(jj,:) = subflt(ii,3:11);
        end
    end
end

% create dnum & snum
dlin = round(linspace(1,Nd,Nd)'); slin = round(linspace(1,Ns,Ns));
dmat = dlin(1:end,ones(1,Ns));    smat = slin(ones(Nd,1),1:end);
dnum = reshape(dmat,[],1);        snum = reshape(smat,[],1);

% newflt=[dnum snum xx yy z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX]
newflt = [ dnum snum x1 y1 z1 z2 len str dip slips ];
[ prjflt,xyzflt ] = GTdef_prjflt3uni(newflt);
