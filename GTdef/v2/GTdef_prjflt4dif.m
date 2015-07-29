function [ newflt,prjflt,xyzctr,xyztop1 ] = GTdef_prjflt4dif(flt,subflt,dipin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_prjflt4dif				        %
%									        %
% Project distributed-slip fault4 geometry and slip information                 %
% onto surface geographic coordinate				                %
%									        %
% INPUT:					  		  	        %
%  flt = [ x1 y1 x2 y2 z1 z2 dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns ] %
%  x1,y1 - one endpoint among the two endpoints of the master fault             %
%  x2,y2 - the other endpoint among the two endpoints 		                %
%          both in the local cartesian coordinate system	                %
%    z1  - vertical burial depth (top of fault)                                 %  
%    z2  - vertical locking depth (bottom of fault)                             %
%    dip - down from Horiz, right looking from the endpoint 1 [0 180]           %
%   rake - Aki-Richards convention                                              %
%    rs  - rake-slip (rake direction +)                                         %
%    ts  - tensile-slip (opening +)                                             %
%  rake0,rakeX - rake is usually fixed, currently dummy parameters              %
%      rs0,ts0 - lower bounds for slips				                %
%      rsX,tsX - upper bounds for slips				                %
%    Nd  - number of rows defining the subfaults along dip 	                %
%    Ns  - number of rows defining the subfaults along strike 		        %
% subflt = [ dnum snum rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]                 %
%   dnum - row number for subfaults					        %
%   snum - column number for subfaults				  	        %
%   rake - rake direction                                                       %
%  rs,ts - subfault slips						        %
%      rake0,rakeX - rake is usually fixed, currently dummy parameters          %
%  rs0,ts0,rsX,tsX - subfault slip bounds			                %
%                                                                               %
% OUTPUT:                                                                       %
% prjflt = [ dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2          %
%            xctr  yctr rake rs ts ]  					        %
%   [xtop1 ytop1], [xbot1 ybot1], [xbot2 ybot2], and [xtop2 ytop2]	        %
%   are the surface projection of four points that confine 		        %
%   the fault interface 						        %
%   They are in a counterclockwise sense looking from the RHS of endpoint       %
%   [xtop1 ytop1] and [xbot1 ybot1] correspond [xx yy] at z1 and z2	        %
%   [xx yy] is surface projection of one fault endpoint 	  	        %
% xyzctr  = [ xx;yy;zz ] (3*nn) for center points                               %
% xyztop1 = [ xx;yy;zz ] (3*nn) for top upper corners                           %
%   zz elevation positive upward                                                %
%									        %
% first created by Lujia Feng Mon May 14 01:34:13 SGT 2012                      %
% output xyzctr & xyztop1 lfeng Tue Jun 12 18:19:51 SGT 2012                    %
% output newflt lfeng Wed Jun 13 15:53:43 SGT 2012                              %
% last modified by Lujia Feng Wed Jun 13 15:59:03 SGT 2012                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 18], error('GTdef_prjflt4dif ERROR: need a 1*18 fault vector as input!');  end

% initialization
% Note: flt is a 1-row vector for the master fault
mx1 = flt(1); my1 = flt(2); mx2 = flt(3); my2 = flt(4); 
mz1 = flt(5); mz2 = flt(6); mdip = flt(7);
mslips = flt(8:16);				% slip block
Nd = round(flt(17)); Ns = round(flt(18));	% number of rows and columns
subflt_num = Nd*Ns;				% subfault num
mstr = GTdef_strike(mx1,my1,mx2,my2);
unit = ones(subflt_num,1);

% exclusively specify dips for rows between z1 and z2 depth range
if ~isempty(dipin)
    [ x1,y1,x2,y2,z1,z2,dip,ddip ] = GTdef_diffdips(mx1,my1,mx2,my2,mz1,mz2,mstr,mdip,Nd,Ns,dipin);
else
    dip = mdip*unit; 			% duplicate dips
    zlin = linspace(mz1,mz2,Nd+1)';
    z1mat = zlin(1:end-1,ones(1,Ns)); z2mat = zlin(2:end,ones(1,Ns));	% depths
    xlin = linspace(mx1,mx2,Ns+1); ylin = linspace(my1,my2,Ns+1); 
    x1mat = xlin(ones(Nd,1),1:end-1); y1mat = ylin(ones(Nd,1),1:end-1);	% endpoint 1
    x2mat = xlin(ones(Nd,1),2:end);   y2mat = ylin(ones(Nd,1),2:end);	% endpoint 2
    x1 = reshape(x1mat,[],1);  y1 = reshape(y1mat,[],1);
    x2 = reshape(x2mat,[],1);  y2 = reshape(y2mat,[],1);
    z1 = reshape(z1mat,[],1);  z2 = reshape(z2mat,[],1);
end

% form subfaults
slips = mslips(unit,:);				% duplicate slips of master-fault

%  subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]
if ~isempty(subflt)
    num = size(subflt,1); mat = [Nd Ns];
    for ii = 1:num
        dnum = round(subflt(ii,1)); snum = round(subflt(ii,2));
        jj = sub2ind(mat,dnum,snum);
	slips(jj,:) = subflt(ii,3:11);
    end
end

dlin = round(linspace(1,Nd,Nd)'); slin = round(linspace(1,Ns,Ns));
dmat = dlin(1:end,ones(1,Ns)); smat = slin(ones(Nd,1),1:end);
dnum = reshape(dmat,[],1);  snum = reshape(smat,[],1);
newflt = [ dnum snum x1 y1 x2 y2 z1 z2 dip slips ];

[ prjflt,xyzctr,xyztop1 ] = GTdef_prjflt4uni(newflt);
