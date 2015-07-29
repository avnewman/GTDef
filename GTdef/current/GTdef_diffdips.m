function [ x1,y1,x2,y2,z1,z2,dip,ddip ] = GTdef_diffdips(mx1,my1,mx2,my2,mz1,mz2,mstr,mdip,Nd,Ns,dipin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          GTdef_diffdips.m				  %
% Calculate dip,x1,y1,x2,y2,z1,z2 of subfaults for bended faults	  %
%								          %
% INPUT:					  		  	  %
%    master fault info							  %
%    dipin - [ dip z1 z2 rows ]						  %
%            additional dip info					  %
%									  %
% OUTPUT: all column vectors					          %
%    x1,y1 - one endpoint among the two endpoints of the master fault     %
%    x2,y2 - the other endpoint among the two endpoints 		  %
%            both in the local cartesian coordinate system	          %
%    z1  - vertical burial depth (top of fault)                           %  
%    z2  - vertical locking depth (bottom of fault)                       %
%    dip - down from Horiz, right looking from the endpoint 1 [0 180]     %
%									  %
% first created by Lujia Feng Mon Dec  7 07:12:25 EST 2009		  %
% last modified by Lujia Feng Wed Dec  1 18:55:32 EST 2010		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zvalue = [ mz1;dipin(:,2);dipin(:,3);mz2 ];		% find all zvalues
num = size(zvalue,1);
for ii = 1:num
    for jj = ii+1:num
        if zvalue(ii)==zvalue(jj);
            zvalue(ii) = 9999; continue
        end
    end
end
zvalue = zvalue(find(zvalue~=9999));			% find non-zero values
zvalue = sort(zvalue);					% ascending sorting
num = size(zvalue,1)-1;					% this is the number of layers with the same dip
numin = size(dipin,1);
rows = zeros(num,1); dipval = mdip*ones(num,1);
for ii = 1:num
    for jj = 1:numin
        if zvalue(ii)==dipin(jj,2)
        	rows(ii) = dipin(jj,4);	dipval(ii) = dipin(jj,1); continue
        end
    end
end
rowsum = sum(rows);  
if rowsum < Nd
    rowzero = find(rows==0) 				% find layers not specified with number of rows yet
    rows(rowzero) = round((Nd-rowsum)/size(rowzero,1));
    rowsum = sum(rows)  
    if rowsum ~= Nd, error('Dip is not specified correctly!'); end
end
coldip = []; zlin = []; ddipsum = 0;
% xx and yy also change
xlin = linspace(mx1,mx2,Ns+1); ylin = linspace(my1,my2,Ns+1); 
rownum = rows(1);
x1mat = xlin(ones(rownum,1),1:end-1); y1mat = ylin(ones(rownum,1),1:end-1);	% endpoint 1
x2mat = xlin(ones(rownum,1),2:end);   y2mat = ylin(ones(rownum,1),2:end);	% endpoint 2
oldx1 = mx1; oldx2 = mx2; oldy1 = my1; oldy2 = my2;
% analyze layer by layer ( one layer means pateches with the same dip )
for ii = 1:num
    rownum = rows(ii); lz1 = zvalue(ii); lz2 = zvalue(ii+1); ldip2 = dipval(ii);
    ddipsum = ddipsum + (lz2-lz1)/rownum/sind(ldip2);
    coldip = [ coldip;ldip2*ones(rownum,1) ];
    zlin0 = linspace(lz1,lz2,rownum+1);
	zlin = [ zlin zlin0(1:end-1)];
    if ii~=1
         ldip1 = dipval(ii-1);
         newx1 = oldx1+(lz1/tand(ldip1)-lz1/tand(ldip2)).*cosd(mstr); oldx1 = newx1;
         newy1 = oldy1-(lz1/tand(ldip1)-lz1/tand(ldip2)).*sind(mstr); oldy1 = newy1;
         newx2 = oldx2+(lz1/tand(ldip1)-lz1/tand(ldip2)).*cosd(mstr); oldx2 = newx2;
         newy2 = oldy2-(lz1/tand(ldip1)-lz1/tand(ldip2)).*sind(mstr); oldy2 = newy2;
	     xlin = linspace(newx1,newx2,Ns+1); ylin = linspace(newy1,newy2,Ns+1); 
	     x1mat = [ x1mat; xlin(ones(rownum,1),1:end-1) ]; y1mat = [ y1mat; ylin(ones(rownum,1),1:end-1)];
	     x2mat = [ x2mat; xlin(ones(rownum,1),2:end) ];   y2mat = [ y2mat; ylin(ones(rownum,1),2:end)];
    end
end
ddip = ddipsum/num;
dip = coldip(:,ones(1,Ns)); dip = reshape(dip,[],1);
zlin = [ zlin mz2 ]';
z1mat = zlin(1:end-1,ones(1,Ns)); z2mat = zlin(2:end,ones(1,Ns));	% depths
x1 = reshape(x1mat,[],1);  y1 = reshape(y1mat,[],1);
x2 = reshape(x2mat,[],1);  y2 = reshape(y2mat,[],1);
z1 = reshape(z1mat,[],1);  z2 = reshape(z2mat,[],1);
