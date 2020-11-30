function [ x1,y1,x2,y2,z1,z2,dip,ddip ] = GTdef_diffdips(mx1,my1,mx2,my2,mz1,mz2,mstr,dipin,Nd,Ns,sweepAngle)
%                                                        1   2   3   4   5   6   7    8     9  10 11

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          GTdef_diffdips.m				  %
% Calculate x1,y1,x2,y2,z1,z2,dip,ddip of subfaults for bended faults	  %
%								          %
% INPUT:					  		  	  %
%    master fault info							  %
%    mx1,mx2,my1,my2,mz1,mz2,mstr,Nd,Ns                                   %
%    dipin - [ dip z1 z2 rows ] additional dip info                       %
%    sweepAngle - angle between normal direction of fault trace and east  %
%                E=0 N=90 W=180 S=270 [deg]                               %
%									  %
% OUTPUT: all column vectors					          %
%    x1,y1 - one endpoint among the two endpoints                         %
%    x2,y2 - the other endpoint among the two endpoints 		  %
%            both in the local cartesian coordinate system	          %
%    z1  - vertical burial depth (top of fault)                           %
%    z2  - vertical locking depth (bottom of fault)                       %
%    dip - down from Horiz, right looking from the endpoint 1 [0 180]     %
%   ddip - average patch size along-dip                                   %
%									  %
% first created by Lujia Feng Mon Dec  7 07:12:25 EST 2009		  %
% added sweepAngle lfeng Wed Nov  5 19:44:17 SGT 2014                     %
% simplified dip input and added more format check lfeng Wed Nov 12 2014  %
% last modified by Lujia Feng Wed Nov 12 17:57:04 SGT 2014                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 10
   sweepAngle = 0;
end

dips = dipin(:,1); % dip values
z1   = dipin(:,2); 
z2   = dipin(:,3);
rows = dipin(:,4); % number of layers for each dip

dipNum = size(dips,1);  % number of dips
rowSum = sum(rows);     % total number of layers

% format check
if size(z1)>size(unique(z1)) | size(z2)>size(unique(z2))
    error('GTdef_diffdips ERROR: depths for dips can not be repeated!');
end

if z1~=sort(z1) | z2~=sort(z2) 
    error('GTdef_diffdips ERROR: dip must be specified with an increasing depth!');
end

if  size(z1)>1 & size(z2)>1 & z1(2:end)~=z1(1:end-1)
    error('GTdef_diffdips ERROR: z1 & z2 do not match!');
end

if rowSum ~= Nd, error('GTdef_diffdips ERROR: dip is not specified correctly!'); end
if mz1 ~= z1(1),   error('GTdef_diffdips ERROR: dip z1 do not match z1 for the main fault!'); end
if mz2 ~= z2(end), error('GTdef_diffdips ERROR: dip z2 do not match z2 for the main fault!'); end

dipall = []; zlin = []; ddipSum = 0;
% xx and yy also change
xlin   = linspace(mx1,mx2,Ns+1); 
ylin   = linspace(my1,my2,Ns+1); 
rownum = rows(1); 
x1mat  = xlin(ones(rownum,1),1:end-1); y1mat = ylin(ones(rownum,1),1:end-1);	% endpoint 1
x2mat  = xlin(ones(rownum,1),2:end);   y2mat = ylin(ones(rownum,1),2:end);	% endpoint 2
oldx1  = mx1; oldx2 = mx2; oldy1 = my1; oldy2 = my2;
% analyze dip by dip (one dip may have sevral patches)
for ii = 1:dipNum
    rownum  = rows(ii); 
    lz1     = z1(ii); 
    lz2     = z2(ii); 
    ldip2   = dips(ii);
    ddipSum = ddipSum + (lz2-lz1)/rownum/sind(ldip2); % sum of dipPatchSize for each dip
    dipall  = [ dipall;ldip2*ones(rownum,1) ];        % dip for each layer
    zlin0   = linspace(lz1,lz2,rownum+1);
    zlin    = [ zlin zlin0(1:end-1)];                 % depth vector
    if ii~=1
         ldip1 = dips(ii-1);
         dwidth = lz1/tand(ldip1)-lz1/tand(ldip2);
         newx1 = oldx1+dwidth.*cosd(mstr); oldx1 = newx1;
         newy1 = oldy1-dwidth.*sind(mstr); oldy1 = newy1;
         newx2 = oldx2+dwidth.*cosd(mstr); oldx2 = newx2;
         newy2 = oldy2-dwidth.*sind(mstr); oldy2 = newy2;
         xlin  = linspace(newx1,newx2,Ns+1); 
         ylin  = linspace(newy1,newy2,Ns+1); 
         x1mat = [ x1mat; xlin(ones(rownum,1),1:end-1) ]; y1mat = [ y1mat; ylin(ones(rownum,1),1:end-1)];
         x2mat = [ x2mat; xlin(ones(rownum,1),2:end) ];   y2mat = [ y2mat; ylin(ones(rownum,1),2:end)];
    end
end

ddip  = ddipSum/dipNum;
dpmat = dipall(:,ones(1,Ns));
zlin  = [ zlin mz2 ]';
z1mat = zlin(1:end-1,ones(1,Ns)); z2mat = zlin(2:end,ones(1,Ns));	% depths

% adjust endpoints according to sweepAngle except for 1st row
if sweepAngle~=0
   for ii=2:Nd
      ldip1  = dpmat(ii-1,:);
      ldip2  = dpmat(ii,:);
      % get the top left corner of each patch
      dwidth    = z1mat(ii-1,:)/tand(ldip1);
      x1TopLeft = x1mat(ii-1,:)+dwidth.*cosd(mstr);
      y1TopLeft = y1mat(ii-1,:)-dwidth.*sind(mstr);
      x2TopRigt = x2mat(ii-1,:)+dwidth.*cosd(mstr);
      y2TopRigt = y2mat(ii-1,:)-dwidth.*sind(mstr);
      % apply sweep angle
      dwidth    = (z2mat(ii-1,:)-z1mat(ii-1,:))/tand(ldip1);
      x3TopLeft = x1TopLeft+dwidth.*cosd(sweepAngle);
      y3TopLeft = y1TopLeft+dwidth.*sind(sweepAngle);
      x4TopRigt = x2TopRigt+dwidth.*cosd(sweepAngle);
      y4TopRigt = y2TopRigt+dwidth.*sind(sweepAngle);
      % get back the interception point at surface
      dwidth      = z1mat(ii,:)/tand(ldip2);
      x1mat(ii,:) = x3TopLeft-dwidth.*cosd(mstr);
      y1mat(ii,:) = y3TopLeft+dwidth.*sind(mstr);
      x2mat(ii,:) = x4TopRigt-dwidth.*cosd(mstr);
      y2mat(ii,:) = y4TopRigt+dwidth.*sind(mstr);
   end
end

% note: x1mat y1mat z1mat etc will be reshaped column-wise!
dip = reshape(dpmat,[],1);
x1  = reshape(x1mat,[],1);  y1 = reshape(y1mat,[],1);
x2  = reshape(x2mat,[],1);  y2 = reshape(y2mat,[],1);
z1  = reshape(z1mat,[],1);  z2 = reshape(z2mat,[],1);
