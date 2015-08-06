function [ ] = GTdef_efm2GTdef(filename,fltname,flttype)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              GTdef_efm2GTdef 			               %
%									       %
% Convert output from Emma's code to input for GTdef.m                         %
%                                                                              %
% INPUT:                                                                       %
% (1) filename - output from Emma's code                                       %
% 1      2      3            4          5         6    7      8   9            %
% MidLon MidLat MidDepth(km) Length(km) Width(km) rake strike dip slip(cm)     %
% Depth positive downward                                                      %
% First along strike, then along dip                                           %
% strike should be all the same, but rake can be different!                    %
%		                                                               %
% (2) fltname - fault name in GTdef input                                      %
% all units in m                                                               %
%                                                                              %
% (3) flttype - fault type [ 1 or 3 ]                                          %
%                                                                              %
% OUTPUT:                                                                      %
% GTdef fault type-1 or fault type-3                                           %
%									       %
% first created by lfeng Tue May 22 14:28:50 SGT 2012                          %
% last modified by lfeng Wed May 23 00:01:34 SGT 2012                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read in 
if ~exist(filename,'file'), error('GTdef_efm2GTdef ERROR: %s does not exist!',filename); end

fin = fopen(filename,'r');
efm_cell = textscan(fin,'%f %f %f %f %f %f %f %f %f','CommentStyle','#');       % Delimiter default white space
efm = cell2mat(efm_cell);
fclose(fin);

% Note: length unit is converted from km to m; slip unit is converted from cm to m
lon  = efm(:,1); lat = efm(:,2); dep = 1e3*efm(:,3); len = 1e3*efm(:,4); width = 1e3*efm(:,5);
rake = efm(:,6); str = efm(:,7); dip = efm(:,8); slip = 1e-2*efm(:,9); 
Ns   = sum(dep==dep(1));
depM = reshape(dep,Ns,[]);
Nd   = size(depM,2);

% strike and dip all the same
Tstr  = str(1);
Trot  = Tstr-360;
Tdip  = dip(1);
if ~all(str==Tstr), error('GTdef_efm2GTdef ERROR: all strike should be the same!'); end 
if ~all(dip==Tdip), error('GTdef_efm2GTdef ERROR: all dip should be the same!'); end 

% get start point
rake1 = rake(1);
x1    = -dep(1)/tand(Tdip);
y1    = -0.5*len(1);
[ lon0,lat0 ] = ckm2LLd(x1,y1,lon(1),lat(1),Trot);

% get total fault length
lenM = reshape(len,Ns,Nd);
Tlen = sum(lenM(:,1));

% get bottom & top depth for each layer
depM = reshape(dep,Ns,Nd); depList = depM(1,:)';
widM = reshape(width,Ns,Nd); widList = widM(1,:)';
z2 = depList + 0.5*widList.*sind(Tdip);
z1 = [ 0; z2(1:end-1) ];

% get rake and slip
slipM = reshape(slip,Ns,Nd);
rakeM = reshape(rake,Ns,Nd);

% write out
[ ~,basename,~ ] = fileparts(filename);
% fault type 1
if flttype==1
    foutname = [ basename '_flt1.in' ];
    fout = fopen(foutname,'w');
    %fault type name lon lat z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns
    fprintf(fout,'#fault type name lon lat z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns\n');
    fprintf(fout,'fault  1 %s  %-14.8f %-12.8f %14.8e %-12.8e %12.8e %-12.8e %-5.2f  %10.5f %-8.5f %-8.5f  %-5.4f %-5.4f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d %d\n',fltname,lon0,lat0,z1(1),z2(end),Tlen,Tstr,Tdip,0,0,0,0,0,0,0,0,0,Nd,Ns);
    
    % subfault name dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX
    fprintf(fout,'\n#subfault name dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX\n');
    for jj=1:Nd
        for ii=1:Ns
            slp = slipM(ii,jj); rk  = rakeM(ii,jj);
            %if slp==0.0, continue; end
	    ss = slp*cosd(rk); ds = slp*sind(rk);
            fprintf(fout,'     subfault %s  %5d %5d  %12.8f  %12.8f 0.0   0.0 0.0   0.0 0.0   0.0 0.0\n',fltname,jj,ii,ss,ds);
        end
    end
end

% fault type 3
if flttype==3
    foutname = [ basename '_flt3.in' ];
    fout = fopen(foutname,'w');
    %fault type name lon lat z1 z2 len str dip rake rs ts rake0 rakeX rs0  rsX ts0 tsX Nd Ns
    fprintf(fout,'#fault type name lon lat z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns\n');
    fprintf(fout,'fault  3 %s  %-14.8f %-12.8f %14.8e %-12.8e %12.8e %-12.8e %-5.2f  %10.5f %-8.5f %-8.5f  %-5.4f %-5.4f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d %d\n',fltname,lon0,lat0,z1(1),z2(end),Tlen,Tstr,Tdip,rake1,0,0,rake1,rake1,0,0,0,0,Nd,Ns);
    
    %subfault name dnum snum rake rs ts rake0 rakeX rs0 rsX ts0 tsX
    fprintf(fout,'\n#subfault name dnum snum rake rs ts rake0 rakeX rs0 rsX ts0 tsX\n');
    for jj=1:Nd
        for ii=1:Ns
            slp = slipM(ii,jj);
            %if slp==0.0, continue; end
    	rk  = rakeM(ii,jj);
            fprintf(fout,'     subfault %s  %5d %5d  %10.2f %12.8f 0.0  %-7.2f %-7.2f  0.0 0.0  0.0  0.0\n',fltname,jj,ii,rk,slp,rk,rk);
        end
    end
end

%dip fault name  dip   z1    z2    rows
fprintf(fout,'\n#dip fault name  dip   z1    z2    rows\n');
dipnum = length(z1);
for ii =1:dipnum
    fprintf(fout,'dip   %s  %8.4f  %-12.4e  %-12.4e  %-d\n',fltname,Tdip,z1(ii),z2(ii),1);
end
fclose(fout);
