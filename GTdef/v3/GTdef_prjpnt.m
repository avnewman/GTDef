function [ ] = GTdef_prjpnt(fout,pnt,flt_type,flt_name,flt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                GTdef_prjpnt				        %
% Project point location onto the fault surface coordinate		        %
% only works for master fault of fault1, fault2, fault3, and fault4	        %
%									        %
% INPUT:					  		  	        %
% fout - output file handle						        %
% flt1 = [xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]       %
% flt2 = [x1 y1 x2 y2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]         %
% flt3 = [xx yy z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns] %
% flt4 = [x1 y1 x2 y2 z1 z2 dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns]   %
%    input coordinate is cartesian					        %
%                                                                               %
% OUTPUT: (output to a file, the format is)				        %
% [ point 3 pnt.name lon lat z Ue Un Uv eUe eUn eUv weight 		        %
%   flt_name Dstr1 Dstr2 Ddip Dvert ]					        %
%    Dstr1 - distance along fault strike from endpoint 1		        %
%            strike direction is positive				        %
%    Dstr2 - distance along fault strike from endpoint 2		        %
%            strike direction is negative				        %
%    Ddip - distance along fault dip from the fault axis		        %
%    Dvert - distance vertically to the fault plane			        %
%                                                                               %
% first created by Lujia Feng Fri May  1 14:15:53 EDT 2009		        %
% use cell array of strings for names lfeng Wed Dec  1 15:19:07 EST 2010        %
% added Dstr2 lfeng Fri Dec 10 14:40:19 EST 2010			        %
% used structure lfeng Thu Feb 23 10:21:11 SGT 2012			        %
% last modified by Lujia Feng Thu Feb 23 10:21:20 SGT 2012		        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datanum = size(flt,2);
if datanum~=18  
    error('GTdef_prjpnt ERROR: need a n*18 fault vector as input!'); 
end

flt_num = size(flt,1);			% fault num
x1 = flt(:,1); y1 = flt(:,2); 

pnt_out = [ pnt.loc pnt.disp pnt.err pnt.wgt ];

if flt_type==1||flt_type==3
    z1 = flt(:,3); z2 = flt(:,4); 
    len = flt(:,5); str = flt(:,6); dip = flt(:,7);
    x2 = x1+len.*sind(str); y2 = y1+len.*cosd(str);	% endpoint 2
else
    x2 = flt(:,3); y2 = flt(:,4); 
    z1 = flt(:,5); z2 = flt(:,6); dip = flt(:,7);
    len = sqrt((x1-x2).^2+(y1-y2).^2);
    str = GTdef_strike(x1,y1,x2,y2);
end

for ii = 1:pnt.num
    xx = []; yy = []; pstr = []; dist = []; strdiff = [];
    Ddip = []; Dstr = []; Dvert = []; pnt_cur = [];
    xx = pnt.crt(ii,1)*ones(flt_num,1);
    yy = pnt.crt(ii,2)*ones(flt_num,1); 
    pstr = GTdef_strike(x1,y1,xx,yy);
    dist = sqrt((x1-xx).^2+(y1-yy).^2);
    strdiff = pstr-str; 
    Ddip = dist.*sind(strdiff);			% Ddip can be negative if it's on the other side of fault plane
    Dstr1 = dist.*cosd(strdiff);		% Dstr can be negative if it's off the fault plane
    Dstr2 = len-Dstr1;
    Dvert = Ddip.*tand(dip);
    pnt_cur = pnt_out(ii,:);
    for jj =1:flt_num
        fprintf(fout,'point 3 %s %-14.8f %-12.8f %-6.4e %10.5f %-10.5f %-10.5f %8.5f %-8.5f %-8.5f %-5.2f %s %12.4e %12.4e %12.4e %12.4e\n',pnt.name{ii},pnt_cur,flt_name{jj},Dstr1(jj),Dstr2(jj),Ddip(jj),Dvert(jj));
    end
end
