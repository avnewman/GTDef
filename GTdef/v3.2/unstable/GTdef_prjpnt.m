function [ ] = GTdef_prjpnt(fout,pnt,fltType,fltName,flt,strin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              GTdef_prjpnt				        %
% Project point location onto one fault system		                        %
% only works for master fault of fault1, fault2, fault3, and fault4	        %
%									        %
% INPUT:					  		  	        %
% fout    - output file handle						        %
% fltType - fault type [1234]                                                   %
% fltName - fault name                                                          %
% flt  - fault parameters                                                       %
% flt1 = [xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]       %
% flt2 = [x1 y1 x2 y2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]         %
% flt3 = [xx yy z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns] %
% flt4 = [x1 y1 x2 y2 z1 z2 dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns]   %
%    note: input coordinate is cartesian                                        %
% strin - strike addon info for the master fault                                %
%        = [ x1 y1 x2 y2 columns sweepAngle ]                                   %
%                                                                               %
% OUTPUT: (output to a file, the format is)				        %
% [ point 3 pnt.name lon lat z Ue Un Uv eUe eUn eUv weight 		        %
%   fltName Dstr1 Dstr2 Ddip Dvert ]					        %
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
% changed to dealing with one fault system lfeng Fri Oct 24 19:19:15 SGT 2014   %
% changed point output lfeng Fri Oct 24 19:28:41 SGT 2014                       %
% changed name to have a length of 20 letters lfeng Tue Jun 14 13:05:14 SGT 2016%
% last modified by Lujia Feng Tue Jun 14 13:05:20 SGT 2016                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6, strin = []; end

dataNum = size(flt,2);
if dataNum~=18  
    error('GTdef_prjpnt ERROR: need a n*18 fault vector as input!'); 
end

x1 = flt(1); y1 = flt(2); 
dip = flt(7);

pntOut = [ pnt.loc pnt.disp pnt.err pnt.wgt ];

if isempty(strin)
   if fltType==1 || fltType==3
      z1  = flt(3); 
      z2  = flt(4); 
      len = flt(5); 
      str = flt(6); 
      x2  = x1+len.*sind(str); 
      y2  = y1+len.*cosd(str);     % endpoint 2
   else
      x2  = flt(3); 
      y2  = flt(4); 
      z1  = flt(5); 
      z2  = flt(6); 
      len = sqrt((x1-x2).^2+(y1-y2).^2);
      str = GTdef_strike(x1,y1,x2,y2);
   end
   for ii = 1:pnt.num
       xx      = pnt.crt(ii,1);
       yy      = pnt.crt(ii,2); 
       pstr    = GTdef_strike(x1,y1,xx,yy);
       dist    = sqrt((x1-xx).^2+(y1-yy).^2);
       strdiff = pstr-str; 
       Ddip    = dist.*sind(strdiff);                      % Ddip can be negative if it's on the other side of fault plane
       Dstr1   = dist.*cosd(strdiff);                      % Dstr can be negative if it's off the fault plane
       Dstr2   = len-Dstr1;
       Dvert   = Ddip.*tand(dip);
       cpnt    = pntOut(ii,:);
       fprintf(fout,'point 3 %-30s %12.6f %12.6f %12.4e %10.5f %10.5f %10.5f %8.5f %8.5f %8.5f %5.2f %-10s %12.4e %12.4e %12.4e %12.4e\n',...
               pnt.name{ii},cpnt,fltName,Dstr1,Dstr2,Ddip,Dvert);
   end
else
   sx1 = strin(:,1); sy1 = strin(:,2);
   sx2 = strin(:,3); sy2 = strin(:,4);
   len = sqrt((sx1-sx2).^2+(sy1-sy2).^2);
   for ii = 1:pnt.num
       xx = pnt.crt(ii,1);
       yy = pnt.crt(ii,2); 
       dd = sqrt((sx1-xx).^2+(sy1-yy).^2);
       [ ddsort,ind ] = sort(dd);                          % find the closest segment
       x1closest  = sx1(ind(1)); y1closest = sy1(ind(1));
       x2closest  = sx2(ind(1)); y2closest = sy2(ind(1));
       lenclosest = sqrt((x1closest-x2closest).^2+(y1closest-y2closest).^2); 
       strclosest = GTdef_strike(x1closest,y1closest,x2closest,y2closest);
       Dstr1_0    = sum(len(1:ind(1)-1));
       Dstr2_0    = sum(len(ind(1)+1:end));

       pstr       = GTdef_strike(x1closest,y1closest,xx,yy);
       dist       = sqrt((x1closest-xx).^2+(y1closest-yy).^2);
       strdiff    = pstr-strclosest; 
       Ddip       = dist.*sind(strdiff);                   % Ddip can be negative if it's on the other side of fault plane
       Dstr1      = Dstr1_0+dist.*cosd(strdiff);           % Dstr can be negative if it's off the fault plane
       Dstr2      = Dstr2_0+lenclosest-Dstr1;
       Dvert      = Ddip.*tand(dip);
       cpnt       = pntOut(ii,:);
       fprintf(fout,'point 3 %-30s %12.6f %12.6f %12.4e %10.5f %10.5f %10.5f %8.5f %8.5f %8.5f %5.2f %-10s %12.4e %12.4e %12.4e %12.4e\n',...
               pnt.name{ii},cpnt,fltName,Dstr1,Dstr2,Ddip,Dvert);
   end
end
