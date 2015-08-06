function [ xsect_name,xsect ] = GTdef_xsection(flt_type,flt_name,flt,dipin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               GTdef_xsection				        %
% Calculate cross section of fault interface				        %
%									        %
% INPUT:					  		  	        %
% flt1 = [xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]       %
% flt2 = [x1 y1 x2 y2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]         %
% flt3 = [xx yy z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns] %
% flt4 = [x1 y1 x2 y2 z1 z2 dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns]   %
% input coordinate is cartesian					                %
% dipin = [ dip z1 z2 rows ]						        %
%                                                                               %
% OUTPUT: 								        %
% [ index dip x1 z1 x2 z2 width rows ]                            	        %	
%  xsect_name - master fault name [cell]				        %
%  index    - position of each segment counting from surface  	  	        %
%  dip      - dip angle [degree] (consider positive only)		        %
%  x1,x2    - distance from fault vertical trace [m]                            %
%  z1,z2    - shallowest and deepest depth of each segment [m]                  %
%  width    - surface apparent width of each segment [m]                        %
%  rows     - num of patches on each segment                                    %
%                                                                               %
% first created by Lujia Feng Thu Dec  2 02:22:28 EST 2010		        %
% corrected a bug when dip is not specified lfeng Tue Jul  2 20:07:37 SGT 2013  %
% last modified by Lujia Feng Tue Jul  2 20:10:52 SGT 2013                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datanum = size(flt,2);
if datanum~=18  
   error('GTdef_xsection ERROR: need a n*18 fault vector as input!'); 
end

if flt_type==1 || flt_type==3
   z1 = flt(:,3); z2 = flt(:,4); dip = flt(:,7);
elseif flt_type==2 || flt_type==4
   z1 = flt(:,5); z2 = flt(:,6); dip = flt(:,7);
else
   error('GTdef_xsection ERROR: fault type is wrong!');
end

if isempty(dipin)
   xsect_name = flt_name;					% flt_name - cell array
   x1 = zeros(size(dip)); num = ones(size(dip));
   width = bsxfun(@rdivide,z2-z1,tand(dip));
   x2 = width;
   xsect = [ 1 dip x1 z1 x2 z2 width num ];
else
   dipnum = size(dipin,1);
   xsect_name = cell(dipnum,1);
   for ii = 1:dipnum, xsect_name{ii} = flt_name; end		% flt_name - strings
   xsect = zeros(dipnum,8);
   x2 = 0;
   dipin = sortrows(dipin,2);		                        % order from surface down
   for ii=1:dipnum
      dip = dipin(ii,1); z1 = dipin(ii,2); z2 = dipin(ii,3); rows = dipin(ii,4); 
      width = (z2-z1)/tand(dip);
      x1 = x2; x2 = x1+width;
      xsect(ii,:) = [ ii dip x1 z1 x2 z2 width rows ];                     
   end
end
