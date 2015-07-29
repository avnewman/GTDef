function [ xsect_name,xsect ] = GTdef_xsection(flt_type,flt_name,flt,dipin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               GTdef_xsection				  %
% Calculate cross section of fault interface				  %
%									  %
% INPUT:					  		  	  %
%  fault1 & fault3							  %
%  flt = [ xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX ]     %
%  fault2 & fault4							  %
%  flt = [ x1 y1 x2 y2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX ]       %
%  input coordinate is cartesian					  %
%  dip is only used for fault3 & fault4					  %
%  dipin = [ dip z1 z2 rows ]						  %
%                                                                         %
% OUTPUT: 								  %
% [ index dip x1 z1 x2 z2 width rows ]                            	  %	
%  xsect_name - master fault name [cell]				  %
%  index    - position of each segment counting from surface  	  	  %
%  dip      - dip angle [degree] (consider positive only)		  %
%  x1,x2    - distance from fault vertical trace [m]                      %
%  z1,z2    - shallowest and deepest depth of each segment [m]            %
%  width    - surface apparent width of each segment [m]                  %
%  rows     - num of patches on each segment                              %
%                                                                         %
% related functions:							  %
% GTdef_project.m							  %
% first created by Lujia Feng Thu Dec  2 02:22:28 EST 2010		  %
% last modified by Lujia Feng Thu Dec  2 05:04:11 EST 2010		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datanum = size(flt,2);
if datanum~=16&&datanum~=18, error('need a n*16 or n*18 fault vector for GTdef_xsection'); end

if flt_type==1||flt_type==3
   z1 = flt(:,3); z2 = flt(:,4); dip = flt(:,7);
elseif flt_type==2||flt_type==4
   z1 = flt(:,5); z2 = flt(:,6); dip = flt(:,7);
else
   error('fault type is wrong!');
end

if isempty(dipin)
   xsect_name = flt_name;					% flt_name - cell array for fault 1 & 2
   x1 = zeros(size(dip)); num = ones(size(dip));
   width = bsxfun(@rdivide,z2-z1,tand(dip));
   x2 = width;
   xsect = [ dip x1 z1 x2 z2 width num ];
elseif flt_type==3||flt_type==4
   dipnum = size(dipin,1);
   xsect_name = cell(dipnum,1);
   for ii = 1:dipnum, xsect_name{ii} = flt_name; end		% flt_name - strings for fault 3 & 4
   xsect = zeros(dipnum,8);
   x2 = 0;
   dipin = sortrows(dipin,2);		% order from surface down
   for ii=1:dipnum
      dip = dipin(ii,1); z1 = dipin(ii,2); z2 = dipin(ii,3); rows = dipin(ii,4); 
      width = (z2-z1)/tand(dip);
      x1 = x2; x2 = x1+width;
      xsect(ii,:) = [ ii dip x1 z1 x2 z2 width rows ];                     
   end
else
   error('fault type 1 & 2 can not have different dips!');
end
