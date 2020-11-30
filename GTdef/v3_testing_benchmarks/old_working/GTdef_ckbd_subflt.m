function [ ] = GTdef_ckbd_subflt(subflt_name,Nd,Ns,dd,ss,p1,p0)

% run example
% GTdef_ckbd_subflt('nicoya',30,40,3,4,[0 -0.0823 0 0 0 0 0 0 0 ],[ 0 0 0 0 0 0 0 0 0 ])
% GTdef_ckbd_subflt('nicoya',30,40,3,4,[ 0 0 0 0 0 0 0 0 0 ],[0 -0.0823 0 0 0 0 0 0 0 ])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          GTdef_ckbd_subflt				  %
% Create the checkerboard input for subfaults				  %
%									  %
% INPUT:								  %
%   subflt_name - main fault name that is used for output		  %
%   Nd - number of patches along dip				          %
%   Ns - number of patches along strike					  %
%   dd - number of patches along strike that are combined together	  %
%   ss - number of patches along dip that are combined together	  	  %
%   p1 - pattern 1 (black,1)						  %
%      = [ ss ds ts ss0 ssX ds0 dsX ts0 tsX ]				  %
%      = [ ss ds ts 0 0 0 0 0 0 ] (for forward)   			  %
%   p0 - pattern 0 (white,0)						  %
%      = [ ss ds ts ss0 ssX ds0 dsX ts0 tsX ]				  %
%      = [ 0 0 0 0 0 0 0 0 0 ] (for forward)   			  	  %
%									  %
%  The example is Nd=4; Ns=6; dd=2; ss=3				  %
%    _________________________						  %
%    |	 |   |   |	 |   |   					  %
%    | 1 | 1 | 1 | 0 | 0 | 0 |                                            %
%    |___|___|___|___|___|___|                                            %
%    |	 |   |   |	 |   |                                            %
%    | 1 | 1 | 1 | 0 | 0 | 0 |                                            %
%    |___|___|___|___|___|___|                                            %
%    |	 |   |   |	 |   |   					  %
%    | 0 | 0 | 0 | 1 | 1 | 1 |                                            %
%    |___|___|___|___|___|___|                                            %
%    |	 |   |   |	 |   |                                            %
%    | 0 | 0 | 0 | 1 | 1 | 1 |                                            %
%    |___|___|___|___|___|___|                                            %
% 									  %
% OUTPUT:								  %
%   is part of the forward checkerboard test input file			  %
%   subfault name dnum snum  ss ds ts ss0 ssX ds0 dsX ts0 tsX		  %
%								 	  %
% first created by Lujia Feng Wed Feb 24 08:56:45 EST 2010		  %
% last modified by Lujia Feng Wed Feb 24 10:09:18 EST 2010		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default output file name
filename = 'GTdef_ckbd_subflt.out';
fout = fopen(filename,'w');

for jj = 1:Ns 
   is = ceil(jj/ss);
   ms = mod(is,2);
   for ii = 1:Nd
      id = ceil(ii/dd);
      md = mod(id,2);
      if md==ms
         fprintf(fout,'     subfault %s  %5d %5d  %10.5f %-8.5f %-8.5f  %-5.2f %-5.2f  %-5.2f %-5.2f  %-5.2f %-5.2f\n',subflt_name,ii,jj,p1);
      else
         fprintf(fout,'     subfault %s  %5d %5d  %10.5f %-8.5f %-8.5f  %-5.2f %-5.2f  %-5.2f %-5.2f  %-5.2f %-5.2f\n',subflt_name,ii,jj,p0);
      end
   end
end

fclose(fout);
