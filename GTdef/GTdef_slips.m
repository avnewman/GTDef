function [ flt1_out,flt2_out,subflt_name,subflt_out ]...
          = GTdef_slips(lb,ub,xx,flt1_num,flt1,flt2_num,flt2,...
                        flt3_name,flt3_num,flt3,...
			flt4_name,flt4_num,flt4,...
			flt5_name,flt5_num,flt5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_slips				        %
% Having the final values for the slips, put back the slips in the same	        %
% format as the input file		  				        %
%									        %
% INPUT:					  		  	        %
% Each fault has three (strike, dip, and tensile) components, so                %
% slip_num = flt_num*3                                                          %
%  xx   - final values for ss,ds,ts 	[slip_num*1]                            %
%  lb   - lower bounds for ss,ds,ts 	[slip_num*1]                            %
%  ub   - upper bounds for ss,ds,ts	[slip_num*1]			        %
%  flt1 - [lon lat z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]     	%
%  flt2 - [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]   	%
%  flt3 - [lon1 lat1 z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]	%
%  flt4 - [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]%
%  flt5 - [ ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns dlen slen ]		  	%
%  subflt - [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		  	%
%									        %
% OUTPUT:                                                                       %
%  flt1_out - [lon lat z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]	%
%             with ss, ds, ts replaced by inverted values			%
%  flt2_out - [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]  %
%             with ss, ds, ts replaced by inverted values			%
% Full info about subfaults							%
%  subflt_name - master-fault names						%
%  subflt_out - [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]			%
%                                                                               %
% first created by Lujia Feng Wed May  6 16:21:17 EDT 2009		        %
% added fault5 lfeng Fri Dec 11 12:57:14 EST 2009				%
% use cell array of strings for names lfeng Wed Dec  1 17:50:08 EST 2010	%
% last modified by Lujia Feng Wed Dec  1 17:50:15 EST 2010			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reset lb/ub from -inf/inf to 0 for slips whose values are fixed
fix_ind = find(isinf(lb));
lb(fix_ind) = 0;
ub(fix_ind) = 0;

first = 0; last = 0;
flt1_out = []; flt2_out = []; flt3_out = []; flt4_out = []; 
subflt_name = {}; subflt_out = []; subflt = [];
%%%%%%%%%% fault 1 %%%%%%%%%%
if flt1_num~=0
    first = last+1;
    last  = last+flt1_num*3;
    slp1  = reshape(xx(first:last),[],3);
    flt1_out = [ flt1(:,1:7) slp1 flt1(:,11:16) ];
end

%%%%%%%%%% fault 2 %%%%%%%%%%
if flt2_num~=0
    first = last+1;
    last  = last+flt2_num*3;
    slp2  = reshape(xx(first:last),[],3);
    flt2_out = [ flt2(:,1:7) slp2 flt2(:,11:16) ];
end

%%%%%%%%%% fault 3 %%%%%%%%%%
if flt3_num~=0
    for ii = 1:flt3_num
        Nd = flt3(ii,17); Ns = flt3(ii,18); flt_num = Nd*Ns;
        first = last+1;   last = last+flt_num*3;
        slp3 = reshape(xx(first:last),[],3);
        lb3 = reshape(lb(first:last),[],3);
        ub3 = reshape(ub(first:last),[],3);
	dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
	subflt = [ dnum snum slp3 lb3(:,1) ub3(:,1) lb3(:,2) ub3(:,2) lb3(:,3) ub3(:,3) ];
	subflt_out = [ subflt_out; subflt ];

	flt_name = flt3_name{ii};
        name = cell(flt_num,1);
        for ii = 1:flt_num, name{ii} = flt_name; end
	subflt_name = [ subflt_name; name ];   
    end
end

%%%%%%%%%% fault 4 %%%%%%%%%%
if flt4_num~=0
    for ii = 1:flt4_num
        Nd = flt4(ii,17); Ns = flt4(ii,18); flt_num = Nd*Ns;
        first = last+1;   last = last+flt_num*3;
        slp4 = reshape(xx(first:last),[],3);
        lb4 = reshape(lb(first:last),[],3);
        ub4 = reshape(ub(first:last),[],3);
	dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
	subflt = [ dnum snum slp4 lb4(:,1) ub4(:,1) lb4(:,2) ub4(:,2) lb4(:,3) ub4(:,3) ];
	subflt_out = [ subflt_out; subflt ];

	flt_name = flt4_name{ii};
        name = cell(flt_num,1);
        for ii = 1:flt_num, name{ii} = flt_name; end
	subflt_name = [ subflt_name; name ];   
    end
end

%%%%%%%%%% fault 5 %%%%%%%%%%
if flt5_num~=0
    for ii = 1:flt5_num
        Nd = flt5(ii,10); Ns = flt5(ii,11); flt_num = Nd*Ns;
        first = last+1;   last = last+flt_num*3;
        slp5 = reshape(xx(first:last),[],3);
        lb5 = reshape(lb(first:last),[],3);
        ub5 = reshape(ub(first:last),[],3);
	dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
	subflt = [ dnum snum slp5 lb5(:,1) ub5(:,1) lb5(:,2) ub5(:,2) lb5(:,3) ub5(:,3) ];
	subflt_out = [ subflt_out; subflt ];

	flt_name = flt5_name(ii,:);
        name = cell(flt_num,1);
        for ii = 1:flt_num, name{ii} = flt_name; end
	subflt_name = [ subflt_name; name ];   
    end
end
