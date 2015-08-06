function [ flt1_out,flt2_out,flt3_out,flt4_out,flt5_out,subflt_out,subflt_name ]...
          = GTdef_slips(lb,ub,xx,flt1,flt2,flt3,flt4,flt5,subflt_in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   GTdef_slips				                 %
% Having the final values for the slips, put back the slips in the same	                 %
% format as the input file		  				                 %
%									                 %
% INPUT:					  		  	                 %
%----------------------------------------------------------------------------------------%
% fault1 & fault2 have three (strike, dip, and tensile) components, so                   %
% slip_num = flt_num*3                                                                   %
%      xx  - final values for ss,ds,ts 	[slip_num*1]                                     %
%      lb  - lower bounds for ss,ds,ts 	[slip_num*1]                                     %
%      ub  - upper bounds for ss,ds,ts	[slip_num*1]			                 %
% flt1.flt - [lon1 lat1 z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]        %
% flt2.flt - [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]      %
%      subflt.flt - [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		  	 %
%----------------------------------------------------------------------------------------%
% fault3 & fault4 have two (rake and tensile) components, so                             %
% slip_num = flt_num*2                                                                   %
%      xx  - final values for rs,ts 	[slip_num*1]                                     %
%      lb  - lower bounds for rs,ts 	[slip_num*1]                                     %
%      ub  - upper bounds for rs,ts	[slip_num*1]			                 %
% flt3.flt - [lon1 lat1 z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns]  %
% flt4.flt - [lon1 lat1 lon2 lat2 z1 z2 dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns]%
%      subflt.flt - [ dnum snum rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]                 %
%----------------------------------------------------------------------------------------%
% flt5.flt - [ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]                                    %
%									                 %
% OUTPUT:                                                                                %
% flt1_out - [lon lat z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]	         %
%             with ss, ds, ts replaced by inverted values			         %
% flt2_out - [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]      %
%             with ss, ds, ts replaced by inverted values			         %
% Full info about subfaults							         %
% subflt_name - master-fault names						         %
% subflt_out  - [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]			         %
%                                                                                        %
% first created by Lujia Feng Wed May  6 16:21:17 EDT 2009		                 %
% added fault5 lfeng Fri Dec 11 12:57:14 EST 2009				         %
% used cell array of strings for names lfeng Wed Dec  1 17:50:08 EST 2010	         %
% used structure lfeng Wed Feb 22 19:45:15 SGT 2012				         %
% merged flt1 & flt3 and flt2 & flt4 lfeng Wed May  9 10:40:55 SGT 2012                  %
% added flt5 lfeng Sat Dec  1 22:00:47 SGT 2012                                          %
% corrected flt3 & flt4 lfeng Mon Mar 11 20:27:58 SGT 2013                               %
% last modified by Lujia Feng Mon Mar 11 20:40:43 SGT 2013                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reset lb/ub from -inf/inf to 0 for slips whose values are fixed
fix_ind = find(isinf(lb));
lb(fix_ind) = 0;
ub(fix_ind) = 0;

first = 0; last = 0;
flt1_out = flt1.flt; flt2_out = flt2.flt; flt3_out = flt3.flt; flt4_out = flt4.flt; flt5_out = flt5.flt;
subflt_name = {}; subflt_out = []; subflt = [];
%%%%%%%%%% fault 1 %%%%%%%%%%
if flt1.num~=0
    comp_num = 3;
    for ii = 1:flt1.num
        Nd = flt1.flt(ii,17); Ns = flt1.flt(ii,18); fltNum = Nd*Ns;
        first = last+1; last = last+fltNum*comp_num;
	if fltNum==1
            slp1 = reshape(xx(first:last),[],comp_num);
            flt1_out(ii,8:10) = slp1;
	else
            slp1 = reshape(xx(first:last),[],comp_num);
            lb1  = reshape(lb(first:last),[],comp_num);
            ub1  = reshape(ub(first:last),[],comp_num);
	    dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	    slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
	    subflt     = [ dnum snum slp1 lb1(:,1) ub1(:,1) lb1(:,2) ub1(:,2) lb1(:,3) ub1(:,3) ];
	    subflt_out = [ subflt_out; subflt ];

	    flt_name = flt1.name{ii};
            name     = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = flt_name; end
	    subflt_name = [ subflt_name; name ];   
	end
    end
end

%%%%%%%%%% fault 2 %%%%%%%%%%
if flt2.num~=0
    comp_num = 3;
    for ii = 1:flt2.num
        Nd = flt2.flt(ii,17); Ns = flt2.flt(ii,18); fltNum = Nd*Ns;
        first = last+1; last = last+fltNum*comp_num;
	if fltNum==1
            slp2 = reshape(xx(first:last),[],comp_num);
            flt2_out(ii,8:10) = slp2;
	else
            slp2 = reshape(xx(first:last),[],comp_num);
            lb2  = reshape(lb(first:last),[],comp_num);
            ub2  = reshape(ub(first:last),[],comp_num);
	    dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	    slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
	    subflt     = [ dnum snum slp2 lb2(:,1) ub2(:,1) lb2(:,2) ub2(:,2) lb2(:,3) ub2(:,3) ];
	    subflt_out = [ subflt_out; subflt ];

	    flt_name = flt2.name{ii};
            name = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = flt_name; end
	    subflt_name = [ subflt_name; name ];   
	end
    end
end

%%%%%%%%%% fault 3 %%%%%%%%%%
if flt3.num~=0
   comp_num = 2;
   for ii = 1:flt3.num
        Nd = flt3.flt(ii,17); Ns = flt3.flt(ii,18); fltNum = Nd*Ns;
        first = last+1;   last = last+fltNum*comp_num;
        if fltNum==1
            slp3 = reshape(xx(first:last),[],comp_num);
            flt3_out(ii,9:10) = slp3;
        else
            slp3 = reshape(xx(first:last),[],comp_num);
            lb3  = reshape(lb(first:last),[],comp_num);
            ub3  = reshape(ub(first:last),[],comp_num);
            dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
            slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
            % rake of master faults
            rake  = flt3.flt(ii,8);  rlin  = ones(fltNum,1).*rake;
            rake0 = flt3.flt(ii,11); rlin0 = ones(fltNum,1).*rake0;
            rakeX = flt3.flt(ii,12); rlinX = ones(fltNum,1).*rakeX;
            subflt = [ dnum snum rlin slp3 rlin0 rlinX lb3(:,1) ub3(:,1) lb3(:,2) ub3(:,2) ];
            % subfault names
            flt_name = flt3.name{ii};
            name = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = flt_name; end
            subflt_name = [ subflt_name; name ];   
    	    % rake of subfaults
    	    sub_ind  = strcmpi(flt_name,subflt_in.name);
            subflt3  = subflt_in.flt(sub_ind,:);
            numT = size(subflt3,1);
            for ii=1:numT
                dd = subflt3(ii,1); ss = subflt3(ii,2);
                ind = subflt(:,1)==dd & subflt(:,2)==ss;
                subflt(ind,[3 6 7])  = subflt3(ii,[3 6 7]);
            end
            subflt_out = [ subflt_out; subflt ];
        end
    end
end

%%%%%%%%%% fault 4 %%%%%%%%%%
if flt4.num~=0
    comp_num = 2;
    for ii = 1:flt4.num
        Nd = flt4.flt(ii,17); Ns = flt4.flt(ii,18); fltNum = Nd*Ns;
        first = last+1;   last = last+fltNum*comp_num;
	if fltNum==1
            slp4 = reshape(xx(first:last),[],comp_num);
            flt4_out(ii,9:10) = slp4;
	else
            slp4 = reshape(xx(first:last),[],comp_num);
            lb4  = reshape(lb(first:last),[],comp_num);
            ub4  = reshape(ub(first:last),[],comp_num);
	    dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	    slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
	    % rake of master faults
	    rake  = flt4.flt(ii,8);  rlin  = ones(fltNum,1)*rake;
	    rake0 = flt4.flt(ii,11); rlin0 = ones(fltNum,1)*rake0;
	    rakeX = flt4.flt(ii,12); rlinX = ones(fltNum,1)*rakeX;
	    subflt = [ dnum snum rlin slp4 rlin0 rlinX lb4(:,1) ub4(:,1) lb4(:,2) ub4(:,2) ];
	    % subfault names
	    flt_name = flt4.name{ii};
            name = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = flt_name; end
	    subflt_name = [ subflt_name; name ];   
    	    % rake of subfaults
    	    sub_ind  = strcmpi(flt_name,subflt_in.name);
	    subflt4  = subflt_in.flt(sub_ind,:);
	    numT = size(subflt4,1);
	    for ii=1:numT
	        dd = subflt4(ii,1); ss = subflt4(ii,2);
		ind = subflt(:,1)==dd & subflt(:,2)==ss;
		subflt(ind,[3 6 7])  = subflt4(ii,[3 6 7]);
	    end
	    subflt_out = [ subflt_out; subflt ];
	end
    end
end

%%%%%%%%%% fault 5 %%%%%%%%%%
if flt5.num~=0
    comp_num = 3;
    for ii = 1:flt5.num
        Nd = flt5.flt(ii,10); Ns = flt5.flt(ii,11); fltNum = Nd*Ns;
        first = last+1; last = last+fltNum*comp_num;
	if fltNum==1
            slp1 = reshape(xx(first:last),[],comp_num);
            flt5_out(ii,1:3) = slp1;
	else
            slp5 = reshape(xx(first:last),[],comp_num);
            lb5  = reshape(lb(first:last),[],comp_num);
            ub5  = reshape(ub(first:last),[],comp_num);
	    dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	    slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
	    subflt     = [ dnum snum slp5 lb5(:,1) ub5(:,1) lb5(:,2) ub5(:,2) lb5(:,3) ub5(:,3) ];
	    subflt_out = [ subflt_out; subflt ];

	    flt_name = flt5.name{ii};
            name     = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = flt_name; end
	    subflt_name = [ subflt_name; name ];   
	end
    end
end
