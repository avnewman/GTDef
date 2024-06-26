function [ flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt ] = GTdef_update_slips(earth,modspace,flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   GTdef_update_slips				             %
% Update the final values for the slips                                                      %
% Put back the slips in the same format as the input file		  		     %
%									                     %
% INPUT:					  		  	                     %
% earth    - earth structure                                                                 %
% modspace - model structure                                                                 %
% flt?     - fault structure	  	                                                     %
% subflt   - subfault structure                                                              %
%									                     %
% OUTPUT:                                                                                    %
% ------------------------------------------------------------------------------------------ %
% (1) add .out to fault structures                                                           %
% flt?.out - with ss, ds, ts replaced by inverted values			             %
%									                     %
% (2) update flt?.Min, flt?.slip, and flt?.rake                                              %
%									                     %
% (3) add .out & .outname to subfault structures                                             %
% subflt.outname - master-fault names						             %
% subflt.out     = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]			     %
%                                                                                            %
% first created by Lujia Feng Wed May  6 16:21:17 EDT 2009		                     %
% added fault5 lfeng Fri Dec 11 12:57:14 EST 2009				             %
% used cell array of strings for names lfeng Wed Dec  1 17:50:08 EST 2010	             %
% used structure lfeng Wed Feb 22 19:45:15 SGT 2012				             %
% merged flt1 & flt3 and flt2 & flt4 lfeng Wed May  9 10:40:55 SGT 2012                      %
% added flt5 lfeng Sat Dec  1 22:00:47 SGT 2012                                              %
% corrected flt3 & flt4 lfeng Mon Mar 11 20:27:58 SGT 2013                                   %
% changed from GTdef_slips.m to GTdef_update_slips.m lfeng Fri Mar 20 14:49:44 SGT 2015      %
% added updating flt?.Min lfeng Fri Mar 20 15:21:26 SGT 2015                                 %
% added calculating flt?.slip & flt?.rake lfeng Wed Mar 25 01:20:22 SGT 2015                 %
% removed calculating flt?.slip & flt?.rake lfeng Wed Mar 25 16:15:48 SGT 2015               %
% added fault5 lfeng Tue Jun 23 17:29:48 SGT 2015                                            %
% added fault6, changed old fault6 to fault7 lfeng Wed Jun  1 18:48:36 SGT 2016              %
% last modified by Lujia Feng Thu Jun  2 10:29:27 SGT 2016                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb    = modspace.lb;
ub    = modspace.ub;
xx    = modspace.xx;
%x0    = modspace.x0;
etype = earth.type;

% reset lb/ub from -inf/inf to 0 for slips whose values are fixed
fix_ind = find(isinf(lb));
lb(fix_ind) = 0;
ub(fix_ind) = 0;

first = 0; last = 0;
subflt.outname = {}; subflt.out = []; newsubflt = [];
%%%%%%%%%% fault 1 %%%%%%%%%%
%             1    2    3  4  5   6   7   8  9  10 11  12  13  14  15  16  17 18
% flt1.flt = [lon1 lat1 z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]
if flt1.num~=0
    compNum = 3;
    for ii = 1:flt1.num
        Nd = flt1.flt(ii,17); Ns = flt1.flt(ii,18); fltNum = Nd*Ns;
        first = last+1; last = last+fltNum*compNum;
        % single fault
        if fltNum==1
            slp1 = reshape(xx(first:last),[],compNum);
            flt1.out(ii,8:10)  = slp1;
        % multiple subfaults
        else
            slp1 = reshape(xx(first:last),[],compNum);
            lb1  = reshape(lb(first:last),[],compNum);
            ub1  = reshape(ub(first:last),[],compNum);
	    dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	    slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
            newsubflt  = [ dnum snum slp1 lb1(:,1) ub1(:,1) lb1(:,2) ub1(:,2) lb1(:,3) ub1(:,3) ];
	    subflt.out = [ subflt.out; newsubflt ];

	    fltName = flt1.name{ii};
            name     = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = fltName; end
	    subflt.outname = [ subflt.outname; name ];   
	end

        % update flt1.slip & flt1.rake
        ss = slp1(:,1); 
        ds = slp1(:,2); 
        slip = sqrt(ss.^2+ds.^2);
        rake = atan2(ds,ss).*180/pi;
        %flt1.slip{ii} = slip;
        %flt1.rake{ii} = rake;

        % update flt1.Min
        %                1   2     3     4   5   6    7     8  9  10
        % Okada   Min = [len width depth dip str east north ss ds ts]     [flt_num*10]
        %                1    2     3    4     5      6     7   8   9
        % layered Min = [slip north east depth length width str dip rake] [flt_num*9]
        if strcmpi(etype,'homogeneous')
            flt1.Min{ii}(:,8:10) = slp1;
        else
            flt1.Min{ii}(:,1)   = slip;
            flt1.Min{ii}(:,end) = rake;
        end
    end
end

%%%%%%%%%% fault 2 %%%%%%%%%%
%             1    2    3    4    5  6  7   8  9  10 11  12  13  14  15  16  17 18
% flt2.flt = [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]
if flt2.num~=0
    compNum = 3;
    for ii = 1:flt2.num
        Nd = flt2.flt(ii,17); Ns = flt2.flt(ii,18); fltNum = Nd*Ns;
        first = last+1; last = last+fltNum*compNum;
        % single fault
	if fltNum==1
            slp2 = reshape(xx(first:last),[],compNum);
            flt2.out(ii,8:10) = slp2;
        % multiple subfaults
	else
            slp2 = reshape(xx(first:last),[],compNum);
            lb2  = reshape(lb(first:last),[],compNum);
            ub2  = reshape(ub(first:last),[],compNum);
	    dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	    slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
	    newsubflt  = [ dnum snum slp2 lb2(:,1) ub2(:,1) lb2(:,2) ub2(:,2) lb2(:,3) ub2(:,3) ];
	    subflt.out = [ subflt.out; newsubflt ];

	    fltName = flt2.name{ii};
            name = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = fltName; end
	    subflt.outname = [ subflt.outname; name ];   
	end

        % update flt2.slip & flt2.rake
        ss = slp2(:,1); 
        ds = slp2(:,2); 
        slip = sqrt(ss.^2+ds.^2);
        rake = atan2(ds,ss).*180/pi;
        %flt2.slip{ii} = slip;
        %flt2.rake{ii} = rake;

        % update flt2.Min
        %                1   2     3     4   5   6    7     8  9  10
        % Okada   Min = [len width depth dip str east north ss ds ts]     [flt_num*10]
        %                1    2     3    4     5      6     7   8   9
        % layered Min = [slip north east depth length width str dip rake] [flt_num*9]
        if strcmpi(etype,'homogeneous')
            flt2.Min{ii}(:,8:10) = slp2;
        else
            flt2.Min{ii}(:,1)   = slip;
            flt2.Min{ii}(:,end) = rake;
        end
    end
end

%%%%%%%%%% fault 3 %%%%%%%%%%
%             1    2    3  4  5   6   7   8    9  10 11    12    13  14  15  16  17 18
% flt3.flt = [lon1 lat1 z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns]
if flt3.num~=0
   compNum = 2;
   for ii = 1:flt3.num
        Nd = flt3.flt(ii,17); Ns = flt3.flt(ii,18); fltNum = Nd*Ns;
        first = last+1;   last = last+fltNum*compNum;
	% single fault
        if fltNum==1
            slp3 = reshape(xx(first:last),[],compNum);
            flt3.out(ii,9:10) = slp3;
        % multiple subfaults
        else
            slp3 = reshape(xx(first:last),[],compNum);
            lb3  = reshape(lb(first:last),[],compNum);
            ub3  = reshape(ub(first:last),[],compNum);
            dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
            slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
            % rake of master faults
            rake  = flt3.flt(ii,8);  rlin  = ones(fltNum,1).*rake;
            rake0 = flt3.flt(ii,11); rlin0 = ones(fltNum,1).*rake0;
            rakeX = flt3.flt(ii,12); rlinX = ones(fltNum,1).*rakeX;
            newsubflt = [ dnum snum rlin slp3 rlin0 rlinX lb3(:,1) ub3(:,1) lb3(:,2) ub3(:,2) ];
            % subfault names
            fltName = flt3.name{ii};
            name = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = fltName; end
            subflt.outname = [ subflt.outname; name ];   
    	    % rake of subfaults
    	    sub_ind  = strcmpi(fltName,subflt.name);
            subflt3  = subflt.flt(sub_ind,:);
            numT = size(subflt3,1);
            for ii=1:numT
                dd  = subflt3(ii,1); ss = subflt3(ii,2);
                ind = newsubflt(:,1)==dd & newsubflt(:,2)==ss;
                newsubflt(ind,[3 6 7])  = subflt3(ii,[3 6 7]);
            end
            subflt.out = [ subflt.out; newsubflt ];
        end

        % update flt3.slip & flt3.rake
        rs   = slp3(:,1); % rake slip
        ts   = slp3(:,2); % tensile slip
        rake = newsubflt(:,3);
        %flt3.slip{ii} = rs;
        %flt3.rake{ii} = rake;

        % update flt3.Min
        %                1   2     3     4   5   6    7     8  9  10
        % Okada   Min = [len width depth dip str east north ss ds ts]     [flt_num*10]
        %                1    2     3    4     5      6     7   8   9
        % layered Min = [slip north east depth length width str dip rake] [flt_num*9]
        if strcmpi(etype,'homogeneous')
            ss = rs.*cosd(rake); 
            ds = rs.*sind(rake);
            flt3.Min{ii}(:,8:10) = [ss ds ts];
        else
            flt3.Min{ii}(:,1) = rs;
            % rake cannot be changed
        end
    end
end

%%%%%%%%%% fault 4 %%%%%%%%%%
%             1    2    3    4    5  6  7   8    9  10 11    12    13  14  15  16  17 18
% flt4.flt = [lon1 lat1 lon2 lat2 z1 z2 dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns]
if flt4.num~=0
    compNum = 2;
    for ii = 1:flt4.num
        Nd = flt4.flt(ii,17); Ns = flt4.flt(ii,18); fltNum = Nd*Ns;
        first = last+1;   last = last+fltNum*compNum;
	% single fault
        if fltNum==1
            slp4 = reshape(xx(first:last),[],compNum);
            flt4.out(ii,9:10) = slp4;
        % multiple subfaults
        else
            slp4 = reshape(xx(first:last),[],compNum);
            lb4  = reshape(lb(first:last),[],compNum);
            ub4  = reshape(ub(first:last),[],compNum);
	    dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	    slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
	    % rake of master faults
	    rake  = flt4.flt(ii,8);  rlin  = ones(fltNum,1)*rake;
	    rake0 = flt4.flt(ii,11); rlin0 = ones(fltNum,1)*rake0;
	    rakeX = flt4.flt(ii,12); rlinX = ones(fltNum,1)*rakeX;
	    newsubflt = [ dnum snum rlin slp4 rlin0 rlinX lb4(:,1) ub4(:,1) lb4(:,2) ub4(:,2) ];
	    % subfault names
	    fltName = flt4.name{ii};
            name = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = fltName; end
	    subflt.outname = [ subflt.outname; name ];   
    	    % rake of subfaults
    	    sub_ind  = strcmpi(fltName,subflt.name);
	    subflt4  = subflt.flt(sub_ind,:);
	    numT = size(subflt4,1);
	    for ii=1:numT
	        dd = subflt4(ii,1); ss = subflt4(ii,2);
		ind = newsubflt(:,1)==dd & newsubflt(:,2)==ss;
		newsubflt(ind,[3 6 7])  = subflt4(ii,[3 6 7]);
	    end
	    subflt.out = [ subflt.out; newsubflt ];
	end

        % update flt4.slip & flt4.rake
        rs   = slp4(:,1); % rake slip
        ts   = slp4(:,2); % tensile slip
        rake = newsubflt(:,3);
        %flt4.slip{ii} = rs;
        %flt4.rake{ii} = rake;

        % update Min
        %                1   2     3     4   5   6    7     8  9  10
        % Okada   Min = [len width depth dip str east north ss ds ts]     [flt_num*10]
        %                1    2     3    4     5      6     7   8   9
        % layered Min = [slip north east depth length width str dip rake] [flt_num*9]
        if strcmpi(etype,'homogeneous')
            ss = rs.*cosd(rake); 
            ds = rs.*sind(rake);
            flt4.Min{ii}(:,8:10) = [ss ds ts];
        else
            flt4.Min{ii}(:,1) = rs;
            % rake cannot be changed
        end
    end
end

%%%%%%%%%% fault 5 %%%%%%%%%%
%             1  2  3  4   5   6   7   8   9   10 11
% flt5.flt = [ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]
if flt5.num~=0
    compNum = 3;
    for ii = 1:flt5.num
        Nd = flt5.out(ii,10); Ns = flt5.out(ii,11); fltNum = Nd*Ns;
        first = last+1; last = last+fltNum*compNum;
        if fltNum==1
            slp5 = reshape(xx(first:last),[],compNum);
            flt5.out(ii,1:3) = slp5;
        else
            slp5 = reshape(xx(first:last),[],compNum);
            lb5  = reshape(lb(first:last),[],compNum);
            ub5  = reshape(ub(first:last),[],compNum);
            dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
            slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
            newsubflt     = [ dnum snum slp5 lb5(:,1) ub5(:,1) lb5(:,2) ub5(:,2) lb5(:,3) ub5(:,3) ];
            subflt.out = [ subflt.out; newsubflt ];

            fltName = flt5.name{ii};
            name     = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = fltName; end
            subflt.outname = [ subflt.outname; name ];   
	end

        % update flt5.slip & flt5.rake
        ss = slp5(:,1); 
        ds = slp5(:,2); 
        slip = sqrt(ss.^2+ds.^2);
        rake = atan2(ds,ss).*180/pi;
        %flt5.slip{ii} = slip;
        %flt5.rake{ii} = rake;

        % update flt5.Min
        %                1   2     3     4   5   6    7     8  9  10
        % Okada   Min = [len width depth dip str east north ss ds ts]     [flt_num*10]
        %                1    2     3    4     5      6     7   8   9
        % layered Min = [slip north east depth length width str dip rake] [flt_num*9]
        if strcmpi(etype,'homogeneous')
            flt5.Min{ii}(:,8:10) = slp5;
        else
            flt5.Min{ii}(:,1)   = slip;
            flt5.Min{ii}(:,end) = rake;
        end
    end
end

%%%%%%%%%% fault 6 %%%%%%%%%%
%             1    2  3  4     5     6   7   8   9   10 11
% flt6.flt = [rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns]
if flt6.num~=0
    compNum = 2;
    for ii = 1:flt6.num
        Nd = flt6.flt(ii,10); Ns = flt6.flt(ii,11); fltNum = Nd*Ns;
        first = last+1; last = last+fltNum*compNum;
	if fltNum==1
            slp6 = reshape(xx(first:last),[],compNum);
            flt6.out(ii,2:3) = slp6;
	else
            slp6 = reshape(xx(first:last),[],compNum);
            lb6  = reshape(lb(first:last),[],compNum);
            ub6  = reshape(ub(first:last),[],compNum);
	    dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	    slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
            % rake of master faults
            rake  = flt6.flt(ii,1);  rlin  = ones(fltNum,1).*rake;
            rake0 = flt6.flt(ii,4);  rlin0 = ones(fltNum,1).*rake0;
            rakeX = flt6.flt(ii,5);  rlinX = ones(fltNum,1).*rakeX;
            newsubflt = [ dnum snum rlin slp6 rlin0 rlinX lb6(:,1) ub6(:,1) lb6(:,2) ub6(:,2) ];
            % subfault names
            fltName = flt6.name{ii};
            name = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = fltName; end
            subflt.outname = [ subflt.outname; name ];   
    	    % rake of subfaults
    	    sub_ind  = strcmpi(fltName,subflt.name);
            subflt6  = subflt.flt(sub_ind,:);
            numT = size(subflt6,1);
            for ii=1:numT
                dd  = subflt6(ii,1); ss = subflt6(ii,2);
                ind = newsubflt(:,1)==dd & newsubflt(:,2)==ss;
                newsubflt(ind,[3 6 7])  = subflt6(ii,[3 6 7]);
            end
            subflt.out = [ subflt.out; newsubflt ];
	end

        % update flt6.slip & flt6.rake
        rs   = slp6(:,1); % rake slip
        ts   = slp6(:,2); % tensile slip
        rake = newsubflt(:,3);
        %flt6.slip{ii} = rs;
        %flt6.rake{ii} = rake;

        % update flt6.Min
        %                1   2     3     4   5   6    7     8  9  10
        % Okada   Min = [len width depth dip str east north ss ds ts]     [flt_num*10]
        %                1    2     3    4     5      6     7   8   9
        % layered Min = [slip north east depth length width str dip rake] [flt_num*9]
        if strcmpi(etype,'homogeneous')
            ss = rs.*cosd(rake); 
            ds = rs.*sind(rake);
            flt6.Min{ii}(:,8:10) = [ss ds ts];
        else
            flt6.Min{ii}(:,1) = rs;
            % rake cannot be changed
        end
    end
end

%%%%%%%%%% fault 7 %%%%%%%%%%
%             1  2  3  4   5   6   7   8   9   10 11
% flt7.flt = [ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]
if flt7.num~=0
    compNum = 3;
    for ii = 1:flt7.num
        Nd = flt7.flt(ii,10); Ns = flt7.flt(ii,11); fltNum = Nd*Ns;
        first = last+1; last = last+fltNum*compNum;
	if fltNum==1
            slp5 = reshape(xx(first:last),[],compNum);
            flt7.out(ii,1:3) = slp5;
	else
            slp5 = reshape(xx(first:last),[],compNum);
            lb5  = reshape(lb(first:last),[],compNum);
            ub5  = reshape(ub(first:last),[],compNum);
	    dlin = [ 1:Nd ]'; dmat = dlin(:,ones(Ns,1)); dnum = reshape(dmat,[],1);
	    slin = [ 1:Ns ];  smat = slin(ones(Nd,1),:); snum = reshape(smat,[],1);
	    newsubflt     = [ dnum snum slp5 lb5(:,1) ub5(:,1) lb5(:,2) ub5(:,2) lb5(:,3) ub5(:,3) ];
	    subflt.out = [ subflt.out; newsubflt ];

	    fltName = flt7.name{ii};
            name     = cell(fltNum,1);
            for ii = 1:fltNum, name{ii} = fltName; end
	    subflt.outname = [ subflt.outname; name ];   
	end
    end
end
