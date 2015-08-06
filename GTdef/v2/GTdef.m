function [] = GTdef(finName,wnum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              	  GTdef.m                                       %
%              Georgia Tech Matlab program for deformation                      %
%                Lujia Feng; Andrew V. Newman; Ting Chen                        %
%                                                                               %
% INPUT:                                                                        %
% finName - input file name                                                     %
% wnum    - num of matlab parallel workers to be used                           %
%   0: do not use parallel computing                                            %
%   8: up to 8 workers that can be specified                                    %
%  -1: use parallel computing                                                   %
%      but number of workers will be determined by the system                   %
%                                                                               %
% v1:                                                                           %
% first created by Lujia Feng Mon Apr 20 14:15:52 EDT 2009                      %
% last modified by Lujia Feng Tue May 19 01:14:49 EDT 2009                      %
% added coordiante type flag 'coord' lfeng Thu Nov 5 17:12:59 EST 2009          %
% added first derivative modes lfeng Tue Dec  1 14:31:10 EST 2009               %
% modified beta (added beta) and roughness for 1st derivatives                  %
% added 'dip' flag for bended faults lfeng Mon Dec  7 01:04:06 EST 2009         %
% added 'freesurface' flag lfeng Wed Dec  9 17:00:58 EST 2009                   %
% added fault type 5 lfeng Fri Dec 11 10:57:18 EST 2009                         %
% changed 'freesurface' to 'surface' flag lfeng Wed Feb 24 12:46:01 EST 2010	%
% changed 'coord' to string flag lfeng Wed Feb 24 13:40:01 EST 2010             %
% allows input file name include multiple "." besides ".in" lfeng Oct 4 2010	%
% added matlabpool lfeng Wed Dec  1 12:12:00 EST 2010                           %
% edited the origin definition lfeng Thu Apr  7 18:45:35 EDT 2011               %
% v2:                                                                           %
% used strucutres for passing parameters lfeng Wed Feb 22 03:39:11 SGT 2012     %
% added layered earth model lfeng Tue Feb 28 03:44:29 SGT 2012                  %
% new definition used for fault1, fault2, fault3 & fault4 lfeng May 8 2012      %
% added polyconic projection lfeng Thu Jun  7 12:23:26 SGT 2012                 %
% last modified lfeng Wed Jun 13 18:09:47 SGT 2012                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% specify matlabpool for parallel computing %%%%%%%%%%%%%%%%%%%%%%%%%%%
if wnum==0		% do not use parallel computing
    if matlabpool('size')>0
       matlabpool close
    end
elseif wnum<0		% use parallel computing, but do not specify num of workers
    if matlabpool('size')==0
       matlabpool
    end
elseif wnum<=8
    if matlabpool('size')==0
       matlabpool('open',int32(wnum));
    end
else
    error('GTdef ERROR: Matlabpool input is wrong!!!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'.......... reading the input file ...........\t');
tic
[ coord,smooth,surf,beta,rigidity,poisson,...
  earth,edgrn,layer,...
  flt1,flt2,flt3,flt4,flt5,...
  bndry,subflt,dip,...
  pnt,bsl,prf,grd,...
  sspnt,ssflt1,ssflt2 ] = GTdef_open(finName);
toc

%basename = strtok(finName,'.');	% noly works for names without "."
cellname = regexp(finName,'\.in','split');
basename = char(cellname(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set origin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the middle point of the region as origin for the local cartesian coordinate
if strcmp(coord,'geo') || strcmp(coord,'geo_polyconic')
    lonlist = []; latlist = [];
    if flt1.num~=0
        lonlist = [ lonlist;flt1.flt(:,1) ]; latlist = [ latlist;flt1.flt(:,2) ];
    end
    if flt2.num~=0
        lonlist = [ lonlist;flt2.flt(:,1) ]; latlist = [ latlist;flt2.flt(:,2) ];
    end
    if flt3.num~=0
        lonlist = [ lonlist;flt3.flt(:,1) ]; latlist = [ latlist;flt3.flt(:,2) ];
    end
    if flt4.num~=0
        lonlist = [ lonlist;flt4.flt(:,1) ]; latlist = [ latlist;flt4.flt(:,2) ];
    end
    lon0 = 0.5*(min(lonlist)+max(lonlist));
    lat0 = 0.5*(min(latlist)+max(latlist));
elseif strcmpi(coord,'local')~=1
    error('GTdef ERROR: Coordinate input is wrong!!!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% layered earth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if layered model is used, green function library is built up here
edgrnfcts = [];
if strcmp(earth,'layered')
fprintf(1,'\n..... calculating point source library ......\t');
tic
    %---------- only need to create green functions once ----------
    %fedgrnName = [ basename '_edgrn.inp' ];
    %GTdef_write_edgrn_input(fedgrnName,edgrn,layer);
    %folderName = 'edgrnfcts';
    %% create green function folder if it does not exist
    %if ~exist(folderName,'dir'), mkdir(folderName); end
    %system(['echo ' fedgrnName ' | /Users/lfeng/matlab/edgrn2.0']);
    %%system(['echo ' fedgrnName ' | ./edgrn2.0']);
    %--------------------------------------------------------------

    % read in point source green functions
    [ edgrn,edgrnfcts ] = GTdef_read_edgrn_output(edgrn);
toc
end

pnt.crt = []; pnt.obs = []; pnt.obs_err = []; pnt.obs_wgt = []; pnt.coef = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% point data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pnt.num~=0
fprintf(1,'\n........... processing point data ...........\t');
tic
    % convert point data from geographic to local cartesian coordinate
    if strcmp(coord,'geo')
       [pxx,pyy] = LL2ckmd(pnt.loc(:,1),pnt.loc(:,2),lon0,lat0,0);
    end
    if strcmp(coord,'geo_polyconic')
       [pxx,pyy] = latlon_to_xy(pnt.loc(:,1),pnt.loc(:,2),lon0,lat0);
    end
    if strcmp(coord,'local')
       pxx = pnt.loc(:,1); pyy = pnt.loc(:,2);
    end
    pzz = pnt.loc(:,3);
    zz_ind = pzz>0; pzz(zz_ind) = 0;                % positive depths are all set to be zero   
    pnt.crt = [pxx'; pyy'; pzz'];                   % cartesian - 3*n matrix [xx;yy;zz]; it is just Xin
    % prepare the point observation data
    pnt.obs = reshape(pnt.disp,[],1);			    % (3*n)*1 observation vector [east;north;vertical]
    pnt.obs_err = reshape(pnt.err,[],1);		    % (3*n)*1 error vector [east;north;vertical]
    pnt.obs_wgt = [pnt.wgt;pnt.wgt;pnt.wgt]; 		% (3*n)*1 weight vector [east;north;vertical]
    pnt.coef = sqrt(pnt.obs_wgt)./pnt.obs_err;
toc
end

bsl.crt = []; bsl.obs = []; bsl.obs_err = []; bsl.obs_wgt = []; bsl.coef = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% baseline data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bsl.num~=0
fprintf(1,'\n.......... processing baseline data .........\t');
tic
    % convert from geographic to local cartesian coordinate
    if strcmp(coord,'geo')
       [bx1,by1] = LL2ckmd(bsl.loc(:,1),bsl.loc(:,2),lon0,lat0,0);
       [bx2,by2] = LL2ckmd(bsl.loc(:,4),bsl.loc(:,5),lon0,lat0,0);
    end
    if strcmp(coord,'geo_polyconic')
       [bx1,by1] = latlon_to_xy(bsl.loc(:,1),bsl.loc(:,2),lon0,lat0);
       [bx2,by2] = latlon_to_xy(bsl.loc(:,4),bsl.loc(:,5),lon0,lat0);
    end
    if strcmp(coord,'local')
       bx1 = bsl.loc(:,1); by1  = bsl.loc(:,2);
       bx2 = bsl.loc(:,4); by2  = bsl.loc(:,5);
    end
    bz1 = bsl.loc(:,3);	 bz2 = bsl.loc(:,6);
    bsl.crt = [bx1'; by1'; bz1';bx2'; by2'; bz2'];      % cartesian - 6*n matrix [bx1;by1;bz1;bx2;by2;bz2]; it is just Bin
    % prepare the baseline observation data
    bsl.obs = reshape(bsl.disp,[],1);                   % (4*n)*1 observation vector [east;north;vertical:length]
    bsl.obs_err = reshape(bsl.err,[],1);                % (4*n)*1 error vector [east;north;vertical:length]
    bsl.obs_wgt = [bsl.wgt;bsl.wgt;bsl.wgt;bsl.wgt];	% (4*n)*1 weight vector [east;north;vertical:length]
    bsl.coef = sqrt(bsl.obs_wgt)./bsl.obs_err;
toc
end

nod.loc = []; nod.crt = []; nod.lon = []; nod.lat = []; nod.name = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% profile & grid data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if prf.num~=0||grd.num~=0
fprintf(1,'\n....... processing profile & grid data ......\t');
tic
    if prf.num~=0
        for ii = 1:prf.num
            plon = []; plat = []; pname = {};
    	[ plon,plat,pname ] = GTdef_profile(prf.prf(ii,:),prf.name{ii});
    	nod.lon = [ nod.lon plon ]; nod.lat = [ nod.lat plat ];
        	nod.name = [ nod.name; pname ];   
        end
    end
    if grd.num~=0
        for ii = 1:grd.num
            glon = []; glat = []; gname = {};
    	[ glon,glat,gname ] = GTdef_grid(grd.grd(ii,:),grd.name{ii});
    	nod.lon = [ nod.lon glon ]; nod.lat = [ nod.lat glat ];
        	nod.name = [ nod.name; gname ];   
        end
    end
    nod_zz = nan(length(nod.lon),1);
    nod.loc = [ nod.lon' nod.lat' nod_zz ];
    if strcmp(coord,'geo')
       [nod_xx,nod_yy] = LL2ckmd(nod.lon,nod.lat,lon0,lat0,0);
    end
    if strcmp(coord,'geo_polyconic')
       [nod_xx,nod_yy] = latlon_to_xy(nod.lon,nod.lat,lon0,lat0);
    end
    if strcmp(coord,'local')
       nod_xx = nod.lon; nod_yy = nod.lat;
    end
    nod_zz = zeros(1,length(nod_xx));
    nod.crt = [ nod_xx;nod_yy;nod_zz ];			% cartesian - 3*n matrix [xx;yy;zz]; it is just Nin
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stress point data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sspnt.num~=0
fprintf(1,'\n....... processing stress point data ........\t');
tic
    % convert point data from geographic to local cartesian coordinate
    if strcmp(coord,'geo')
       [sspxx,sspyy] = LL2ckmd(sspnt.loc(:,1),sspnt.loc(:,2),lon0,lat0,0);
    end
    if strcmp(coord,'geo_polyconic')
       [sspxx,sspyy] = latlon_to_xy(sspnt.loc(:,1),sspnt.loc(:,2),lon0,lat0);
    end
    if strcmp(coord,'local')
       sspxx = sspnt.loc(:,1); sspyy = sspnt.loc(:,2);
    end
    sspzz = sspnt.loc(:,3);
    zz_ind = sspzz>0;  sspzz(zz_ind) = 0;               % positive depths are all set to be zero
    sspnt.crt = [sspxx'; sspyy'; sspzz'];               % cartesian - 3*n matrix; it is just Xin
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stress fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ssflt1.fltnum~=0
fprintf(1,'\n....... processing stress fault type-1 ......\t');
tic
    if strcmp(coord,'geo')
       [x1,y1] = LL2ckmd(ssflt1.flt(:,1),ssflt1.flt(:,2),lon0,lat0,0);
    end
    if strcmp(coord,'geo_polyconic')
       [x1,y1] = latlon_to_xy(ssflt1.flt(:,1),ssflt1.flt(:,2),lon0,lat0);
    end
    if strcmp(coord,'local')
       x1 = ssflt1.flt(:,1); y1 = ssflt1.flt(:,2);
    end
    ssflt1.flt = [x1 y1 ssflt1.flt(:,3:end)];
    [ ssflt1 ] = GTdef_stressfault1(ssflt1,dip);
    if strcmp(coord,'geo')
       [lon1,lat1] = ckm2LLd(ssflt1.crt(1,:),ssflt1.crt(2,:),lon0,lat0,0);
       ssflt1.loc = [ lon1; lat1; ssflt1.crt(3,:) ]';
    end
    if strcmp(coord,'geo_polyconic')
       [lon1,lat1] = xy_to_latlon(ssflt1.crt(1,:),ssflt1.crt(2,:),lon0,lat0);
       ssflt1.loc = [ lon1; lat1; ssflt1.crt(3,:) ]';
    end
    if strcmp(coord,'local')
       ssflt1.loc = ssflt1.crt';
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stress fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ssflt2.fltnum~=0
fprintf(1,'\n....... processing stress fault type-2 ......\t');
tic
    if strcmp(coord,'geo')
       [x2_1,y2_1] = LL2ckmd(ssflt2.flt(:,1),ssflt2.flt(:,2),lon0,lat0,0);
       [x2_2,y2_2] = LL2ckmd(ssflt2.flt(:,3),ssflt2.flt(:,4),lon0,lat0,0);
    end
    if strcmp(coord,'geo_polyconic')
       [x2_1,y2_1] = latlon_to_xy(ssflt2.flt(:,1),flt2.flt(:,2),lon0,lat0);
       [x2_2,y2_2] = latlon_to_xy(ssflt2.flt(:,3),flt2.flt(:,4),lon0,lat0);
    end
    if strcmp(coord,'local')
       x2_1 = ssflt2.flt(:,1); y2_1 = ssflt2.flt(:,2);
       x2_2 = ssflt2.flt(:,3); y2_2 = ssflt2.flt(:,4);
    end
    ssflt2.flt = [x2_1 y2_1 x2_2 y2_2 ssflt2.flt(:,5:end)];
    [ ssflt2 ] = GTdef_stressfault2(ssflt2,dip);
    if strcmp(coord,'geo')
       [lon2,lat2] = ckm2LLd(ssflt2.crt(1,:),ssflt2.crt(2,:),lon0,lat0,0);
       ssflt2.loc = [ lon2; lat2; ssflt2.crt(3,:) ]';
    end
    if strcmp(coord,'geo_polyconic')
       [lon2,lat2] = xy_to_latlon(ssflt2.crt(1,:),ssflt2.crt(2,:),lon0,lat0);
       ssflt2.loc = [ lon2; lat2; ssflt2.crt(3,:) ]';
    end
    if strcmp(coord,'local')
       ssflt2.loc = ssflt2.crt';
    end
toc
end

% form everything that is needed for x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)
Xgrn = []; Bgrn = []; Ngrn = []; sm = []; Aeq = []; beq = []; lb = []; ub = []; x0 = [];
% sm_abs for calculate absolute 1st derivative (strain)
sm_abs = [];
% stress calculation
Min0 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt1.num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    if strcmp(coord,'geo')
       [x1,y1] = LL2ckmd(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
    end
    if strcmp(coord,'geo_polyconic')
       [x1,y1] = latlon_to_xy(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0);
    end
    if strcmp(coord,'local')
       x1 = flt1.flt(:,1); y1 = flt1.flt(:,2);
    end
    for ii = 1:flt1.num
        cflt_name = flt1.name{ii};
    	% find subfaults for the master fault
    	sub_ind = strcmpi(cflt_name,subflt.name);
        % find dips for the master fault
    	dip_ind = strcmpi(cflt_name,dip.name);
        [ Xgrn1,Bgrn1,Ngrn1,sm1,sm1_abs,Aeq1,beq1,lb1,ub1,x01,Min01 ] = ...
          GTdef_fault1dif([x1(ii) y1(ii) flt1.flt(ii,3:end)],...
	                 subflt.flt(sub_ind,:),dip.dip(dip_ind,:),pnt.crt,bsl.crt,nod.crt,...
	                 earth,rigidity,poisson,edgrn,edgrnfcts,smooth,surf);
    	[ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = GTdef_addall(Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0,...
	                                                Xgrn1,Bgrn1,Ngrn1,sm1,sm1_abs,Aeq1,beq1,lb1,ub1,x01);
        Min0 = [ Min0; Min01 ];
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt2.num~=0
fprintf(1,'\n.......... processing fault type-2 ..........\t');
tic
    if strcmp(coord,'geo')
       [x2_1,y2_1] = LL2ckmd(flt2.flt(:,1),flt2.flt(:,2),lon0,lat0,0);
       [x2_2,y2_2] = LL2ckmd(flt2.flt(:,3),flt2.flt(:,4),lon0,lat0,0);
    end
    if strcmp(coord,'geo_polyconic')
       [x2_1,y2_1] = latlon_to_xy(flt2.flt(:,1),flt2.flt(:,2),lon0,lat0);
       [x2_2,y2_2] = latlon_to_xy(flt2.flt(:,3),flt2.flt(:,4),lon0,lat0);
    end
    if strcmp(coord,'local')
       x2_1 = flt2.flt(:,1); y2_1 = flt2.flt(:,2);
       x2_2 = flt2.flt(:,3); y2_2 = flt2.flt(:,4);
    end
    for ii = 1:flt2.num
        cflt_name = flt2.name{ii};
    	% find subfaults for the master fault
    	sub_ind = strcmpi(cflt_name,subflt.name);
	% find dips for the master fault
    	dip_ind = strcmpi(cflt_name,dip.name);
        [ Xgrn2,Bgrn2,Ngrn2,sm2,sm2_abs,Aeq2,beq2,lb2,ub2,x02,Min02 ] = ...
	GTdef_fault2dif([x2_1(ii) y2_1(ii) x2_2(ii) y2_2(ii) flt2.flt(ii,5:end)],...
	                subflt.flt(sub_ind,:),dip.dip(dip_ind,:),pnt.crt,bsl.crt,nod.crt,...
		        earth,rigidity,poisson,edgrn,edgrnfcts,smooth,surf);
    	[ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = GTdef_addall(Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0,...
	                                                Xgrn2,Bgrn2,Ngrn2,sm2,sm2_abs,Aeq2,beq2,lb2,ub2,x02);
        Min0 = [ Min0; Min02 ];
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt3.num~=0
fprintf(1,'\n.......... processing fault type-3 ..........\t');
tic
    if strcmp(coord,'geo')
       [x3,y3] = LL2ckmd(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0,0);
    end
    if strcmp(coord,'geo_polyconic')
       [x3,y3] = latlon_to_xy(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0);
    end
    if  strcmp(coord,'local')
       x3 = flt3.flt(:,1); y3 = flt3.flt(:,2);
    end
    for ii = 1:flt3.num
        cflt_name = flt3.name{ii};
    	% find subfaults for the master fault
    	sub_ind = strcmpi(cflt_name,subflt.name);
	% find dips for the master fault
    	dip_ind = strcmpi(cflt_name,dip.name);
        [ Xgrn3,Bgrn3,Ngrn3,sm3,sm3_abs,Aeq3,beq3,lb3,ub3,x03,Min03 ] = ...
	GTdef_fault3dif([x3(ii) y3(ii) flt3.flt(ii,3:end)],...
	                subflt.flt(sub_ind,:),dip.dip(dip_ind,:),pnt.crt,bsl.crt,nod.crt,...
	                earth,rigidity,poisson,edgrn,edgrnfcts,smooth,surf);
    	[ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = GTdef_addall(Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0,...
	                                                Xgrn3,Bgrn3,Ngrn3,sm3,sm3_abs,Aeq3,beq3,lb3,ub3,x03);
        Min0 = [ Min0; Min03 ];
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault4 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt4.num~=0
fprintf(1,'\n.......... processing fault type-4 ..........\t');
tic
    if strcmp(coord,'geo')
       [x4_1,y4_1] = LL2ckmd(flt4.flt(:,1),flt4.flt(:,2),lon0,lat0,0);
       [x4_2,y4_2] = LL2ckmd(flt4.flt(:,3),flt4.flt(:,4),lon0,lat0,0);
    end
    if strcmp(coord,'geo')
       [x4_1,y4_1] = latlon_to_xy(flt4.flt(:,1),flt4.flt(:,2),lon0,lat0);
       [x4_2,y4_2] = latlon_to_xy(flt4.flt(:,3),flt4.flt(:,4),lon0,lat0);
    end
    if strcmp(coord,'local')
       x4_1 = flt4.flt(:,1); y4_1 = flt4.flt(:,2);
       x4_2 = flt4.flt(:,3); y4_2 = flt4.flt(:,4);
    end
    for ii = 1:flt4.num
        cflt_name = flt4.name{ii};
    	% find subfaults for the master fault
    	sub_ind = strcmpi(cflt_name,subflt.name);
	% find dips for the master fault
    	dip_ind = strcmpi(cflt_name,dip.name);
        [ Xgrn4,Bgrn4,Ngrn4,sm4,sm4_abs,Aeq4,beq4,lb4,ub4,x04,Min04 ] = ...
	GTdef_fault4dif([x4_1(ii) y4_1(ii) x4_2(ii) y4_2(ii) flt4.flt(ii,5:end)],...
	                subflt.flt(sub_ind,:),dip.dip(dip_ind,:),pnt.crt,bsl.crt,nod.crt,...
		        earth,rigidity,poisson,edgrn,edgrnfcts,smooth,surf);
    	[ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = GTdef_addall(Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0,...
	                                                Xgrn4,Bgrn4,Ngrn4,sm4,sm4_abs,Aeq4,beq4,lb4,ub4,x04);
        Min0 = [ Min0; Min04 ];
    end
toc
end
Min0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault5 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt5.num~=0
fprintf(1,'\n.......... processing fault type-5 ..........\t');
tic
        if strcmp(coord,'geo')
           [bx1,by1] = LL2ckmd(bndry.bd(:,3), bndry.bd(:,4), lon0,lat0,0);	% upleft point
           [bx2,by2] = LL2ckmd(bndry.bd(:,6), bndry.bd(:,7), lon0,lat0,0);	% lower left point
           [bx3,by3] = LL2ckmd(bndry.bd(:,9), bndry.bd(:,10),lon0,lat0,0);	% lower right point
           [bx4,by4] = LL2ckmd(bndry.bd(:,12),bndry.bd(:,13),lon0,lat0,0);	% upleft point
           bndry.crt = [ bx1 by1 bndry.bd(:,5) bx2 by2 bndry.bd(:,8) bx3 by3 bndry.bd(:,11) bx4 by4 bndry.bd(:,14) ];
        end
        if strcmp(coord,'geo_polyconic')
           [bx1,by1] = latlon_to_xy(bndry.bd(:,3), bndry.bd(:,4), lon0,lat0);	% upleft point
           [bx2,by2] = latlon_to_xy(bndry.bd(:,6), bndry.bd(:,7), lon0,lat0);	% lower left point
           [bx3,by3] = latlon_to_xy(bndry.bd(:,9), bndry.bd(:,10),lon0,lat0);	% lower right point
           [bx4,by4] = latlon_to_xy(bndry.bd(:,12),bndry.bd(:,13),lon0,lat0);	% upleft point
           bndry.crt = [ bx1 by1 bndry.bd(:,5) bx2 by2 bndry.bd(:,8) bx3 by3 bndry.bd(:,11) bx4 by4 bndry.bd(:,14) ];
        end
        if strcmp(coord,'local')
    	   bndry.crt = bndry.bd;
        end
        for ii = 1:flt5.num
           cflt_name = flt5.name{ii};
           % find subfaults for the master fault
           sub_ind = strcmpi(cflt_name,subflt.name);
           % find the boundary for fault 5
           bnd_ind = strcmpi(cflt_name,bndry.name);
           [ Xgrn5,Bgrn5,Ngrn5,sm5,sm5_abs,Aeq5,beq5,lb5,ub5,x05 ] = ...;
    	   GTdef_fault5(flt5.flt,subflt.flt(sub_ind,:),bndry.bd(bnd_ind,:),pnt.crt,bsl.crt,nod.crt,rigidity,poisson,smooth,surf);
           [ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = ...
    	   GTdef_addall(Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0,...
    	                Xgrn5,Bgrn5,Ngrn5,sm5,sm5_abs,Aeq5,beq5,lb5,ub5,x05);
        end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% forward only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lb==-Inf
fprintf(1,'\n........... doing forward modeling ..........\t');
tic
    [ sm_use ] = GTdef_condense(sm);
    [ sm_abs ] = GTdef_condense(sm_abs);
    foutName = [ basename '_fwd.out' ];
    % mod_info = [ data_num slip_num ndf rss rms wrrs wrms chi2 rchi2 r_1d r_2d ];
    [ mod_info,pnt.out,bsl.out,nod.out ] = GTdef_forward(Xgrn,Bgrn,Ngrn,sm,sm_abs,lb,ub,x0,pnt,bsl,nod,smooth);
    if ~isempty(sspnt) || ~isempty(ssflt1) || ~isempty(ssflt2)
        fstressName = [ basename '_stress.out' ];
        [ sspnt,ssflt1,ssflt2 ] = GTdef_calc_stress(sspnt,ssflt1,ssflt2,Min0,earth,rigidity,poisson,edgrn,layer,edgrnfcts);
        GTdef_output_stress(fstressName,sspnt,ssflt1,ssflt2);
    end
    % forward models do not change slips
    flt1.out = flt1.flt; flt2.out = flt2.flt; flt3.out = flt3.flt; flt4.out = flt4.flt; 
    subflt.out = subflt.flt; subflt.outname = subflt.name;
    GTdef_output(foutName,coord,'none','none',0,rigidity,poisson,earth,edgrn,layer,...
    		 flt1,flt2,flt3,flt4,flt5,bndry,subflt,dip,pnt,bsl,prf,grd,nod,mod_info);
toc
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% inversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n............. doing inversion .............\t');
    % condense the smoothing matrix by removing rows of all zeros
    [ sm_use ] = GTdef_condense(sm);
    [ sm_abs ] = GTdef_condense(sm_abs);

    % do inversion for each beta
    betaNum = length(beta);
    for ii = 1:betaNum
        bt = beta(ii);
        fprintf(1,'\n............. beta = %16.5f .............\t',bt);
        tic
            if strcmp(smooth,'2d')
                kp = sqrt(bt);
                if kp>=1
                    foutName = strcat(basename,'_kp',num2str(kp,'%-.0f'),'.out');
                else
                    foutName = strcat(basename,'_kp',num2str(kp,'%-.5f'),'.out');
                end
            else
                if bt>=1
                    foutName = strcat(basename,'_bt',num2str(bt,'%-.0f'),'.out');
                else
                    foutName = strcat(basename,'_bt',num2str(bt,'%-.5f'),'.out');
                end
            end
            [ xx ] = GTdef_invert(Xgrn,Bgrn,sm_use,Aeq,beq,lb,ub,x0,pnt,bsl,bt);
            % faults info
            [ flt1.out,flt2.out,flt3.out,flt4.out,subflt.out,subflt.outname ] = ...
	    GTdef_slips(lb,ub,xx,flt1,flt2,flt3,flt4,flt5,subflt);
            [ mod_info(ii,:),pnt.out,bsl.out,nod.out ] = ...
	    GTdef_forward(Xgrn,Bgrn,Ngrn,sm_use,sm_abs,lb,ub,xx,pnt,bsl,nod,smooth);
            % output results
            GTdef_output(foutName,coord,smooth,surf,bt,rigidity,poisson,earth,edgrn,layer,...
                         flt1,flt2,flt3,flt4,flt5,bndry,subflt,dip,pnt,bsl,prf,grd,nod,mod_info(ii,:));
        toc
    end
    fsumName = [ basename '_inv.out' ];
    GTdef_summary(fsumName,beta,mod_info);
end

% close up matlabpool for parallel computing
if matlabpool('size')>0
   matlabpool close
end
