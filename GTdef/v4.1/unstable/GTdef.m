function [] = GTdef(finName,wnum)

%   Please follow documentation for instructions at
%      https://avnewman.github.io/GTDef/documentation/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    GTdef.m                                           %
%                  Georgia Tech Matlab program for deformation                         %
%                    Lujia Feng; Andrew V. Newman; Ting Chen                           %
%                                                                                      %
% INPUT:                                                                               %
% finName - input file name                                                            %
% wnum    - num of matlab parallel workers to be used                                  %
%   #: number of workers specified                                                     %
%   0: do not use parallel computing                                                   %
%  -1: use parallel computing                                                          %
%      but number of workers will be determined by the system                          %
%    : if option is excluded, will use open pool, if available                         %
%                                                                                      %
% v1:                                                                                  %
% first created by Lujia Feng Mon Apr 20 14:15:52 EDT 2009                             %
% added coordiante type flag 'coord' lfeng Thu Nov 5 17:12:59 EST 2009                 %
% added first derivative modes lfeng Tue Dec  1 14:31:10 EST 2009                      %
% modified beta (added beta) and roughness for 1st derivatives                         %
% added 'dip' flag for bended faults lfeng Mon Dec  7 01:04:06 EST 2009                %
% added 'freesurface' flag lfeng Wed Dec  9 17:00:58 EST 2009                          %
% added fault type 5 lfeng Fri Dec 11 10:57:18 EST 2009                                %
% changed 'freesurface' to 'surface' flag lfeng Wed Feb 24 12:46:01 EST 2010           %
% changed 'coord' to string flag lfeng Wed Feb 24 13:40:01 EST 2010                    %
% allows input file name include multiple "." besides ".in" lfeng Oct 4 2010           %
% added matlabpool lfeng Wed Dec  1 12:12:00 EST 2010                                  %
% edited the origin definition lfeng Thu Apr  7 18:45:35 EDT 2011                      %
% v2:                                                                                  %
% used strucutres for passing parameters lfeng Wed Feb 22 03:39:11 SGT 2012            %
% added layered earth model lfeng Tue Feb 28 03:44:29 SGT 2012                         %
% new definition used for fault1, fault2, fault3 & fault4 lfeng May 8 2012             %
% added polyconic projection lfeng Thu Jun  7 12:23:26 SGT 2012                        %
% added fault5 lfeng Fri Nov 30 14:52:01 SGT 2012                                      %
% added saving greensfns lfeng Mon Aug  5 14:17:44 SGT 2013                            %
% added 'strike' flag for curved faults lfeng Fri Oct 24 14:47:22 SGT 2014             %
% added addon to combine 'dip' & 'strike' lfeng Fri Oct 24 15:06:09 SGT 2014           %
% added sweepAngle lfeng Wed Nov  5 19:34:02 SGT 2014                                  %
% modified greensfns output lfeng Fri Nov 14 16:16:28 SGT 2014                         %
% corrected addon.crt error with Paul Morgan lfeng Fri Dec 12 16:55:38 SGT 2014        %
% v3:                                                                                  %
% added modspace (model space) structure for lsqlin lfeng Thu Mar 19 17:32:29 SGT 2015 %
% added Min and xyzctr to fault structures lfeng Thu Mar 19 20:18:30 SGT 2015          %
% added earth structure lfeng Fri Mar 20 20:37:32 SGT 2015                             %
% added xyzflt lfeng Tue Mar 24 11:21:45 SGT 2015                                      %
% added Min, SSgrn, DSgrn, and TSgrn to xyzflt lfeng Thu Mar 26 15:54:56 SGT 2015      %
% added modspace.sdropflag lfeng Thu Mar 26 17:31:07 SGT 2015                          %
% added output individual Xgrn lfeng Fri Jun 12 12:21:25 SGT 2015                      %
% added fault5 lfeng with Paul Morgan Wed Jun 17 17:36:52 SGT 2015                     %
% fixed fault4 bug lfeng with Paul Morgan Wed Jun 24 12:31:10 SGT 2015                 %
% fixed profile & grid lfeng Mon Jul 27 16:14:16 SGT 2015                              %
% 3D geometry for Okada models can be imported using fault5 lfeng Tue Aug  4 SGT 2015  %
% fixed greensfns for fault5 lfeng Wed Aug  5 15:57:54 SGT 2015                        %
% added InSAR los lfeng Tue Nov  3 10:47:25 SGT 2015                                   %
% added Matlab equivalent of edgrn lfeng Thu Feb  4 14:51:26 SGT 2016                  %
% replaced matlabpool with parpool for newer Matlab AVN Tue Apr 19 15:27:25 EDT 2016   %
% added external geometry to fault3 for fault5, rename old fault6 to fault7 lfeng 2016 %
% added output resolution matrix information (see GTdef_input) anewman May 10 2016     %
% added optional .mat file output (see GTdef_input) anewman May 18 17:32:55 UTC 2016   %
% corrected fault type 7 bugs lfeng Fri Jun 10 01:04:06 SGT 2016                       %
% added project option for runs. AVN: Thu Jun 18 20:25:58 EDT 2020                     %
% last modified Andrew Newman Thu Jun 18 20:25:58 EDT 2020                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% specify matlabpool for parallel computing %%%%%%%%%%%%%%%%%%%%%%%%%%%
% updated for use by parpool. AVN 5/3/16
% depending on your local setup, parpool may be opened automatically due to
% parfor loops within the code.
if exist('wnum')
   if wnum>0
     localpool = parpool(int32(wnum));
   elseif wnum<0
     localpool = parpool;
   else
     fprintf(1,'GTdef WARNING: parpool is not used.');
   end
else
   fprintf(1,'GTdef WARNING: Will use parpool if already running.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n.......... reading the input file ...........\t');
tic
[ modspace,earth,...
  flt1,flt2,flt3,flt4,flt5,flt6,flt7,...
  subflt,addon,...
  pnt,los,bsl,prf,grd,...
  sspnt,ssflt1,ssflt2 ] = GTdef_open(finName);
toc

lat0=0; lon0=0;   % needs to be defined even if using local coords. AVN 7/13/17
%basename = strtok(finName,'.');	% noly works for names without "."
%cellname = regexp(finname,'\.(in|out)','split');
%basename = char(cellname(1));
[ ~,basename,~ ] = fileparts(finName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set origin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the middle point of the region as origin for the local cartesian coordinate
if strcmp(modspace.coord,'geo') || strcmp(modspace.coord,'geo_polyconic')
    if isempty(modspace.origin)
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
    else
        lon0 = modspace.origin(1); lat0 = modspace.origin(2);
    end
elseif strcmpi(modspace.coord,'local')~=1
    error('GTdef ERROR: coordinate input is wrong!!!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% layered earth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(earth.type,'layered')
fprintf(1,'\n..... calculating point source library ......\t');
tic
   folderName = 'edgrnfcts';
   fssName = [ folderName '/izmhs.ss' ];
   fdsName = [ folderName '/izmhs.ds' ];
   fclName = [ folderName '/izmhs.cl' ];
   if exist(fssName,'file') && exist(fdsName,'file') && exist(fclName,'file')
      % read in point source green functions
      [ earth.edgrn,earth.edgrnfcts ] = GTdef_read_edgrn_output(earth.edgrn);
   else
      %---------- option 1: using Fortran code (edcmp/edgrn) to create green functions ----------
      % only need to run once
      %fedgrnName = [ basename '_edgrn.inp' ];
      %GTdef_write_edgrn_input(fedgrnName,edgrn,layer);
      %folderName = 'edgrnfcts';
      %% create green function folder if it does not exist
      %if ~exist(folderName,'dir'), mkdir(folderName); end
      %system(['echo ' fedgrnName ' | /Users/lfeng/matlab/edgrn2.0']);
      %%system(['echo ' fedgrnName ' | ./edgrn2.0']);

      %---------- option 2: using Matlab code to create green functions (recommended) ----------
      [ earth ] = GTdef_edgrn(earth);
   end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% point data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pnt.num~=0
fprintf(1,'\n........... processing point data ...........\t');
tic
    % convert point data from geographic to local cartesian coordinate
    switch modspace.coord
       case 'geo'
          [pxx,pyy] = LL2ckmd(pnt.loc(:,1),pnt.loc(:,2),lon0,lat0,0);
       case 'geo_polyconic'
          [pxx,pyy] = latlon_to_xy(pnt.loc(:,1),pnt.loc(:,2),lon0,lat0,0);
       case 'local'
          pxx = pnt.loc(:,1); pyy = pnt.loc(:,2);
    end
    pzz = pnt.loc(:,3);
    zz_ind = pzz>0; pzz(zz_ind) = 0;                % positive depths are all set to be zero
    pnt.crt = [pxx pyy pzz];                        % cartesian - n*3 matrix [xx yy zz]; it is just Xin'
    % prepare the point observation data
    pnt.obs = reshape(pnt.disp,[],1);               % (3*n)*1 observation vector [east;north;vertical]
    pnt.obs_err = reshape(pnt.err,[],1);            % (3*n)*1 error vector [east;north;vertical]
    pnt.obs_wgt = [pnt.wgt;pnt.wgt;pnt.wgt];        % (3*n)*1 weight vector [east;north;vertical]
    %pnt.coef = pnt.obs_wgt./pnt.obs_err.^2;        % (3*n)*1 coefficient vector  # v. used in fork by A. Williamson
    pnt.coef = sqrt(pnt.obs_wgt)./pnt.obs_err;      % (3*n)*1 coefficient vector
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% InSAR los data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if los.num~=0
fprintf(1,'\n........... processing point data ...........\t');
tic
    % convert InSAR los data from geographic to local cartesian coordinate
    switch modspace.coord
       case 'geo'
          [pxx,pyy] = LL2ckmd(los.loc(:,1),los.loc(:,2),lon0,lat0,0);
       case 'geo_polyconic'
          [pxx,pyy] = latlon_to_xy(los.loc(:,1),los.loc(:,2),lon0,lat0,0);
       case 'local'
          pxx = los.loc(:,1); pyy = los.loc(:,2);
    end
    pzz = los.loc(:,3);
    zz_ind = pzz>0; pzz(zz_ind) = 0;                % positive depths are all set to be zero
    los.crt = [ pxx pyy pzz los.dir ];              % cartesian - n*3 matrix [xx yy zz] & add los dir from ground to satellite
    % prepare the los observation data
    los.obs     = los.disp;                         % (1*n)*1 observation vector [los]
    los.obs_err = los.err;                          % (1*n)*1 error vector [los]
    los.obs_wgt = los.wgt;                          % (1*n)*1 weight vector [los]
    %los.coef    = los.obs_wgt./los.obs_err.^2;     % (1*n)*1 coefficient vector
    los.coef    = sqrt(los.obs_wgt)./los.obs_err;   % (1*n)*1 coefficient vector
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% baseline data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bsl.num~=0
fprintf(1,'\n.......... processing baseline data .........\t');
tic
    % convert from geographic to local cartesian coordinate
    switch modspace.coord
       case 'geo'
          [bx1,by1] = LL2ckmd(bsl.loc(:,1),bsl.loc(:,2),lon0,lat0,0);
          [bx2,by2] = LL2ckmd(bsl.loc(:,4),bsl.loc(:,5),lon0,lat0,0);
       case 'geo_polyconic'
          [bx1,by1] = latlon_to_xy(bsl.loc(:,1),bsl.loc(:,2),lon0,lat0,0);
          [bx2,by2] = latlon_to_xy(bsl.loc(:,4),bsl.loc(:,5),lon0,lat0,0);
       case 'local'
          bx1 = bsl.loc(:,1); by1 = bsl.loc(:,2);
          bx2 = bsl.loc(:,4); by2 = bsl.loc(:,5);
    end
    bz1 = bsl.loc(:,3);	 bz2 = bsl.loc(:,6);
    bsl.crt = [bx1 by1 bz1 bx2 by2 bz2];                % cartesian - n*6 matrix [bx1 by1 bz1 bx2 by2 bz2]; it is just Bin'
    % prepare the baseline observation data
    bsl.obs = reshape(bsl.disp,[],1);                   % (4*n)*1 observation vector [east;north;vertical:length]
    bsl.obs_err = reshape(bsl.err,[],1);                % (4*n)*1 error vector [east;north;vertical:length]
    bsl.obs_wgt = [bsl.wgt;bsl.wgt;bsl.wgt;bsl.wgt];	% (4*n)*1 weight vector [east;north;vertical:length]
    %bsl.coef = bsl.obs_wgt./bsl.obs_err.^2;            % (4*n)*1 coefficient vector
    bsl.coef = sqrt(bsl.obs_wgt)./bsl.obs_err;          % (4*n)*1 coefficient vector
toc
end

nod.loc = []; nod.crt = []; nod.name = {};
nod_lon = []; nod_lat = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% profile & grid data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if prf.num~=0||grd.num~=0
fprintf(1,'\n....... processing profile & grid data ......\t');
tic
    if prf.num~=0
        for ii = 1:prf.num
            [ plon,plat,pname ] = GTdef_profile(prf.prf(ii,:),prf.name{ii});
            nod_lon  = [ nod_lon; plon ];
            nod_lat  = [ nod_lat; plat ];
            nod.name = [ nod.name; pname ];
        end
    end
    if grd.num~=0
        for ii = 1:grd.num
    	    [ glon,glat,gname ] = GTdef_grid(grd.grd(ii,:),grd.name{ii});
    	    nod_lon  = [ nod_lon; glon ];
	    nod_lat  = [ nod_lat; glat ];
            nod.name = [ nod.name; gname ];
        end
    end
    nod_zz = nan(size(nod_lon));
    nod.loc = [ nod_lon nod_lat nod_zz ];
    switch modspace.coord
       case 'geo'
          [nod_xx,nod_yy] = LL2ckmd(nod_lon,nod_lat,lon0,lat0,0);
       case 'geo_polyconic'
          [nod_xx,nod_yy] = latlon_to_xy(nod_lon,nod_lat,lon0,lat0,0);
       case 'local'
          nod_xx = nod_lon; nod_yy = nod_lat;
    end
    nod_zz = zeros(size(nod_xx));
    nod.crt = [ nod_xx nod_yy nod_zz ];			% cartesian - n*3 matrix [xx yy zz]; it is just Nin'
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stress point data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sspnt.num~=0
fprintf(1,'\n....... processing stress point data ........\t');
tic
    % convert point data from geographic to local cartesian coordinate
    switch modspace.coord
       case 'geo'
          [sspxx,sspyy] = LL2ckmd(sspnt.loc(:,1),sspnt.loc(:,2),lon0,lat0,0);
       case 'geo_polyconic'
          [sspxx,sspyy] = latlon_to_xy(sspnt.loc(:,1),sspnt.loc(:,2),lon0,lat0,0);
       case 'local'
          sspxx = sspnt.loc(:,1); sspyy = sspnt.loc(:,2);
    end
    sspzz = sspnt.loc(:,3);
    zz_ind = sspzz>0;  sspzz(zz_ind) = 0;               % positive depths are all set to be zero
    sspnt.crt = [sspxx'; sspyy'; sspzz'];               % cartesian - 3*n matrix; it is just Xin
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% strike variations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if addon.strnum~=0
fprintf(1,'\n....... processing strike variations .......\t');
tic
    % convert strike controlling points from geographic to local cartesian coordinate
    switch modspace.coord
       case 'geo'
          [sx1,sy1] = LL2ckmd(addon.str(:,1),addon.str(:,2),lon0,lat0,0);
          [sx2,sy2] = LL2ckmd(addon.str(:,3),addon.str(:,4),lon0,lat0,0);
       case 'geo_polyconic'
          [sx1,sy1] = latlon_to_xy(addon.str(:,1),addon.str(:,2),lon0,lat0,0);
          [sx2,sy2] = latlon_to_xy(addon.str(:,3),addon.str(:,4),lon0,lat0,0);
       case 'local'
          sx1 = addon.str(:,1); sy1 = addon.str(:,2);
          sx2 = addon.str(:,3); sy2 = addon.str(:,4);
    end
    addon.crt = [ sx1 sy1 sx2 sy2 addon.str(:,[5 6]) ];
toc
else
    addon.crt = addon.str;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stress fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ssflt1.fltnum~=0
fprintf(1,'\n....... processing stress fault type-1 ......\t');
tic
    switch modspace.coord
       case 'geo'
          [x1,y1] = LL2ckmd(ssflt1.flt(:,1),ssflt1.flt(:,2),lon0,lat0,0);
       case 'geo_polyconic'
          [x1,y1] = latlon_to_xy(ssflt1.flt(:,1),ssflt1.flt(:,2),lon0,lat0,0);
       case 'local'
          x1 = ssflt1.flt(:,1); y1 = ssflt1.flt(:,2);
    end
    ssflt1.flt = [x1 y1 ssflt1.flt(:,3:end)];
    [ ssflt1 ] = GTdef_stressfault1(ssflt1,addon);
    switch modspace.coord
       case 'geo'
          [lon1,lat1] = ckm2LLd(ssflt1.crt(1,:),ssflt1.crt(2,:),lon0,lat0,0);
          ssflt1.loc  = [ lon1; lat1; ssflt1.crt(3,:) ]';
       case 'geo_polyconic'
          [lon1,lat1] = xy_to_latlon(ssflt1.crt(1,:),ssflt1.crt(2,:),lon0,lat0,0);
          ssflt1.loc  = [ lon1; lat1; ssflt1.crt(3,:) ]';
       case 'local'
          ssflt1.loc  = ssflt1.crt';
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stress fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ssflt2.fltnum~=0
fprintf(1,'\n....... processing stress fault type-2 ......\t');
tic
    switch modspace.coord
       case 'geo'
          [x2_1,y2_1] = LL2ckmd(ssflt2.flt(:,1),ssflt2.flt(:,2),lon0,lat0,0);
          [x2_2,y2_2] = LL2ckmd(ssflt2.flt(:,3),ssflt2.flt(:,4),lon0,lat0,0);
       case 'geo_polyconic'
          [x2_1,y2_1] = latlon_to_xy(ssflt2.flt(:,1),flt2.flt(:,2),lon0,lat0,0);
          [x2_2,y2_2] = latlon_to_xy(ssflt2.flt(:,3),flt2.flt(:,4),lon0,lat0,0);
       case 'local'
          x2_1 = ssflt2.flt(:,1); y2_1 = ssflt2.flt(:,2);
          x2_2 = ssflt2.flt(:,3); y2_2 = ssflt2.flt(:,4);
    end
    ssflt2.flt = [x2_1 y2_1 x2_2 y2_2 ssflt2.flt(:,5:end)];
    [ ssflt2 ] = GTdef_stressfault2(ssflt2,addon);
    switch modspace.coord
       case 'geo'
          [lon2,lat2] = ckm2LLd(ssflt2.crt(1,:),ssflt2.crt(2,:),lon0,lat0,0);
          ssflt2.loc = [ lon2; lat2; ssflt2.crt(3,:) ]';
       case 'geo_polyconic'
          [lon2,lat2] = xy_to_latlon(ssflt2.crt(1,:),ssflt2.crt(2,:),lon0,lat0,0);
          ssflt2.loc = [ lon2; lat2; ssflt2.crt(3,:) ]';
       case 'local'
          ssflt2.loc = ssflt2.crt';
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt1.num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    switch modspace.coord
        case 'geo'
           [x1,y1] = LL2ckmd(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
        case 'geo_polyconic'
           [x1,y1] = latlon_to_xy(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
        case 'local'
           x1 = flt1.flt(:,1); y1 = flt1.flt(:,2);
    end
    newflt1 = [ x1 y1 flt1.flt(:,3:end) ];
    for ii = 1:flt1.num
        cfname = flt1.name{ii};
        cflt   = newflt1(ii,:);
    	% find subfaults for the master fault
    	subInd = strcmpi(cfname,subflt.name);
        % find dips for the master fault
    	dipInd = strcmpi(cfname,addon.dipname);
        [ modspace,xyzflt,Xgrn1 ] = GTdef_fault1dif(modspace,cflt,subflt.flt(subInd,:),addon.dip(dipInd,:),...
	                                            pnt.crt,los.crt,bsl.crt,nod.crt,earth);
        flt1.xyzflt{ii} = xyzflt;
        % save green's functions
        if strcmpi(modspace.grnflag,'on')
           [ ~,prjflt1,~ ] = GTdef_prjflt1dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:));
           GTdef_save_greensfns(cfname,pnt,prjflt1,Xgrn1,modspace.coord,lon0,lat0,0);
        end
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt2.num~=0
fprintf(1,'\n.......... processing fault type-2 ..........\t');
tic
    switch modspace.coord
        case 'geo'
           [x2_1,y2_1] = LL2ckmd(flt2.flt(:,1),flt2.flt(:,2),lon0,lat0,0);
           [x2_2,y2_2] = LL2ckmd(flt2.flt(:,3),flt2.flt(:,4),lon0,lat0,0);
        case 'geo_polyconic'
           [x2_1,y2_1] = latlon_to_xy(flt2.flt(:,1),flt2.flt(:,2),lon0,lat0,0);
           [x2_2,y2_2] = latlon_to_xy(flt2.flt(:,3),flt2.flt(:,4),lon0,lat0,0);
        case 'local'
           x2_1 = flt2.flt(:,1); y2_1 = flt2.flt(:,2);
           x2_2 = flt2.flt(:,3); y2_2 = flt2.flt(:,4);
    end
    newflt2 = [ x2_1 y2_1 x2_2 y2_2 flt2.flt(:,5:end) ];
    for ii = 1:flt2.num
        cfname = flt2.name{ii};
        cflt   = newflt2(ii,:);
    	% find subfaults for the master fault
    	subInd = strcmpi(cfname,subflt.name);
        % find dips for the master fault
    	dipInd = strcmpi(cfname,addon.dipname);
        % find strikes for the master fault
    	strInd = strcmpi(cfname,addon.strname);
        [ modspace,xyzflt,Xgrn2 ] = GTdef_fault2dif(modspace,cflt,subflt.flt(subInd,:),addon.dip(dipInd,:),addon.crt(strInd,:),...
                                                    pnt.crt,los.crt,bsl.crt,nod.crt,earth);
        flt2.xyzflt{ii} = xyzflt;
        % save green's functions
        if strcmpi(modspace.grnflag,'on')
            [ ~,prjflt2,~ ] = GTdef_prjflt2dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:),addon.crt(strInd,:));
            GTdef_save_greensfns(cfname,pnt,prjflt2,Xgrn2,modspace.coord,lon0,lat0,0);
        end
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt3.num~=0
fprintf(1,'\n.......... processing fault type-3 ..........\t');
tic
    switch modspace.coord
        case 'geo'
           [x3,y3] = LL2ckmd(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0,0);
        case 'geo_polyconic'
           [x3,y3] = latlon_to_xy(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0,0);
        case 'local'
           x3 = flt3.flt(:,1); y3 = flt3.flt(:,2);
    end
    newflt3 = [x3 y3 flt3.flt(:,3:end)];
    for ii = 1:flt3.num
        cfname = flt3.name{ii};
        cflt   = newflt3(ii,:);
    	% find subfaults for the master fault
    	subInd = strcmpi(cfname,subflt.name);
        % find dips for the master fault
    	dipInd = strcmpi(cfname,addon.dipname);
        [ modspace,xyzflt,Xgrn3 ] = GTdef_fault3dif(modspace,cflt,subflt.flt(subInd,:),addon.dip(dipInd,:),...
                                                    pnt.crt,los.crt,bsl.crt,nod.crt,earth);
        flt3.xyzflt{ii} = xyzflt;
        % save green's functions
        if strcmpi(modspace.grnflag,'on')
            [ ~,prjflt3,~ ] = GTdef_prjflt3dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:));
            GTdef_save_greensfns(cfname,pnt,prjflt3,Xgrn3,modspace.coord,lon0,lat0,0);
        end
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault4 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt4.num~=0
fprintf(1,'\n.......... processing fault type-4 ..........\t');
tic
    switch modspace.coord
       case 'geo'
          [x4_1,y4_1] = LL2ckmd(flt4.flt(:,1),flt4.flt(:,2),lon0,lat0,0);
          [x4_2,y4_2] = LL2ckmd(flt4.flt(:,3),flt4.flt(:,4),lon0,lat0,0);
       case 'geo_polyconic'
          [x4_1,y4_1] = latlon_to_xy(flt4.flt(:,1),flt4.flt(:,2),lon0,lat0,0);
          [x4_2,y4_2] = latlon_to_xy(flt4.flt(:,3),flt4.flt(:,4),lon0,lat0,0);
       case 'local'
          x4_1 = flt4.flt(:,1); y4_1 = flt4.flt(:,2);
          x4_2 = flt4.flt(:,3); y4_2 = flt4.flt(:,4);
    end
    newflt4 = [x4_1 y4_1 x4_2 y4_2 flt4.flt(:,5:end)];
    for ii = 1:flt4.num
        cfname = flt4.name{ii};
        cflt   = newflt4(ii,:);
    	% find subfaults for the master fault
    	subInd = strcmpi(cfname,subflt.name);
        % find dips for the master fault
    	dipInd = strcmpi(cfname,addon.dipname);
        % find strikes for the master fault
    	strInd = strcmpi(cfname,addon.strname);
        [ modspace,xyzflt,Xgrn4 ] = GTdef_fault4dif(modspace,cflt,subflt.flt(subInd,:),addon.dip(dipInd,:),addon.crt(strInd,:),...
                                                    pnt.crt,los.crt,bsl.crt,nod.crt,earth);
        flt4.xyzflt{ii} = xyzflt;
        % save green's functions
        if strcmpi(modspace.grnflag,'on')
            [ ~,prjflt4,~ ] = GTdef_prjflt4dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:),addon.crt(strInd,:));
            GTdef_save_greensfns(cfname,pnt,prjflt4,Xgrn4,modspace.coord,lon0,lat0,0);
        end
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault5 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt5.num~=0
fprintf(1,'\n.......... processing fault type-5 ..........\t');
tic
    for ii = 1:flt5.num
        cfname = flt5.name{ii};
        cflt   = flt5.flt(ii,:);
        % find geometry file for the fault
        geoname = flt5.geoname{ii};
        % find column names for the fault
        colname = flt5.colname{ii};
        % find subfaults for the master fault
        subInd = strcmpi(cfname,subflt.name);

        [ modspace,xyzflt,Xgrn5,newflt ] = GTdef_fault5(modspace,geoname,colname,cflt,subflt.flt(subInd,:),pnt.crt,los.crt,bsl.crt,nod.crt,earth);

        flt5.out(ii,:)  = newflt; % update Nd & Ns if not provided
        flt5.xyzflt{ii} = xyzflt;
        % save green's functions
        if strcmpi(modspace.grnflag,'on')
            [ ~,prjflt5,~ ] = GTdef_prjflt5(modspace,geoname,colname,cflt,subflt.flt(subInd,:));
            GTdef_save_greensfns(cfname,pnt,prjflt5,Xgrn5,modspace.coord,lon0,lat0,0);
        end
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault6 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt6.num~=0
fprintf(1,'\n.......... processing fault type-6 ..........\t');
tic
    for ii = 1:flt6.num
        cfname = flt6.name{ii};
        cflt   = flt6.flt(ii,:);
        % find geometry file for the fault
        geoname = flt6.geoname{ii};
        % find column names for the fault
        colname = flt6.colname{ii};
        % find subfaults for the master fault
        subInd = strcmpi(cfname,subflt.name);

        [ modspace,xyzflt,Xgrn6,newflt ] = GTdef_fault6(modspace,geoname,colname,cflt,subflt.flt(subInd,:),pnt.crt,los.crt,bsl.crt,nod.crt,earth);

        flt6.out(ii,:)  = newflt; % update Nd & Ns if not provided
        flt6.xyzflt{ii} = xyzflt;
        % save green's functions
        if strcmpi(modspace.grnflag,'on')
            [ ~,prjflt6,~ ] = GTdef_prjflt6(modspace,geoname,colname,cflt,subflt.flt(subInd,:));
            GTdef_save_greensfns(cfname,pnt,prjflt6,Xgrn6,modspace.coord,lon0,lat0,0);
        end
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault7 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt7.num~=0
fprintf(1,'\n.......... processing fault type-7 ..........\t');
tic
    for ii = 1:flt7.num
       cfname = flt7.name{ii};
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       % find green functions for the fault
       cgrname = flt7.grname{ii};
       % read greens functions
       [ siteList,siteloc,vertices,grnList,grnfns ] = PyLith_read_greensfns(cgrname);
       % trim greens functions
       [ siteloc,grnList,grnfns ] = PyLith_trim_greensfns(pnt.name,siteList,siteloc,grnList,grnfns);
       [ Xgrn7,Lgrn7,Bgrn7,Ngrn7,sm7,sm7_abs,Aeq7,beq7,lb7,ub7,x07 ] = ...
       GTdef_fault7(flt7.flt(ii,:),subflt.flt(subInd,:),vertices,grnfns,modspace.smooth,modspace.surf);
       [ modspace ] = GTdef_addall(modspace,Xgrn7,Lgrn7,Bgrn7,Ngrn7,sm7,sm7_abs,Aeq7,beq7,lb7,ub7,x07);
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% forward only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if modspace.lb==-Inf
    % set up forward parameters
    modspace.xx = modspace.x0;

    % forward calculation
    fprintf(1,'\n........... doing forward modeling ..........\t');
    tic
    [ modspace,pnt,los,bsl,nod ] = GTdef_forward(modspace,pnt,los,bsl,nod);
    toc
    if strcmp(modspace.sdropflag,'on')
        % calculate stress drop
        fprintf(1,'\n.......... calculating stress drop ..........\t');
        tic
        [ flt1,flt2,flt3,flt4,flt5 ] = GTdef_calc_stressdrop(earth,flt1,flt2,flt3,flt4,flt5);
        toc
    end

    if (strcmpi(modspace.proj,'on'))  % write out porjections
      ldir=pwd;
      myFiles=dir(fullfile(ldir,strcat(basename,'_fwd.out'))); % added to inlcude any fwd model
      GTdef_project(myFiles.name)
    end


    % stress calculation
%    if sspnt.num~=0 || ssflt1.fltnum~=0 || ssflt2.fltnum~=0
%        fstressName = [ basename '_stress.out' ];
%        [ sspnt,ssflt1,ssflt2 ] = GTdef_calc_stress(sspnt,ssflt1,ssflt2,Min0,earth);
%        GTdef_output_stress(fstressName,sspnt,ssflt1,ssflt2);
%    end
    foutName = [ basename '_fwd.out' ];
    GTdef_output(foutName,earth,modspace,0,flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt,addon,pnt,los,bsl,prf,grd,nod);
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% inversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n............. doing inversion .............\t');
    % do inversion for each beta
    betaNum = length(modspace.beta);
    for ii = 1:betaNum
        bt = modspace.beta(ii);
        fprintf(1,'\n............. beta = %10d .............\t',bt);
        tic
        if strcmp(modspace.smooth,'2d')
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
        % inversion
        [ modspace ] = GTdef_invert(modspace,pnt,los,bsl,bt);
        % forward calculation
        [ modspace,pnt,los,bsl,nod ] = GTdef_forward(modspace,pnt,los,bsl,nod);
        % update fault slips
        [ flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt ] = GTdef_update_slips(earth,modspace,flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt);
        if strcmp(modspace.sdropflag,'on')
            % calculate stress drop
            [ flt1,flt2,flt3,flt4,flt5,flt6 ] = GTdef_calc_stressdrop(earth,flt1,flt2,flt3,flt4,flt5,flt6);
        end
        % output results
        GTdef_output(foutName,earth,modspace,bt,flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt,addon,pnt,los,bsl,prf,grd,nod);
	% output resolution matrix
        if ~isempty(modspace.res)
            GTdef_resolution(foutName,modspace,flt1,flt2,flt3,flt4,flt5,flt6,subflt,addon,pnt,los,bsl);
        end
        toc
    end

    fsumName = [ basename '_inv.out' ];
    GTdef_summary(fsumName,modspace);

    if (strcmpi(modspace.mat,'on'))
        fmatName = [ basename '.mat' ];
	save(fmatName);
    end

    if (strcmpi(modspace.proj,'on'))  % write out porjections
      ldir=pwd;
      myFiles=dir(fullfile(ldir,strcat(basename,'_kp*.out')));
      for k=1:length(myFiles)
        GTdef_project(myFiles(k).name)
      end
    end

end

% close up parpool for parallel computing
if exist('wnum') && wnum~=0
   delete(localpool);
end
