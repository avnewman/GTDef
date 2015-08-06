function [ ] = GTdef_GTdef2faults(finName,pntpos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         GTdef_GTdef2faults                           %
%                                                                      %
% Convert GTdef format to lists of subfaults                           %
% only fault2 has been verified!                                       %
%                                                                      %
% INPUT:                                                               %
% finName - GTdef input/output file                                    %
% pntpos  - position for points                                        %
%    'center' or 'centre' - center points                              %
%    'top'    - top points                                             %
%                                                                      %
% OUTPUT:                                                              %
% a file ended with '_faults.out'                                      %
%                                                                      %
% first created by Lujia Feng Mon Jul  9 11:26:06 SGT 2012             %
% added origin lfeng Thu Dec  5 21:55:50 SGT 2013                      %
% last modified by Lujia Feng Thu Dec  5 21:56:53 SGT 2013             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'.......... reading the input file ...........\t');
tic
[ coord,origin,smooth,surf,beta,grnflag,...
  rigidity,poisson,...
  earth,edgrn,layer,...
  flt1,flt2,flt3,flt4,flt5,...
  subflt,dip,...
  pnt,bsl,prf,grd,...
  sspnt,ssflt1,ssflt2 ] = GTdef_open(finName);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set origin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the middle point of the region as origin for the local cartesian coordinate
if strcmpi(coord,'geo') || strcmpi(coord,'geo_polyconic')
    if isempty(origin)
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
        lon0 = origin(1); lat0 = origin(2);
    end
elseif strcmpi(coord,'local')~=1
    error('GTdef_GTdef2faults ERROR: Coordinate input is wrong!!!');
end

[ ~,basename,~ ] = fileparts(finName);
foutName = [ basename '_faults.out' ];
fout = fopen(foutName,'w');
fprintf(fout,'#------------------------------------------------------------------------------------------------\n');
fprintf(fout,'# (1)no (2)dnum (3)snum (4)lon (5)lat (6)z (7)length (8)width (9)strike (10)dip (11)rake (12)slip\n');
fprintf(fout,'#                         [deg]  [deg]  [m]    [m]       [m]     [deg]     [deg]   [deg]     [m]\n');
fprintf(fout,'#------------------------------------------------------------------------------------------------\n');

newflt = []; xyztop = []; xyzctr = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt1.num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    % convert to cartesian coordinate
    if strcmpi(coord,'geo')
        [x1,y1] = LL2ckmd(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
        [x1,y1] = latlon_to_xy(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
    end
    if strcmpi(coord,'local')
        x1 = flt1.flt(:,1); y1 = flt1.flt(:,2);
    end
    tmpflt1 = [x1 y1 flt1.flt(:,3:end)];
    for ii = 1:flt1.num
        Nd = flt1.flt(ii,17); Ns = flt1.flt(ii,18); flt_num = Nd*Ns;
        cflt_name = flt1.name{ii};
        % find the subfaults for the master fault
        sub_ind = strcmpi(cflt_name,subflt.name);
        % find dips for the master fault
        dip_ind = strcmpi(cflt_name,dip.name);
        [ newflt1,~,xyzctr1,xyztop1 ] = GTdef_prjflt1dif(tmpflt1(ii,:),subflt.flt(sub_ind,:),dip.dip(dip_ind,:));
        newflt = [ newflt; newflt1 ];
        xyztop = [ xyztop  xyztop1 ];
        xyzctr = [ xyzctr  xyzctr1 ];
    end
    if strcmpi(pntpos,'center') || strcmpi(pntpos,'centre')
        xx =  xyzctr(1,:)';  yy = xyzctr(2,:)'; zz = -xyzctr(3,:)';
    elseif strcmpi(pntpos,'top')
        xx =  xyztop(1,:)';  yy = xyztop(2,:)'; zz = -xyztop(3,:)';
    else
        error('GTdef_GTdef2faults ERROR: outputing points is wrong!!!');
    end
    % convert to geographic coordinate
    if strcmpi(coord,'geo')
        [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
        [lon,lat] = xy_to_latlon(xx,yy,lon0,lat0,0);
    end
    %            (1)  (2)  (3) (4) (5) (6) (7) (8) (9) (10-18)
    % newflt = [ dnum snum x1  y1  z1  z2  len str dip slips ];
    % slips  = [ ss ds ts ss0 ssX ds0 dsX ts0 tsX ];
    dnum  = newflt(:,1);  snum = newflt(:,2);
    x1    = newflt(:,3);  y1   = newflt(:,4);
    z1    = newflt(:,5);  z2   = newflt(:,6);
    len   = newflt(:,7);  str  = newflt(:,8);
    dp    = newflt(:,9);
    ss    = newflt(:,10); ds   = newflt(:,11);
    width = (z2-z1)./sind(dp);
    slip  = sqrt(ss.^2+ds.^2);
    rake  = atan2(ds,ss).*180/pi;
    subflt_num = size(newflt,1);
    % (1)no (2)dnum (3)snum (4)lon (5)lat (6)z (7)length (8)width (9)strike (10)dip (11)rake (12)slip
    for ii =1:subflt_num
        fprintf(fout,'%3d %3d %3d %12.5f %10.5f  %10.5e   %11.5e %12.5e  %8.5f %7.2f %7.2f %12.6f\n',...
	        ii,dnum(ii),snum(ii),lon(ii),lat(ii),zz(ii),len(ii),width(ii),str(ii),dp(ii),rake(ii),slip(ii));
    end
toc
end

newflt = []; xyztop = []; xyzctr = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt2.num~=0
fprintf(1,'\n.......... processing fault type-2 ..........\t');
tic
    % convert to cartesian coordinate
    if strcmpi(coord,'geo')
       [x2_1,y2_1] = LL2ckmd(flt2.flt(:,1),flt2.flt(:,2),lon0,lat0,0);
       [x2_2,y2_2] = LL2ckmd(flt2.flt(:,3),flt2.flt(:,4),lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
       [x2_1,y2_1] = latlon_to_xy(flt2.flt(:,1),flt2.flt(:,2),lon0,lat0,0);
       [x2_2,y2_2] = latlon_to_xy(flt2.flt(:,3),flt2.flt(:,4),lon0,lat0,0);
    end
    if strcmpi(coord,'local')
       x2_1 = flt2.flt(:,1); y2_1 = flt2.flt(:,2);
       x2_2 = flt2.flt(:,3); y2_2 = flt2.flt(:,4);
    end
    tmpflt2 = [ x2_1 y2_1 x2_2 y2_2 flt2.flt(:,5:end) ];
    for ii = 1:flt2.num
        Nd = flt2.flt(ii,17); Ns = flt2.flt(ii,18); flt_num = Nd*Ns;
        cflt_name = flt2.name{ii};
        % find the subfaults for the master fault
        sub_ind = strcmpi(cflt_name,subflt.name);
        % find dips for the master fault
        dip_ind = strcmpi(cflt_name,dip.name);
        [ newflt2,~,xyzctr2,xyztop2 ] = GTdef_prjflt2dif(tmpflt2(ii,:),subflt.flt(sub_ind,:),dip.dip(dip_ind,:));
        newflt = [ newflt; newflt2 ];
        xyztop = [ xyztop  xyztop2 ];
        xyzctr = [ xyzctr  xyzctr2 ];
    end
    if strcmpi(pntpos,'center') || strcmpi(pntpos,'centre')
        xx =  xyzctr(1,:)';  yy = xyzctr(2,:)'; zz = -xyzctr(3,:)';
    elseif strcmpi(pntpos,'top')
        xx =  xyztop(1,:)';  yy = xyztop(2,:)'; zz = -xyztop(3,:)';
    else
        error('GTdef_GTdef2faults ERROR: outputing points is wrong!!!');
    end
    % convert to geographic coordinate
    if strcmpi(coord,'geo')
        [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
        [lon,lat] = xy_to_latlon(xx,yy,lon0,lat0,0);
    end
    %            (1)  (2)  (3) (4) (5) (6) (7) (8) (9) (10-18)
    % newflt = [ dnum snum x1  y1  x2  y2  z1  z2  dip slips ];
    % slips  = [ ss ds ts ss0 ssX ds0 dsX ts0 tsX ];
    dnum  = newflt(:,1);  snum = newflt(:,2);
    x1    = newflt(:,3);  y1   = newflt(:,4);
    x2    = newflt(:,5);  y2   = newflt(:,6);
    z1    = newflt(:,7);  z2   = newflt(:,8);
    dp    = newflt(:,9);
    ss    = newflt(:,10); ds   = newflt(:,11);
    str   = 90-180*atan2(y2-y1,x2-x1)/pi;   % degree CW from N [0-360]
    len   = sqrt((x2-x1).^2+(y2-y1).^2);
    width = (z2-z1)./sind(dp);
    slip  = sqrt(ss.^2+ds.^2);
    rake  = atan2(ds,ss).*180/pi;
    subflt_num = size(newflt,1);
    % (1)no (2)dnum (3)snum (4)lon (5)lat (6)z (7)length (8)width (9)strike (10)dip (11)rake (12)slip
    for ii =1:subflt_num
        fprintf(fout,'%3d %3d %3d %12.5f %10.5f  %10.5e   %11.5e %12.5e  %8.5f %7.2f %7.2f %12.6f\n',...
	        ii,dnum(ii),snum(ii),lon(ii),lat(ii),zz(ii),len(ii),width(ii),str(ii),dp(ii),rake(ii),slip(ii));
    end
toc
end

newflt = []; xyztop = []; xyzctr = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt3.num~=0
fprintf(1,'\n.......... processing fault type-3 ..........\t');
tic
    % convert to cartesian coordinate
    if strcmpi(coord,'geo')
        [x3,y3] = LL2ckmd(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
        [x3,y3] = latlon_to_xy(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0,0);
    end
    if strcmpi(coord,'local')
        x3 = flt3.flt(:,1); y3 = flt3.flt(:,2);
    end
    tmpflt3 = [x3 y3 flt3.flt(:,3:end)];
    for ii = 1:flt3.num
        Nd = flt3.flt(ii,17); Ns = flt3.flt(ii,18); flt_num = Nd*Ns;
        cflt_name = flt3.name{ii};
        % find the subfaults for the master fault
        sub_ind = strcmpi(cflt_name,subflt.name);
        % find dips for the master fault
        dip_ind = strcmpi(cflt_name,dip.name);
        [ newflt3,~,xyzctr3,xyztop3 ] = GTdef_prjflt3dif(tmpflt3(ii,:),subflt.flt(sub_ind,:),dip.dip(dip_ind,:));
        newflt = [ newflt; newflt3 ];
        xyztop = [ xyztop  xyztop3 ];
        xyzctr = [ xyzctr  xyzctr3 ];
    end
    if strcmpi(pntpos,'center') || strcmpi(pntpos,'centre')
        xx =  xyzctr(1,:)';  yy = xyzctr(2,:)'; zz = -xyzctr(3,:)';
    elseif strcmpi(pntpos,'top')
        xx =  xyztop(1,:)';  yy = xyztop(2,:)'; zz = -xyztop(3,:)';
    else
        error('GTdef_GTdef2faults ERROR: outputing points is wrong!!!');
    end
    % convert to geographic coordinate
    if strcmpi(coord,'geo')
        [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
        [lon,lat] = xy_to_latlon(xx,yy,lon0,lat0,0);
    end
    %            (1)  (2)  (3) (4) (5) (6) (7) (8) (9) (10-18)
    % newflt = [ dnum snum x1  y1  z1  z2  len str dip slips];
    % slips  = [ rake rs ts rake0 rakeX rs0 rsX ts0 tsX ];
    dnum  = newflt(:,1);  snum = newflt(:,2);
    x1    = newflt(:,3);  y1   = newflt(:,4);
    z1    = newflt(:,5);  z2   = newflt(:,6);
    len   = newflt(:,7);  str  = newflt(:,8);
    dp    = newflt(:,9);
    rake  = newflt(:,10); slip = newflt(:,11);
    width = (z2-z1)./sind(dp);
    subflt_num = size(newflt,1);
    % (1)no (2)dnum (3)snum (4)lon (5)lat (6)z (7)length (8)width (9)strike (10)dip (11)rake (12)slip
    for ii =1:subflt_num
        fprintf(fout,'%3d %3d %3d %12.5f %10.5f  %10.5e   %11.5e %12.5e  %8.5f %7.2f %7.2f %12.6f\n',...
	        ii,dnum(ii),snum(ii),lon(ii),lat(ii),zz(ii),len(ii),width(ii),str(ii),dp(ii),rake(ii),slip(ii));
    end
toc
end

newflt = []; xyztop = []; xyzctr = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault4 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt4.num~=0
fprintf(1,'\n.......... processing fault type-4 ..........\t');
tic
    % convert to cartesian coordinate
    if strcmpi(coord,'geo')
       [x4_1,y4_1] = LL2ckmd(flt4.flt(:,1),flt4.flt(:,2),lon0,lat0,0);
       [x4_2,y4_2] = LL2ckmd(flt4.flt(:,3),flt4.flt(:,4),lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
       [x4_1,y4_1] = latlon_to_xy(flt4.flt(:,1),flt4.flt(:,2),lon0,lat0,0);
       [x4_2,y4_2] = latlon_to_xy(flt4.flt(:,3),flt4.flt(:,4),lon0,lat0,0);
    end
    if strcmpi(coord,'local')
       x4_1 = flt4.flt(:,1); y4_1 = flt4.flt(:,2);
       x4_2 = flt4.flt(:,3); y4_2 = flt4.flt(:,4);
    end
    tmpflt4 = [ x4_1 y4_1 x4_2 y4_2 flt4.flt(:,5:end) ];
    for ii = 1:flt4.num
        Nd = flt4.flt(ii,17); Ns = flt4.flt(ii,18); flt_num = Nd*Ns;
        cflt_name = flt4.name{ii};
        % find the subfaults for the master fault
        sub_ind = strcmpi(cflt_name,subflt.name);
        % find dips for the master fault
        dip_ind = strcmpi(cflt_name,dip.name);
        [ newflt4,~,xyzctr4,xyztop4 ] = GTdef_prjflt4dif(tmpflt4(ii,:),subflt.flt(sub_ind,:),dip.dip(dip_ind,:));
        newflt = [ newflt; newflt4 ];
        xyztop = [ xyztop  xyztop4 ];
        xyzctr = [ xyzctr  xyzctr4 ];
    end
    if strcmpi(pntpos,'center') || strcmpi(pntpos,'centre')
        xx =  xyzctr(1,:)';  yy = xyzctr(2,:)'; zz = -xyzctr(3,:)';
    elseif strcmpi(pntpos,'top')
        xx =  xyztop(1,:)';  yy = xyztop(2,:)'; zz = -xyztop(3,:)';
    else
        error('GTdef_GTdef2faults ERROR: outputing points is wrong!!!');
    end
    % convert to geographic coordinate
    if strcmpi(coord,'geo')
        [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
        [lon,lat] = xy_to_latlon(xx,yy,lon0,lat0,0);
    end
    %            (1)  (2)  (3) (4) (5) (6) (7) (8) (9) (10-18)
    % newflt = [ dnum snum x1  y1  x2  y2  z1  z2  dip slips ];
    % slips  = [ rake rs ts rake0 rakeX rs0 rsX ts0 tsX ];
    dnum  = newflt(:,1);  snum = newflt(:,2);
    x1    = newflt(:,3);  y1   = newflt(:,4);
    x2    = newflt(:,5);  y2   = newflt(:,6);
    z1    = newflt(:,7);  z2   = newflt(:,8);
    dp    = newflt(:,9);
    rake  = newflt(:,10); slip = newflt(:,11);
    str   = 90-180*atan2(y2-y1,x2-x1)/pi;   % degree CW from N [0-360]
    len   = sqrt((x2-x1).^2+(y2-y1).^2);
    width = (z2-z1)./sind(dp);
    subflt_num = size(newflt,1);
    % (1)no (2)dnum (3)snum (4)lon (5)lat (6)z (7)length (8)width (9)strike (10)dip (11)rake (12)slip
    for ii =1:subflt_num
        fprintf(fout,'%3d %3d %3d %12.5f %10.5f  %10.5e   %11.5e %12.5e  %8.5f %7.2f %7.2f %12.6f\n',...
	        ii,dnum(ii),snum(ii),lon(ii),lat(ii),zz(ii),len(ii),width(ii),str(ii),dp(ii),rake(ii),slip(ii));
    end
toc
end

fclose(fout);
