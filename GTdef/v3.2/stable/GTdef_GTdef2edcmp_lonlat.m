function [ ] = GTdef_GTdef2edcmp_lonlat(finName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         GTdef_GTdef2edcmp_lonlat                        %
%                                                                         %
% Convert GTdef format to edcmp format                                    %
% used to verify layered GTdef code                                       %
% only works for fault type-3                                             %
%                                                                         %
% INPUT:                                                                  %
% finName - GTdef input/output file                                       %
%                                                                         %
% OUTPUT:                                                                 %
% edcmp.inp file                                                          %
%                                                                         %
% first created by Lujia Feng Tue Jun 12 15:54:16 SGT 2012                %
% added origin lfeng Thu Dec  5 21:51:14 SGT 2013                         %
% modified outputs of GTdef_prjflt3dif.m lfeng Wed Apr 27 SGT 2016        %
% last modified by Lujia Feng 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'.......... reading the input file ...........\t');
tic
[ coord,origin,~,~,~,~,~,~,~,~,~,...
  flt1,flt2,flt3,flt4,~,~,~,...
  ~,subflt,dip,pnt,los,~,~,~,~,~,~ ] = GTdef_open(finName);
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
    error('GTdef_GTdef2edcmp_lonlat ERROR: coordinate input is wrong!!!');
end

[ ~,basename,~ ] = fileparts(finName);
foutName = [ basename '_edcmp.inp' ];
fout = fopen(foutName,'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% point data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n......... processing the point data .........\t');
tic
if pnt.num~=0
    fprintf(fout,'#------------------------------------------------------------------------------\n');
    fprintf(fout,'#     OBSERVATION ARRAY\n');
    fprintf(fout,'#------------------------------------------------------------------------------\n');
    fprintf(fout,'  0\n');                              % 0 - irregular positions
    fprintf(fout,'  %d\n ',pnt.num);
    for ii=1:pnt.num-1
        fprintf(fout,' (%f,%f),',pnt.loc(ii,[2 1]));  % switch x & y (GTdef:x-east;y-north; edcmp:x-north;y-east)
    end
    fprintf(fout,' (%f,%f)\n',pnt.loc(end,[2 1]));
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% specify output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'#     OUTPUTS\n');
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'# izmhs.disp  - displacement vectors\n');
fprintf(fout,'# izmhs.strn  - strain tenosrs\n');
fprintf(fout,'# izmhs.strss - stress tensors\n');
fprintf(fout,'# izmhs.tilt  - vertical tilts\n');
fprintf(fout,'# Note: name length < 80 characters\n');
fprintf(fout,'#       directories must exist and be ended by / (unix) or \\ (dos)\n');
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'  ''./edcmpout/''\n');
fprintf(fout,'       1                1                1                1\n');
fprintf(fout,'  ''izmhs.disp''   ''izmhs.strn''   ''izmhs.strss''  ''izmhs.tilt''\n');


newflt34 = []; newflt3 = [];
xyztop34 = []; xyztop3 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt3.num~=0
fprintf(1,'\n.......... processing fault type-3 ..........\t');
tic
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
       [ newflt3,~,xyzflt3 ] = GTdef_prjflt3dif(tmpflt3(ii,:),subflt.flt(sub_ind,:),dip.dip(dip_ind,:));
       xyztop3  = xyzflt3.xyztop1;
       newflt34 = [ newflt34; newflt3 ];
       xyztop34 = [ xyztop34 xyztop3 ];
    end
toc
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output fault data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(newflt34)
    %              (1)  (2)  (3) (4) (5) (6) (7) (8) (9) (10-18)
    % newflt34 = [ dnum snum x1  y1  z1  z2  len str dip slips ];
    xx =  xyztop34(1,:)';  yy = xyztop34(2,:)';
    zz = -xyztop34(3,:)';
    if strcmpi(coord,'geo')
        [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
        [lon,lat] = xy_to_latlon(xx,yy,lon0,lat0,0);
    end
    % switch xx & yy; top side down
    xyztop34 = [ lat lon zz ];
    width = (newflt34(:,6)-newflt34(:,5))./sind(newflt34(:,9));
    subflt_num = size(newflt34,1);
    fprintf(fout,'#------------------------------------------------------------------------------\n');
    fprintf(fout,'#     RECTANGLAR DISLOCATION SOURCES\n');
    fprintf(fout,'#------------------------------------------------------------------------------\n');
    fprintf(fout,'  %d\n',subflt_num);
    fprintf(fout,'#         coord. origin: (%fN, %fE)\n',lat0,lon0);
    fprintf(fout,'#-------------------------------------------------------------------------------\n');
    fprintf(fout,'# no  Slip   xs        ys       zs        length    width   strike   dip  rake\n');
    fprintf(fout,'#-------------------------------------------------------------------------------\n');
    for ii =1:subflt_num
        fprintf(fout,'  %-3d %12.6f  %-12.5f %-11.5f %-12.5e  %-11.5e %-12.5e %-8.5f  %-5.1f %-5.1f\n',...
	        ii,newflt34(ii,11),xyztop34(ii,:),newflt34(ii,7),width(ii),newflt34(ii,[8 9 10]));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% earth model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fout,'#===============================================================================\n');
fprintf(fout,'#     CHOICE OF EARTH MODEL\n');
fprintf(fout,'#===============================================================================\n');
fprintf(fout,'# 1 - layered model or 0 vs homogeneous model\n');
fprintf(fout,'# if 1\n');
fprintf(fout,'# grndir      - directory of the Greens functions = edgrnfcts\n');
fprintf(fout,'# grnfiles(3) - three greens functions files = izmhs.ss; izmhs.ds; izmhs.cl\n');
fprintf(fout,'# if 2\n');
fprintf(fout,'# zrec [m], lambda [Pa], mu [Pa]\n');
fprintf(fout,'# observation depth, the two Lame constants parameters of the homogeneous\n');
fprintf(fout,'#===============================================================================\n');

% This is for layered model
fprintf(fout,'  1\n');
fprintf(fout,'  ''./edgrnfcts/''  ''izmhs.ss''  ''izmhs.ds''  ''izmhs.cl''\n');

% % This is for Okada half space
% fprintf(fout, '  0\n');
% fprintf(fout, '  0.00d+00  %f %f\n', lambda, mu);

fprintf(fout,'#================================end of input===================================\n');
fclose(fout);
