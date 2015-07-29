function [] = GTdef_project(fin_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       	GTdef_project.m                                 % 
%									        %
% (1) Project fault geometry and slip information onto surface geographic       %
% coordinate								        %
% (2) Project points onto fault coordiantes				        %
%									        %
% INPUT									        %
%   fin_name - input file name                                                  %
% OUTPUT								        %
% (1)  create an output file that contains			                %
%  [ flt_name dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2         %
%    xctr  yctr ss ds ts ]  					  	        %
% PARAMETERS								        %
%   flt_name - name of fault					  	        %
%   dnum - row number for subfaults					        %
%       dnum is 0 for single fault					        %
%   snum - column number for subfaults				  	        %
%       snum is 0 for sinlge fualt					        %
%   [xtop1 ytop1], [xbot1 ybot1], [xbot2 ybot2], and [xtop2 ytop2]	        %
%   are the surface vertical projection of four points that confine 	        %
%   the fault interface 						        %
%   They are in a counterclockwise sense looking from the RHS of endpoint       %
%   [xtop1 ytop1] and [xbot1 ybot1] correspond [xx yy] at z1 and z2	        %
%   [x1 y1] and [x2 y2] are surface vertical projection of fault endpoints      %
% (2) another file that contains points info				        %
%									        %
% first created by lfeng Fri Dec  4 19:51:19 EST 2009			        %
% added flag 'dip' by lfeng Mon Dec  7 01:05:28 EST 2009		        %
% added fault 5 by lfeng Sat Dec 12 00:10:15 EST 2009			        %
% changed 'coord' to string flag lfeng Wed Feb 24 13:40:01 EST 2010	        %
% corrected errors for fault1&2 lfeng Sat Jun 19 11:27:37 EDT 2010 Greece       %
% fixed a bug related to no point data lfeng Fri Oct  1 18:15:16 EDT 2010       %
% allows input file name include multiple "." besides ".in" lfeng Oct 5 2010	%
% added ".out" as input file name lfeng Mon Oct 25 18:29:53 EDT 2010		%
% added cross-section file lfeng Thu Dec  2 05:14:18 EST 2010			%
% added Dstr2 lfeng Fri Dec 10 14:42:02 EST 2010				%
% made origion only dependent on faults lfeng Sun Apr 10 17:11:15 EDT 2011	%
% used structure lfeng Thu Feb 23 10:12:16 SGT 2012				%
% merged fault1 & fault3 and fault2 & fault4 lfeng Sun May 13 23:41:32 SGT 2012 %
% added polyconic projection lfeng Thu Jun  7 13:27:43 SGT 2012                 %
% last modified by lfeng Wed Jun 13 16:00:57 SGT 2012                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'.......... reading the input file ...........\t');
tic
[ coord,~,~,~,~,~,~,~,~,...
  flt1,flt2,flt3,flt4,flt5,...
  ~,subflt,dip,pnt,~,~,~,~,~,~ ] = GTdef_open(fin_name);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set origin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the middle point of the region as origin for the local cartesian coordinate
if strcmpi(coord,'geo') || strcmpi(coord,'geo_polyconic')
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
    error('GTdef_project ERROR: Coordinate input is wrong!!!');
end

%basename = strtok(fin_name,'.');	% noly works for names without "."
cellname = regexp(fin_name,'\.(in|out)','split');
basename = char(cellname(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% point data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pnt.crt = [];
if pnt.num~=0
fprintf(1,'\n......... processing the point data .........\t');
tic
    % convert point data from geographic to local cartesian coordinate
    if strcmpi(coord,'geo')
       [pxx,pyy] = LL2ckmd(pnt.loc(:,1),pnt.loc(:,2),lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
       [pxx,pyy] = latlon_to_xy(pnt.loc(:,1),pnt.loc(:,2),lon0,lat0);
    end
    if strcmpi(coord,'local')
       pxx = pnt.loc(:,1); pyy = pnt.loc(:,2);
    end
    pzz = zeros(length(pxx),1);
    pnt.crt = [pxx pyy pzz];                   	% cartesian - 3*n matrix [xx;yy;zz]; it is just Xin

    % output point file
    fout_pname = strcat(basename,'_point.out');
    fpnt = fopen(fout_pname,'w');
    fprintf(fpnt,'#(1)point (2)3 (3)name (4)lon (5)lat (6)z (7)Ue (8)Un (9)Uv (10)eUe (11)eUn (12)eUv (13)weight (14)fault name (15)Dstr1  (16)Dstr2 (17)Ddip (18)Dvert\n'); 
toc
end


xsect = []; xsect_name = {};
prjflt12 = []; flt_name12 = {}; 
prjflt34 = []; flt_name34 = {};
prjflt1  = []; prjflt2 = []; prjflt3 = []; prjflt4 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt1.num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    if strcmpi(coord,'geo')
       [x1,y1] = LL2ckmd(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
       [x1,y1] = latlon_to_xy(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0);
    end
    if strcmpi(coord,'local')
       x1 = flt1.flt(:,1); y1 = flt1.flt(:,2);
    end
    newflt1 = [ x1 y1 flt1.flt(:,3:end) ];
    if ~isempty(pnt.crt)
       GTdef_prjpnt(fpnt,pnt,1,flt1.name,newflt1);
    end
    for ii = 1:flt1.num
       Nd = flt1.flt(ii,17); Ns = flt1.flt(ii,18); flt_num = Nd*Ns;
       cflt_name = flt1.name{ii};
       % find the subfaults for the master fault
       sub_ind = strcmpi(cflt_name,subflt.name);
       % find dips for the master fault
       dip_ind = strcmpi(cflt_name,dip.name);
       % cross-section projection
       [ xsect1_name,xsect1 ] = GTdef_xsection(1,cflt_name,newflt1(ii,:),dip.dip(dip_ind,:));
       xsect = [ xsect; xsect1 ];
       xsect_name = [ xsect_name; xsect1_name ];   
       % surface projection
       [ ~,prjflt1,~,~ ] = GTdef_prjflt1dif(newflt1(ii,:),subflt.flt(sub_ind,:),dip.dip(dip_ind,:));
       prjflt12 = [ prjflt12; prjflt1 ];
       name1 = cell(flt_num,1);
       for ii = 1:flt_num, name1{ii} = cflt_name; end
       flt_name12 = [ flt_name12; name1 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt2.num~=0
fprintf(1,'\n.......... processing fault type-2 ..........\t');
tic
    if strcmpi(coord,'geo')
       [x2_1,y2_1] = LL2ckmd(flt2.flt(:,1),flt2.flt(:,2),lon0,lat0,0);
       [x2_2,y2_2] = LL2ckmd(flt2.flt(:,3),flt2.flt(:,4),lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
       [x2_1,y2_1] = latlon_to_xy(flt2.flt(:,1),flt2.flt(:,2),lon0,lat0);
       [x2_2,y2_2] = latlon_to_xy(flt2.flt(:,3),flt2.flt(:,4),lon0,lat0);
    end
    if strcmpi(coord,'local')
       x2_1 = flt2.flt(:,1); y2_1 = flt2.flt(:,2);
       x2_2 = flt2.flt(:,3); y2_2 = flt2.flt(:,4);
    end
    newflt2 = [ x2_1 y2_1 x2_2 y2_2 flt2.flt(:,5:end) ];
    if ~isempty(pnt.crt)
       GTdef_prjpnt(fpnt,pnt,2,flt2.name,newflt2);
    end
    for ii = 1:flt2.num
       Nd = flt2.flt(ii,17); Ns = flt2.flt(ii,18); flt_num = Nd*Ns;
       cflt_name = flt2.name{ii};
       % find the subfaults for the master fault
       sub_ind = strcmpi(cflt_name,subflt.name);
       % find dips for the master fault
       dip_ind = strcmpi(cflt_name,dip.name);
       [ xsect2_name,xsect2 ] = GTdef_xsection(2,cflt_name,newflt2(ii,:),dip.dip(dip_ind,:));
       xsect = [ xsect; xsect2 ];
       xsect_name = [ xsect_name; xsect2_name ];   
       [ ~,prjflt2,~,~ ] = GTdef_prjflt2dif(newflt2(ii,:),subflt.flt(sub_ind,:),dip.dip(dip_ind,:));
       prjflt12 = [ prjflt12; prjflt2 ];
       name2 = cell(flt_num,1);
       for ii = 1:flt_num, name2{ii} = cflt_name; end
       flt_name12 = [ flt_name12; name2 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt3.num~=0
fprintf(1,'\n.......... processing fault type-3 ..........\t');
tic
    if strcmpi(coord,'geo')
       [x3,y3] = LL2ckmd(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
       [x3,y3] = latlon_to_xy(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0);
    end
    if strcmpi(coord,'local')
       x3 = flt3.flt(:,1); y3 = flt3.flt(:,2);
    end
    newflt3 = [x3 y3 flt3.flt(:,3:end)];
    if ~isempty(pnt.crt)
       GTdef_prjpnt(fpnt,pnt,3,flt3.name,newflt3);
    end
    for ii = 1:flt3.num
       Nd = flt3.flt(ii,17); Ns = flt3.flt(ii,18); flt_num = Nd*Ns;
       cflt_name = flt3.name{ii};
       % find the subfaults for the master fault
       sub_ind = strcmpi(cflt_name,subflt.name);
       % find dips for the master fault
       dip_ind = strcmpi(cflt_name,dip.name);
       [ xsect3_name,xsect3 ] = GTdef_xsection(3,cflt_name,newflt3(ii,:),dip.dip(dip_ind,:));
       xsect = [ xsect; xsect3 ];
       xsect_name = [ xsect_name; xsect3_name ];   
       [ ~,prjflt3,~,~ ] = GTdef_prjflt3dif(newflt3(ii,:),subflt.flt(sub_ind,:),dip.dip(dip_ind,:));
       prjflt34 = [ prjflt34; prjflt3 ];
       name3 = cell(flt_num,1);
       for ii = 1:flt_num, name3{ii} = cflt_name; end
       flt_name34 = [ flt_name34; name3 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault4 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt4.num~=0
fprintf(1,'\n.......... processing fault type-4 ..........\t');
tic
    if strcmpi(coord,'geo')
       [x4_1,y4_1] = LL2ckmd(flt4.flt(:,1),flt4.flt(:,2),lon0,lat0,0);
       [x4_2,y4_2] = LL2ckmd(flt4.flt(:,3),flt4.flt(:,4),lon0,lat0,0);
    end
    if strcmpi(coord,'geo_polyconic')
       [x4_1,y4_1] = latlon_to_xy(flt4.flt(:,1),flt4.flt(:,2),lon0,lat0);
       [x4_2,y4_2] = latlon_to_xy(flt4.flt(:,3),flt4.flt(:,4),lon0,lat0);
    end
    if strcmpi(coord,'local')
       x4_1 = flt4.flt(:,1); y4_1 = flt4.flt(:,2);
       x4_2 = flt4.flt(:,3); y4_2 = flt4.flt(:,4);
    end
    newflt4 = [x4_1 y4_1 x4_2 y4_2 flt4.flt(:,5:end)];
    if ~isempty(pnt.crt)
       GTdef_prjpnt(fpnt,pnt,4,flt4.name,newflt4);
    end
    for ii = 1:flt4.num
       Nd = flt4.flt(ii,17); Ns = flt4.flt(ii,18); flt_num = Nd*Ns;
       cflt_name = flt4.name{ii};
       % find the subfaults for the master fault
       sub_ind = strcmpi(cflt_name,subflt.name);
       % find dips for the master fault
       dip_ind = strcmpi(cflt_name,dip.name);
       [ xsect4_name,xsect4 ] = GTdef_xsection(4,cflt_name,newflt4(ii,:),dip.dip(dip_ind,:));
       xsect = [ xsect; xsect4 ];
       xsect_name = [ xsect_name; xsect4_name ];   
       [ ~,prjflt4,~,~ ] = GTdef_prjflt4dif(newflt4(ii,:),subflt.flt(sub_ind,:),dip.dip(dip_ind,:));
       prjflt34 = [ prjflt34; prjflt4 ];
       name4 = cell(flt_num,1);
       for ii = 1:flt_num, name4{ii} = cflt_name; end
       flt_name34 = [ flt_name34; name4 ];   
    end
toc
end
if ~isempty(pnt.crt), fclose(fpnt); end

% surface file
fout_fname = strcat(basename,'_surface.out');
fout = fopen(fout_fname,'w');

if ~isempty(prjflt12)
    %              (1)  (2)  (3)    (4)   (5)  (6)   (7)   (8)   (9)    (10) (11) (12) (13) (14) (15)
    % prjflt12 = [ dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2 xctr yctr  ss   ds   ts ]
    if strcmpi(coord,'geo')
        xx = prjflt12(:,[3 5 7 9 11]);  yy = prjflt12(:,[4 6 8 10 12]);
        [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,0);
        newprjflt12 = [ prjflt12(:,1:2) lon(:,1) lat(:,1) lon(:,2) lat(:,2) lon(:,3) lat(:,3) lon(:,4) lat(:,4) lon(:,5) lat(:,5) prjflt12(:,13:15) ];
    end
    if strcmpi(coord,'geo_polyconic')
        xx = prjflt12(:,[3 5 7 9 11]); yy = prjflt12(:,[4 6 8 10 12]);
        lon = zeros(size(xx));         lat = zeros(size(yy));
        for ii=1:5
            [lon(:,ii),lat(:,ii)] = xy_to_latlon(xx(:,ii),yy(:,ii),lon0,lat0);
        end        
        newprjflt12 = [ prjflt12(:,1:2) lon(:,1) lat(:,1) lon(:,2) lat(:,2) lon(:,3) lat(:,3) lon(:,4) lat(:,4) lon(:,5) lat(:,5) prjflt12(:,13:15) ];
    end
    if strcmpi(coord,'local')
        newprjflt12 = prjflt12;
    end
    fprintf(fout,'#(1)name (2)dnum (3)snum (4)xtop1 (5)ytop1 (6)xbot1 (7)ybot1 (8)xbot2 (9)ybot2 (10)xtop2 (11)ytop2 (12)xctr (13)yctr (14)ss[m] (15)ds[m] (16)ts[m]\n'); 
    [ row,col ] = size(newprjflt12);
    for ii =1:row
        name = flt_name12{ii};
        flt  = newprjflt12(ii,:);
        fprintf(fout,'%s %-3d %-3d %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-8.5f %-8.5f %-8.5f\n',name,flt);
    end
end

if ~isempty(prjflt34)
    %              (1)  (2)  (3)    (4)   (5)  (6)   (7)   (8)   (9)    (10) (11) (12) (13) (14) (15)
    % prjflt34 = [ dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2 xctr yctr rake  rs   ts ]
    if strcmpi(coord,'geo')
        xx = prjflt34(:,[3 5 7 9 11]);  yy = prjflt34(:,[4 6 8 10 12]);
        [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,0);
        newprjflt34 = [ prjflt34(:,1:2) lon(:,1) lat(:,1) lon(:,2) lat(:,2) lon(:,3) lat(:,3) lon(:,4) lat(:,4) lon(:,5) lat(:,5) prjflt34(:,13:15) ];
    end
    if strcmpi(coord,'geo_polyconic')
        xx = prjflt34(:,[3 5 7 9 11]);  yy = prjflt34(:,[4 6 8 10 12]);
        lon = zeros(size(xx));         lat = zeros(size(yy));
        for ii=1:5
            [lon(:,ii),lat(:,ii)] = xy_to_latlon(xx(:,ii),yy(:,ii),lon0,lat0);
        end         
        newprjflt34 = [ prjflt34(:,1:2) lon(:,1) lat(:,1) lon(:,2) lat(:,2) lon(:,3) lat(:,3) lon(:,4) lat(:,4) lon(:,5) lat(:,5) prjflt34(:,13:15) ];
    end
    if strcmpi(coord,'local')
        newprjflt34 = prjflt34;
    end
    fprintf(fout,'#(1)name (2)dnum (3)snum (4)xtop1 (5)ytop1 (6)xbot1 (7)ybot1 (8)xbot2 (9)ybot2 (10)xtop2 (11)ytop2 (12)xctr (13)yctr (14)rake[deg] (15)rs[m] (16)ts[m]\n'); 
    [ row,col ] = size(newprjflt34);
    for ii =1:row
        name = flt_name34{ii};
        flt  = newprjflt34(ii,:);
        fprintf(fout,'%s %-3d %-3d %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-8.2f %-8.5f %-8.5f\n',name,flt);
    end
end
fclose(fout);

% cross-section file
fout_fname = strcat(basename,'_xsection.out');
fout = fopen(fout_fname,'w');

fprintf(fout,'# (1)fault name (2)index (3)dip (4)x1[m] (5)z1[m] (6)x2[m] (7)z2[m] (8)width[m] (9)rows\n'); 
row = size(xsect,1);

for ii =1:row
    name = xsect_name{ii};
    sct = xsect(ii,:);
    fprintf(fout,'%s  %-5d  %6.2f    %12.4e   %-14.6e   %-12.4e %-14.6e    %-14.6e  %-d\n',name,sct);
end
fclose(fout);
