function [] = GTdef_project(finName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       	GTdef_project.m                                 % 
%									        %
% (1) Project fault geometry and slip information onto surface geographic       %
% coordinate								        %
%     multiple faults are allowed!                                              %
% (2) Project points onto fault coordiantes				        %
%									        %
% INPUT									        %
%   finName - input file name                                                   %
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
% first created by Lujia Feng Fri Dec  4 19:51:19 EST 2009			%
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
% added origin lfeng Thu Dec  5 21:49:25 SGT 2013                               %
% added addon to combine 'dip' & 'strike' lfeng Fri Oct 24 15:06:09 SGT 2014    %
% added sweepAngle lfeng Wed Nov 12 10:57:14 SGT 2014                           %
% changed output from *_surface.out to *_patches.out lfeng Wed Nov 12 2014      %
% fixed typos for fault3 fault4 with Paul M. lfeng Fri Dec 12 11:03:32 SGT 2014 %
% added modspace structure lfeng Tue Mar 24 13:12:58 SGT 2015                   %
% addef fault5 lfeng Tue Jun 23 18:47:29 SGT 2015                               %
% last modified by Lujia Feng Tue Jun 23 19:48:49 SGT 2015                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'.......... reading the input file ...........\t');
tic
[ modspace,~,...
  flt1,flt2,flt3,flt4,flt5,flt6,...
  subflt,addon,...
  pnt,los,~,~,~,...
  ~,~,~ ] = GTdef_open(finName);
toc

coord  = modspace.coord;
origin = modspace.origin;

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
   error('GTdef_project ERROR: coordinate input is wrong!!!');
end

[ ~,basename,~ ] = fileparts(finName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% strike variations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addon.crt = [];
if addon.strnum~=0
fprintf(1,'\n....... processing strike variations .......\t');
tic
    % convert strike controlling points from geographic to local cartesian coordinate
    switch coord
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% point data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pnt.crt = [];
if pnt.num~=0
fprintf(1,'\n......... processing the point data .........\t');
tic
    % convert point data from geographic to local cartesian coordinate
    switch coord
       case 'geo'
          [pxx,pyy] = LL2ckmd(pnt.loc(:,1),pnt.loc(:,2),lon0,lat0,0);
       case 'geo_polyconic'
          [pxx,pyy] = latlon_to_xy(pnt.loc(:,1),pnt.loc(:,2),lon0,lat0,0);
       case 'local'
          pxx = pnt.loc(:,1); pyy = pnt.loc(:,2);
    end
    pzz = zeros(length(pxx),1);
    pnt.crt = [pxx pyy pzz];                   	% cartesian - 3*n matrix [xx;yy;zz]; it is just Xin

    % output point file
    fpntName = strcat(basename,'_point.out');
    fpnt = fopen(fpntName,'w');
    fprintf(fpnt,'#(1)point (2)3 (3)name (4)lon (5)lat (6)z (7)Ue (8)Un (9)Uv (10)eUe (11)eUn (12)eUv (13)weight (14)fault name (15)Dstr1  (16)Dstr2 (17)Ddip (18)Dvert\n'); 
toc
end


xsect = []; xsectName = {};
prjflt12 = []; flt12Name = {}; 
prjflt34 = []; flt34Name = {};
prjflt1  = []; prjflt2 = []; prjflt3 = []; prjflt4 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt1.num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    switch coord
       case 'geo'
          [x1,y1] = LL2ckmd(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
       case 'geo_polyconic'
          [x1,y1] = latlon_to_xy(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
       case 'local'
          x1 = flt1.flt(:,1); y1 = flt1.flt(:,2);
    end
    newflt1 = [ x1 y1 flt1.flt(:,3:end) ];
    for ii = 1:flt1.num
       Nd = flt1.flt(ii,17); Ns = flt1.flt(ii,18); fltNum = Nd*Ns;
       cfname = flt1.name{ii};
       cflt   = newflt1(ii,:);
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       % find dips for the master fault
       dipInd = strcmpi(cfname,addon.dipname);

       % point projection
       if ~isempty(pnt.crt)
          GTdef_prjpnt(fpnt,pnt,1,cfname,cflt);
       end
       % cross-section projection
       [ xsect1Name,xsect1 ] = GTdef_xsection(1,cfname,cflt,addon.dip(dipInd,:));
       xsect     = [ xsect; xsect1 ];
       xsectName = [ xsectName; xsect1Name ];   
       % surface projection
       [ ~,prjflt1,~ ] = GTdef_prjflt1dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:));
       prjflt12 = [ prjflt12; prjflt1 ];
       name1    = cell(fltNum,1);
       for ii = 1:fltNum, name1{ii} = cfname; end
       flt12Name = [ flt12Name; name1 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt2.num~=0
fprintf(1,'\n.......... processing fault type-2 ..........\t');
tic
    switch coord
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
       Nd = flt2.flt(ii,17); Ns = flt2.flt(ii,18); fltNum = Nd*Ns;
       cfname = flt2.name{ii};
       cflt   = newflt2(ii,:);
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       % find dips for the master fault
       dipInd = strcmpi(cfname,addon.dipname);
       % find strikes for the master fault
       strInd = strcmpi(cfname,addon.strname);

       % point projection
       if ~isempty(pnt.crt)
          GTdef_prjpnt(fpnt,pnt,2,cfname,cflt,addon.crt(strInd,:));
       end
       % cross-section projection
       [ xsect2Name,xsect2 ] = GTdef_xsection(2,cfname,cflt,addon.dip(dipInd,:));
       xsect     = [ xsect; xsect2 ];
       xsectName = [ xsectName; xsect2Name ];   
       % surface projection
       [ ~,prjflt2,~ ] = GTdef_prjflt2dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:),addon.crt(strInd,:));
       prjflt12 = [ prjflt12; prjflt2 ];
       name2    = cell(fltNum,1);
       for ii = 1:fltNum, name2{ii} = cfname; end
       flt12Name = [ flt12Name; name2 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt3.num~=0
fprintf(1,'\n.......... processing fault type-3 ..........\t');
tic
    % convert to cartesian coordinate
    switch coord
       case 'geo'
          [x3,y3] = LL2ckmd(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0,0);
       case 'geo_polyconic'
          [x3,y3] = latlon_to_xy(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0,0);
       case 'local'
          x3 = flt3.flt(:,1); y3 = flt3.flt(:,2);
    end
    newflt3 = [x3 y3 flt3.flt(:,3:end)];
    for ii = 1:flt3.num
       Nd = flt3.flt(ii,17); Ns = flt3.flt(ii,18); fltNum = Nd*Ns;
       cfname = flt3.name{ii};
       cflt   = newflt3(ii,:);
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       % find dips for the master fault
       dipInd = strcmpi(cfname,addon.dipname);

       % point projection
       if ~isempty(pnt.crt)
          GTdef_prjpnt(fpnt,pnt,3,cfname,cflt);
       end
       % cross-section projection
       [ xsect3Name,xsect3 ] = GTdef_xsection(3,cfname,cflt,addon.dip(dipInd,:));
       xsect     = [ xsect; xsect3 ];
       xsectName = [ xsectName; xsect3Name ];   
       % surface projection
       [ ~,prjflt3,~ ] = GTdef_prjflt3dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:));
       prjflt34 = [ prjflt34; prjflt3 ];
       name3 = cell(fltNum,1);
       for ii = 1:fltNum, name3{ii} = cfname; end
       flt34Name = [ flt34Name; name3 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault4 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt4.num~=0
fprintf(1,'\n.......... processing fault type-4 ..........\t');
tic
    switch coord
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
       Nd = flt4.flt(ii,17); Ns = flt4.flt(ii,18); fltNum = Nd*Ns;
       cfname = flt4.name{ii};
       cflt   = newflt4(ii,:);
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       % find dips for the master fault
       dipInd = strcmpi(cfname,addon.dipname);
       % find strikes for the master fault
       strInd = strcmpi(cfname,addon.strname);

       % point projection
       if ~isempty(pnt.crt)
          GTdef_prjpnt(fpnt,pnt,4,cfname,cflt,addon.crt(strInd,:));
       end
       % cross-section projection
       [ xsect4Name,xsect4 ] = GTdef_xsection(4,cfname,cflt,addon.dip(dipInd,:));
       xsect     = [ xsect; xsect4 ];
       xsectName = [ xsectName; xsect4Name ];   
       % surface projection
       [ ~,prjflt4,~ ] = GTdef_prjflt4dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:),addon.crt(strInd,:));
       prjflt34 = [ prjflt34; prjflt4 ];
       name4 = cell(fltNum,1);
       for ii = 1:fltNum, name4{ii} = cfname; end
       flt34Name = [ flt34Name; name4 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault5 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt5.num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    for ii = 1:flt5.num
       cfname  = flt5.name{ii};
       cflt    = flt5.flt(ii,:);
       % find geometry file for the fault
       geoname = flt5.geoname{ii};
       % find column names for the fault
       colname = flt5.colname{ii};
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       
       % point & cross-section projections do ont apply to fault5

       % surface projection
       [ ~,prjflt5,~ ] = GTdef_prjflt5(modspace,geoname,colname,cflt,subflt.flt(subInd,:));
       prjflt12 = [ prjflt12; prjflt5 ];
       fltNum   = size(prjflt5,1);
       name5    = cell(fltNum,1);
       for ii = 1:fltNum, name5{ii} = cfname; end
       flt12Name = [ flt12Name; name5 ];   
    end
toc
end

if ~isempty(pnt.crt)  
    fclose(fpnt); 
    fprintf(1,'\nGTdef_project output %s\n',fpntName); 
end

% surface file
foutName = strcat(basename,'_patches.out');
fout     = fopen(foutName,'w');

if ~isempty(prjflt12)
    %              (1)  (2)  (3)   (4)   (5)   (6)   (7)   (8)   (9)   (10)  (11)  (12)  (13)  (14)  (15) (16) (17) (18) (19) (20)
    % prjflt12 = [ dnum snum xtop1 ytop1 ztop1 xbot1 ybot1 zbot1 xbot2 ybot2 zbot2 xtop2 ytop2 ztop2 xctr yctr zctr ss   ds   ts ]
    if strcmpi(coord,'geo')
        xx = prjflt12(:,[3 6 9 12 15]);  yy = prjflt12(:,[4 7 10 13 16]);
        [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,0);
        newprjflt12 = [ prjflt12(:,1:2) lon(:,1) lat(:,1) prjflt12(:,5) lon(:,2) lat(:,2) prjflt12(:,8) lon(:,3) lat(:,3) prjflt12(:,11) ...
	                lon(:,4) lat(:,4) prjflt12(:,14) lon(:,5) lat(:,5) prjflt12(:,17) prjflt12(:,18:20) ];
    end
    if strcmpi(coord,'geo_polyconic')
        xx = prjflt12(:,[3 6 9 12 15]);  yy = prjflt12(:,[4 7 10 13 16]);
        lon = zeros(size(xx));         lat = zeros(size(yy));
        for ii=1:5
            [lon(:,ii),lat(:,ii)] = xy_to_latlon(xx(:,ii),yy(:,ii),lon0,lat0,0);
        end        
        newprjflt12 = [ prjflt12(:,1:2) lon(:,1) lat(:,1) prjflt12(:,5) lon(:,2) lat(:,2) prjflt12(:,8) lon(:,3) lat(:,3) prjflt12(:,11) ...
	                lon(:,4) lat(:,4) prjflt12(:,14) lon(:,5) lat(:,5) prjflt12(:,17) prjflt12(:,18:20) ];
    end
    if strcmpi(coord,'local')
        newprjflt12 = prjflt12;
    end
    fprintf(fout,'#(1)name (2)dnum (3)snum (4)xtop1 (5)ytop1 (6)ztop1 (7)xbot1 (8)ybot1 (9)zbot1 (10)xbot2 (11)ybot2 (12)zbot2 (13)xtop2 (14)ytop2 (15)ztop2 (16)xctr (17)yctr (18)zcrt (19)ss[m] (20)ds[m] (21)ts[m]\n'); 
    [ row,col ] = size(newprjflt12);
    for ii =1:row
        name = flt12Name{ii};
        flt  = newprjflt12(ii,:);
        fprintf(fout,'%-10s %4d %4d %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %10.5f %10.5f %10.5f\n',name,flt);
    end
end

if ~isempty(prjflt34)
    %              (1)  (2)  (3)   (4)   (5)   (6)   (7)   (8)   (9)   (10)  (11)  (12)  (13)  (14)  (15) (16) (17) (18) (19) (20)
    % prjflt34 = [ dnum snum xtop1 ytop1 ztop1 xbot1 ybot1 zbot1 xbot2 ybot2 zbot2 xtop2 ytop2 ztop2 xctr yctr zctr rake rs   ts ]
    if strcmpi(coord,'geo')
        xx = prjflt34(:,[3 6 9 12 15]);  yy = prjflt34(:,[4 7 10 13 16]);
        [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,0);
        newprjflt34 = [ prjflt34(:,1:2) lon(:,1) lat(:,1) prjflt34(:,5) lon(:,2) lat(:,2) prjflt34(:,8) lon(:,3) lat(:,3) prjflt34(:,11) ...
	                lon(:,4) lat(:,4) prjflt34(:,14) lon(:,5) lat(:,5) prjflt34(:,17) prjflt34(:,18:20) ];
    end
    if strcmpi(coord,'geo_polyconic')
        xx  = prjflt34(:,[3 6 9 12 15]); yy  = prjflt34(:,[4 7 10 13 16]);
        lon = zeros(size(xx));           lat = zeros(size(yy));
        for ii=1:5
            [lon(:,ii),lat(:,ii)] = xy_to_latlon(xx(:,ii),yy(:,ii),lon0,lat0,0);
        end         
        newprjflt34 = [ prjflt34(:,1:2) lon(:,1) lat(:,1) prjflt34(:,5) lon(:,2) lat(:,2) prjflt34(:,8) lon(:,3) lat(:,3) prjflt34(:,11) ...
	                lon(:,4) lat(:,4) prjflt34(:,14) lon(:,5) lat(:,5) prjflt34(:,17) prjflt34(:,18:20) ];
    end
    if strcmpi(coord,'local')
        newprjflt34 = prjflt34;
    end
    fprintf(fout,'#(1)name (2)dnum (3)snum (4)xtop1 (5)ytop1 (6)ztop1 (7)xbot1 (8)ybot1 (9)zbot1 (10)xbot2 (11)ybot2 (12)zbot2 (13)xtop2 (14)ytop2 (15)ztop2 (16)xctr (17)yctr (18)zcrt (19)rake[deg] (20)rs[m] (21)ts[m]\n'); 
    [ row,col ] = size(newprjflt34);
    for ii =1:row
        name = flt34Name{ii};
        flt  = newprjflt34(ii,:);
        fprintf(fout,'%-10s %4d %4d %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %10.5f %10.5f %10.5f\n',name,flt);
    end
end
fclose(fout);
fprintf(1,'\nGTdef_project output %s\n',foutName); 

% cross-section file
if ~isempty(xsect)
    foutName = strcat(basename,'_xsection.out');
    fout     = fopen(foutName,'w');
    
    fprintf(fout,'# (1)fault name (2)index (3)dip (4)x1[m] (5)z1[m] (6)x2[m] (7)z2[m] (8)width[m] (9)rows\n'); 
    row = size(xsect,1);
    
    for ii =1:row
        name = xsectName{ii};
        sct = xsect(ii,:);
        fprintf(fout,'%-10s  %5d  %6.2f    %12.4e   %14.6e   %12.4e %14.6e    %14.6e  %d\n',name,sct);
    end
    fclose(fout);
    fprintf(1,'\nGTdef_project output %s\n',foutName); 
end
