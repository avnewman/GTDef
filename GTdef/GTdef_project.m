function [] = GTdef_project(fin_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       GTdef_project.m                                  	% 
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
% last modified by lfeng Sun Apr 10 17:12:54 EDT 2011				%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'.......... reading the input file ...........\t');
tic
 [ coord,smooth,fsurf,beta,rigidity,poisson,... 
   flt1_name,flt1_num,flt1,flt2_name,flt2_num,flt2,...
   flt3_name,flt3_num,flt3,flt4_name,flt4_num,flt4,... 
   ~,~,~,~,~,...
   subflt_name,subflt,dip_name,dip,...
   pnt_name,pnt_num,pnt_loc,pnt_disp,pnt_err,pnt_wgt,... 
   bsl_name,bsl_num,bsl_loc,bsl_disp,bsl_err,bsl_wgt,...
   prf_name,prf_num,prf, grd_name,grd_num,grd ] = GTdef_open(fin_name);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set origin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the middle point of the region as origin for the local cartesian coordinate
if strcmpi(coord,'geo')
   lonlist = []; latlist = [];
   if flt1_num~=0
       lonlist = [ lonlist;flt1(:,1) ]; latlist = [ latlist;flt1(:,2) ];
   end
   if flt2_num~=0
       lonlist = [ lonlist;flt2(:,1) ]; latlist = [ latlist;flt2(:,2) ];
   end
   if flt3_num~=0
       lonlist = [ lonlist;flt3(:,1) ]; latlist = [ latlist;flt3(:,2) ];
   end
   if flt4_num~=0
       lonlist = [ lonlist;flt4(:,1) ]; latlist = [ latlist;flt4(:,2) ];
   end
   lon0 = 0.5*(min(lonlist)+max(lonlist));
   lat0 = 0.5*(min(latlist)+max(latlist));
elseif strcmpi(coord,'local')~=1
    error('Coordinate input is wrong!!!');
end

%basename = strtok(fin_name,'.');	% noly works for names without "."
cellname = regexp(fin_name,'\.(in|out)','split');
basename = char(cellname(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% point data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n......... processing the point data .........\t');
tic
pnt_crt = [];
if pnt_num~=0
    % convert point data from geographic to local cartesian coordinate
    if strcmpi(coord,'geo')
       [pxx,pyy] = LL2ckmd(pnt_loc(:,1),pnt_loc(:,2),lon0,lat0,0);
    end
    if strcmpi(coord,'local')
       pxx = pnt_loc(:,1); pyy = pnt_loc(:,2);
    end
    pzz = zeros(length(pxx),1);
    pnt_crt = [pxx pyy pzz];                   	% cartesian - 3*n matrix [xx;yy;zz]; it is just Xin

    % output point file
    fout_pname = strcat(basename,'_point.out');
    fpnt = fopen(fout_pname,'w');
    fprintf(fpnt,'#(1)point (2)3 (3)name (4)lon (5)lat (6)z (7)Ue (8)Un (9)Uv (10)eUe (11)eUn (12)eUv (13)weight (14)fault name (15)Dstr1  (16)Dstr2 (17)Ddip (18)Dvert\n'); 
end
toc


xsect = []; xsect_name = {};
prjflt = []; flt_name = {}; prjflt1 = []; prjflt2 = []; prjflt3 = []; prjflt4 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt1_num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    if strcmpi(coord,'geo')
       [x1,y1] = LL2ckmd(flt1(:,1),flt1(:,2),lon0,lat0,0);
    end
    if strcmpi(coord,'local')
       x1 = flt1(:,1); y1 = flt1(:,2);
    end
    newflt1 = [ x1 y1 flt1(:,3:end) ];
    [ xsect1_name,xsect1 ] = GTdef_xsection(1,flt1_name,newflt1,[]);
    xsect = [ xsect; xsect1 ];
    xsect_name = [ xsect_name; xsect1_name ];   
    if ~isempty(pnt_crt)
       GTdef_prjpnt(fpnt,pnt_name,pnt_loc,pnt_crt,pnt_disp,pnt_err,pnt_wgt,1,flt1_name,newflt1);
    end
    newflt1 = [ zeros(flt1_num,2) x1 y1 flt1(:,3:end) ];
    [ prjflt1 ] = GTdef_prjfault1(newflt1);
    prjflt = [ prjflt; prjflt1 ];
    flt_name = [ flt_name; flt1_name ];   
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt2_num~=0
fprintf(1,'\n.......... processing fault type-2 ..........\t');
tic
    if strcmpi(coord,'geo')
       [x2_1,y2_1] = LL2ckmd(flt2(:,1),flt2(:,2),lon0,lat0,0);
       [x2_2,y2_2] = LL2ckmd(flt2(:,3),flt2(:,4),lon0,lat0,0);
    end
    if strcmpi(coord,'local')
       x2_1 = flt2(:,1); y2_1 = flt2(:,2);
       x2_2 = flt2(:,3); y2_2 = flt2(:,4);
    end
    newflt2 = [ x2_1 y2_1 x2_2 y2_2 flt2(:,5:end) ];
    [ xsect2_name,xsect2 ] = GTdef_xsection(2,flt2_name,newflt2,[]);
    xsect = [ xsect; xsect2 ];
    xsect_name = [ xsect_name; xsect2_name ];   
    if ~isempty(pnt_crt)
       GTdef_prjpnt(fpnt,pnt_name,pnt_loc,pnt_crt,pnt_disp,pnt_err,pnt_wgt,2,flt2_name,newflt2);
    end
    newflt2 = [ zeros(flt2_num,2) x2_1 y2_1 x2_2 y2_2 flt2(:,5:end) ];
    [ prjflt2 ] = GTdef_prjfault2(newflt2);
    prjflt = [ prjflt; prjflt2 ];
    flt_name = [ flt_name; flt2_name ];   
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt3_num~=0
fprintf(1,'\n.......... processing fault type-3 ..........\t');
tic
    if strcmpi(coord,'geo')
       [x3,y3] = LL2ckmd(flt3(:,1),flt3(:,2),lon0,lat0,0);
    end
    if strcmpi(coord,'local')
       x3 = flt3(:,1); y3 = flt3(:,2);
    end
    newflt3 = [x3 y3 flt3(:,3:end)];
    if ~isempty(pnt_crt)
       GTdef_prjpnt(fpnt,pnt_name,pnt_loc,pnt_crt,pnt_disp,pnt_err,pnt_wgt,3,flt3_name,newflt3);
    end
    for ii = 1:flt3_num
       Nd = flt3(ii,17); Ns = flt3(ii,18); flt_num = Nd*Ns;
       cflt_name = flt3_name{ii};
       % find the subfaults for the master fault
       sub_ind = strcmpi(cflt_name,subflt_name);
       % find dips for the master fault
       dip_ind = strcmpi(cflt_name,dip_name);
       [ xsect3_name,xsect3 ] = GTdef_xsection(3,cflt_name,newflt3(ii,:),dip(dip_ind,:));
       xsect = [ xsect; xsect3 ];
       xsect_name = [ xsect_name; xsect3_name ];   
       [ prjflt3 ] = GTdef_prjfault3(newflt3(ii,:),subflt(sub_ind,:),dip(dip_ind,:));
       prjflt = [ prjflt; prjflt3 ];
       name3 = cell(flt_num,1);
       for ii = 1:flt_num, name3{ii} = cflt_name; end
       flt_name = [ flt_name; name3 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault4 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt4_num~=0
fprintf(1,'\n.......... processing fault type-4 ..........\t');
tic
    if strcmpi(coord,'geo')
       [x4_1,y4_1] = LL2ckmd(flt4(:,1),flt4(:,2),lon0,lat0,0);
       [x4_2,y4_2] = LL2ckmd(flt4(:,3),flt4(:,4),lon0,lat0,0);
    end
    if strcmpi(coord,'local')
       x4_1 = flt4(:,1); y4_1 = flt4(:,2);
       x4_2 = flt4(:,3); y4_2 = flt4(:,4);
    end
    newflt4 = [x4_1 y4_1 x4_2 y4_2 flt4(:,5:end)];
    if ~isempty(pnt_crt)
       GTdef_prjpnt(fpnt,pnt_name,pnt_loc,pnt_crt,pnt_disp,pnt_err,pnt_wgt,4,flt4_name,newflt4);
    end
    for ii = 1:flt4_num
       Nd = flt4(ii,17); Ns = flt4(ii,18); flt_num = Nd*Ns;
       cflt_name = flt4_name{ii};
       % find the subfaults for the master fault
       sub_ind = strcmpi(cflt_name,subflt_name);
       % find dips for the master fault
       dip_ind = strcmpi(cflt_name,dip_name);
       [ xsect4_name,xsect4 ] = GTdef_xsection(4,cflt_name,newflt4(ii,:),dip(dip_ind,:));
       xsect = [ xsect; xsect4 ];
       xsect_name = [ xsect_name; xsect4_name ];   
       [ prjflt4 ] = GTdef_prjfault4(newflt4(ii,:),subflt(sub_ind,:),dip(dip_ind,:));
       prjflt = [ prjflt; prjflt4 ];
       name4 = cell(flt_num,1);
       for ii = 1:flt_num, name4{ii} = cflt_name; end
       flt_name = [ flt_name; name4 ];   
    end
toc
end
if ~isempty(pnt_crt), fclose(fpnt); end

%            (1)  (2)  (3)    (4)   (5)  (6)   (7)   (8)   (9)    (10) (11) (12) (13) (14) (15)
% prjflt = [ dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2 xctr yctr  ss   ds   ts ]
if strcmpi(coord,'geo')
    xx = prjflt(:,[3 5 7 9 11]);  yy = prjflt(:,[4 6 8 10 12]);
    [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,0);
    newprjflt = [ prjflt(:,1:2) lon(:,1) lat(:,1) lon(:,2) lat(:,2) lon(:,3) lat(:,3) lon(:,4) lat(:,4) lon(:,5) lat(:,5) prjflt(:,13:15) ];
end
if strcmpi(coord,'local')
    newprjflt = prjflt;
end

% surface file
fout_fname = strcat(basename,'_surface.out');
fout = fopen(fout_fname,'w');

fprintf(fout,'#(1)name (2)dnum (3)snum (4)xtop1 (5)ytop1 (6)xbot1 (7)ybot1 (8)xbot2 (9)ybot2 (10)xtop2 (11)ytop2 (12)xctr (13)yctr (14)ss[m] (15)ds[m] (16)ts[m]\n'); 
[ row,col ] = size(newprjflt);

for ii =1:row
    name = flt_name{ii};
    flt  = newprjflt(ii,:);
    fprintf(fout,'%s %-3d %-3d %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-8.5f %-8.5f %-8.5f\n',name,flt);
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
