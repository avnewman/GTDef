function [] = GTdef_project_EW(fin_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         GTdef_project_EW.m                              %
%									  %
% Project fault geometry and slip information onto surface geographic 	  %
% coordinate ASSUMING FAULT STRIKE EW					  %
% Note: this is used for GMT grdrotater, not real fault orientation	  %
%									  %
% INPUT									  %
%   fin_name - input file name                                            %
% OUTPUT								  %
% (1)  create an output file that contains			          %
%  [ flt_name dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2   %
%    xctr  yctr ss ds ts ]  					  	  %
% PARAMETERS								  %
%   flt_name - name of fault					  	  %
%   dnum - row number for subfaults					  %
%       dnum is 0 for single fault					  %
%   snum - column number for subfaults				  	  %
%       snum is 0 for sinlge fualt					  %
%   [xtop1 ytop1], [xbot1 ybot1], [xbot2 ybot2], and [xtop2 ytop2]	  %
%   are the surface vertical projection of four points that confine 	  %
%   the fault interface 						  %
%   They are in a counterclockwise sense looking from the RHS of endpoint %
%   [xtop1 ytop1] and [xbot1 ybot1] correspond [xx yy] at z1 and z2	  %
%   [x1 y1] and [x2 y2] are surface vertical projection of fault endpoints%
%									  %
% based on GTdef_project.m lfeng Tue Jul 27 02:07:11 EDT 2010		  %
% last modified by lfeng Wed Dec  1 16:22:56 EST 2010			  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
   if pnt_num~=0
       lonlist = [ lonlist;pnt_loc(:,1) ]; latlist = [ latlist;pnt_loc(:,2) ];
   end
   if bsl_num~=0
       lonlist = [ lonlist;bsl_loc(:,1);bsl_loc(:,4) ]; latlist = [ latlist;bsl_loc(:,2);bsl_loc(:,5) ];
   end
   if prf_num~=0
       lonlist = [ lonlist;prf(:,1);prf(:,3) ]; latlist = [ latlist;prf(:,2);prf(:,4) ];
   end
   if grd_num~=0
       lonlist = [ lonlist;grd(:,3);grd(:,5) ]; latlist = [ latlist;grd(:,4);grd(:,6) ];
   end
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
    newflt1 = [ zeros(flt1_num,2) x1 y1 flt1(:,3:end) ];
    [ prjflt1 ] = GTdef_prjfault1_EW(newflt1);
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
    newflt2 = [ zeros(flt2_num,2) x2_1 y2_1 x2_2 y2_2 flt2(:,5:end) ];
    [ prjflt2 ] = GTdef_prjfault2_EW(newflt2);
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
    for ii = 1:flt3_num
       Nd = flt3(ii,17); Ns = flt3(ii,18); flt_num = Nd*Ns;
       cflt_name = flt3_name{ii};
       % find the subfaults for the master fault
       sub_ind = strcmpi(cflt_name,subflt_name);
       % find dips for the master fault
       dip_ind = strcmpi(cflt_name,dip_name);
       [ prjflt3 ] = GTdef_prjfault3_EW(newflt3(ii,:),subflt(sub_ind,:),dip(dip_ind,:));
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
    for ii = 1:flt4_num
       Nd = flt4(ii,17); Ns = flt4(ii,18); flt_num = Nd*Ns;
       cflt_name = flt4_name{ii};
       % find the subfaults for the master fault
       sub_ind = strcmpi(cflt_name,subflt_name);
       % find dips for the master fault
       dip_ind = strcmpi(cflt_name,dip_name);
       [ prjflt4 ] = GTdef_prjfault4_EW(newflt4(ii,:),subflt(sub_ind,:),dip(dip_ind,:));
       prjflt = [ prjflt; prjflt4 ];
       name4 = cell(flt_num,1);
       for ii = 1:flt_num, name4{ii} = cflt_name; end
       flt_name = [ flt_name; name4 ];   
    end
toc
end

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

fout_fname = strcat(basename,'_surface_EW.out');
fout = fopen(fout_fname,'w');

fprintf(fout,'#name(1) dnum(2) snum(3) xtop1(4) ytop1(5) xbot1(6) ybot1(7) xbot2(8) ybot2(9) xtop2(10) ytop2(11) xctr(12) yctr(13) ss[m](14) ds[m](15) ts[m](16)\n'); 
[ row,col ] = size(newprjflt);

for ii =1:row
    name = flt_name{ii};
    flt  = newprjflt(ii,:);
    fprintf(fout,'%s %-3d %-3d %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-12.5f %-11.5f %-8.5f %-8.5f %-8.5f\n',name,flt);
end
fclose(fout);
