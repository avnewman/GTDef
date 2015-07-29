function [] = GTdef(fin_name,wnum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              GTdef.m                                    	%
% 	       Georgia Tech Matlab program for deformation	  	  	%
%		 Lujia Feng; Andrew V. Newman; Ting Chen   		  	%
%									  	%
% INPUT:									%
% fin_name - input file name                                            	%
% wnum     - num of matlab parallel workers to be used				%
%   0: do not use parallel computing						%
%   8: up to 8 workers that can be specified                                    %
%  -1: use parallel computing                                                   %
%      but number of workers will be determined by the system                   %
%									  	%
% v1:									  	%
%  first created by Lujia Feng Mon Apr 20 14:15:52 EDT 2009		  	%
%  last modified by Lujia Feng Tue May 19 01:14:49 EDT 2009		  	%
% v1.1:									  	%
%  added coordiante type flag 'coord' lfeng Thu Nov 5 17:12:59 EST 2009   	%
% v1.2:									  	%
% added first derivative modes lfeng Tue Dec  1 14:31:10 EST 2009	  	%
% modified beta (added beta) and roughness for 1st derivatives 	  	  	%
% added 'dip' flag for bended faults lfeng Mon Dec  7 01:04:06 EST 2009	  	%
% added 'freesurface' flag lfeng Wed Dec  9 17:00:58 EST 2009		  	%
% added fault type 5 lfeng Fri Dec 11 10:57:18 EST 2009			  	%
% changed 'freesurface' to 'surface' flag lfeng Wed Feb 24 12:46:01 EST 2010	%
% changed 'coord' to string flag lfeng Wed Feb 24 13:40:01 EST 2010		%
% allows input file name include multiple "." besides ".in" lfeng Oct 4 2010	%
% added matlabpool lfeng Wed Dec  1 12:12:00 EST 2010				%
% edited the origin definition lfeng Thu Apr  7 18:45:35 EDT 2011		%
% last modified lfeng Thu Apr  7 18:45:48 EDT 2011				%
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
   error('Matlabpool input is wrong!!!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'.......... reading the input file ...........\t');
tic
 [ coord,smooth,surf,beta,rigidity,poisson,... 
   flt1_name,flt1_num,flt1, flt2_name,flt2_num,flt2,...
   flt3_name,flt3_num,flt3, flt4_name,flt4_num,flt4,...
   ~,~,~,~,~,...
   subflt_name,subflt,dip_name,dip,...
   pnt_name,pnt_num,pnt_loc,pnt_disp,pnt_err,pnt_wgt,... 
   bsl_name,bsl_num,bsl_loc,bsl_disp,bsl_err,bsl_wgt,...
   prf_name,prf_num,prf,grd_name,grd_num,grd ] = GTdef_open(fin_name);
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set origin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the middle point of the region as origin for the local cartesian coordinate
if strcmp(coord,'geo')
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
elseif strcmp(coord,'local')~=1
    error('Coordinate input is wrong!!!');
end

pnt_crt = []; pnt_obs = []; pnt_obs_err = []; pnt_obs_wgt = []; pnt_coef = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% point data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n......... processing the point data .........\t');
tic
if pnt_num~=0
    % convert point data from geographic to local cartesian coordinate
    if strcmp(coord,'geo')
       [pxx,pyy] = LL2ckmd(pnt_loc(:,1),pnt_loc(:,2),lon0,lat0,0);
    end
    if strcmp(coord,'local')
       pxx = pnt_loc(:,1); pyy = pnt_loc(:,2);
    end
    pzz = zeros(1,length(pxx));
    pnt_crt = [pxx'; pyy'; pzz];                   	% cartesian - 3*n matrix [xx;yy;zz]; it is just Xin
    % prepare the point observation data
    pnt_obs = reshape(pnt_disp,[],1);			% (3*n)*1 observation vector [east;north;vertical]
    pnt_obs_err = reshape(pnt_err,[],1);		% (3*n)*1 error vector [east;north;vertical]
    pnt_obs_wgt = [pnt_wgt;pnt_wgt;pnt_wgt]; 		% (3*n)*1 weight vector [east;north;vertical]
    pnt_coef = sqrt(pnt_obs_wgt)./pnt_obs_err;
end
toc

bsl_crt = []; bsl_obs = []; bsl_obs_err = []; bsl_obs_wgt = []; bsl_coef = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% baseline data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n........ processing the baseline data .......\t');
tic
if bsl_num~=0
    % convert from geographic to local cartesian coordinate
    if strcmp(coord,'geo')
       [bx1,by1] = LL2ckmd(bsl_loc(:,1),bsl_loc(:,2),lon0,lat0,0);
       [bx2,by2] = LL2ckmd(bsl_loc(:,4),bsl_loc(:,5),lon0,lat0,0);
    end
    if strcmp(coord,'local')
       bx1 = bsl_loc(:,1); by1  = bsl_loc(:,2);
       bx2 = bsl_loc(:,4); by2  = bsl_loc(:,5);
    end
    bz1 = bsl_loc(:,3);	 bz2 = bsl_loc(:,6);
    bsl_crt = [bx1'; by1'; bz1';bx2'; by2'; bz2'];      % cartesian - 6*n matrix [bx1;by1;bz1;bx2;by2;bz2]; it is just Bin
    % prepare the baseline observation data
    bsl_obs = reshape(bsl_disp,[],1);			% (4*n)*1 observation vector [east;north;vertical:length]
    bsl_obs_err = reshape(bsl_err,[],1);		% (4*n)*1 error vector [east;north;vertical:length]
    bsl_obs_wgt = [bsl_wgt;bsl_wgt;bsl_wgt;bsl_wgt];	% (4*n)*1 weight vector [east;north;vertical:length]
    bsl_coef = sqrt(bsl_obs_wgt)./bsl_obs_err;
end
toc


nod_loc = []; nod_crt = []; nod_lon = []; nod_lat = []; nod_name = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% profile & grid data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n..... processing the profile & grid data ....\t');
tic
if prf_num~=0
    for ii = 1:prf_num
        plon = []; plat = []; pname = {};
	[ plon,plat,pname ] = GTdef_profile(prf(ii,:),prf_name{ii});
	nod_lon = [ nod_lon plon ]; nod_lat = [ nod_lat plat ];
    	nod_name = [ nod_name; pname ];   
    end
end
if grd_num~=0
    for ii = 1:grd_num
        glon = []; glat = []; gname = {};
	[ glon,glat,gname ] = GTdef_grid(grd(ii,:),grd_name{ii});
	nod_lon = [ nod_lon glon ]; nod_lat = [ nod_lat glat ];
    	nod_name = [ nod_name; gname ];   
    end
end
if prf_num~=0||grd_num~=0
    nod_zz = nan(length(nod_lon),1);
    nod_loc = [ nod_lon' nod_lat' nod_zz ];
    if strcmp(coord,'geo')
       [nod_xx,nod_yy] = LL2ckmd(nod_lon,nod_lat,lon0,lat0,0);
    end
    if strcmp(coord,'local')
       nod_xx = nod_lon; nod_yy = nod_lat;
    end
    nod_zz = zeros(1,length(nod_xx));
    nod_crt = [ nod_xx;nod_yy;nod_zz ];			% cartesian - 3*n matrix [xx;yy;zz]; it is just Nin
end
toc


% form everything that is needed for x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)
Xgrn = []; Bgrn = []; Ngrn = []; sm = []; Aeq = []; beq = []; lb = []; ub = []; x0 = [];
% sm_abs for calculate absolute 1st derivative (strain)
sm_abs = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt1_num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    if strcmp(coord,'geo')
       [x1,y1] = LL2ckmd(flt1(:,1),flt1(:,2),lon0,lat0,0);
    end
    if strcmp(coord,'local')
       x1 = flt1(:,1); y1 = flt1(:,2);
    end
    [ Xgrn1,Bgrn1,Ngrn1,sm1,Aeq1,beq1,lb1,ub1,x01 ] = GTdef_fault1([x1 y1 flt1(:,3:end)],pnt_crt,bsl_crt,nod_crt,rigidity,poisson);
    [ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = GTdef_addall(Xgrn,Bgrn,Ngrn,sm,[],Aeq,beq,lb,ub,x0,...
                                                          Xgrn1,Bgrn1,Ngrn1,sm1,[],Aeq1,beq1,lb1,ub1,x01);
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt2_num~=0
fprintf(1,'\n.......... processing fault type-2 ..........\t');
tic
    if strcmp(coord,'geo')
       [x2_1,y2_1] = LL2ckmd(flt2(:,1),flt2(:,2),lon0,lat0,0);
       [x2_2,y2_2] = LL2ckmd(flt2(:,3),flt2(:,4),lon0,lat0,0);
    end
    if  strcmp(coord,'local')
       x2_1 = flt2(:,1); y2_1 = flt2(:,2);
       x2_2 = flt2(:,3); y2_2 = flt2(:,4);
    end
    [ Xgrn2,Bgrn2,Ngrn2,sm2,Aeq2,beq2,lb2,ub2,x02 ] = GTdef_fault2([x2_1 y2_1 x2_2 y2_2 flt2(:,5:end)],pnt_crt,bsl_crt,nod_crt,rigidity,poisson);
    [ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = GTdef_addall(Xgrn,Bgrn,Ngrn,sm,[],Aeq,beq,lb,ub,x0,...
                                                          Xgrn2,Bgrn2,Ngrn2,sm2,[],Aeq2,beq2,lb2,ub2,x02);
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt3_num~=0
fprintf(1,'\n.......... processing fault type-3 ..........\t');
tic
    if strcmp(coord,'geo')
       [x3,y3] = LL2ckmd(flt3(:,1),flt3(:,2),lon0,lat0,0);
    end
    if  strcmp(coord,'local')
       x3 = flt3(:,1); y3 = flt3(:,2);
    end
    for ii = 1:flt3_num
        cflt_name = flt3_name{ii};
    	% find subfaults for the master fault
    	sub_ind = strcmpi(cflt_name,subflt_name);
	% find dips for the master fault
    	dip_ind = strcmpi(cflt_name,dip_name);
	Xgrn3 = []; Bgrn3 = []; Ngrn3 = []; Aeq3 = []; beq3 = []; lb3 = []; ub3 = []; x03 = [];		% reset empty
        [ Xgrn3,Bgrn3,Ngrn3,sm3,sm3_abs,Aeq3,beq3,lb3,ub3,x03 ] = ...
	GTdef_fault3([x3(ii) y3(ii) flt3(ii,3:end)],subflt(sub_ind,:),dip(dip_ind,:),pnt_crt,bsl_crt,nod_crt,rigidity,poisson,smooth,surf);
    	[ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = ...
	GTdef_addall(Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0,...
	             Xgrn3,Bgrn3,Ngrn3,sm3,sm3_abs,Aeq3,beq3,lb3,ub3,x03);
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault4 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt4_num~=0
fprintf(1,'\n.......... processing fault type-4 ..........\t');
tic
    if strcmp(coord,'geo')
       [x4_1,y4_1] = LL2ckmd(flt4(:,1),flt4(:,2),lon0,lat0,0);
       [x4_2,y4_2] = LL2ckmd(flt4(:,3),flt4(:,4),lon0,lat0,0);
    end
    if strcmp(coord,'local')
       x4_1 = flt4(:,1); y4_1 = flt4(:,2);
       x4_2 = flt4(:,3); y4_2 = flt4(:,4);
    end
    for ii = 1:flt4_num
        cflt_name = flt4_name{ii};
    	% find subfaults for the master fault
    	sub_ind = strcmpi(cflt_name,subflt_name);
	% find dips for the master fault
    	dip_ind = strcmpi(cflt_name,dip_name);
	Xgrn4 = []; Bgrn4 = []; Ngrn4 = []; Aeq4 = []; beq4 = []; lb4 = []; ub4 = []; x04 = [];		% reset empty
        [ Xgrn4,Bgrn4,Ngrn4,sm4,sm4_abs,Aeq4,beq4,lb4,ub4,x04 ] = ...
	GTdef_fault4([x4_1(ii) y4_1(ii) x4_2(ii) y4_2(ii) flt4(ii,5:end)],subflt(sub_ind,:),dip(dip_ind,:),pnt_crt,bsl_crt,nod_crt,rigidity,poisson,smooth,surf);
    	[ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = ...
	GTdef_addall(Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0,...
	             Xgrn4,Bgrn4,Ngrn4,sm4,sm4_abs,Aeq4,beq4,lb4,ub4,x04);
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault5 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flt5_num = 0;
if flt5_num~=0
    fprintf(1,'\n.......... processing fault type-5 ..........\t');
    tic
        if strcmp(coord,'geo')
           [bx1,by1] = LL2ckmd(bndry(:,3),bndry(:,4),lon0,lat0,0);		% upleft point
           [bx2,by2] = LL2ckmd(bndry(:,6),bndry(:,7),lon0,lat0,0);		% lower left point
           [bx3,by3] = LL2ckmd(bndry(:,9),bndry(:,10),lon0,lat0,0);		% lower right point
           [bx4,by4] = LL2ckmd(bndry(:,12),bndry(:,13),lon0,lat0,0);	% upleft point
           bndry_crt = [ bx1 by1 bndry(:,5) bx2 by2 bndry(:,8) bx3 by3 bndry(:,11) bx4 by4 bndry(:,14) ];
        end
        if strcmp(coord,'local')
    	bndry_crt = bndry;
        end
        for ii = 1:flt5_num
           cflt_name = flt5_name{ii};
        	% find subfaults for the master fault
        	sub_ind = strcmpi(cflt_name,subflt_name);
        	% find the boundary for fault 5
        	bnd_ind = strcmpi(cflt_name,bndry_name);
    	Xgrn5 = []; Bgrn5 = []; Ngrn5 = []; Aeq5 = []; beq5 = []; lb5 = []; ub5 = []; x05 = [];		% reset empty
            [ Xgrn5,Bgrn5,Ngrn5,sm5,sm5_abs,Aeq5,beq5,lb5,ub5,x05 ] = ...;
    	GTdef_fault5(flt5,subflt(sub_ind,:),bndry(bnd_ind,:),pnt_crt,bsl_crt,nod_crt,rigidity,poisson,smooth,surf);
        	[ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = ...
    	GTdef_addall(Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0,...
    	             Xgrn5,Bgrn5,Ngrn5,sm5,sm5_abs,Aeq5,beq5,lb5,ub5,x05);
        end
    toc
end


%basename = strtok(fin_name,'.');	% noly works for names without "."
cellname = regexp(fin_name,'\.in','split');
basename = char(cellname(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% forward only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lb==-Inf
fprintf(1,'\n.......... doing forward modeling ..........\t');
tic
    [ sm_use ] = GTdef_condense(sm);
    [ sm_abs ] = GTdef_condense(sm_abs);
    fout_name = [ basename '_fwd.out' ];
    % mod_info = [ data_num slip_num ndf rss rms wrrs wrms chi2 rchi2 r_1d r_2d ];
    [ mod_info,pnt_out,bsl_out,nod_out ]...
       = GTdef_forward(Xgrn,Bgrn,Ngrn,sm,sm_abs,lb,ub,x0,pnt_loc,pnt_obs,pnt_obs_err,pnt_wgt,bsl_loc,bsl_obs,bsl_obs_err,bsl_wgt,nod_loc,smooth);
    GTdef_output(fout_name,'none','none',0,rigidity,poisson,flt1_name,flt1_num,flt1,flt2_name,flt2_num,flt2,flt3_name,flt3_num,flt3,...
	  	 flt4_name,flt4_num,flt4,{},0,[],{},[],...
		 subflt_name,subflt,dip_name,dip,pnt_name,pnt_num,pnt_out,bsl_name,bsl_num,bsl_out,...
           	 prf_name,prf_num,prf,grd_name,grd_num,grd,nod_name,nod_out,mod_info);
toc
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% inversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n............. doing inversion .............\t');
    % condense the smoothing matrix by removing rows of all zeros
    [ sm_use ] = GTdef_condense(sm);
    [ sm_abs ] = GTdef_condense(sm_abs);

    % do inversion for each beta
    beta_num = length(beta);
    for ii = 1:beta_num
    bt = beta(ii);
    fprintf(1,'\n............. beta = %16.5f .............\t',bt);
    tic
	if strcmp(smooth,'2d')
	    kp = sqrt(bt);
            if kp>=1
                fout_name = strcat(basename,'_kp',num2str(kp,'%-.0f'),'.out');
	    else
                fout_name = strcat(basename,'_kp',num2str(kp,'%-.5f'),'.out');
	    end
	else
            if bt>=1
                fout_name = strcat(basename,'_bt',num2str(bt,'%-.0f'),'.out');
	    else
                fout_name = strcat(basename,'_bt',num2str(bt,'%-.5f'),'.out');
	    end
	end
        [ xx ] = GTdef_invert(Xgrn,Bgrn,sm_use,Aeq,beq,lb,ub,x0,pnt_obs,pnt_coef,bsl_obs,bsl_coef,bt);
	% faults info
        [ flt1_out,flt2_out,subflt_name,subflt_out ]...
          = GTdef_slips(lb,ub,xx,flt1_num,flt1,flt2_num,flt2,flt3_name,flt3_num,flt3,flt4_name,flt4_num,flt4,{},0,[]);
        [ mod_info(ii,:),pnt_out,bsl_out,nod_out ]...
          = GTdef_forward(Xgrn,Bgrn,Ngrn,sm_use,sm_abs,lb,ub,xx,pnt_loc,pnt_obs,pnt_obs_err,pnt_wgt,bsl_loc,bsl_obs,bsl_obs_err,bsl_wgt,nod_loc,smooth);
	% output results
        GTdef_output(fout_name,smooth,surf,bt,rigidity,poisson,flt1_name,flt1_num,flt1_out,flt2_name,flt2_num,flt2_out,flt3_name,flt3_num,flt3,...
              	     flt4_name,flt4_num,flt4,{},0,[],{},[],...
		     subflt_name,subflt_out,dip_name,dip,pnt_name,pnt_num,pnt_out,bsl_name,bsl_num,bsl_out,...
               	     prf_name,prf_num,prf,grd_name,grd_num,grd,nod_name,nod_out,mod_info(ii,:));
    toc
    end
    fout_sum = [ basename '_inv.out' ];
    GTdef_summary(fout_sum,beta,mod_info);
end

% close up matlabpool for parallel computing
if matlabpool('size')>0
   matlabpool close
end
