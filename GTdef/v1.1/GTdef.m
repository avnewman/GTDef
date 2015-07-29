function [] = GTdef(fin_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              GTdef.m                                    %
% 	       Georgia Tech Matlab program for deformation	  	  %
%		 Ting Chen; Andrew V. Newman; Lujia Feng		  %
%									  %
% INPUT									  %
%   fin_name - input file name                                            %
%									  %
% v1:									  %
%  first created by Lujia Feng Mon Apr 20 14:15:52 EDT 2009		  %
%  last modified by Lujia Feng Tue May 19 01:14:49 EDT 2009		  %
% v1.1:									  %
%  added coordiante type flag 'coord' lfeng Thu Nov 5 17:12:59 EST 2009   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'.......... reading the input file ...........\t');
tic
 [ coord,kappa,rigidity,poisson,... 
   flt1_name,flt1_num,flt1, flt2_name,flt2_num,flt2,...
   flt3_name,flt3_num,flt3, flt4_name,flt4_num,flt4,... 
   subflt_name,subflt_num,subflt,...
   pnt_name,pnt_num,pnt_loc,pnt_disp,pnt_err,pnt_wgt,... 
   bsl_name,bsl_num,bsl_loc,bsl_disp,bsl_err,bsl_wgt,...
   prf_name,prf_num,prf, grd_name,grd_num,grd ] = GTdef_open(fin_name);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set origin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the middle point of the region as origin for the local cartesian coordinate
if coord==1
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
end


pnt_crt = []; pnt_obs = []; pnt_obs_err = []; pnt_obs_wgt = []; pnt_coef = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% point data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n......... processing the point data .........\t');
tic
if pnt_num~=0
    % convert point data from geographic to local cartesian coordinate
    if coord==1
       [pxx,pyy] = LL2ckmd(pnt_loc(:,1),pnt_loc(:,2),lon0,lat0,0);
    end
    if coord==2
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
    if coord==1
       [bx1,by1] = LL2ckmd(bsl_loc(:,1),bsl_loc(:,2),lon0,lat0,0);
       [bx2,by2] = LL2ckmd(bsl_loc(:,4),bsl_loc(:,5),lon0,lat0,0);
    end
    if coord==2
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


nod_loc = []; nod_crt = []; nod_lon = []; nod_lat = []; nod_name = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% profile & grid data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n..... processing the profile & grid data ....\t');
tic
if prf_num~=0
    for ii = 1:prf_num
        plon = []; plat = []; pname = '';
	[ plon,plat,pname ] = GTdef_profile(prf(ii,:),prf_name(ii,:));
	nod_lon = [ nod_lon plon ]; nod_lat = [ nod_lat plat ];
    	nod_name = strvcat(nod_name,pname);   
    end
end
if grd_num~=0
    for ii = 1:grd_num
        glon = []; glat = []; gname = '';
	[ glon,glat,gname ] = GTdef_grid(grd(ii,:),grd_name(ii,:));
	nod_lon = [ nod_lon glon ]; nod_lat = [ nod_lat glat ];
    	nod_name = strvcat(nod_name,gname);   
    end
end
if prf_num~=0||grd_num~=0
    nod_zz = nan(length(nod_lon),1);
    nod_loc = [ nod_lon' nod_lat' nod_zz ];
    if coord==1
       [nod_xx,nod_yy] = LL2ckmd(nod_lon,nod_lat,lon0,lat0,0);
    end
    if coord==2
       nod_xx = nod_lon; nod_yy = nod_lat;
    end
    nod_zz = zeros(1,length(nod_xx));
    nod_crt = [ nod_xx;nod_yy;nod_zz ];			% cartesian - 3*n matrix [xx;yy;zz]; it is just Nin
end
toc


% form everything that is needed for x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)
Xgrn = []; Bgrn = []; Ngrn = []; sm = []; Aeq = []; beq = []; lb = []; ub = []; x0 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt1_num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    if coord==1
       [x1,y1] = LL2ckmd(flt1(:,1),flt1(:,2),lon0,lat0,0);
    end
    if coord==2
       x1 = flt1(:,1); y1 = flt1(:,2);
    end
    [ Xgrn1,Bgrn1,Ngrn1,sm1,Aeq1,beq1,lb1,ub1,x01 ] = GTdef_fault1([x1 y1 flt1(:,3:end)],pnt_crt,bsl_crt,nod_crt,rigidity,poisson);
    [ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0 ] = GTdef_addall(Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0,...
                                                          Xgrn1,Bgrn1,Ngrn1,sm1,Aeq1,beq1,lb1,ub1,x01);
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt2_num~=0
fprintf(1,'\n.......... processing fault type-2 ..........\t');
tic
    if coord==1
       [x2,y2] = LL2ckmd(flt2(:,1),flt2(:,2),lon0,lat0,0);
    end
    if coord==2
       x2 = flt2(:,1); y2 = flt2(:,2);
    end
    [ Xgrn2,Bgrn2,Ngrn2,sm2,Aeq2,beq2,lb2,ub2,x02 ] = GTdef_fault2([x2 y2 flt2(:,3:end)],pnt_crt,bsl_crt,nod_crt,rigidity,poisson);
    [ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0 ] = GTdef_addall(Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0,...
                                                          Xgrn2,Bgrn2,Ngrn2,sm2,Aeq2,beq2,lb2,ub2,x02);
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt3_num~=0
fprintf(1,'\n.......... processing fault type-3 ..........\t');
tic
    if coord==1
       [x3,y3] = LL2ckmd(flt3(:,1),flt3(:,2),lon0,lat0,0);
    end
    if coord==2
       x3 = flt3(:,1); y3 = flt3(:,2);
    end
    for ii = 1:flt3_num
    	% find the subfaults for the master fault
    	sub_ind = strmatch(flt3_name(ii,:),subflt_name,'exact');
	Xgrn3 = []; Bgrn3 = []; Ngrn3 = []; Aeq3 = []; beq3 = []; lb3 = []; ub3 = []; x03 = [];		% reset empty
        [ Xgrn3,Bgrn3,Ngrn3,sm3,Aeq3,beq3,lb3,ub3,x03 ] = GTdef_fault3([x3(ii) y3(ii) flt3(ii,3:end)],subflt(sub_ind,:),pnt_crt,bsl_crt,nod_crt,rigidity,poisson);
    	[ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0 ] = GTdef_addall(Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0,...
							      Xgrn3,Bgrn3,Ngrn3,sm3,Aeq3,beq3,lb3,ub3,x03);
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault4 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt4_num~=0
fprintf(1,'\n.......... processing fault type-4 ..........\t');
tic
    if coord==1
       [x4,y4] = LL2ckmd(flt4(:,1),flt4(:,2),lon0,lat0,0);
    end
    if coord==2
       x4 = flt4(:,1); y4 = flt4(:,2);
    end
    for ii = 1:flt4_num
    	% find the subfaults for the master fault
    	sub_ind = strmatch(flt4_name(ii,:),subflt_name,'exact');
	Xgrn4 = []; Bgrn4 = []; Ngrn4 = []; Aeq4 = []; beq4 = []; lb4 = []; ub4 = []; x04 = [];		% reset empty
        [ Xgrn4,Bgrn4,Ngrn4,sm4,Aeq4,beq4,lb4,ub4,x04 ] = GTdef_fault4([x4(ii) y4(ii) flt4(ii,3:end)],subflt(sub_ind,:),pnt_crt,bsl_crt,nod_crt,rigidity,poisson);
    	[ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0 ] = GTdef_addall(Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0,...
							      Xgrn4,Bgrn4,Ngrn4,sm4,Aeq4,beq4,lb4,ub4,x04);
    end
toc
end


basename = strtok(fin_name,'.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% forward only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lb==-Inf
fprintf(1,'\n.......... doing forward modeling ..........\t');
tic
    fout_name = strcat(basename,'_fwd.out');
    % mod_info = [ data_num slip_num ndf rss rms wrrs wrms chi2 rchi2 roughness ];
    [ mod_info,pnt_out,bsl_out,nod_out ]...
       = GTdef_forward(Xgrn,Bgrn,Ngrn,[],lb,ub,x0,pnt_loc,pnt_obs,pnt_obs_err,pnt_wgt,bsl_loc,bsl_obs,bsl_obs_err,bsl_wgt,nod_loc);
    GTdef_output(fout_name,0,rigidity,poisson,flt1_name,flt1_num,flt1,flt2_name,flt2_num,flt2,flt3_name,flt3_num,flt3,...
	  	 flt4_name,flt4_num,flt4,subflt_name,subflt,pnt_name,pnt_num,pnt_out,bsl_name,bsl_num,bsl_out,...
           	 prf_name,prf_num,prf,grd_name,grd_num,grd,nod_name,nod_out,mod_info);
toc
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% inversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n............. doing inversion .............\t');
    % condense the smoothing matrix by removing rows of all zeros
    [ sm_row sm_col ] = size(sm);
    zero_row = zeros(1,sm_col);
    ind = [];
    for ii = 1:sm_row
        if ~isequal(sm(ii,:),zero_row)
            ind = [ ind ii ];
        end
    end
    sm_use = sm(ind,:);

    % do inversion for each kappa
    kappa_num = length(kappa);
    for ii = 1:kappa_num
    kp = kappa(ii);
    fprintf(1,'\n............. kappa = %10.5f .............\t',kp);
    tic
        if kp>=1
            fout_name = strcat(basename,'_kp',num2str(kp,'%-.0f'),'.out');
	else
            fout_name = strcat(basename,'_kp',num2str(kp,'%-.5f'),'.out');
	end
        [ xx ] = GTdef_invert(Xgrn,Bgrn,sm_use,Aeq,beq,lb,ub,x0,pnt_obs,pnt_coef,bsl_obs,bsl_coef,kp);
	% faults info
        [ flt1_out,flt2_out,subflt_name,subflt_out ]...
          = GTdef_slips(lb,ub,xx,flt1_num,flt1,flt2_num,flt2,flt3_name,flt3_num,flt3,flt4_name,flt4_num,flt4);
        [ mod_info(ii,:),pnt_out,bsl_out,nod_out ]...
           = GTdef_forward(Xgrn,Bgrn,Ngrn,sm_use,lb,ub,xx,pnt_loc,pnt_obs,pnt_obs_err,pnt_wgt,bsl_loc,bsl_obs,bsl_obs_err,bsl_wgt,nod_loc);
	% output results
        GTdef_output(fout_name,kp,rigidity,poisson,flt1_name,flt1_num,flt1_out,flt2_name,flt2_num,flt2_out,flt3_name,flt3_num,flt3,...
              	     flt4_name,flt4_num,flt4,subflt_name,subflt_out,pnt_name,pnt_num,pnt_out,bsl_name,bsl_num,bsl_out,...
               	     prf_name,prf_num,prf,grd_name,grd_num,grd,nod_name,nod_out,mod_info(ii,:));
    toc
    end
    fout_sum = strcat(basename,'_inv.out');
    GTdef_summary(fout_sum,kappa,mod_info);
end
