function [ ] = GTdef_ckbd_mask(model_name,ckbdin_name,ckbdout_name,slip,percent,badnum,Nd,Ns,dd,ss)

% e.g. GTdef_ckbd_mask('NIC_30_40_bt16900.out','NIC_30_40_ckbd_fwd_surface.out','NIC_30_40_ckbdinv_bt10_surface.out','dd',0.2,3,30,40,3,4);
% e.g. GTdef_ckbd_mask('NIC_30_40_bt16900.out','NIC_30_40_ckbd2_fwd_surface.out','NIC_30_40_ckbdinv2_bt10_surface.out','dd',0.2,3,30,40,3,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           GPS_ckbd_mask.m				        %
% model_name - a regular output file, not a surface file			%  
% The other three input files are surface files					%
% slip    - slip type								%
%    ds - dip slip; ss - strike slip; ts - tensile slip                         %
% percent - if the difference larger than this percentage, 			%
%           this data point is ignored						%
% badnum  - number of bad sufault in each block					%
% Nd,Ns,dd,ss are the parameters that create the checkerboard array		%
%										%
% Note: subflt order is column-wise, can not be random				%
%										%
% first created by lfeng Tue Dec  7 08:20:54 EST 2010				%
% print out as a new output file lfeng Thu Apr 14 11:57:08 EDT 2011		%
% used structure lfeng Thu Feb 23 14:27:39 SGT 2012				%
% last modified by lfeng Thu Feb 23 14:44:03 SGT 2012				%
% no test has been done after last modification					%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ coord,origin,smooth,surf,beta,grnflag,...
  rigidity,poisson,...
  earth,edgrn,layer,...
  flt1,flt2,flt3,flt4,flt5,...
  bndry,subflt,dip,...
  pnt,bsl,prf,grd ] = GTdef_open(model_name);

fckbdin  = fopen(ckbdin_name,'r');
fckbdout = fopen(ckbdout_name,'r');

% keep column 2 3 & 15 (dip-slip)
ckbdin_cell  = textscan(fckbdin,'%*s %d %d %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %*f','CommentStyle','#');
ckbdout_cell = textscan(fckbdout,'%*s %d %d %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %*f','CommentStyle','#');
% column vectors
ckbdin_dnum = ckbdin_cell{1}; ckbdin_snum = ckbdin_cell{2};  ckbdin = ckbdin_cell{3};
ckbdout_dnum = ckbdout_cell{1}; ckbdout_snum = ckbdout_cell{2};  ckbdout = ckbdout_cell{3};

fclose(fckbdin); fclose(fckbdout); 

%  subflt - [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		  	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mask out %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(slip,'ss')		% identify the column for slip
    colnum = 3;
elseif strcmpi(slip,'ds')
    colnum = 4;
elseif strcmpi(slip,'ts')
    colnum = 5;
else
    error('GTdef_ckbd_mask ERROR: Slip is not correctly specified!');
end
slipmax = max(abs(ckbdin));	% for checkerboard, maximum slips are the same

diff = bsxfun(@minus,ckbdin,ckbdout);
diff =abs(diff)/slipmax;

Bdnum = ceil(Nd/dd);		% block num each column
Bsnum = ceil(Ns/ss);		% block num each row
Bdiff = reshape(diff,Nd,Ns);
Bmodel = reshape(subflt.flt(:,colnum),Nd,Ns);

flag = zeros(Bdnum,Bsnum);
for ii=0:Bdnum-1
   for jj=0:Bsnum-1
       % exclude 4 corners
       if (ii==0&&jj==0) || (ii==0&&jj==Bsnum-1) || (ii==Bdnum-1&&jj==0) || (ii==Bdnum-1&&jj==Bsnum-1) 
          Bmodel(ii*dd+1:ii*dd+dd,jj*ss+1:jj*ss+ss) = NaN;
	  flag(ii+1,jj+1) = 1;	% mark as bad
	  continue
       end
       ind = Bdiff(ii*dd+1:ii*dd+dd,jj*ss+1:jj*ss+ss)>percent;
       num = sum(sum(ind));
       if num>badnum
          Bmodel(ii*dd+1:ii*dd+dd,jj*ss+1:jj*ss+ss) = NaN;
	  flag(ii+1,jj+1) = 1;	% mark as bad
       end
   end
end

% identify isolated good blocks
for ii=1:Bdnum
   for jj=1:Bsnum
       if flag(ii,jj)==1, continue; end
       total = 0;
       for nn=-1:1
           for mm=-1:1
	      if (ii+nn)<=0, continue; end
	      if (ii+nn)>Bdnum, continue; end
	      if (jj+mm)<=0, continue; end
	      if (jj+mm)>Bsnum, continue; end
	      if flag(ii+nn,jj+mm)==1, total = total+1; end
	   end
       end
       if total>4
          Bmodel((ii-1)*dd+1:(ii-1)*dd+dd,(jj-1)*ss+1:(jj-1)*ss+ss) = NaN;
       end
   end
end
model = reshape(Bmodel,[],1);
subflt.flt(:,colnum) = model;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ ~,basename,~ ] = fileparts(model_name);
fout_name = [ basename '_mask.out' ];

pnt.num = 0; bsl.num = 0; prf.num = 0; grd.num = 0;

GTdef_output(fout_name,coord,'none','none',0,rigidity,poisson,earth,edgrn,layer,...
    	     flt1,flt2,flt3,flt4,flt5,bndry,subflt,dip,pnt,bsl,prf,grd,nod,mod_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output magnitude estimate for Nicoya %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output magnitude estimate
% for Nicoya using 5x5 km^2 patches
ind = ~isnan(model);
% number of useful patches
patchnum = sum(ind)
% total slip
sliptotal = -sum(model(ind))
rigidity
M0 = rigidity*5e3*5e3*sliptotal
% from N-m to dyn-cm
M0 = M0*1e7;
% one year
Wm = 2/3*(log10(M0)-16.1)
% 60 years
M0_60 = M0*60
Wm = 2/3*(log10(M0_60)-16.1)
