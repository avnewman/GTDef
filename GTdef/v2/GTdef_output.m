function [] = GTdef_output(filename,...
			   coord,smooth,surf,beta,rigidity,poisson,...
		           earth,edgrn,layer,...
           		   flt1,flt2,flt3,flt4,flt5,...   
			   subflt,dip,...
          		   pnt,bsl,prf,grd,nod,mod_info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_output.m				                 %
% 		  function to output the model results				         %
%										         %
% INPUT:									         %
% (1) Controlling Parameters: 							         %
%  coord        (string)                {geo}					         %
%  smooth	(string)		{2d}					         %
%  surf         (string)		{free}		   			         %
%  beta 	(scalar)							         %
%										         %
% (2) Earth Structure:								         %
% Either of the two types of earth structure can be used.			         %
%  earth = homogeneous								         %
%     		rigidity		(scalar)		{30e9 Pa}                %
%		poisson			(scalar)		{0.25}		         %
%  earth = layered		        					         %
%		edgrn.nl        	(scalar)                                         %
% 		edgrn.obsz     		(scalar)                                         %
% 		edgrn.nr	        (scalar)                                         %
% 		edgrn.minr,edgrn.maxr   (scalar)                                         %
% 		edgrn.nz                (scalar)                                         %
% 		edgrn.minz,edgrn.maxz   (scalar)                                         %
%		edgrn.srate		(scalar)				         %	
%    		layer - [ id depth vp vs ro ]	(nn*5)				         %
%										         %
% (3) Faults:									         % 
% each fault has a structure to store corresponding data			         %
% flt? structure: flt?.name flt?.num flt?.flt					         %
% subflt structure: subflt.name subflt.num subflt.flt subflt.out subflt.outname	         %
% flt1.flt - [lon1 lat1 z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]        %
% flt2.flt - [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]      %
%      subflt.flt - [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		  	 %
% flt3.flt - [lon1 lat1 z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns ] %
% flt4.flt - [lon1 lat1 lon2 lat2 z1 z2 dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns]%
%      subflt.flt - [ dnum snum rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]                 %
% dip structure: dip.name dip.num & dip.dip					         %
% dip.dip  - [ dip z1 z2 rows ]						                 %
%                                                                                        %
% OUTPUT: flt1.out flt2.out subflt.out subflt.outname				         %
%										         %
% (4) Data:									         % 
%  mod_info = [ data_num slip_num ndf rss rms wrrs wrms chi2 rchi2            	         %
%               r_1d r_2d strain ]				  		         %
%     slip_num - number of free slips					  	         %
%     data_num - number of data points (not including nan)                 	         %
%     ndf      - number of degrees of freedom                              	         %
%  Note: since if we introduce smoothing, slips are not independent	  	         %  
%  we don't really know the real ndf.					  	         %
%     rss      - residual sum of squares [m^2]                                           %
%     rms      - root mean square of rss [m]                                  	         %
%     wrss     - weighted residual sum of squares [m^2]			  	         %
%     wrms     = sqrt(wrss/data_num) [m]					         %
%     chi2     - chi-square                                                	         %
%     rchi2    - reduced chi-square                                        	         %
%     r_1d     - the average 1st derivative sum of each patch [cm/km]		         %
%                depends on the finite difference method used		  	         %
%     r_2d     - the average 2nd derivative sum of each patch		                 %
%                (eq. 5 in Jonsson_etal_BSSA_2002) [cm/km^2]		                 %
%                 r represents roughness					         %
%    strain   - average strain of each patch [cm/km]				         %
%                (absolute value of 1st derivaitves)			  	         %
%  pnt.out  - [lon lat zz Ue Un Uv eUe eUn eUv weight]                    	         %
%  bsl.out  - [lon1 lat1 z1 lon2 lat2 z2 Ue Un Uv Ul eUe eUn eUv eUl wgt] 	         %
%  nod.out  - [lon lat zz Ue Un Uv eUe eUn eUv weight]			  	         %
%										         %
% OUTPUT: an output file                                                                 %
%                                                                                        %
% first created by Lujia Feng Wed May  6 20:58:39 EDT 2009			         %
% added beta lfeng Wed Dec  2 23:42:51 EST 2009					         %
% added 1st derivative r_1d lfeng Thu Dec  3 01:27:48 EST 2009			         %
% added flag 'dip' by lfeng Mon Dec  7 01:05:28 EST 2009		  	         %
% added fault5 by lfeng Fri Dec 11 13:00:11 EST 2009				         %
% corrected the wrong 'freesurface' flag and changed it to 'surface' flag	         %
%    lfeng Wed Feb 24 13:26:18 EST 2010					                 %
% added the units in the output	lfeng Wed Jul 21 17:06:07 EDT 2010		         %
% use cell array of strings for names lfeng Wed Dec  1 17:36:42 EST 2010	         %
% use 4 digits after . for slip constraints lfeng Thu Dec  9 03:52:06 EST 2010           %
% test existence before output lfeng Thu Apr 14 12:58:16 EDT 2011		         %
% used structure & added layered output lfeng Wed Feb 22 14:31:31 SGT 2012	         %
% merged fault1 & fault3 and fault2 & fault4 lfeng Thu May 10 17:03:16 SGT 2012          %
% modified fault5 lfeng Sat Dec  1 00:32:44 SGT 2012                                     %
% last modified by Lujia Feng Sat Dec  1 22:52:48 SGT 2012                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fout = fopen(filename,'w');
                                                                                
%%%%%%%%%% model info %%%%%%%%%%
if ~isempty(mod_info)
    fprintf(fout,'#data_num\t%d\n#slip_num\t%d\n#ndf     \t%d\n#rss      \t%-12.5e [m^2]\n#rms      \t%-12.5e [m]\n#wrss   \t%-12.5e [m^2]\n#wrms   \t%-12.5e [m]\n#chi2    \t%-12.5e\n#rchi2   \t%-12.5e\n#r_1d   \t%-12.5e [cm/km]\n#r_2d   \t%-12.5e [cm/km^2]\n#strain   \t%-12.5e [cm/km]\n',mod_info);
    fprintf(fout,'\n');
end

%%%%%%%%%% smooth parameters %%%%%%%%%%
if ~isempty(beta)
    kappa = sqrt(beta);
    fprintf(fout,'kappa   \t%-12.5f\nbeta     \t%-12.5f\n\n',kappa,beta);
end

%%%%%%%%%% earth model %%%%%%%%%%
if strcmpi(earth,'homogeneous')||strcmpi(earth,'homo')
    fprintf(fout,'earth\thomogeneous\t%-10.2e\t%-6.4f\n',rigidity,poisson);
    fprintf(fout,'\n');
elseif strcmpi(earth,'layered')
    fprintf(fout,'earth\tlayered   %.0f %-8.4e  %.0f %-8.4e %-8.4e  %.0f %-8.4e %-8.4e   %.4f\n',...
            edgrn.nl,edgrn.obsz,edgrn.nr,edgrn.minr,edgrn.maxr,edgrn.nz,edgrn.minz,edgrn.maxz,edgrn.srate);
    fprintf(fout,'\n');
end

%%%%%%%%%% layers %%%%%%%%%%
if ~isempty(layer)
    fprintf(fout,'layer   %.0f   %8.4e  %-8.4e  %-8.4e  %-8.4e\n',layer');
    fprintf(fout,'\n');
end

%%%%%%%%%% smooth parameters %%%%%%%%%%
if ~isempty(coord)
    fprintf(fout,'coord   \t%s\n',coord);
end
if ~isempty(smooth)
    fprintf(fout,'smooth   \t%s\n',smooth);
end
if ~isempty(surf)
    fprintf(fout,'surface  \t%s\n',surf);
end
fprintf(fout,'\n');

%%%%%%%%%%  fault 1  %%%%%%%%%%
% note: output flt1.out
for ii =1:flt1.num
    flt_name = flt1.name{ii};
    fprintf(fout,'fault 1 %s  %-14.8f %-12.8f %-6.4e %12.4e %-12.4e %-5.2f %-5.2f  %10.5f %-8.5f %-8.5f  %-5.4f %-5.4f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d %d\n',flt_name,flt1.out(ii,:));
    ind = strcmpi(flt_name,subflt.outname);
    num = sum(ind);
    subflt1 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %10.5f %10.5f %10.5f  %8.4f %8.4f  %8.4f %8.4f  %8.4f %8.4f\n',flt_name,subflt1(jj,:));
    end
end

%%%%%%%%%%  fault 2  %%%%%%%%%%
% note: output flt2.out
for ii =1:flt2.num
    flt_name = flt2.name{ii};
    fprintf(fout,'fault 2 %s  %-14.8f %-12.8f %14.8f %-12.8f %12.4e %-12.4e %-5.2f  %10.5f %-8.5f %-8.5f  %-5.4f %-5.4f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d %d\n',flt_name,flt2.out(ii,:));
    ind = strcmpi(flt_name,subflt.outname);
    num = sum(ind);
    subflt2 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %10.5f %10.5f %10.5f  %8.4f %8.4f  %8.4f %8.4f  %8.4f %8.4f\n',flt_name,subflt2(jj,:));
    end
end

%%%%%%%%%%  fault 3  %%%%%%%%%%
for ii =1:flt3.num
    flt_name = flt3.name{ii};
    fprintf(fout,'fault 3 %s  %-14.8f %-12.8f %-6.4e %12.4e %-12.4e %-5.2f %-5.2f  %10.2f %-8.5f %-8.5f  %-7.2f %-7.2f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d %d\n',flt_name,flt3.out(ii,:));
    ind = strcmpi(flt_name,subflt.outname);
    num = sum(ind);
    subflt3 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %10.2f %12.5f %12.5f  %7.2f %7.2f  %8.4f %8.4f  %8.4f %8.4f\n',flt_name,subflt3(jj,:));
    end
end

%%%%%%%%%%  fault 4  %%%%%%%%%%
for ii =1:flt4.num
    flt_name = flt4.name{ii};
    fprintf(fout,'fault 4 %s  %-14.8f %-12.8f %14.8f %-12.8f %12.4e %-12.4e %-5.2f  %10.2f %-8.5f %-8.5f  %-7.2f %-7.2f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d %d\n',flt_name,flt4.out(ii,:));
    ind = strcmpi(flt_name,subflt.outname);
    num = sum(ind);
    subflt4 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %10.2f %12.5f %12.5f  %7.2f %7.2f  %8.4f %8.4f  %8.4f %8.4f\n',flt_name,subflt4(jj,:));
    end
end

%%%%%%%%%%  fault 5  %%%%%%%%%%
for ii =1:flt5.num
    flt_name = flt5.name{ii};
    grn_name = flt5.grname{ii};
    fprintf(fout,'fault 5   %s   %s   %10.5f %10.5f %10.5f  %5.4f %5.4f  %5.4f %5.4f  %5.4f %5.4f    %d  %d \n',flt_name,grn_name,flt5.out(ii,:));
    ind = strcmpi(flt_name,subflt.outname);
    num = sum(ind);
    subflt5 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %12.5f %12.5f %12.5f  %8.4f %8.4f  %8.4f %8.4f  %8.4f %8.4f\n',flt_name,subflt5(jj,:));
    end
end

%%%%%%%%%%   dip   %%%%%%%%%%
for ii =1:dip.num
    fprintf(fout,'dip   %s  %8.4f  %-12.4e  %-12.4e  %-d\n',dip.name{ii},dip.dip(ii,:));
end

%%%%%%%%%%  point  %%%%%%%%%%
if strcmp(coord,'local')
    for ii =1:pnt.num
        fprintf(fout,'point 3 %s\t%14.5e %14.5e %14.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n', pnt.name{ii},pnt.out(ii,:));
    end
else
    for ii =1:pnt.num
        fprintf(fout,'point 3 %s\t%14.8f %12.8f %12.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n', pnt.name{ii},pnt.out(ii,:));
    end
end

%%%%%%%%%%  baseline  %%%%%%%%%%
if strcmp(coord,'local')
    for ii =1:bsl.num
        fprintf(fout,'baseline 3 %s\t%14.5e %14.5e %14.5e  %14.8f %12.8f %6.4e  %10.5f %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f %10.5f  %5.2f\n', bsl.name{ii},bsl.out(ii,:));
    end
else
    for ii =1:bsl.num
        fprintf(fout,'baseline 3 %s\t%14.8f %12.8f %12.5e  %14.8f %12.8f %6.4e  %10.5f %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f %10.5f  %5.2f\n', bsl.name{ii},bsl.out(ii,:));

    end
end

%%%%%%%%%%  profile  %%%%%%%%%%
for ii =1:prf.num
    name = prf.name{ii};
    name_len = length(name);
    fprintf(fout,'#profile %s  %-14.8f %-12.8f  %14.8f %-12.8f    %d\n',name,prf.prf(ii,:));
    ind = strncmpi(name,nod.name,name_len);
    num = sum(ind); 
    cnod_name = nod.name(ind);
    cnod_out  = nod.out(ind,:);
    if strcmp(coord,'local')
        for jj = 1:num
            fprintf(fout,'point 3 %s\t%14.5e %14.5e %14.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n',cnod_name{jj},cnod_out(jj,:));
        end
    else
        for jj = 1:num
            fprintf(fout,'point 3 %s\t%14.8f %12.8f %12.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n',cnod_name{jj},cnod_out(jj,:));
        end
    end
end

%%%%%%%%%%  grid  %%%%%%%%%%
for ii =1:grd.num
    name = grd.name{ii};
    name_len = length(name);
    fprintf(fout,'#grid %s %-5.2f %-5.2f    %12.8f %-12.8f  %10.4f %-8.4f    %d  %d\n',name,grd.grd(ii,:));
    ind = strncmpi(name,nod.name,name_len);
    num = sum(ind); 
    cnod_name = nod.name(ind);
    cnod_out  = nod.out(ind,:);
    if strcmp(coord,'local')
        for jj = 1:num
        	fprintf(fout,'point 3 %s\t%14.5e %14.5e %14.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n',cnod_name{jj},cnod_out(jj,:));
        end
    else
        for jj = 1:num
        	fprintf(fout,'point 3 %s\t%14.8f %12.8f %12.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n',cnod_name{jj},cnod_out(jj,:));
    
        end
    end
end
fclose(fout);
