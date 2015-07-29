function [] = GTdef_output(filename,kappa,rigidity,poisson,...
           		   flt1_name,flt1_num,flt1,...
	  		   flt2_name,flt2_num,flt2,...
          		   flt3_name,flt3_num,flt3,...
	  		   flt4_name,flt4_num,flt4,...
	  		   subflt_name,subflt,...
          		   pnt_name,pnt_num,pnt_out,...
          		   bsl_name,bsl_num,bsl_out,...
           		   prf_name,prf_num,prf,...
	   		   grd_name,grd_num,grd,...
	  		   nod_name,nod_out,mod_info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_output.m				        %
% 		  function to output the model results				%
%										%
% INPUT:									%
% (1) Parameters: 								%
%  kappa 	scalar								%
%  rigidity	scalar			{30e9 Pa}               		%
%  poisson	scalar			{0.25}					%
%										%
% (2) Faults:									% 
%  flt1 - [lon lat z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]     	%
%  flt2 - [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]  	%
%  flt3 - [lon1 lat1 z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]	%
%  flt4 - [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]%
%  subflt_out - [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		  	%
%										%
% (3) Data:									% 
%  mod_info = [ data_num slip_num ndf rss rms wrrs wrms chi2 rchi2 roughness ]	%
%     slip_num - number of free slips					  	%
%     data_num - number of data points (not including nan)                 	%
%     ndf      - number of degrees of freedom                              	%
%     rss      - residual sum of squares                                   	%
%     rms      - root mean square of rss                                   	%
%     chi2     - chi-square                                                	%
%     rchi2    - reduced chi-square                                        	%
%  pnt_out  - [lon lat zz Ue Un Uv eUe eUn eUv weight]                    	%
%  bsl_out  - [lon1 lat1 z1 lon2 lat2 z2 Ue Un Uv Ul eUe eUn eUv eUl wgt] 	%
%  nod_out  - [lon lat zz Ue Un Uv eUe eUn eUv weight]			  	%
%										%
% OUTPUT: an output file                                                        %
%                                                                               %
% first created by Lujia Feng Wed May  6 20:58:39 EDT 2009			%
% last modified by Lujia Feng Tue May 12 16:33:16 EDT 2009			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout = fopen(filename,'w');
                                                                                
%%%%%%%%%% model info %%%%%%%%%%
fprintf(fout,'#data_num\t%d\n#slip_num\t%d\n#ndf     \t%d\n#rss     \t%-12.5e\n#rms     \t%-12.5e\n#wrss    \t%-12.5e\n#wrms    \t%-12.5e\n#chi2    \t%-12.5e\n#rchi2   \t%-12.5e\n#roughness\t%-12.5e\n',mod_info);

%%%%%%%%%% parameters %%%%%%%%%%
fprintf(fout,'kappa   \t%-12.5f\nrigidity\t%-10.2e\npoisson \t%-6.4f\n',kappa,rigidity,poisson);

%%%%%%%%%%  fault 1  %%%%%%%%%%
for ii =1:flt1_num
    fprintf(fout,'fault 1 %s  %-14.8f %-12.8f %-6.4e %12.4e %-12.4e %-5.2f %-5.2f  %10.5f %-8.5f %-8.5f  %-5.2f %-5.2f  %-5.2f %-5.2f  %-5.2f %-5.2f\n',flt1_name(ii,:),flt1(ii,:));
end

%%%%%%%%%%  fault 2  %%%%%%%%%%
for ii =1:flt2_num
    fprintf(fout,'fault 2 %s  %-14.8f %-12.8f %14.8f %-12.8f %12.4e %-12.4e %-5.2f  %10.5f %-8.5f %-8.5f  %-5.2f %-5.2f  %-5.2f %-5.2f  %-5.2f %-5.2f\n',flt2_name(ii,:),flt2(ii,:));
end

%%%%%%%%%%  fault 3  %%%%%%%%%%
for ii =1:flt3_num
    flt_name = flt3_name(ii,:);
    fprintf(fout,'fault 3 %s  %-14.8f %-12.8f %-6.4e %12.4e %-12.4e %-5.2f %-5.2f  %10.5f %-8.5f %-8.5f  %-5.2f %-5.2f  %-5.2f %-5.2f  %-5.2f %-5.2f    %d %d\n',flt3_name(ii,:),flt3(ii,:));
    ind = strmatch(flt_name,subflt_name,'exact');
    subflt_num = length(ind);
    for jj = 1:subflt_num
        fprintf(fout,'     subfault %s  %-d %-d  %10.5f %-8.5f %-8.5f  %-5.2f %-5.2f  %-5.2f %-5.2f  %-5.2f %-5.2f\n',flt_name,subflt(ind(jj),:));
    end
end

%%%%%%%%%%  fault 4  %%%%%%%%%%
for ii =1:flt4_num
    flt_name = flt4_name(ii,:);
    fprintf(fout,'fault 4 %s  %-14.8f %-12.8f %14.8f %-12.8f %12.4e %-12.4e %-5.2f  %10.5f %-8.5f %-8.5f  %-5.2f %-5.2f  %-5.2f %-5.2f  %-5.2f %-5.2f    %d %d\n',flt_name,flt4(ii,:));
    ind = strmatch(flt_name,subflt_name,'exact');
    subflt_num = length(ind);
    for jj = 1:subflt_num
        fprintf(fout,'     subfault %s  %-d %-d  %10.5f %-8.5f %-8.5f  %-5.2f %-5.2f  %-5.2f %-5.2f  %-5.2f %-5.2f\n',flt_name,subflt(ind(jj),:));
    end
end

%%%%%%%%%%  point  %%%%%%%%%%
for ii =1:pnt_num
    fprintf(fout,'point 3 %s     %-14.8f %-12.8f %-6.4e  %10.5f %-8.5f %-8.5f  %8.5f %-8.5f %-8.5f  %-5.2f\n', pnt_name(ii,:),pnt_out(ii,:));
end

%%%%%%%%%%  baseline  %%%%%%%%%%
for ii =1:bsl_num
    fprintf(fout,'baseline 3 %s     %-14.8f %-12.8f %-6.4e  %14.8f %-12.8f %-6.4e  %10.5f %-8.5f %-8.5f %-8.5f  %8.5f %-8.5f %-8.5f %-8.5f  %5.2f\n', bsl_name(ii,:),bsl_out(ii,:));
end

%%%%%%%%%%  profile  %%%%%%%%%%
for ii =1:prf_num
    name = prf_name(ii,:);
    fprintf(fout,'#profile %s  %-14.8f %-12.8f  %14.8f %-12.8f    %d\n',name,prf(ii,:));
    ind = strmatch(name,nod_name);
    nod_num = length(ind); 
    for jj = 1:nod_num
    	fprintf(fout,'point 3 %s     %-14.8f %-12.8f %-6.4e  %10.5f %-8.5f %-8.5f  %8.5f %-8.5f %-8.5f  %-5.2f\n', nod_name(ind(jj),:),nod_out(ind(jj),:));
    end
end

%%%%%%%%%%  grid  %%%%%%%%%%
for ii =1:grd_num
    name = grd_name(ii,:);
    fprintf(fout,'#grid %s %-5.2f %-5.2f  %14.8f %-12.8f  %10.4f %-8.4f    %d  %d\n',name,grd(ii,:));
    ind = strmatch(name,nod_name); 
    nod_num = length(ind);
    for jj = 1:nod_num
    	fprintf(fout,'point 3 %s     %-14.8f %-12.8f %-6.4e  %10.5f %-8.5f %-8.5f  %8.5f %-8.5f %-8.5f  %-5.2f\n',nod_name(ind(jj),:),nod_out(ind(jj),:));
    end
end
fclose(fout);
