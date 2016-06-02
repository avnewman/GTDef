function [ flt1,flt2,flt3,flt4,flt5,flt6 ] = GTdef_calc_stressdrop(earth,flt1,flt2,flt3,flt4,flt5,flt6) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          GTdef_cal_stressdrop                           %
% calculate stress drop for each rupture                                  %
%                                                                         %
% INPUT:                                                                  %
% earth structure                                                         %
% flt?  structure 							  %
%                                                                         %
% INTERMEDIATE:                                                           %
% disp   = [ Ux Uy Uz ] (3 row vectors)					  %
% strain = [exx eyy ezz eyz exz exy] (6 row vectors)                      %
% stress = [sxx syy szz syz sxz sxy] (6 row vectors)                      %
% tilt (2 row vectors) 							  %
%                                                                         % 
% OUTPUT:                                                                 %
% Fault Structures: flt?                                                  %
% flt?.sdrop - stress drop for each main fault                            %
% slip-averaged stress drop                                               %
% averaged by the absolute value of slips                                 %
%                                                                         % 
%         /                                                               %
%         | tau * s dA                                                    %
%         /                                                               %
%   dtau = --------------                                                 %
%             /                                                           %
%             | s dA                                                      %
%             /                                                           %
%                                                                         % 
% where tau is the amplitude of shear stress on the fault, s is the       %
% coseismic slip on the fault                                             %
% ----------------------------------------------------------------------- % 
% It has been have verified by unicycle utility                           %
% addpath ~lfeng/Mercurial/unicycle/matlab                                %
% import unicycle.*                                                       %
% sourceFilename = '.flt'                                                 %
% G    = 30e3                                                             %
% nu   = 0.25                                                             %
% dtau = unicycle.utils.computeStressDrop(sourceFilename,G,nu)            %
% ----------------------------------------------------------------------- % 
%                                                                         % 
% Reference:                                                              %
% Noda, H., N. Lapusta, and H. Kanamori (2013),                           %
% Comparison of average stress drop measures for ruptures                 %
% with heterogeneous stress change and implications for                   %
% earthquake physics, Geophys. J. Int., 193(3),                           %
% 1691–1712, doi:10.1093/gji/ggt074.                                      %
%                                                                         %
% first created by Lujia Feng Fri Mar 20 10:44:10 SGT 2015                %
% added Min, SSgrn, DSgrn, and TSgrn to xyzflt lfeng Thu Mar 26 SGT 2015  %
% added fault5 lfeng Tue Jun 23 17:20:20 SGT 2015                         %
% added fault6 lfeng Thu Jun  2 10:31:29 SGT 2016                         %
% last modified by Lujia Feng Thu Jun  2 10:36:20 SGT 2016                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

etype = earth.type;

% fault type 1
if flt1.num~=0
   for ii=1:flt1.num
      %                1   2     3     4   5   6    7     8  9  10
      % Okada   Min = [len width depth dip str east north ss ds ts]     [flt_num*10]
      %                1    2     3    4     5      6     7   8   9
      % layered Min = [slip north east depth length width str dip rake] [flt_num*9]
      Min  = flt1.xyzflt{ii}.Min;
      Xin  = flt1.xyzflt{ii}.xyzctr;
      nuv  = flt1.xyzflt{ii}.nuv;
      suv  = flt1.xyzflt{ii}.suv;
      duv  = flt1.xyzflt{ii}.duv;
      
      if strcmpi(etype,'homogeneous')
         ss = Min(:,8);
         ds = Min(:,9);
      else
         slip = Min(:,1);
         rake = Min(:,9);
         ss   = slip.*cosd(rake); 
         ds   = slip.*sind(rake);
      end

      patchnum = size(Min,1);
      disp   = 0;
      strain = 0;
      stress = 0;
      tilt   = 0;

      [ disp,strain,stress,tilt ] = GTdef_calc(earth,Min,Xin);
      patchnum = size(stress,1);
      slipSum   = 0;
      sliptsSum = 0;
      for jj=1:patchnum
          % slip vector
          slipV = suv(jj,:).*ss(jj,:) + duv(jj,:).*ds(jj,:);

          % 6 independent stress components = [sxx syy szz syz sxz sxy]
          ss6 = stress(jj,:);

          % full stress tensor = [ sxx sxy sxz ]
          %                      [ syx syy syz ]
          %                      [ szx szy szz ]
          ssT = diag(ss6(1:3))+diag(ss6([6 4]),-1)+diag(ss6([6 4]),1)+diag(ss6(5),-2)+diag(ss6(5),2);       

	  % normal direction to patches = [ nxx nyy nzz ]
          nnV = nuv(jj,:); 

          % full traction vector on patches [ txx tyy tzz ] = [ nxx nyy nzz ]*[ sxx sxy sxz ]
          % t = n*s                                                           [ syx syy syz ]
          %                                                                   [ szx szy szz ]
          ttV = nnV*ssT;
          % matrix notation is the same as the equation below
          % ttV = [ sum(ssT(:,1).*nuv(jj,:)') sum(ssT(:,2).*nuv(jj,:)') sum(ssT(:,3).*nuv(jj,:)') ];

          % normal traction vector on patches
          tnV = sum(ttV.*nnV)*nnV;

          % shear traction vector on patches
          tsV = ttV - tnV;

          % sum up
          slipSum   = slipSum   + sqrt(sum(slipV.*slipV));
          %sliptsSum = sliptsSum + sum(slipV.*tsV); wrong
          sliptsSum = sliptsSum + sqrt(sum(slipV.*slipV))*sqrt(sum(tsV.*tsV)); 
      end
      % dtau = sum(slip.*ts)/sum(slip);
      flt1.sdrop(ii) = sliptsSum/slipSum;
   end
end

% fault type 2
if flt2.num~=0
   for ii=1:flt2.num
      Min  = flt2.xyzflt{ii}.Min;
      Xin  = flt2.xyzflt{ii}.xyzctr;
      nuv  = flt2.xyzflt{ii}.nuv;
      suv  = flt2.xyzflt{ii}.suv;
      duv  = flt2.xyzflt{ii}.duv;

      if strcmpi(etype,'homogeneous')
         ss = Min(:,8);
         ds = Min(:,9);
      else
         slip = Min(:,1);
         rake = Min(:,9);
         ss   = slip.*cosd(rake); 
         ds   = slip.*sind(rake);
      end

      [ disp,strain,stress,tilt ] = GTdef_calc(earth,Min,Xin);
      patchnum = size(stress,1);
      slipSum   = 0;
      sliptsSum = 0;
      for jj=1:patchnum
          slipV = suv(jj,:).*ss(jj,:) + duv(jj,:).*ds(jj,:);
          ss6 = stress(jj,:);
          ssT = diag(ss6(1:3))+diag(ss6([6 4]),-1)+diag(ss6([6 4]),1)+diag(ss6(5),-2)+diag(ss6(5),2);       
          nnV = nuv(jj,:); 
          ttV = nnV*ssT;
          tnV = sum(ttV.*nnV)*nnV;
          tsV = ttV - tnV;
          slipSum   = slipSum   + sqrt(sum(slipV.*slipV));
          sliptsSum = sliptsSum + sqrt(sum(slipV.*slipV))*sqrt(sum(tsV.*tsV)); 
      end
      flt2.sdrop(ii) = sliptsSum/slipSum;
   end
end

% fault type 3
if flt3.num~=0
   for ii=1:flt3.num
      Min  = flt3.xyzflt{ii}.Min;
      Xin  = flt3.xyzflt{ii}.xyzctr;
      nuv  = flt3.xyzflt{ii}.nuv;
      suv  = flt3.xyzflt{ii}.suv;
      duv  = flt3.xyzflt{ii}.duv;

      if strcmpi(etype,'homogeneous')
         ss = Min(:,8);
         ds = Min(:,9);
      else
         slip = Min(:,1);
         rake = Min(:,9);
         ss   = slip.*cosd(rake); 
         ds   = slip.*sind(rake);
      end

      [ disp,strain,stress,tilt ] = GTdef_calc(earth,Min,Xin);
      patchnum = size(stress,1);
      slipSum   = 0;
      sliptsSum = 0;
      for jj=1:patchnum
          slipV = suv(jj,:).*ss(jj,:) + duv(jj,:).*ds(jj,:);
          ss6 = stress(jj,:);
          ssT = diag(ss6(1:3))+diag(ss6([6 4]),-1)+diag(ss6([6 4]),1)+diag(ss6(5),-2)+diag(ss6(5),2);       
          nnV = nuv(jj,:); 
          ttV = nnV*ssT;
          tnV = sum(ttV.*nnV)*nnV;
          tsV = ttV - tnV;
          slipSum   = slipSum   + sqrt(sum(slipV.*slipV));
          sliptsSum = sliptsSum + sqrt(sum(slipV.*slipV))*sqrt(sum(tsV.*tsV)); 
      end
      flt3.sdrop(ii) = sliptsSum/slipSum;
   end
end

% fault type 4
if flt4.num~=0
   for ii=1:flt4.num
      Min  = flt4.xyzflt{ii}.Min;
      Xin  = flt4.xyzflt{ii}.xyzctr;
      nuv  = flt4.xyzflt{ii}.nuv;
      suv  = flt4.xyzflt{ii}.suv;
      duv  = flt4.xyzflt{ii}.duv;

      if strcmpi(etype,'homogeneous')
         ss = Min(:,8);
         ds = Min(:,9);
      else
         slip = Min(:,1);
         rake = Min(:,9);
         ss   = slip.*cosd(rake); 
         ds   = slip.*sind(rake);
      end

      [ disp,strain,stress,tilt ] = GTdef_calc(earth,Min,Xin);
      patchnum = size(stress,1);
      slipSum   = 0;
      sliptsSum = 0;
      for jj=1:patchnum
          slipV = suv(jj,:).*ss(jj,:) + duv(jj,:).*ds(jj,:);
          ss6 = stress(jj,:);
          ssT = diag(ss6(1:3))+diag(ss6([6 4]),-1)+diag(ss6([6 4]),1)+diag(ss6(5),-2)+diag(ss6(5),2);       
          nnV = nuv(jj,:); 
          ttV = nnV*ssT;
          tnV = sum(ttV.*nnV)*nnV;
          tsV = ttV - tnV;
          slipSum   = slipSum   + sqrt(sum(slipV.*slipV));
          sliptsSum = sliptsSum + sqrt(sum(slipV.*slipV))*sqrt(sum(tsV.*tsV)); 
      end
      flt4.sdrop(ii) = sliptsSum/slipSum;
   end
end

% fault type 5 is the same as fault type 1
if flt5.num~=0
   for ii=1:flt5.num
      %                1   2     3     4   5   6    7     8  9  10
      % Okada   Min = [len width depth dip str east north ss ds ts]     [flt_num*10]
      %                1    2     3    4     5      6     7   8   9
      % layered Min = [slip north east depth length width str dip rake] [flt_num*9]
      Min  = flt5.xyzflt{ii}.Min;
      Xin  = flt5.xyzflt{ii}.xyzctr;
      nuv  = flt5.xyzflt{ii}.nuv;
      suv  = flt5.xyzflt{ii}.suv;
      duv  = flt5.xyzflt{ii}.duv;
      
      if strcmpi(etype,'homogeneous')
         ss = Min(:,8);
         ds = Min(:,9);
      else
         slip = Min(:,1);
         rake = Min(:,9);
         ss   = slip.*cosd(rake); 
         ds   = slip.*sind(rake);
      end

      patchnum = size(Min,1);
      disp   = 0;
      strain = 0;
      stress = 0;
      tilt   = 0;

      [ disp,strain,stress,tilt ] = GTdef_calc(earth,Min,Xin);
      patchnum = size(stress,1);
      slipSum   = 0;
      sliptsSum = 0;
      for jj=1:patchnum
          slipV = suv(jj,:).*ss(jj,:) + duv(jj,:).*ds(jj,:);
          ss6 = stress(jj,:);
          ssT = diag(ss6(1:3))+diag(ss6([6 4]),-1)+diag(ss6([6 4]),1)+diag(ss6(5),-2)+diag(ss6(5),2);       
          nnV = nuv(jj,:); 
          ttV = nnV*ssT;
          tnV = sum(ttV.*nnV)*nnV;
          tsV = ttV - tnV;
          slipSum   = slipSum   + sqrt(sum(slipV.*slipV));
          sliptsSum = sliptsSum + sqrt(sum(slipV.*slipV))*sqrt(sum(tsV.*tsV)); 
      end
      flt5.sdrop(ii) = sliptsSum/slipSum;
   end
end

% fault type 6 is the same as fault type 3
if flt6.num~=0
   for ii=1:flt6.num
      Min  = flt6.xyzflt{ii}.Min;
      Xin  = flt6.xyzflt{ii}.xyzctr;
      nuv  = flt6.xyzflt{ii}.nuv;
      suv  = flt6.xyzflt{ii}.suv;
      duv  = flt6.xyzflt{ii}.duv;

      %                1   2     3     4   5   6    7     8  9  10
      % Okada   Min = [len width depth dip str east north ss ds ts]     [flt_num*10]
      %                1    2     3    4     5      6     7   8   9
      % layered Min = [slip north east depth length width str dip rake] [flt_num*9]
      if strcmpi(etype,'homogeneous')
         ss = Min(:,8);
         ds = Min(:,9);
      else
         slip = Min(:,1);
         rake = Min(:,9);
         ss   = slip.*cosd(rake); 
         ds   = slip.*sind(rake);
      end

      [ disp,strain,stress,tilt ] = GTdef_calc(earth,Min,Xin);
      patchnum = size(stress,1);
      slipSum   = 0;
      sliptsSum = 0;
      for jj=1:patchnum
          slipV = suv(jj,:).*ss(jj,:) + duv(jj,:).*ds(jj,:);
          ss6 = stress(jj,:);
          ssT = diag(ss6(1:3))+diag(ss6([6 4]),-1)+diag(ss6([6 4]),1)+diag(ss6(5),-2)+diag(ss6(5),2);       
          nnV = nuv(jj,:); 
          ttV = nnV*ssT;
          tnV = sum(ttV.*nnV)*nnV;
          tsV = ttV - tnV;
          slipSum   = slipSum   + sqrt(sum(slipV.*slipV));
          sliptsSum = sliptsSum + sqrt(sum(slipV.*slipV))*sqrt(sum(tsV.*tsV)); 
      end
      flt6.sdrop(ii) = sliptsSum/slipSum;
   end
end
