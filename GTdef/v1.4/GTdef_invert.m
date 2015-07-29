function [ xx ] = GTdef_invert(Xgrn,Bgrn,sm,Aeq,beq,lb,ub,x0,...
            		   pnt_obs,pnt_coef,bsl_obs,bsl_coef,beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_invert				  %
% Invert the slip values using						  %
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)					  %
% Here we have no inequalities, so set A=[];b=[]			  %
%									  %
% INPUT:					  		  	  %
% Each fault has three (strike, dip, and tensile) components, so          %
% slip_num = flt_num*3                                                    %
%  Xgrn - displacements [east;north;vertical] for different sites   	  %
%         from unit slips [(3*nn)*slip_num] 				  %
%         (nn is the  number of sites)  				  %
%  Bgrn - length changes [east;north;vertical;absolute] for 	  	  %
%         different baselines from unit slips [(4*nn)*slip_num] 	  %
%         (nn is the  number of baselines)  				  %
%  sm  - condensed smoothing matrix with rows of all zeros removed	  %
%  Aeq - left-hand side matrix for linear equalities  [slip_num*slip_num] %
%  beq - right-hand side vector for linear equalities [slip_num*1]        %
%  x0  - initial values for ss,ds,ts 	[slip_num*1]                      %
%  lb  - lower bounds for ss,ds,ts 	[slip_num*1]                      %
%  ub  - upper bounds for ss,ds,ts	[slip_num*1]			  %
%									  %
% Point:    pnt_obs  - [east;north;vert]	     (3nn*1)		  %
%   	    pnt_coef - [east;north;vert]	     (3nn*1)    	  %
% Baseline: bsl_obs  - [east;north;vert;length]      (4nn*1)    	  %
%   	    bsl_coef - [east;north;vert;length]      (4nn*1)    	  %
%  beta - smoothing parameter			     scalar		  %
%           beta = kappa^2						  %
%									  %
% OUTPUT:                                                                 %
%  xx - final values for ss,ds,ts from inversion     [slip_num*1]	  %
%                                                                         %
% first created by Lujia Feng Tue May  5 20:21:47 EDT 2009		  %
% last modified by Lujia Feng Mon May 11 16:40:22 EDT 2009		  %
% used beta = beta^2 instead of beta lfeng Wed Dec  2 23:35:12 EST 2009   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = []; d = [];
% add error & weight to the equations
if ~isempty(Xgrn)
    eq_num = size(Xgrn,1);
    for ii = 1:eq_num
        Xgrn(ii,:) = Xgrn(ii,:).*pnt_coef(ii);
	pnt_obs(ii) = pnt_obs(ii).*pnt_coef(ii);
    end
    ind = find(~isnan(pnt_obs));		% exclude nan values
    C = [ C;Xgrn(ind,:) ]; d = [ d;pnt_obs(ind) ];
end
if ~isempty(Bgrn)
    eq_num = size(Bgrn,1);
    for ii = 1:eq_num
        Bgrn(ii,:) = Bgrn(ii,:).*bsl_coef(ii,1);
	bsl_obs(ii,1) = bsl_obs(ii,1).*bsl_coef(ii,1);
    end
    ind = find(~isnan(bsl_obs));		% exclude nan values
    C = [ C;Bgrn(ind,:) ]; d = [ d;bsl_obs(ind) ];
end

% add beta to the smoothing matrix
sm = sm.*beta;

% combine the green matrices and the smoothing matrix
d_sm = zeros(size(sm,1),1);
C = [ C;sm ]; d = [ d;d_sm ];

ind = find(isinf(lb));
Aeq_red = Aeq(ind,:); beq_red = beq(ind);

options = optimset('MaxIter',1000,'TolFun',1e-30);
[xx,resnorm] = lsqlin(C,d,[],[],Aeq_red,beq_red,lb,ub,x0,options);
%[xx,resnorm,res,exitflag,output] = lsqlin(C,d,[],[],Aeq_red,beq_red,lb,ub,x0);
fprintf(1,'resnorm = %-12.5e\n',resnorm);
