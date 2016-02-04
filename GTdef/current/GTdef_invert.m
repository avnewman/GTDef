function [ modspace ] = GTdef_invert(modspace,pnt,los,bsl,beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_invert				  %
% Invert the slip values using						  %
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)					  %
% Here we have no inequalities, so set A=[];b=[]			  %
%									  %
% INPUT:					  		  	  %
% modspace - model structure                                              %
% pnt      - point structure	  	                                  %
% los      - los structure                                                %
% bsl      - baseline structure                                           %
% beta     - current beta value                                           %
%									  %
% OUTPUT:                                                                 %
% add .xx to modspace                                                     %
% modspace.xx - final values for ss,ds,ts from inversion     [slip_num*1] %
%                                                                         %
% first created by Lujia Feng Tue May  5 20:21:47 EDT 2009		  %
% used beta = beta^2 instead of beta lfeng Wed Dec  2 23:35:12 EST 2009   %
% used structure lfeng Wed Feb 22 19:36:43 SGT 2012			  %
% added modspace structure lfeng Thu Mar 19 17:32:29 SGT 2015             %
% added condensing sm & sm_abs lfeng Fri Mar 20 20:28:15 SGT 2015         %
% added InSAR los & Lgrn lfeng Tue Nov  3 23:17:37 SGT 2015               %
% last modified by Lujia Feng Tue Nov  3 23:57:38 SGT 2015                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xgrn = modspace.Xgrn;  
Lgrn = modspace.Lgrn;  
Bgrn = modspace.Bgrn;
Aeq  = modspace.Aeq;
beq  = modspace.beq;
lb   = modspace.lb;
ub   = modspace.ub;
x0   = modspace.x0;

% condense the smoothing matrix by removing rows of all zeros
sm     = GTdef_condense(modspace.sm);
sm_abs = GTdef_condense(modspace.sm_abs);
modspace.sm     = sm;
modspace.sm_abs = sm_abs;

C = []; d = [];
% add error & weight to the equations
if ~isempty(Xgrn)
    eqNum = size(Xgrn,1); % equation number
    for ii = 1:eqNum
        Xgrn(ii,:)  = Xgrn(ii,:).*pnt.coef(ii);
        pnt.obs(ii) = pnt.obs(ii).*pnt.coef(ii);
    end
    ind = find(~isnan(pnt.obs));		% exclude nan values
    C = [ C;Xgrn(ind,:) ]; d = [ d;pnt.obs(ind) ];
end

if ~isempty(Lgrn)
    eqNum = size(Lgrn,1);
    for ii = 1:eqNum
        Lgrn(ii,:)  = Lgrn(ii,:).*los.coef(ii);
        los.obs(ii) = los.obs(ii).*los.coef(ii);
    end
    ind = find(~isnan(los.obs));		% exclude nan values
    C = [ C;Lgrn(ind,:) ]; d = [ d;los.obs(ind) ];
end

if ~isempty(Bgrn)
    eqNum = size(Bgrn,1);
    for ii = 1:eqNum
        Bgrn(ii,:)    = Bgrn(ii,:).*bsl.coef(ii,1);
	bsl.obs(ii,1) = bsl.obs(ii,1).*bsl.coef(ii,1);
    end
    ind = find(~isnan(bsl.obs));		% exclude nan values
    C = [ C;Bgrn(ind,:) ]; d = [ d;bsl.obs(ind) ];
end

% add beta to the smoothing matrix
sm = sm.*beta;

% combine the green matrices and the smoothing matrix
d_sm = zeros(size(sm,1),1);
C = [ C;sm ]; d = [ d;d_sm ];

ind = find(isinf(lb));
Aeq_red = Aeq(ind,:); beq_red = beq(ind);

options = optimset('MaxIter',2000,'TolFun',1e-30);
[xx,resnorm] = lsqlin(C,d,[],[],Aeq_red,beq_red,lb,ub,x0,options);
%[xx,resnorm,res,exitflag,output] = lsqlin(C,d,[],[],Aeq_red,beq_red,lb,ub,x0);
fprintf(1,'resnorm = %-12.5e\n',resnorm);

modspace.xx = xx;
