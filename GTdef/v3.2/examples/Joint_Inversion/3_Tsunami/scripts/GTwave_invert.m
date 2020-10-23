function [modspace] = GTwave_invert(modspace,fault,kappa,beta)

C=modspace.C;
d=modspace.d;


num_patch=length(C);
lb=zeros(1,num_patch);
ub=ones(1,num_patch);  
x0=1*ones(1,num_patch);

slip_num = num_patch; % number of patches	
xunit = x0; xunit(xunit~=0) = 1;		% set unit slips for nonzero initial slips
Aineq = zeros(slip_num*2,slip_num); bineq = zeros(slip_num*2,1);   % useful only for fault type 3 & 4
Aeq   = zeros(slip_num,slip_num);  beq   = zeros(slip_num,1);


Nd=fault.Nd; Ns=fault.Ns;
ds=fault.len/Ns;
dd=ds;  % FIX THIS
[ sm ] = GTdef_sm2d_free(dd,ds,Nd,Ns);


[ sm ] = GTwave_condense(sm);
sm = sm.*beta;
d_sm = zeros(size(sm,1),1);  % append zeros to data vector
%%

ind_fixed = find(lb==-Inf);	% index for fixed slips 
sm(ind_fixed) = 0;		% don't do smoothing for them

nrealdata=size(C,1);
modspace.nrealdata=nrealdata;
C = [ C;sm ]; d = [ d';d_sm ];  % and add smothing matrix to GF matrix

modspace.lb=lb;  modspace.ub=ub; modspace.x0=x0;
modspace.C=C; modspace.d=d;

%% WAVEFORM INVERSION

options = optimset('MaxIter',2000,'TolFun',1e-30);
d=double(d);
C=full(C);
%[xx,resnorm] = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0,options);
[xx,resnorm] = lsqlin(C,d,[],[],[],[],lb,ub,x0,options);

modspace.xx=xx;
modspace.resnorm=resnorm;
modspace.sm=sm;
end

