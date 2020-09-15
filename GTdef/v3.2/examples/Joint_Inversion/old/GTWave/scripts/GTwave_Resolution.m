function [modspace] = GTwave_Resolution(modspace,fault)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
lb=modspace.lb;
nrealdata=modspace.nrealdata;
C=modspace.C;       d=modspace.d;
x0=modspace.x0;     sm=modspace.sm;
 Ns=fault.Ns; ds_len=fault.len/Ns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 ix=find(isfinite(lb'));
 id=[[1:nrealdata]';ix+nrealdata];
 C0=C(id,:);                        % operator matrix (includes smoothing)
 d=d(1:nrealdata);
 x0=x0(ix);
 
 % Develop weighted and damped over-determined generalized inverse matrix, Gg
 % for the model resolution matrix  (See Menke, 2012 -- Chapter 4)
 G=C0(1:nrealdata,:);                 % strip off the smoothing part of the matrix
     
 e2I=sm.^2;                          % weighted smoothing matrix
 Gg=inv(G'*G+e2I)*G';                % overdetermined.

 R=Gg*G ;                            % Model Resolution Matrix 
 N=G*Gg ;                            % Data Resolution Matrix
 
 R_ds=diag(R);
 R_ss=zeros(length(R_ds),1);   % TEMPORARY FIX WHEN SS AND TS ARE ADDED TO INVERSION
 R_ts=zeros(length(R_ds),1);
 R_out=[R_ss,R_ds,R_ts];
 
 % resolution spread parameter
 Li=ds_len;  % km PATCH SIZE for a 30x30 subfault size
 r_i=Li./sqrt(R_out(:,2)); % resolution spread parameter (in km)
 
 modspace.R=R;
 modspace.ri=r_i;
 modspace.R_ds=R_ds;    modspace.R_ss=R_ss;    modspace.R_ts=R_ts;
 

end

