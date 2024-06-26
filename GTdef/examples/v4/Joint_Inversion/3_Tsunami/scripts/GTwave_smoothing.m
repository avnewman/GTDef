function [sm] = GTwave_smoothing(dd,ds,Nd,Ns)
% THIS IS FROM L.Feng's GTdef_sm2d_free via GTdef
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                            GTdef_sm2d_free                              %
 % Calculate the second-order derivatives of the slips of the subfaults    %
 % over the fault surface using 3-point central finite-difference          %
 % approximation                                                           %
 % REFERENCE:                                                              %
 % Jonsson, et al. (2002), BSSA, 92(4), 1377-1389 eq (A1)                  %
 % Two-dimensional, 2nd-order, finite-difference central approximation sum %
 % (ith row and jth colunm)                                                %
 %  Si-1,j - 2Si,j + Si+1,j       Si,j-1 - 2Si,j + Si,j+1                  %
 % -------------------------  +  -------------------------                 %
 %          dd*dd                         ds*ds                            %
 %                                                                         %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = Ns*Nd;                     % total patch/slip number
  mx = 1/ds^2;                    % the patches on left and right
  my = 1/dd^2;                    % the patches on top and bottom
  mxy = -2*mx-2*my;               % the patch in the center
  mxy_free = -2*mx-my;            % the free-surface patch in the center
  
  if nn==1, error('GTdef_sm2d_free ERROR: only one patch. No need to smooth!'); end
  
  % more layers
  d0 = mxy*ones(Nd,Ns);           % diagonals = 0
  d0(1,:) = mxy_free;
  d0 = reshape(d0,[],1);
  dn = mx*ones(nn,1);             % diagonal +/- Nd
  % sub-diagonal = diagonal - 1
  t1 = my*ones(Nd-1,Ns);  
  t2 = zeros(1,Ns);
  t3 = [ t1;t2 ];
  d10 = reshape(t3,[],1);
  % super-diagonal = diagonal + 1
  d01 = d10(end:-1:1);
  B = [ dn d10 d0 d01 dn ];       % diagonal columns
  ind = [ -Nd -1 0 1 Nd ];        % index for diagonal columns
  sm = spdiags(B,ind,nn,nn);
end

