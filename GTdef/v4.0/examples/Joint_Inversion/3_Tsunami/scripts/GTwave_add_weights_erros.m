function [modspace,tsu] = GTwave_add_weights_erros(modspace,tsu)

Tgrn=modspace.Tgrn;
%d=modspace.d;
%C=modspace.C;
C=[];
d=[];

if ~isempty(Tgrn)
    eqNum = size(Tgrn,1); % equation number
    for ii = 1:eqNum
        mod_Tgrn(ii,:)  = Tgrn(ii,:).*tsu.coef(ii);
        mod_tsu.obs(ii) = tsu.obs(ii).*tsu.coef(ii);
    end
    ind = find(~isnan(tsu.obs));                % exclude nan values
    
    C = [ C;mod_Tgrn(ind,:) ]; 
    d = [ d;mod_tsu.obs(ind) ];
end

modspace.C=C;   modspace.d=d; 
end

