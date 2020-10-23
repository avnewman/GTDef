function [dart,tsu,modspace] = GTwave_build_gf_matrix(dart,tsu,modspace)

cd GFs
modspace.Tgrn=[];
Tgrn=[];%modspace.Tgrn;

k=1;
while 1
M=dir(sprintf('%s%03.f%s','win_SS_',k,'_*.txt'));
    if ~isempty(M)
        
        
        r=[];
        for i = 1:length(M)
            tmp=load(M(i).name);
            r=[r;tmp(:,1)];
         end
         Tgrn=[Tgrn,r];
    k=k+1;
    else
        break
    end
end
modspace.Tgrn=Tgrn;


cd ..
end

