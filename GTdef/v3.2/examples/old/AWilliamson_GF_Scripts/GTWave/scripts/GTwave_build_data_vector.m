function [dart,tsu] = GTwave_build_data_vector(dart,tsu)

cd Observed
tsu.obs=[];
tsu.obs_wgt=[];
tsu.obs_err=[];
ID=dart.id;
for i = 1: length(ID)
FNAME= strcat('win_Dart',num2str(dart.id(i)),'.txt');
tmp=load(FNAME);
tsu.obs=[tsu.obs ; tmp(:,1)];
tsu.obs_err=[tsu.obs_err ; ones(length(tmp(:,1)),1)*tsu.err(i)];
tsu.obs_wgt=[tsu.obs_wgt ; ones(length(tmp(:,1)),1)*tsu.wgt(i)];
end
tsu.coef=sqrt(tsu.obs_wgt)./((tsu.obs_err)/100);

cd ..
end

