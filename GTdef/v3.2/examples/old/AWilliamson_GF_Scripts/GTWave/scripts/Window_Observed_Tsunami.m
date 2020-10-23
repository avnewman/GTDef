function [Wave] = Window_Observed_Tsunami(starttime,endtime, filename)
% Window an observed DART gauge to between starttime and endtime (both in
% seconds)


Wave=load(filename);
t=Wave(:,2);

loc=find(t >= starttime & t <= endtime);

WIN_t=t(loc);  WIN_Observed=Wave(loc,1);
Wave=[double(WIN_Observed),double(WIN_t),];
file_out= strcat('win_', filename);

save(file_out,'Wave','-ascii')

end

