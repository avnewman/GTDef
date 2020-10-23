function [dart] = GTwave_window_greens_fun(dart)
WINDOW=dart.window;
cd GFs
k=1;
while 1
M=dir(sprintf('%s%03.f%s','SS_',k,'_*.txt'));
    if ~isempty(M)
        for i = 1:length(M)
            Window_Observed_Tsunami(WINDOW(i,1),WINDOW(i,2), M(i).name); 
        end
    k=k+1;
    else
        break
    end
end
 cd ..
end

