function [dart] = GTwave_window_observed(dart)
WINDOW=dart.window;
cd Observed
dart.win_length=[];
dart.win_t=[];
for i = 1: length(WINDOW)
    
    FNAME= strcat('Dart',num2str(dart.id(i)),'.txt');
    [Wave]=Window_Observed_Tsunami(WINDOW(i,1),WINDOW(i,2), FNAME);
    dart.win_length=[dart.win_length; dart.id(i),length(Wave(:,1))];
    dart.win_t=[dart.win_t ; dart.id(i), Wave(1,2), (Wave(2,2)-Wave(1,2))];
end
cd ..
end

