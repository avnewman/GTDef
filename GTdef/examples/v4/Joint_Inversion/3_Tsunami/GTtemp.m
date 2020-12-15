clc; clear all;

filename='Windowing_DART_input.txt';
[dart] = GTwave_open(filename);


%% Window observed data
tic
fprintf('%s\n','......Windowing Observed Tsunami')
GTwave_window_observed(dart);
fprintf('%s\n','......Windowing Observed Tsunami COMPLETE')
toc
%% Window Green's Function Data
tic
fprintf('%s\n','......Windowing Greens Functions')
GTwave_window_greens_fun(dart)
fprintf('%s\n','......Windowing Greens Functions COMPLETE')
toc

%% DATA VECTOR (DART)
% lines up the detieded, windowed tsunami data
