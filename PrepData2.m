% PrepData2.m
% PrepData loads the selected file, and prepares the data so that it is useful
% PrepData2 is built from code fragments scavenged from Carlotta's code, so that I m sure I'm doing the same thing
function [PID time f fs Valve ntrials] = PrepData2(filename)
% load file
load(filename)

% some parameters
deltat = 10^(-4); % sampling rate of data
win = 0.03;
sliding = 0.003;
memory = round(1/sliding); % confusingly, memory will be the time step in the processed data
shift_input = round(0.1/sliding);

% bin and average
[time PIDs sr f fs] = bin_traces(PID, stA, deltat, win, sliding);

% process PID
PID = mean(squeeze(PIDs(1,:,:)), 1);      
sdPID = std(squeeze(PIDs(1,:,:)), 1);  
clear PIDs

% find out when the trace starts and stops
[~, Valve] = sliding_average(squeeze(stim_signal(1,1,:)), deltat, win, sliding);
Valve(Valve <= 0.5) = 0;
Valve(Valve>0) = 1;
t_start=find(Valve,1,'first') + 5*memory; % we will throw away data before this; i.e, 5 seconds after stimulus starts
t_end=find(Valve,1,'last') - 5*memory; % we will throw away data after this; i.e, 5 seconds before stimulus ends

ntrialsPID = size(PID,2); % number of trials of PID
ntrials = size(stA,2);    % number of trials of neuron

% crop traces
time = time(t_start:t_end);   
PID = PID(t_start:t_end); 
sdPID = sdPID(t_start:t_end);    
f = f(t_start:t_end);
fs = fs(t_start:t_end);
Valve = Valve(t_start:t_end);
clear t_start t_end
