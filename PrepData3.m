% PrepData3.m
% PrepData loads the selected file, and prepares the data so that it is useful
% PrepData3 is built from my own code, as opposed to PrepData2, which uses Carlotta's code. Specifically, the firing rate computation is significantly different. 
function [PID time f Valve uncropped] = PrepData3(filename,varargin)
% load file
load(filename)

full_trace = 0;
if nargin  == 2
	full_trace = varargin{1};
end


% some parameters
deltat = 10^(-4); % sampling rate of data
win = 0.03;
sliding = 0.003;
memory = round(1/sliding); % confusingly, memory will be the time step in the processed data

% extract f
time=deltat:deltat:length(PID)*deltat;
f = mean(spiketimes2f(squeeze(stA),time,0.01,sliding));

% build time vector
time = sliding:sliding:sliding*length(f);

% filter PID
PID = mean(squeeze(PID(1,:,:)), 1); 
PID = filtfilt(ones(1,sliding/deltat)/(sliding/deltat),1,PID);     
% subsample
s = round(sliding/deltat);
PID = PID(s:s:end);

% process stim_signal
stim_signal = squeeze(stim_signal);
if isvector(stim_signal)
else
	stim_signal = mean(stim_signal);
end
Valve = stim_signal(s:s:end);
Valve(Valve <= 0.5) = 0;
Valve(Valve>0) = 1;

% find out when the trace starts and stops
if ~full_trace
	t_start=find(Valve,1,'first') + 5*memory; % we will throw away data before this; i.e, 5 seconds after stimulus starts
	t_end=find(Valve,1,'last') - 5*memory; % we will throw away data after this
else
	t_start = 1;
	t_end = length(Valve);
end

uncropped.PID = PID;
uncropped.Valve =Valve;
uncropped.f = f;
uncropped.time = time;
uncropped.t_start = t_start - 5*memory;
uncropped.t_end = t_end + 5*memory;

% crop traces
time = time(t_start:t_end);   
PID = PID(t_start:t_end);    
f = f(t_start:t_end);
Valve = Valve(t_start:t_end);
