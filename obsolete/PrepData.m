% PrepData.m
% PrepData loads the selected file, and prepares the data so that it is useful
% PrepData
function [PID t f fs stim ntrials on_durations off_durations] = PrepData(filename)
load(filename)
win = 3e-2;
sliding = 3e-3;
ntrials = size(stA,2);
[t, PID, ~, f, fs] = bin_traces(PID, stA, 1e-4, win, sliding);

% bin stim
stim = squeeze(stim_signal);
stim_signal = squeeze(stim_signal);
if ~isvector(stim)
	% if there are multiple trials, average them
	stim = mean(stim);
	stim_signal = mean(stim_signal);
end
stim(stim>0.5) = 1;
stim(stim<0.5) = 0;
d = diff(stim);

ons = find(d>0.5);
offs = find(d<-0.5);
on_durations = (offs-ons)*1e-4;
ons(1) = [];
offs(end) = [];
off_durations = (ons-offs)*1e-4;

f = squeeze(f);
fs = squeeze(fs);
PID = (squeeze(mean(PID,2)));



stim=stim(1:floor(sliding/1e-4):end);
stim(end-9:end)=[]; % god knows why


% now chop off extra bits of the trace
start = round(find(stim_signal>0.5,1,'first')*1e-4+5);
disp('Keeping data from:')
disp(start)
start=start/sliding;
disp('Keeping data till:')
stop = round(find(stim_signal>0.5,1,'last')*1e-4-5);
disp(stop)
stop = stop/sliding;
PID = PID(start:stop);
PID = PID - min(PID);
PID = PID/max(PID);

f = f(start:stop);
fs = fs(start:stop);
t = t(start:stop) - min(t);
stim = stim(start:stop);

% conform to Carlotta's standards
f = f';
fs = fs';
