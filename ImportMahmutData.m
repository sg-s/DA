% Convert data from what Mahmut churns out into a usable format and append it to the global database. 

% specify where to load
i = 6;

% specify neuron
data(i).neuron = 'ab3';
data(i).odor = 'ethyl acetate';

% inherit time
data(i).time = data(i-1).time;

% load the PID
PID = squeeze(ab3.PID)';
t = 1e-4:1e-4:70;

PID2= [];
for j = 1:width(PID)
	PID2(:,j) = interp1(t,PID(:,j),data(i).time);

	% remove baseline
	PID2(:,j) = PID2(:,j) - mean(PID2(1:1000,j));


end

% nuke NaNs
PID2(isnan(PID2)) = 0;

data(i).PID = PID2;

% load spike times as is
data(i).spiketimes = squeeze(ab3.tSPKA);
data(i).spiketimeB = squeeze(ab3.tSPKB);

% bin and save firing rates
data(i).ORN = spiketimes2f(squeeze(ab3.tSPKA),t,3e-3);
data(i).ORNB = spiketimes2f(squeeze(ab3.tSPKB),t,3e-3);