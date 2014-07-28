% Convert data from what Mahmut churns out into a usable format and append it to the global database. 

% specify where to load
i = 6;

% specify neuron
data(i).neuron = 'ab2';
data(i).odor = 'ethyl acetate';

% inherit time
data(i).time = data(i-1).time;

% load the PID
PID = squeeze(ab2.PID)';
t = 1e-4:1e-4:70;

PID2= [];
for j = 1:width(PID)
	PID2(:,j) = interp1(t,PID(:,j),data(i).time);
end

data(i).PID = PID2;

% load spike times as is
data(i).spiketimes = ab2.tSPKA;
data(i).spiketimeB = ab2.tSPKB;

% bin and save firing rates
data(i).ORN = spiketimes2f(squeeze(ab2.tSPKA),t,3e-3);
data(i).ORNB = spiketimes2f(squeeze(ab2.tSPKB),t,3e-3);