%% findRhoForHistoryLengths
% finds rho for gain vs. mean stim plots for various history lengths
% this version works only when ons and offs are defined
%
function rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths)

rho = NaN(length(history_lengths),1);

for i = 1:length(history_lengths)
	% filter the stimulus
	shat = computeSmoothedStimulus(stim,history_lengths(i));

	% find the mean stimulus in the windows
	mean_stim = findMeanInWindows(ons,offs,shat);

	rho(i) = spear(mean_stim,gain);
end
