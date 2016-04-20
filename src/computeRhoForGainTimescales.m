%% computeRhoForGainTimescales.m
% computes the Spearman correlation for a sequence of whiffs, for various history lengths
% usage:
% rho = computeRhoForGainTimescales(gain,t,S)
% where
% gain is a vector of gains for various whiffs
% and t is a vector as long as gain indicating the whiff times (in the same units as S)
% and S is a stimulus vector 

function rho = computeRhoForGainTimescales(gain,t,S,tau_gain)

rho = NaN*tau_gain;

for i = 1:length(tau_gain)
	this_tau = round(tau_gain(i));

	% filter the stimulus
	shat = computeSmoothedStimulus(S,this_tau);

	% find the stimulus at whiff times
	rho(i) = spear(vectorise(shat(t)),vectorise(gain));

end


