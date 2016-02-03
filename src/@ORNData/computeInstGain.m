% computeInstGain.m
% computes inst. gain for class ORNdata
% 
% created by Srinivas Gorur-Shandilya at 6:33 , 23 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [o] = computeInstGain(o)

% check if data is there
if ~isempty(o.firing_rate) && ~isempty(o.firing_projected) && ~isempty(o.stimulus)
	% OK, we can find the firing rate inst. gain. findInstGain uses a hashed cache, so blindly call it, hoping that it intelligently uses the cache 

	% respect constrains on where to calculate this 
	if isempty(o.use_these_trials)
		utt = 1:o.n_trials;
	else
		utt = o.use_these_trials;		
	end
	if isempty(o.use_this_segment)
		uts = 1:length(o.stimulus);
	else
		uts = o.use_this_segment;
	end

	s = nanmean(o.stimulus(uts,utt),2);
	r = nanmean(o.firing_rate(uts,utt),2);
	p = nanmean(o.firing_projected(uts,utt),2);
	[~,o.inst_gain_firing,o.inst_gain_firing_err] = findInstGain(s,r,p,o.timescale_inst_gain);

else
	error('156hd not coded')
end