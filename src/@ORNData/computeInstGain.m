% computeInstGain.m
% computes inst. gain for class ORNdata
% 
% created by Srinivas Gorur-Shandilya at 6:33 , 23 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [obj] = computeInstGain(obj)

% check if data is there
if ~isempty(obj.firing_rate) && ~isempty(obj.firing_projected) && ~isempty(obj.stimulus)
	% OK, we can find the firing rate inst. gain. findInstGain uses a hashed cache, so blindly call it, hoping that it intelligently uses the cache 
	s = nanmean(obj.stimulus,2);
	r = nanmean(obj.firing_rate,2);
	p = nanmean(obj.firing_projected,2);
	[~,obj.inst_gain_firing,obj.inst_gain_firing_err] = findInstGain(s,r,p,obj.timescale_inst_gain);

else
	error('156hd not coded')
end