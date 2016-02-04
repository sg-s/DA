% computeInstGain.m
% computes inst. gain for class ORNdata
% 
% created by Srinivas Gorur-Shandilya at 6:33 , 23 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [o] = computeInstGain(o,useLN)

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

if nargin < 2
	useLN = false;
end

s = nanmean(o.stimulus(uts,utt),2);

% check if data is there
if ~isempty(o.firing_rate) && ~isempty(o.firing_projected) && ~isempty(o.stimulus)
	% OK, we can find the firing rate inst. gain. findInstGain uses a hashed cache, so blindly call it, hoping that it intelligently uses the cache 
	r = nanmean(o.firing_rate(uts,utt),2);
	p = nanmean(o.firing_projected(uts,utt),2);
	if useLN
		temp1 = p(:); temp2 = r(:);
		rm_this = isnan(temp1) | isnan(temp2);
		temp1(rm_this) = []; temp2(rm_this) = [];
		ft = fittype('hillFit(x,A,k,n,x_offset)');
		ff = fit(temp1,temp2,ft,'StartPoint',[max(temp2) mean(temp1) 1 0],'Upper',[1e3 max(temp1) 10 Inf],'Lower',[min(temp2)/2 0 0 -Inf]);
		p = ff(p);
	end

	[~,o.inst_gain_firing,o.inst_gain_firing_err] = findInstGain(s,r,p,o.timescale_inst_gain);

end
if ~isempty(o.LFP) && ~isempty(o.LFP_projected) && ~isempty(o.stimulus)
	r = nanmean(o.LFP(uts,utt),2);
	p = nanmean(o.LFP_projected(uts,utt),2);
	useLN = true;
	if useLN
		temp1 = p(:); temp2 = r(:);
		rm_this = isnan(temp1) | isnan(temp2);
		temp1(rm_this) = []; temp2(rm_this) = [];
		ff = fit(temp1,temp2,'poly3');
		p = ff(p);
	end

	[~,o.inst_gain_LFP,o.inst_gain_LFP_err] = findInstGain(s,r,p,o.timescale_inst_gain);
end