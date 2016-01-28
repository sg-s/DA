% plotInstGainVsMeanStim.m
% this is a subfunction called by ORNData/plot
% it makes a plot of inst. gain vs. the mean stimulus in some preceding history length
% this is meant to be called by the overloaded plot, but you can also call it directly
% 
% created by Srinivas Gorur-Shandilya at 1:48 , 28 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function plot_handles = plotInstGainVsMeanStim(plot_here,o,mean_stim_bins,normalise_gain,normalise_preceding_stim,history_length,inst_gain_min_r2)


stimulus = nanmean(o.stimulus,2);
inst_gain = o.inst_gain_firing;

% filter the stimulus using the history length
hl = floor(history_length/o.dt);
shat = filter(ones(hl,1),hl,stimulus);
shat(1:2*hl) = NaN;  % we can't estimate shat here 
inst_gain(1:2*hl) = NaN;

rm_this = o.inst_gain_firing_err<inst_gain_min_r2 | isnan(inst_gain) | isnan(shat) | inst_gain < 0;
shat(rm_this) = [];
inst_gain(rm_this) = [];

if normalise_preceding_stim
	shat = shat - min(shat);
	shat = shat/nanmax(shat);
end


if normalise_gain
	% fix the mean of the gain to be 1
	inst_gain = inst_gain/mean(inst_gain);
end

axes(plot_here)
plot_handles = plotPieceWiseLinear(shat,inst_gain,'nbins',mean_stim_bins);



