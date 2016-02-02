% ORNDATA/PLOT.m
% this method makes plots for class ORNData, overloading the built-in plot function
% usage:
%
% plot(obj, what)
% 
% this is the simplest usage of ORNData/plot. obj is a a ORNData object, and what is a string that is one of the following:
% "gain_analyis", "LN", "Filters"
%
% these strings can also have the following modifiers, allowing you specify precisely what you want to plot. e.g.:
% "gain_analysis.clark", "gain_analysis.binned", gain_analysis.inst"
% 
% Here is a incomplete list of modifiers:
% .clark 		Clark et. al style gain analysis
% .binned 		binned gain analysis
% .inst 		inst. gain analysis
% .firing 		use firing rates
% .LFP 			use LFP traces 
%
%
% plot(obj, where, what)
% 
% you can also specify where you want to plot something by supplying axes handles as the 2nd argument. make sure you provide enough handles to axes for the type of plot you want to make. Some plots (like "gain_analysis" or "LN" need two axes)
% 
% plot(obj, where, what,'Name',value...)
%
% You can override defaults, or provide more options using extensible Name-value syntax. Here is an incomplete list of some attributes you can set/modify:
% attribute 				{default}, other options
% plot_type    				{trial-wise}, sem, std, mean
% history_length 			{5}
% inst_gain_min_r2 			{.8};
% plot_here 
% nbins 					{50}
% mean_stim_bins 			{15}
% normalise_gain 			{true}, false
% normalise_preceding_stim 	{true}, false
% data_filter
% 
% created by Srinivas Gorur-Shandilya at 6:26 , 23 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function plot_handles = plot(varargin)

% defaults
plot_type = 'trial-wise'; % or could be 'sem'
history_length = .5; % in seconds
inst_gain_min_r2 = .8;
plot_here = [];
nbins = 50;
mean_stim_bins = 15;
normalise_gain = true;
normalise_preceding_stim = true;

hl_min = .2;
hl_max = 10;
history_lengths = [logspace(log10(hl_min),log10(.5),15) logspace(log10(.5),log10(hl_max),15)];
history_lengths = unique(history_lengths);


% now figure out what to do based on the string
if ~isempty(strfind(plot_what,'gain_analysis'))
	% OK, we're doing a again analysis. all gain analysis involve two subplots, so get them ready
	if isempty(plot_here)
		figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
		clear plot_here
		plot_here(1) = subplot(1,2,1); hold on
		plot_here(2) = subplot(1,2,2); hold on
	end

	% we default to a binned analysis if not specified 
	if ~isempty(strfind(plot_what,'binned'))
	elseif ~isempty(strfind(plot_what,'clark'))
		error('not coded 85-clark')
		return
	elseif ~isempty(strfind(plot_what,'inst'))
		plot_handles(1).h = plotInstGainVsMeanStim(plot_here(1),o,mean_stim_bins,normalise_gain,normalise_preceding_stim,history_length,inst_gain_min_r2);
		return
	else
	end
	
	if ~isempty(strfind(plot_what,'LFP'))	
		% make a plot of binned gain vs. the mean stim in that bin
		plot_handles(1) = plotBinnedGainVsMeanStim(plot_here(1),o,mean_stim_bins,normalise_gain,normalise_preceding_stim,history_length,true);

		% make a plot of rho vs. history lengths
		plot_handles(2) = plotSpearmanRhoVsHistoryLengths(plot_here(2),o,mean_stim_bins,normalise_gain,normalise_preceding_stim,history_lengths,true);
	else
		% make a plot of binned gain vs. the mean stim in that bin
		plot_handles(1) = plotBinnedGainVsMeanStim(plot_here(1),o,mean_stim_bins,normalise_gain,normalise_preceding_stim,history_length,false);

		% make a plot of rho vs. history lengths
		plot_handles(2) = plotSpearmanRhoVsHistoryLengths(plot_here(2),o,mean_stim_bins,normalise_gain,normalise_preceding_stim,history_lengths,false);
	end

	% cosmetic fixes
	ylabel(plot_here(1),'Gain (norm)')
	if normalise_preceding_stim
		xlabel(plot_here(1),['Mean stim in ' char(10) 'preceding ' oval(history_length) 's (norm)'])
	else
		xlabel(plot_here(1),['Mean stim in ' char(10) 'preceding ' oval(history_length) 's (V)'])
	end
	xlabel(plot_here(2),'History Length (s)')
	ylabel(plot_here(2),'\rho')

	set(plot_here(2),'XScale','log','YLim',[-1 1])





% case 'clark_gain_analysis'
% 	if isempty(plot_here)
% 		figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% 		plot_here(1) = subplot(1,2,1); hold on
% 		plot_here(2) = subplot(1,2,2); hold on
% 	end

% 	time = o.dt*(1:length(o.stimulus));
% 	stimulus = nanmean(o.stimulus,2);
% 	prediction = nanmean(o.firing_projected,2);
% 	response = nanmean(o.firing_rate,2);

% 	% throw out first 5 seconds
% 	time = time(5e3:end);
% 	stimulus = stimulus(5e3:end);
% 	response = response(5e3:end);
% 	prediction = prediction(5e3:end);

% 	% fix the gain to be exactly 1
% 	x = prediction(:);
% 	y = response(:);
% 	rm_this = isnan(x) | isnan(y);
% 	x(rm_this) = [];
% 	y(rm_this) = [];
% 	temp = fit(x,y,'poly1');
% 	prediction = prediction*temp.p1;

% 	% remove trend in stimulus
% 	temp = fit(time(:),stimulus(:),'poly2');
% 	stimulus = stimulus - temp(time) + mean(stimulus);

% 	ph(3) = plot_here(1); ph(4) = plot_here(2);
% 	gainAnalysisWrapper('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph,'engine',@gainAnalysis,'use_cache',true,'example_history_length',history_lengths(15));




% case 'inst_gain_firing -vs_stim'
% 	if isempty(plot_here)
% 		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[600 600]); hold on
% 		plot_here = gca;
% 	end
% 	ploty = o.inst_gain_firing;
% 	ploty(o.inst_gain_firing_err<inst_gain_min_r2) = NaN;
% 	ploty(ploty<0) = NaN;
% 	s = nanmean(o.stimulus,2);
% 	% filter the stimulus
% 	plotx = filter(ones(history_length,1),history_length,s);
% 	rm_this = isnan(plotx) | isnan(ploty);
% 	% censor some transient data
% 	rm_this(1:o.filter_length*2) = true;
% 	rm_this(1:history_length*2) = true;
% 	plotx(rm_this) = []; ploty(rm_this) = [];

% 	l = labelByPercentile(plotx,mean_stim_bins);

% 	if mean_stim_bins < 10
% 		for i = 1:mean_stim_bins
% 			errorbar(plot_here,mean(plotx(l==i)),mean(ploty(l==i)),sem(ploty(l==i)),'k');
% 		end
% 	else
% 		x = NaN(mean_stim_bins,1);
% 		y = x; ye = x;
% 		for i = 1:mean_stim_bins
% 			x(i) = mean(plotx(l==i));
% 			y(i) = mean(ploty(l==i));
% 			ye(i) = sem(ploty(l==i));
% 		end
% 		errorShade(x,y,ye);
% 	end


% 	xlabel(['Mean Stimulus in preceding ' oval(history_length) 'ms (V)'])
% 	ylabel('Inst. Gain (Hz/V)')
% case 'inst_gain_firing'
% 	error
% 	ploty = o.inst_gain_firing;
% otherwise 
% 	error('I dont know what to plot. The allowed values are: "binned_gain_analysis","LN_firing","inst_gain_analysis",""')
% end 

	
