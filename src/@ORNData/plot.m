% plot.m
% this method makes plots for class ORNData, overloading the built-in plot function
% usage:
% plot(where, obj, what)
% plot(obj, what)
% 
% created by Srinivas Gorur-Shandilya at 6:26 , 23 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [plot_handles] = plot(varargin)

% defaults
plot_type = 'trial-wise'; % or could be 'sem'
history_length = 500; % in units of the data
inst_gain_min_r2 = .8;
plot_here = [];
nbins = 50;
mean_stim_bins = 15;
normalise_gain = true;

hl_min = .2;
hl_max = 10;
history_lengths = [logspace(log10(hl_min),log10(.5),15) logspace(log10(.5),log10(10),15)];
history_lengths = unique(history_lengths);


% figure out where to plot
temp = varargin{1};
if strcmp(class(temp),'matlab.graphics.axis.Axes')
	% user has supplied a handlle to a axis; plot here
	plot_here = varargin{1};
	varargin{1} = [];
else
end

% defensive programming
assert(length(varargin)>1,'Not enough input arguments.')
o = varargin{1};
assert(isa(o,'ORNData'),'2nd argument should be an object from the ORNData class')
plot_what = varargin{2};
% handle any additional arguments (they have to be options in name value syntax )
varargin{1,2} = [];
if length(varargin) 
	if iseven(length(varargin))
		for ii = 1:2:length(varargin)-1
			temp = varargin{ii};
        	if ischar(temp)
        		eval(strcat(temp,'=varargin{ii+1};'));
        	end
    	end
	else
    	error('Inputs need to be name value pairs')
    end
end


switch plot_what

case 'binned_gain_analysis'
	if isempty(plot_here)
		figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
		clear plot_here
		plot_here(1) = subplot(1,2,1); hold on
		plot_here(2) = subplot(1,2,2); hold on
	end

	time = o.dt*(1:length(o.stimulus));
	stimulus = nanmean(o.stimulus,2);
	prediction = nanmean(o.firing_projected,2);
	response = nanmean(o.firing_rate,2);

	% throw out first 5 seconds
	time = time(5e3:end);
	stimulus = stimulus(5e3:end);
	response = response(5e3:end);
	prediction = prediction(5e3:end);

	if normalise_gain
		% fix the gain to be exactly 1
		x = prediction(:);
		y = response(:);
		rm_this = isnan(x) | isnan(y);
		x(rm_this) = [];
		y(rm_this) = [];
		temp = fit(x,y,'poly1');
		prediction = prediction*temp.p1;
	end

	% remove trend in stimulus
	temp = fit(time(:),stimulus(:),'poly2');
	stimulus = stimulus - temp(time) + mean(stimulus);

	rho = NaN*history_lengths;


	for j = 1:length(history_lengths)
		hl = floor(history_lengths(j)/o.dt);
		shat = filter(ones(hl,1),hl,stimulus);
		shat(1:hl) = [];  % we can't estimate shat here 
		this_prediction = prediction; this_prediction(1:hl) = [];
		this_response = response; this_response(1:hl) = [];


		l = labelByPercentile(shat,mean_stim_bins);
		x = NaN(mean_stim_bins,1);
		y = NaN(mean_stim_bins,1);
		ye = NaN(mean_stim_bins,1);

		for bi = 1:mean_stim_bins
			x(bi) = nanmean(shat(l==bi));
			p = (this_prediction(l==bi));
			r = (this_response(l==bi));
			rm_this = isnan(p) | isnan(r);
			p(rm_this) = []; r(rm_this) = [];
			[temp] = pca([p r]);
			y(bi) = temp(2,1)/temp(1,1);
			
		end
		rho(j) = spear(x,y);

		if j == 15
			plot(plot_here(1),x,y,'k+')
		end
	end
	plot(plot_here(2),history_lengths,rho,'k+')
	set(plot_here(2),'XScale','log')

case 'clark_gain_analysis'
	if isempty(plot_here)
		figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
		plot_here(1) = subplot(1,2,1); hold on
		plot_here(2) = subplot(1,2,2); hold on
	end

	time = o.dt*(1:length(o.stimulus));
	stimulus = nanmean(o.stimulus,2);
	prediction = nanmean(o.firing_projected,2);
	response = nanmean(o.firing_rate,2);

	% throw out first 5 seconds
	time = time(5e3:end);
	stimulus = stimulus(5e3:end);
	response = response(5e3:end);
	prediction = prediction(5e3:end);

	% fix the gain to be exactly 1
	x = prediction(:);
	y = response(:);
	rm_this = isnan(x) | isnan(y);
	x(rm_this) = [];
	y(rm_this) = [];
	temp = fit(x,y,'poly1');
	prediction = prediction*temp.p1;

	% remove trend in stimulus
	temp = fit(time(:),stimulus(:),'poly2');
	stimulus = stimulus - temp(time) + mean(stimulus);

	ph(3) = plot_here(1); ph(4) = plot_here(2);
	gainAnalysisWrapper('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph,'engine',@gainAnalysis,'use_cache',true,'example_history_length',history_lengths(15));

case 'K_firing'
	if isempty(plot_here)
		figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[500 500]); hold on
		plot_here = gca;
		hold on
	end
	% show the filters
	plot_handles = plotFilters(plot_here,o.filtertime_firing,o.K_firing,plot_type);

case 'LN_firing'
	if isempty(plot_here)
		figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
		plot_here(1) = subplot(1,2,1); hold on
		plot_here(2) = subplot(1,2,2); hold on
	end

	% show the filters
	plot_handles(1).h = plotFilters(plot_here(1),o.filtertime_firing,o.K_firing,plot_type);

	% show the nonlinearity
	plot_handles(2).h = plotNonlinearity(plot_here(2),o.firing_projected,o.firing_rate,plot_type,nbins);
	xlabel(plot_here(2),'Projected Stimulus (V)')
	ylabel(plot_here(2),'Response (Hz)')

case 'firing_rate'
	if isempty(plot_here)
		figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
		plot_here = gca;
	end
	time = o.dt*(1:length(o.stimulus));
	ax = plotyy(time,nanmean(o.firing_rate,2),time,nanmean(o.firing_projected,2));
	set(ax(1),'YLim',[0 max(max(o.firing_rate))])
	set(ax(2),'YLim',[nanmean(min(o.firing_projected)) nanmean(max(o.firing_projected))])

case 'inst_gain_firing -dist'
	temp = o.inst_gain_firing;
	temp(o.inst_gain_firing_err<inst_gain_min_r2) = [];
	temp(temp<0) = [];
	[ploty,x] = histcounts(temp,100);
	ploty = ploty/sum(ploty);
	plotx = diff(x) + x(1:end-1);
	xlabel('Instantaneous Gain (Hz/V')
	ylabel('Probability')
case 'inst_gain_firing -vs_stim'
	if isempty(plot_here)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[600 600]); hold on
		plot_here = gca;
	end
	ploty = o.inst_gain_firing;
	ploty(o.inst_gain_firing_err<inst_gain_min_r2) = NaN;
	ploty(ploty<0) = NaN;
	s = nanmean(o.stimulus,2);
	% filter the stimulus
	plotx = filter(ones(history_length,1),history_length,s);
	rm_this = isnan(plotx) | isnan(ploty);
	% censor some transient data
	rm_this(1:o.filter_length*2) = true;
	rm_this(1:history_length*2) = true;
	plotx(rm_this) = []; ploty(rm_this) = [];

	l = labelByPercentile(plotx,mean_stim_bins);

	if mean_stim_bins < 10
		for i = 1:mean_stim_bins
			errorbar(plot_here,mean(plotx(l==i)),mean(ploty(l==i)),sem(ploty(l==i)),'k');
		end
	else
		x = NaN(mean_stim_bins,1);
		y = x; ye = x;
		for i = 1:mean_stim_bins
			x(i) = mean(plotx(l==i));
			y(i) = mean(ploty(l==i));
			ye(i) = sem(ploty(l==i));
		end
		errorShade(x,y,ye);
	end


	xlabel(['Mean Stimulus in preceding ' oval(history_length) 'ms (V)'])
	ylabel('Inst. Gain (Hz/V)')
case 'inst_gain_firing'
	error
	ploty = o.inst_gain_firing;
otherwise 
	error('I dont know what to plot')
end 

	


function [plot_handles] = plotFilters(plot_here,filtertime,K,plot_type)
	switch plot_type
	case 'trial-wise'
		plot_handles = plot(plot_here,filtertime,K);
	case 'sem'
		plot_handles = errorShade(plot_here,filtertime,nanmean(K,2),sem(K'));
	end	
	xlabel(plot_here,'Filtertime (s)')
	ylabel(plot_here,'Filter')
end

function [plot_handles] = plotNonlinearity(plot_here,pred,resp,plot_type,nbins)
	switch plot_type
	case 'trial-wise'
		for i = 1:width(pred)
			[~,data] = plotPieceWiseLinear(pred(:,i),resp(:,i),'nbins',nbins,'make_plot',false);
			plot_handles(i) = plot(plot_here,data.x,data.y);
		end
	case 'sem'
		for i = 1:width(pred)
			[~,data(i)] = plotPieceWiseLinear(pred(:,i),resp(:,i),'nbins',nbins,'make_plot',false);
		end
		plot_handles = errorShade(plot_here,nanmean([data.x],2),nanmean([data.y],2),nanmean([data.ye],2));
	end	
	
end

end
