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


% figure out where to plot
temp = varargin{1};
if strcmp(class(temp),'matlab.graphics.axis.Axes')
	% user has supplied a handlle to a axis; plot here
	plot_here = varargin{1};
	varargin{1} = [];
else
	% plot where you can
	plot_here = gca;
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

hold on

switch plot_what
case 'K_firing'
	plotx = o.filtertime_firing; 
	ploty = o.K_firing;
	xlabel('Filter Lag (s)')
	ylabel('Filter')
case 'firing_rate'
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
	ploty = o.inst_gain_firing;
	ploty(o.inst_gain_firing_err<inst_gain_min_r2) = NaN;
	ploty(ploty<0) = NaN;
	s = nanmean(o.stimulus,2);
	% filter the stimulus
	plotx = filter(ones(history_length,1),history_length,s);
	rm_this = isnan(plotx) | isnan(ploty);
	rm_this(1:o.filter_length*2) = true;
	plotx(rm_this) = []; ploty(rm_this) = [];

	xlabel(['Mean Stimulus in preceding ' oval(history_length) 'ms (V)'])
	ylabel('Inst. Gain (Hz/V)')
case 'inst_gain_firing'
	error
	ploty = o.inst_gain_firing;
otherwise 
	error('I dont know what to plot')
end 
	
switch plot_type
case 'trial-wise'
	plot_handles = plot(plot_here,plotx,ploty);
case 'sem'
	plot_handles = errorShade(plot_here,plotx,nanmean(ploty,2),sem(ploty'));
end	
