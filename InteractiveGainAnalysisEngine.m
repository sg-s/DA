% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = InteractiveGainAnalysisEngine(data,p,plothere)


% get history lengths
example_history_length = p.t_h;
history_lengths = data.history_lengths;

if ~ismember(example_history_length,history_lengths)
	% pick the closest multiple of 3ms
	example_history_length= (floor(example_history_length/(3e-3))*3e-3);
	history_lengths = sort([history_lengths example_history_length]);
else
	% prevent this from growing out of control
	history_lengths = data.history_lengths;
end


[~,~,~,~,~,example_plot]=GainAnalysis4(data,history_lengths,example_history_length,plothere,NaN*history_lengths,p.frac);

% draw a line to indicate where we are on the history length plot
yy = get(plothere(4),'YLim');
plot(plothere(4),[example_history_length example_history_length],yy,'k.-')
set(plothere(4),'XLim',[3e-3 max(history_lengths)]);

% get the time we are looking at in the time series
t = p.t;
yy = get(plothere(1),'YLim');
plot(plothere(1),[t t],yy,'b-')

yy = get(plothere(2),'YLim');
plot(plothere(2),[t t],yy,'b-')

% find the points in the example plot closest to the point of time we are at. 
fp_high = example_plot.fp_high;
fp_low = example_plot.fp_low;
f_high = example_plot.f_high;
f_low = example_plot.f_low;

if min(abs(example_plot.t_low - t)) > min(abs(example_plot.t_high - t))
	% the closest point is a high stimulus point
	[~,loc]=min(abs(example_plot.t_high - t));
	scatter(plothere(3),fp_high(loc),f_high(loc),64,'filled');
else
	% the closest point is a low stimulus point
	[~,loc]=min(abs(example_plot.t_low - t));
	scatter(plothere(3),fp_low(loc),f_low(loc),64,'filled');
end