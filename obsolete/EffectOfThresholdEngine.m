% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = EffectOfThresholdEngine(data,p,plothere)


% clear axes
for i = 1:length(plothere)
	cla(plothere(i));
end

% plot the original filter
plot(plothere(1),data.filtertime,data.K,'k')

% threshold the prediction
ThresholdedPred = data.LinearFit;
ThresholdedPred(ThresholdedPred<p.x0) = 0;

% show the thresholded prediction vs. the original
plot(plothere(2),data.time,data.LinearFit,'k');
plot(plothere(2),data.time,ThresholdedPred,'r');
legend(plothere(2),{'Linear Pred','After Threshold'})


% also calculate the filter and show it
K2 = FindBestFilter(data.stimulus,ThresholdedPred,[],'filter_length=199;');
ycap = convolve(data.time,data.stimulus,K2,data.filtertime) + mean(ThresholdedPred);
plot(plothere(1),data.filtertime,K2,'r');

% do gain analysis using ycap vs. LinearFit

% assemble data
s = 1; % when we start for the gain analysis
z = length(data.LinearFit); % where we end 
x.response = ThresholdedPred;
x.prediction = ycap;
x.stimulus = data.stimulus;
x.time = data.time(:);
x.filter_length = 200; 

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

ga_plot = zeros(4,1);
ga_plot(3:4) = plothere(3:4);
[~,~,~,~,~,example_plot]=GainAnalysis4(x,history_lengths,example_history_length,ga_plot,NaN(length(history_lengths),2),p.frac);

% draw a line to indicate where we are on the history length plot
yy = get(plothere(4),'YLim');
plot(plothere(4),[example_history_length example_history_length],yy,'k.-')
set(plothere(4),'XLim',[3e-3 max(history_lengths)]);
