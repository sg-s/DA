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

end


GainAnalysis3(data,history_lengths,example_history_length,plothere,NaN*history_lengths);

% draw a line to indicate where we are on the history length plot
yy = get(plothere(4),'YLim');
plot(plothere(4),[example_history_length example_history_length],yy,'k.-')
